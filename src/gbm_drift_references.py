"""Independent numerical minimax diagnostics for a compact physical-drift set.

The interval engine is not imported. Piecewise-affine truncated-lognormal
moments from the existing independent static-hedge reference define all risks.
SciPy proposes an active set, then mpmath refines its KKT equations at 160
significant digits. These values diagnose the interval certificate; they are
not interval proofs or exact feasibility claims.
"""
from fractions import Fraction as F
import itertools

import mpmath as mp
import numpy as np
from scipy.optimize import minimize

import gbm_hedge_references as static


def number(value):
    return static.number(value)


def models(contract, sigma, nodes):
    return [static.parameters(contract, sigma, F(node)) for node in nodes]


def risk(model, holdings):
    return static.fixed_risk(model, holdings)


def finite_grid_optimum(contract, sigma, nodes, kappa, digits=160):
    """Refine an independently proposed finite-grid minimax active set.

    Node risks are divided by a common scale and the funding inequality by
    spot solely to condition the numerical reference. A 1e-90 scaled residual
    threshold is reported, never treated as mathematical exactness.
    """
    nodes=tuple(F(x) for x in nodes)
    with mp.workdps(digits):
        laws=models(contract,sigma,nodes)
        coeff=[static.coefficients(model) for model in laws]
        scale=max(mp.mpf(1),max(x[2] for x in coeff))
        data=[(M/scale,h/scale,c/scale) for M,h,c in coeff]
        costs=(mp.mpf(1),laws[0]['ck']/laws[0]['s'])
        budget=number(kappa)*laws[0]['premium']/laws[0]['s']
        normals=((-mp.mpf(1),mp.mpf(0)),(mp.mpf(0),-mp.mpf(1)),
                 (mp.mpf(1),mp.mpf(0)),(mp.mpf(0),mp.mpf(1)),costs)
        values=(mp.mpf(0),mp.mpf(0),mp.mpf(1),mp.mpf(1),budget)
        def val(d,x):
            M,h,c=d
            return c-2*sum(h[i]*x[i] for i in range(2))+sum(x[i]*M[i,j]*x[j] for i in range(2) for j in range(2))
        if number(kappa)==0:
            theta=(mp.mpf(0),mp.mpf(0))
            return {'holdings':theta,'risk':max(risk(m,theta) for m in laws),
                    'refinement':'zero budget fixes holdings exactly', 'kkt_residual':'0',
                    'node_risks':tuple(risk(m,theta) for m in laws),'success':True}
        floats=[(np.array(M.tolist(),dtype=float),np.array(h,dtype=float).reshape(2),float(c)) for M,h,c in data]
        def fval(d,x):
            M,h,c=d;return c-2*h@x[:2]+x[:2]@M@x[:2]
        def gjac(d,x):
            M,h,_=d;return np.r_[2*h-2*M@x[:2],1.]
        cons=[{'type':'ineq','fun':lambda x,d=d:x[2]-fval(d,x),
               'jac':lambda x,d=d:gjac(d,x)} for d in floats]
        cons.append({'type':'ineq','fun':lambda x:float(budget)-np.array(costs,dtype=float)@x[:2],
                     'jac':lambda x:np.r_[-np.array(costs,dtype=float),0.]})
        initial=np.array((0.,0.,max(fval(d,np.zeros(3)) for d in floats)))
        sol=minimize(lambda x:x[2],initial,jac=lambda x:np.array((0.,0.,1.)),
                     bounds=((0,1),(0,1),(None,None)),constraints=cons,method='SLSQP',
                     options={'ftol':1e-13,'maxiter':2000})
        candidate=sol.x
        qr=[fval(d,candidate) for d in floats]
        # Keep only materially active risks, then try all small active subsets
        # near the proposed optimum. A 2-D epigraph has three primal variables.
        active_nodes=[i for i,v in enumerate(qr) if max(qr)-v<=2e-7*max(1,abs(max(qr)))]
        active_faces=[i for i,(a,b) in enumerate(zip(normals,values))
                      if abs(float(b)-sum(float(a[j])*candidate[j] for j in range(2)))<=2e-7]
        refined=[]
        for qcount in range(1,min(3,len(active_nodes))+1):
            for selected in itertools.combinations(active_nodes,qcount):
                for fcount in range(0,min(3-qcount,len(active_faces))+1):
                    for faces in itertools.combinations(active_faces,fcount):
                        # Unknowns: theta[2], epigraph t, q risk weights and f
                        # financial multipliers. Equalities: gradient[2], q
                        # active risks, f active faces and sum weights=1.
                        def equations(*unknown):
                            x=unknown[:2];t=unknown[2];w=unknown[3:3+qcount];lam=unknown[3+qcount:]
                            result=[]
                            for j in range(2):
                                result.append(sum(w[n]*2*(sum(data[i][0][j,k]*x[k] for k in range(2))-data[i][1][j])
                                                  for n,i in enumerate(selected))
                                              +sum(lam[n]*normals[i][j] for n,i in enumerate(faces)))
                            result.extend(val(data[i],x)-t for i in selected)
                            result.extend(sum(normals[i][j]*x[j] for j in range(2))-values[i] for i in faces)
                            result.append(sum(w)-1)
                            return tuple(result)
                        start=tuple(mp.mpf(float(v)) for v in candidate)+tuple(mp.mpf(1)/qcount for _ in selected)+tuple(mp.mpf('.1') for _ in faces)
                        try:
                            point=tuple(mp.findroot(equations,start,tol=mp.mpf('1e-130'),maxsteps=80,solver='mdnewton'))
                        except (ValueError,ZeroDivisionError,TypeError):
                            continue
                        x=point[:2];t=point[2];weights=point[3:3+qcount];lambdas=point[3+qcount:]
                        residual=max(abs(v) for v in equations(*point))
                        allowance=mp.mpf('1e-90')
                        feasible=all(sum(a[j]*x[j] for j in range(2))-b<=allowance for a,b in zip(normals,values))
                        dominance=all(val(d,x)-t<=allowance for d in data)
                        nonnegative=all(w>=-allowance for w in weights+lambdas)
                        if residual<=allowance and feasible and dominance and nonnegative:
                            refined.append((max(risk(model,x) for model in laws),x,selected,faces,weights,lambdas,residual))
        if not refined:
            raise ArithmeticError(f'Independent 160-digit active-set refinement failed: {sol.message}; nodes {active_nodes}; faces {active_faces}.')
        objective,theta,selected,faces,weights,lambdas,residual=min(refined,key=lambda x:x[0])
        return {'holdings':theta,'risk':objective,'active_nodes':selected,'active_faces':faces,
                'risk_weights':weights,'financial_multipliers':lambdas,
                'kkt_residual':mp.nstr(residual,digits),'node_risks':tuple(risk(m,theta) for m in laws),
                'scipy_success':bool(sol.success),'scipy_message':str(sol.message),
                'refinement':'160-digit independently formed active-set KKT equations; 1e-90 diagnostic tolerance',
                'success':True}


def mixture_dual(contract,sigma,nodes,weights,kappa,multipliers):
    """Independent high-precision complete-square value of a fixed mixture."""
    laws=models(contract,sigma,nodes)
    matrix=mp.matrix(2);cross=mp.matrix(2,1);constant=mp.mpf(0)
    for w,model in zip(weights,laws):
        M,h,c=static.coefficients(model)
        matrix+=number(w)*M;cross+=number(w)*h;constant+=number(w)*c
    A=((-1,0),(0,-1),(1,0),(0,1),(laws[0]['s'],laws[0]['ck']))
    b=(0,0,1,1,number(kappa)*laws[0]['premium'])
    lam=tuple(number(x) for x in multipliers)
    z=mp.matrix([cross[j]-sum(lam[i]*A[i][j] for i in range(5))/2 for j in range(2)])
    return constant-sum(lam[i]*b[i] for i in range(5))-(z.T*(matrix**-1)*z)[0]
