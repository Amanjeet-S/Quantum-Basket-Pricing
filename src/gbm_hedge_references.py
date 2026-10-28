"""Independent high-precision static-hedge diagnostics, not interval proofs.

The reference constructs piecewise-affine claim/gain products directly rather
than the engine's three-gain Gram representation. Inputs are exact rationals.
"""
from fractions import Fraction as F
from itertools import combinations
import mpmath as mp


def number(value):
    if isinstance(value,mp.mpf):return value
    value=F(value)
    return mp.mpf(value.numerator)/value.denominator


def parameters(contract,sigma,mu):
    s,k,h,t,r,q=map(number,(contract.spot,contract.strike,contract.cap,
                           contract.maturity,contract.rate,contract.dividend))
    sig,drift=number(sigma),number(mu)
    j=k+h;sd=sig*mp.sqrt(t);mean=mp.log(s)+(drift-q-sig*sig/2)*t
    discount,stock=mp.exp(-r*t),mp.exp((q-r)*t)
    phi=lambda z:mp.erfc(-z/mp.sqrt(2))/2
    def call(strike):
        d=(mp.log(s/strike)+(r-q+sig*sig/2)*t)/sd
        return s*mp.exp(-q*t)*phi(d)-strike*discount*phi(d-sd)
    ck,cj=call(k),call(j)
    return {'s':s,'k':k,'j':j,'h':h,'t':t,'r':r,'q':q,'sigma':sig,'mu':drift,
            'sd':sd,'mean':mean,'discount':discount,'stock_factor':stock,
            'ck':ck,'cj':cj,'premium':ck-cj,'Phi':phi}


def pieces(model):
    s,k,j,d,a,ck,v=(model[x] for x in ('s','k','j','discount','stock_factor','ck','premium'))
    h=model['h'];mean,sd=model['mean'],model['sd']
    def above(power,boundary=None):
        value=mp.exp(power*mean+mp.mpf(power*power)*sd*sd/2)
        return value if boundary is None else value*model['Phi']((mean+power*sd*sd-mp.log(boundary))/sd)
    total=[above(p) for p in range(3)];lo=[above(p,k) for p in range(3)];hi=[above(p,j) for p in range(3)]
    moments=[tuple(total[i]-lo[i] for i in range(3)),
             tuple(lo[i]-hi[i] for i in range(3)),tuple(hi)]
    gs=(a,-s)
    gains=[(gs,(mp.mpf(0),-ck)),(gs,(d,-d*k-ck)),(gs,(d,-d*k-ck))]
    liability=[(mp.mpf(0),-v),(d,-d*k-v),(mp.mpf(0),d*h-v)]
    return moments,gains,liability


def inner(first,second,moments):
    a,b=first;c,d=second;m0,m1,m2=moments
    return a*c*m2+(a*d+b*c)*m1+b*d*m0


def coefficients(model):
    moments,gains,liability=pieces(model)
    matrix=mp.matrix(2);cross=mp.matrix(2,1);constant=mp.mpf(0)
    for m,g,y in zip(moments,gains,liability):
        for i in range(2):
            cross[i]+=inner(g[i],y,m)
            for j in range(2):matrix[i,j]+=inner(g[i],g[j],m)
        constant+=inner(y,y,m)
    return matrix,cross,constant


def fixed_risk(model,holdings):
    theta=tuple(number(x) for x in holdings)
    moments,gains,liability=pieces(model)
    result=mp.mpf(0)
    for m,g,y in zip(moments,gains,liability):
        residual=tuple(y[i]-theta[0]*g[0][i]-theta[1]*g[1][i] for i in range(2))
        result+=inner(residual,residual,m)
    return result


def solve_two(matrix,rhs):
    determinant=matrix[0][0]*matrix[1][1]-matrix[0][1]*matrix[1][0]
    if determinant==0:raise ZeroDivisionError('Exactly singular numerical reference face.')
    return ((matrix[1][1]*rhs[0]-matrix[0][1]*rhs[1])/determinant,
            (matrix[0][0]*rhs[1]-matrix[1][0]*rhs[0])/determinant)


def optimum(model,kappa):
    """160-digit independently written finite-face QP diagnostic.

    A fixed 1e-130 feasibility allowance resolves arithmetic in equality
    candidates; this is a numerical reference, never a feasibility certificate.
    """
    M,h,c=coefficients(model)
    mat=tuple(tuple(M[i,j] for j in range(2)) for i in range(2))
    dot=lambda a,b:sum(a[i]*b[i] for i in range(2))
    mv=lambda x:tuple(dot(row,x) for row in mat)
    costs=(model['s'],model['ck']);B=number(kappa)*model['premium']
    A=((-1,0),(0,-1),(1,0),(0,1),costs);b=(0,0,1,1,B)
    candidates=[solve_two(mat,tuple(h))]
    for normal,value in zip(A,b):
        centre=tuple(mp.mpf(normal[i])*value/dot(normal,normal) for i in range(2))
        tangent=(-mp.mpf(normal[1]),mp.mpf(normal[0]))
        factor=(dot(tangent,h)-dot(tangent,mv(centre)))/dot(tangent,mv(tangent))
        candidates.append(tuple(centre[i]+factor*tangent[i] for i in range(2)))
    for i,j in combinations(range(5),2):
        try:candidates.append(solve_two((A[i],A[j]),(b[i],b[j])))
        except ZeroDivisionError:pass
    epsilon=mp.mpf('1e-130')
    feasible=[x for x in candidates if all(dot(row,x)-value<=epsilon*max(1,abs(value)) for row,value in zip(A,b))]
    if B==0:feasible=[(mp.mpf(0),mp.mpf(0))]
    if not feasible:raise ArithmeticError('No numerical reference feasible point.')
    # The expanded moment quadratic is safe at the declared high precision;
    # returned risk is recomputed as piecewise squared residuals.
    theta=min(feasible,key=lambda x:c-2*dot(h,x)+dot(x,mv(x)))
    return {'holdings':theta,'risk':fixed_risk(model,theta),'candidates':len(candidates),
            'feasibility_allowance':'1e-130 scaled; reference only'}


def positive_risk(model,holdings):
    """Independent positive Gaussian integration with density/size scaling.

    Normal density is factored at its largest point in each payoff region.
    Each affine residual is also rescaled before integration, preventing tiny
    tail values from being accepted solely on an absolute stopping threshold.
    """
    theta=tuple(number(x) for x in holdings)
    _,gains,liability=pieces(model)
    mean,sd,k,j=model['mean'],model['sd'],model['k'],model['j']
    zk,zj=(mp.log(k)-mean)/sd,(mp.log(j)-mean)/sd
    limits=((-mp.inf,zk),(zk,zj),(zj,mp.inf))
    total=mp.mpf(0)
    for (low,high),g,y in zip(limits,gains,liability):
        slope=y[0]-theta[0]*g[0][0]-theta[1]*g[1][0]
        intercept=y[1]-theta[0]*g[0][1]-theta[1]*g[1][1]
        origin=low if low>0 else high if high<0 else mp.mpf(0)
        scale=abs(slope)*mp.exp(mean+sd*origin)+abs(intercept)
        if scale==0:continue
        def integrand(z):
            exponent=-(z*z-origin*origin)/4
            value=(slope*mp.exp(mean+sd*z+exponent)+intercept*mp.exp(exponent))/scale
            return value*value
        cuts=[low,high]
        for x in (-16,-12,-8,-4,-2,-1,0,1,2,4,8,12,16):
            if low<x<high:cuts.append(mp.mpf(x))
        if low>0:
            for distance in ('0.01','0.05','0.1','0.25','0.5','1','2','4'):
                value=low+mp.mpf(distance)
                if value<high:cuts.append(value)
        elif high<0:
            for distance in ('0.01','0.05','0.1','0.25','0.5','1','2','4'):
                value=high-mp.mpf(distance)
                if value>low:cuts.append(value)
        cuts=sorted(set(cuts))
        normalised=mp.quad(integrand,cuts)
        total+=scale*scale*mp.exp(-origin*origin/2)/mp.sqrt(2*mp.pi)*normalised
    return total
