"""Own-formula auxiliary quadrature and decoder checks for four cases."""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json
from flint import arb,ctx
ROOT=Path(__file__).resolve().parents[1];RESULTS=ROOT/'results'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def A(x):x=F(x);return arb(x.numerator)/x.denominator
def endpoint(x,up):
    m,e=map(int,(x.upper() if up else x.lower()).man_exp())
    return F(m*(1<<e)) if e>=0 else F(m,1<<(-e))
def bounds(x):return endpoint(x,False),endpoint(x,True)
def hi(v):return F(v['upper_exact_fraction'])

def main():
    path=RESULTS/'four_survivor_quantum_resources.json';data=json.loads(path.read_text());prep={p['model']['model_id']:p for p in data['preparations']};ctx.prec=1024;root=arb(6).sqrt();ideals={'z_minus':(4-root)/10,'z_plus':(4+root)/10,'w_minus':arb(1)/2+1/(3*root),'w_plus':arb(1)/2-1/(3*root)}
    constant_checks=0;budget_checks=0;decoder_checks=0
    # The monic degree-two polynomial is z^2-(4/5)z+1/10.
    p=[F(1,10),F(-4,5),F(1)];squared=[F(0)]*5
    for i,a in enumerate(p):
        for j,b in enumerate(p):squared[i+j]+=a*b
    norm=sum(a*F(2,(j+1)*(j+2)) for j,a in enumerate(squared))
    if norm!=F(1,300):raise AssertionError('Beta Gauss monic squared norm differs')
    for row in data['records']:
        preparation=prep[row['model_id']];cfg=row['base_graph_certificate']['configuration'];norms=preparation['payoff_norms'];U=F(cfg['U_star']);M1=F(norms['M1_upper']);M2=F(norms['M2_upper']);M3=F(norms['C_r_upper']['3']);M4=F(norms['C_r_upper']['4']);M6=F(norms['C_r_upper']['6']);mom={k:hi(v) for k,v in preparation['finite_moments']['centred_U_moments'].items()};B=F(cfg['monetary_signal_normaliser'])
        if F(cfg['M1'])!=M1 or F(cfg['M2'])!=M2:raise AssertionError('Decoder norms differ')
        expected_B=min(F(1),M1*U) if row['kind']=='positive_r1' else M2*U*U/2
        if B!=expected_B:raise AssertionError('Monetary signal normaliser differs')
        expected_decode=1 if row['kind']=='positive_r1' else 2
        if F(cfg['monetary_probability_decode_factor'])!=expected_decode:raise AssertionError('Native versus signed decoder differs')
        decoder_checks+=1
        if row['complete_deterministic_budget'] is None:continue
        aux=row['complete_deterministic_budget']['auxiliary_components'];kind=row['kind']
        if kind in ('positive_r1','raw_direct_r2'):
            if F(aux['bias_upper'])!=0 or not aux['auxiliary_free']:raise AssertionError('Auxiliary-free allowance differs')
        elif kind=='boundary_r2':
            if aux['beta_bits']!=8:raise AssertionError('Beta precision changed')
            h=F(1,256);expected=h*h*(M3*mom['3']/12+M4*mom['4']/48)
            if F(aux['bias_upper'])!=expected or aux['Beta_cell_first_moment']!='-h^3/6' or aux['Beta_cell_second_moment_sum']!='h^2/12':raise AssertionError('Beta cell quadrature allowance differs')
        else:
            constants=aux['constants'];node=weight=F(0)
            for name,ideal in ideals.items():
                r=constants[name];l,u=bounds(ideal)
                if not F(r['ideal_lower'])<=l<=u<=F(r['ideal_upper']):raise AssertionError('Gauss constant interval differs')
                x=F(r['actual']);maximum=max(abs(x-l),abs(x-u))
                if maximum>F(r['error_upper']):raise AssertionError('Gauss constant error too small')
                if name!='w_plus':
                    a=l*(1<<28)+F(1,2);b=u*(1<<28)+F(1,2)
                    ix=a.numerator//a.denominator;iy=b.numerator//b.denominator
                    if ix!=iy or F(ix,1<<28)!=x:raise AssertionError('Gauss constant dyadic encoding differs')
                constant_checks+=1
            if F(constants['w_minus']['actual'])+F(constants['w_plus']['actual'])!=1:raise AssertionError('Actual quadrature weights not normalised')
            for suffix in ('minus','plus'):
                node+=F(constants['w_'+suffix]['ideal_upper'])*F(constants['z_'+suffix]['error_upper']);weight+=F(constants['w_'+suffix]['error_upper'])
            trunc=M6*mom['6']/14400;fixed=M3*mom['3']*node/2+M2*mom['2']*weight/2
            if F(aux['Gauss_remainder_upper'])!=trunc or F(aux['fixed_parameter_bias_upper'])!=fixed or F(aux['bias_upper'])!=trunc+fixed:raise AssertionError('Gauss quadrature allowance differs')
        budget_checks+=1
    out=Path(__file__).resolve().parent/'four_quantum_auxiliary_audit.json'
    if out.exists():raise FileExistsError(out)
    record={'schema':'qbp-four-quantum-own-formula-auxiliary-audit-v1','status':'passed','auditor_source_sha256':sha(Path(__file__)),'producer_result_sha256':sha(path),'ideal_constant_precision_bits':1024,'producer_resource_and_certificate_APIs_imported':False,'Beta_Gauss_monic_squared_norm':'1/300','auxiliary_allowance_formula_checks':budget_checks,'normaliser_and_decode_factor_checks':decoder_checks,'independent_Gauss_constant_interval_checks':constant_checks,'scope':'Own-formula Beta midpoint and two-node Gauss allowances, exact signal normalisers, fixed parameter error and signed versus native positive decoding. Finite U moments remain retained producer evidence subject to the separate moment certificates.'}
    out.write_text(json.dumps(record,indent=2)+'\n');print('bounded auxiliary audit passed',flush=True)

if __name__=='__main__':main()
