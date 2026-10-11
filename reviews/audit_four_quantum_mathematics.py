"""Bounded own-formula review of four quantum financial envelopes and priors.

No producer, financial-prior or resource API is imported. This complements
the full literal-resource audit; it does not repeat production compilation.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json,time
from flint import arb,ctx
ROOT=Path(__file__).resolve().parents[1];RESULTS=ROOT/'results'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def digest(v):return hashlib.sha256(json.dumps(v,sort_keys=True,separators=(',',':')).encode()).hexdigest()
def A(x):x=F(x);return arb(x.numerator)/x.denominator
def endpoint(x,up):
    m,e=map(int,(x.upper() if up else x.lower()).man_exp())
    return F(m*(1<<e)) if e>=0 else F(m,1<<(-e))
def bound(x):return endpoint(x,False),endpoint(x,True)
def Phi(x):return (1+(x/arb(2).sqrt()).erf())/2
def phi(x):return (-x*x/2).exp()/(2*arb.pi()).sqrt()
def lo(v):return F(v['lower_exact_fraction'])
def hi(v):return F(v['upper_exact_fraction'])
def enclosed(record,value):
    l,h=bound(value)
    if not lo(record)<=l<=h<=hi(record):raise AssertionError('Fresh 1024-bit value not enclosed')

def cells(data):
    ctx.prec=1024;L=A(6);lower=Phi(-L);mass=Phi(L)-lower;cellmass=mass/1024
    bounds=[-L]+[arb(2).sqrt()*(2*(lower+mass*j/1024)-1).erfinv() for j in range(1,1024)]+[L]
    b1=F(0);b2=arb(0)
    for j in range(1024):
        m=arb(2).sqrt()*(2*(lower+mass*A(F(2*j+1,2048)))-1).erfinv()
        l,h=bound(m);l=l*(1<<28)+F(1,2);h=h*(1<<28)+F(1,2)
        c0=l.numerator//l.denominator;c1=h.numerator//h.denominator
        if c0!=c1 or c0!=data['actual_financial_node_integers'][j]:raise AssertionError('Quantum actual-node mismatch')
        c=F(c0,1<<28);a,b=bounds[j:j+2]
        if not bound(a)[1]<=c<=bound(b)[0]:raise AssertionError('Actual node is outside cell')
        first=phi(a)-phi(b)-A(c)*cellmass
        second=(1+A(c*c))*cellmass+(a-2*A(c))*phi(a)-(b-2*A(c))*phi(b)
        b1+=max(map(abs,bound(first)));b2+=second/2
    enclosed(data['transport_profile']['mass_interval'],mass)
    enclosed(data['transport_profile']['rho2_mass_interval'],b2)
    if b1>hi(data['transport_profile']['rho1_mass_interval']):raise AssertionError('rho1 domination fails')
    return mass,b1,b2

def geometry(model):
    w=list(map(F,model['weights']));r=[list(map(F,row)) for row in model['loadings']];mu=list(map(F,model['log_drifts']));lam=F(model['lambda'])
    if min(w)<=0 or sum(w)!=1:raise AssertionError('Basket weight premise fails')
    p=model['outer_dimension'];mg=sum(x*y for x,y in zip(w,mu));rg=[sum(w[i]*r[i][j] for i in range(len(w))) for j in range(p)]
    d=[m-mg for m in mu];s=[[a-b for a,b in zip(row,rg)] for row in r]
    for m,row in zip(mu,r):
        if m+(lam*lam+sum(v*v for v in row))/2!=F(model['risk_free_rate'])-F(model['dividend_yield']):raise AssertionError('Risk-neutral drift premise fails')
    return w,r,mu,lam,p,mg,rg,d,s

def independent_coefficients_and_moment(preparation):
    ctx.prec=1024;w,r,mu,lam,p,mg,rg,d,s=geometry(preparation['model']);variance=lam*lam+sum(v*v for v in rg);sigma=A(variance).sqrt();ell=A(F(101,100)).log()
    def price(m):
        m=A(m);x=m/sigma;xh=(m-ell)/sigma;growth=100*(m+A(variance/2)).exp()
        return (growth*Phi(x+sigma)-100*Phi(x))-(growth*Phi(xh+sigma)-101*Phi(xh))
    def delta(m):
        m=A(m);growth=100*(m+A(variance/2)).exp()
        return growth*(Phi(m/sigma+sigma)-Phi((m-ell)/sigma+sigma))
    c0=price(mg);c1=-delta(mg)
    for wi,di,si in zip(w,d,s):c1+=A(wi)*A(di+sum(v*v for v in si)/2).exp()*delta(mg+sum(a*b for a,b in zip(rg,si)))
    coef=preparation['coefficient_certificate']
    enclosed(coef['c0'],c0);enclosed(coef['c1'],c1);enclosed(coef['coefficient'],c0+c1)
    raw1=sum((A(wi)*A(di+sum(v*v for v in si)/2).exp() for wi,di,si in zip(w,d,s)),arb(0));raw2=arb(0)
    for i in range(len(w)):
        for k in range(len(w)):
            v=[a+b for a,b in zip(s[i],s[k])]
            raw2+=A(w[i]*w[k])*A(d[i]+d[k]+sum(a*a for a in v)/2).exp()
    U2=raw2-2*raw1+1;enclosed(preparation['continuous_second_moment']['continuous_U2_interval'],U2)
    return 4

def check_direct(preparation,mass,b1,b2):
    ctx.prec=1024;w,r,mu,lam,p,mg,rg,d,s=geometry(preparation['model']);norm=preparation['payoff_norms'];M1=F(norm['M1_upper']);M2=F(norm['M2_upper']);stored=preparation['R1_direct_global_financial_transfer'];value=arb(0)
    for j,row in enumerate(stored['coordinates']):
        aj=max(abs(a[j]) for a in r);D1=M1*(aj+abs(rg[j]));D2=(M1+M2)*(aj*aj+rg[j]*rg[j])
        if F(row['D1_upper'])!=D1 or F(row['D2_upper'])!=D2:raise AssertionError('Direct R1 full-loading derivative mismatch')
        value+=A(D1*b1)+A(D2)*b2
    value=mass**(p-1)*value+(1-mass**p)
    if bound(value)[1]>F(stored['financial_transfer_upper']):raise AssertionError('Direct R1 allowance underestimates fresh 1024-bit bound')
    if F(stored['other_coordinate_mass_upper'])<bound(mass**(p-1))[1]:raise AssertionError('Hybrid mass factor too small')
    if F(stored['tail_upper'])<bound(1-mass**p)[1]:raise AssertionError('Full R1 bounded tail too small')
    enclosed(stored['mass_interval'],mass**p)
    chosen=min(F(stored['financial_transfer_upper']),F(preparation['R1_older_financial_transfer']['financial_transfer_upper']))
    if F(preparation['R1_selected_financial_transfer']['financial_transfer_upper'])!=chosen:raise AssertionError('Deterministic transfer minimum differs')
    return p

def check_previous_prior(row,preparation):
    cert=row['base_graph_certificate'];cfg=cert['configuration'];previous=row['correction_prior'];B=F(cfg['monetary_signal_normaliser']);eta=F(cert['probability_error_upper']);mom=preparation['finite_moments']['centred_U_moments']
    if row['kind']=='positive_r1':
        mean=min(F(1),hi(preparation['finite_moments']['correlated_underlying_gap_interval']),F(cfg['M1'])*hi(mom['1']))
        expected=min(F(1),mean/B+eta)
        if F(previous['ideal_mean_upper'])!=mean or F(previous['probability_upper'])!=expected or not previous['native_positive']:raise AssertionError('Previous positive prior differs')
    else:
        native=cfg['output_encoding']=='positive';expected=min(F(1),hi(mom['2'])/F(cfg['U_star'])**2+(1 if native else 2)*eta)
        if F(previous['probability_upper'])!=expected or previous['native_positive']!=native:raise AssertionError('Previous curvature prior differs')
    return 1

def check_informed_prior(row,preparation,key,financial):
    pr=row[key];previous=row['correction_prior'];cert=row['base_graph_certificate'];B=F(cert['configuration']['monetary_signal_normaliser']);eta=F(cert['probability_error_upper']);M2=F(preparation['payoff_norms']['M2_upper']);c1=preparation['coefficient_certificate']['c1'];moment=preparation['continuous_second_moment'];remainder=M2*hi(moment['continuous_U2_interval'])/2;beta=F(financial['financial_transfer_upper']);mass=financial['mass_interval']
    cl=lo(c1)-remainder;cu=hi(c1)+remainder
    fl=max(F(0),cl-beta)/hi(mass);fu=min(B,(cu+beta)/lo(mass),F(previous['ideal_mean_upper']))
    pl=max(F(0),fl/B-eta);pu=min(F(1),F(previous['probability_upper']),fu/B+eta)
    for name,value in {'continuous_R2_mean_absolute_upper':remainder,'continuous_R1_mean_lower':cl,'continuous_R1_mean_upper':cu,'finite_ideal_mean_lower':fl,'finite_ideal_mean_upper':fu,'probability_lower':pl,'probability_upper':pu,'financial_transfer_upper':beta}.items():
        if F(pr[name])!=value:raise AssertionError('Coefficient-informed '+name+' differs')
    if not 0<=pl<=pu<=1 or pr['continuous_moment_record_sha256']!=digest(moment):raise AssertionError('Informed prior intersection or provenance differs')
    return 1

def check_budget(row,preparation,key,financial):
    bud=row[key];cert=row['base_graph_certificate'];cfg=cert['configuration'];B=F(cfg['monetary_signal_normaliser']);decode=F(cfg['monetary_probability_decode_factor']);base=preparation['coefficient_certificate']['c0' if row['kind']=='positive_r1' else 'coefficient'];mass=financial['mass_interval'];cr=(hi(base)-lo(base))/2;mr=(hi(mass)-lo(mass))/2;arithmetic=decode*B*F(cert['probability_error_upper']);aux=F(bud['auxiliary_components']['bias_upper'])
    if row['kind'] in ('positive_r1','raw_direct_r2') and aux!=0:raise AssertionError('Auxiliary-free method charged inconsistently')
    expected=cr+F(financial['financial_transfer_upper'])+hi(mass)*(aux+arithmetic)+mr*B
    if bud['coefficient_interval']!=base or bud['mass_interval']!=mass or F(bud['total_deterministic_undiscounted_price_error_upper'])!=expected:raise AssertionError('Complete deterministic budget differs')
    return 1

def main():
    started=time.perf_counter();path=RESULTS/'four_survivor_quantum_resources.json';data=json.loads(path.read_text());protocol_path=RESULTS/data['protocol_file'];protocol=json.loads(protocol_path.read_text())
    if sha(protocol_path)!=data['protocol_sha256']:raise AssertionError('Protocol provenance differs')
    for name,value in data['dependency_source_sha256'].items():
        if sha(ROOT/name)!=value:raise AssertionError('Source changed '+name)
    for name,value in data['input_sha256'].items():
        if sha(RESULTS/name)!=value:raise AssertionError('Input changed '+name)
    mass,b1,b2=cells(data);preparations={v['model']['model_id']:v for v in data['preparations']}
    derivative_rows=coefficient_checks=prior_checks=budget_checks=0
    for prep in preparations.values():
        derivative_rows+=check_direct(prep,mass,b1,b2);coefficient_checks+=independent_coefficients_and_moment(prep)
    for row in data['records']:
        prep=preparations[row['model_id']];financial=prep['R1_selected_financial_transfer' if row['kind']=='positive_r1' else 'R2_financial_transfer']
        if row['complete_deterministic_budget'] is not None:budget_checks+=check_budget(row,prep,'complete_deterministic_budget',financial)
        if 'correction_prior' in row:prior_checks+=check_previous_prior(row,prep)
        if row['kind']=='positive_r1':
            prior_checks+=check_informed_prior(row,prep,'coefficient_informed_prior',financial)
            prior_checks+=check_informed_prior(row,prep,'older_coefficient_informed_prior',prep['R1_older_financial_transfer'])
            budget_checks+=check_budget(row,prep,'older_financial_deterministic_budget',prep['R1_older_financial_transfer'])
    wins=[]
    for selected in data['selected_prices']:
        methods={k:v['selected']['result']['T_count'] for k,v in selected['methods'].items()}
        if selected['quantum_resource_winner']!=min(methods,key=methods.get) or selected['quantum_resource_winner']!='positive_r1':raise AssertionError('Quantum ranking differs')
        wins.append({'model_id':selected['model_id'],'logical_T_counts':methods,'winner':'positive_r1'})
    record={'schema':'qbp-four-quantum-bounded-independent-mathematical-review-v1','status':'passed','reviewer_source_sha256':sha(Path(__file__)),'producer_result_sha256':sha(path),'protocol_sha256':sha(protocol_path),'own_formula_precision_bits':1024,'producer_and_resource_APIs_imported':False,'actual_node_checks':1024,'full_loading_derivative_rows':derivative_rows,'fresh_coefficient_and_continuous_moment_interval_checks':coefficient_checks,'prior_formula_checks':prior_checks,'complete_deterministic_budget_checks':budget_checks,'selected_method_checks':wins,'scope':'All full-loading R1 derivatives, raw mass factors, bounded full tails, continuous coefficients and U2 moments, prior intersections and deterministic budgets checked. Literal primitive emission, scalar-table packing coverage and QFT/confidence counts belong to the separate full resource audit; no hardware execution, classical lower bound or theoretical novelty established.','review_seconds_diagnostic':time.perf_counter()-started}
    out=Path(__file__).resolve().parent/'four_quantum_mathematics_audit.json'
    if out.exists():raise FileExistsError(out)
    out.write_text(json.dumps(record,indent=2)+'\n');print('bounded quantum mathematical audit passed',flush=True)

if __name__=='__main__':main()
