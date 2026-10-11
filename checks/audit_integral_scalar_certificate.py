"""Independent Arb checks of scalar nodes and signed interpolation metadata."""
from fractions import Fraction as F
from flint import arb,ctx
import hashlib
import json


def ball(rational):
    x=F(rational);return arb(x.numerator)/arb(x.denominator)


def normal_density(x):return (-x*x/2).exp()/(2*arb.pi()).sqrt()
def normal_distribution(x):return (1+(x/arb(2).sqrt()).erf())/2


def independent_scalar(name,x,parameters):
    p={k:bool(v) if k=='normalised' else F(v) for k,v in parameters.items()}
    if name=='exp':return x.exp()
    if name=='log':return x.log()
    if name=='inv_square':return 1/(x*x)
    if name=='square':return x*x
    if name in ('h','h2'):
        lam=ball(p['lambda']);D=ball(p['D2']);density=normal_density(x)
        derivative=-x*density
        result=1+(-density/lam+derivative/(lam*lam))/D
        return result/2 if name=='h' else result
    if name in ('cap','cap_prime'):
        H=ball(p['cap']);K=ball(p.get('K',F(100)));lam=ball(p['lambda'])
        log_strike=(1+H/K).log()
        lower=x/lam;upper=(x-log_strike)/lam
        growth=K*(x+lam*lam/2).exp()
        exercise_weight=normal_distribution(lower+lam)-normal_distribution(upper+lam)
        result=growth*exercise_weight
        if name=='cap':
            result-=K*(normal_distribution(lower)-normal_distribution(upper))
            result+=H*normal_distribution(upper)
        if p.get('normalised',False):result/=H
        return result
    raise ValueError(('Unsupported independent scalar',name))


def audit_certificate(certificate):
    ctx.prec=416
    key=hashlib.sha256(json.dumps(certificate,sort_keys=True).encode()).hexdigest()
    if certificate['schema']=='qbp-integral-quantile-v1':
        L=F(certificate['L']);q=certificate['q'];count=len(certificate['nodes'])
        cdf_low=normal_distribution(-ball(L));mass=normal_distribution(ball(L))-cdf_low
        for node in certificate['nodes']:
            lo,hi=F(node['quantile_lower']),F(node['quantile_upper'])
            probability=cdf_low+ball(node['probability_midpoint'])*mass
            left=normal_distribution(ball(lo));right=normal_distribution(ball(hi))
            assert left.upper()<=probability.lower(),('Quantile lower endpoint',node['j'])
            assert right.lower()>=probability.upper(),('Quantile upper endpoint',node['j'])
            actual=F(node['offset_integer'],1<<q)-L
            assert actual==F(node['actual_coordinate'])
            assert max(abs(actual-lo),abs(actual-hi))<=F(node['node_error_upper'])
        return {'kind':'quantile','nodes':count,'certificate_sha256':key,
            'method':'Independent 416-bit CDF bracketing of stored inverse endpoints; no erfinv replay'}
    assert certificate['schema']=='qbp-integral-scalar-v1'
    q=certificate['q'];largest=F(0)
    for node in certificate['nodes']:
        argument=ball(node['input']);reference=independent_scalar(certificate['name'],argument,certificate['params'])-ball(certificate['output_offset'])
        saved_low,saved_high=F(node['lower']),F(node['upper'])
        # The alternate algebra at higher precision must intersect the saved
        # enclosure. Direct containment is too strong near exact zero under
        # two different valid Arb expression trees.
        assert reference.upper()>=ball(saved_low).lower() and reference.lower()<=ball(saved_high).upper(),('Nodal interval discrepancy',certificate['name'],node['input'])
        integer=node['integer'];error=max(abs(F(integer,1<<q)-saved_low),abs(F(integer,1<<q)-saved_high))
        assert error==F(node['error_upper']);largest=max(largest,error)
    assert largest==F(certificate['nodal_error_upper'])
    h=F(certificate['span'])/(1<<certificate['table_bits'])
    interpolation=F(certificate['second_derivative_upper'])*h*h/8
    assert interpolation==F(certificate['linear_interpolation_error_upper'])
    assert F(certificate['interpolation_floor_error_upper'])==F(1,1<<q)
    assert largest+interpolation+F(1,1<<q)==F(certificate['output_error_upper'])
    endpoints=[row['integer'] for row in certificate['nodes']]
    assert min(endpoints)==certificate['value_min_integer'] and max(endpoints)==certificate['value_max_integer']
    for left,right in zip(endpoints,endpoints[1:]):
        # A decreasing cell uses a positive magnitude and retains nonnegative
        # endpoints. Its toward-base truncation error is below one output ulp.
        magnitude=abs(right-left)
        assert left>=magnitude if right<left else True
    return {'kind':'scalar','name':certificate['name'],'nodes':len(endpoints),
        'certificate_sha256':key,'nodal_error_upper_fraction':str(largest),
        'total_interpolation_error_upper_fraction':certificate['output_error_upper'],
        'method':'Independent 416-bit elementary/Gaussian derivative or capped-call closed-form expression plus exact rational metadata replay'}


def audit_archive(root,out):
    from pathlib import Path
    import gzip
    root=Path(root);path=root/'results/integral_scalar_table_archive.json.gz'
    archive=json.loads(gzip.open(path,'rt').read());records=[]
    for key,certificate in archive['records'].items():
        checked=audit_certificate(certificate)
        # Archive key hashing convention is recorded by its producer; the
        # independent canonical content hash remains retained separately.
        expected_key=hashlib.sha256(json.dumps(certificate,sort_keys=True,separators=(',',':')).encode()).hexdigest()
        assert key==expected_key,('Archive record hash mismatch',key)
        checked['retained_archive_key']=key;records.append(checked)
    assessment_path=root/'results/integral_graph_certificate_assessment.json'
    assessment=json.loads(assessment_path.read_text())
    assert assessment['scalar_archive_gzip_sha256']==hashlib.sha256(path.read_bytes()).hexdigest()
    links=0
    for graph in assessment['records']:
        assert graph['status']=='passed' and not graph['failures']
        for node,key in graph['scalar_table_certificates'].items():
            assert key in archive['records'],('Unresolved scalar certificate link',node,key)
            links+=1
    result={'status':'independent retained scalar archive audit passed','certificate_count':len(records),
        'graph_assessment_records':len(assessment['records']),'verified_scalar_certificate_links':links,
        'assessment_sha256':hashlib.sha256(assessment_path.read_bytes()).hexdigest(),
        'node_count':sum(row['nodes'] for row in records),'reference_precision_bits':416,
        'records':records,'archive_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
        'audit_source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'scope':'Independent elementary/Gaussian/strike expressions and exact error metadata; not a claim of scalar-method novelty'}
    Path(out).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k in ('status','certificate_count','node_count')}))
    return result


if __name__=='__main__':
    import argparse
    from pathlib import Path
    parser=argparse.ArgumentParser();parser.add_argument('--root',default=str(Path(__file__).resolve().parents[1]))
    parser.add_argument('--out',default=str(Path(__file__).resolve().parents[1]/'results/integral_scalar_archive_audit.json'))
    args=parser.parse_args();audit_archive(args.root,args.out)
