"""Independent exact-statistics and 300-digit confidence endpoint replay.

No project confidence helper is imported. Finite Bernoulli probabilities are
summed as integer numerators over their exact common denominator. The finite
examples support the implementation audit, rather than replace the proof.
"""
from decimal import Decimal, localcontext
from fractions import Fraction as F
import argparse
import hashlib
import json
import math
from pathlib import Path


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def dec(x):
    x=F(x);return Decimal(x.numerator)/Decimal(x.denominator)

def reference(rho,variance,mean,statistics,av,am):
    rho,variance,mean,av,am=map(F,(rho,variance,mean,av,am))
    R=int(statistics['replicates']);m=R//2;D=int(statistics['batch_denominator'])
    observation=F(int(statistics['integer_paired_difference_square_sum']),2*m*D*D)
    sample=F(int(statistics['integer_batch_sum']),R*D)
    with localcontext() as ctx:
        ctx.prec=300
        lv=dec(1/av).ln();lm=dec(2/am).ln();M=dec(rho*rho/2);t=lv/m
        pair=dec(observation)+M*t+(M*M*t*t+2*M*t*dec(observation)).sqrt()
        U=min(dec(variance),pair);b=dec(rho)*lm/(3*R)
        width=b+(b*b+2*U*lm/R).sqrt()
        return {'log_variance':lv,'log_mean':lm,'pair':pair,'variance':U,'width':width,
                'lower':max(Decimal(0),dec(sample)-width),
                'upper':min(dec(mean),dec(sample)+width)}


def exact_bernoulli_coverage(log_variance_upper=None,log_mean_upper=None):
    """Enumerate all (00, mixed, 11) pair types, aggregating orientations.

    The generic deterministic priors are rho=1, Var<=1/4, E B<=1.
    The interval implementation is replayed independently at 300 digits,
    with the actual retained outward logarithms and separately implemented
    160-bit rational square roots.
    All probabilities use p=a/32 and denominator 32**R.
    """
    def up_root(x):
        x=F(x);scaled=x*(1<<320)
        n=(scaled.numerator+scaled.denominator-1)//scaled.denominator
        k=math.isqrt(n)
        if k*k<n:k+=1
        return F(k,1<<160)
    av=am=F(1,2200)
    with localcontext() as ctx:
        ctx.prec=300
        # The retained actual log values have already passed independent
        # 300-digit outward checks. A standalone call uses its own enclosures.
        logv=F(log_variance_upper) if log_variance_upper is not None else F(Decimal(2200).ln()+Decimal('1e-120'))
        logm=F(log_mean_upper) if log_mean_upper is not None else F(Decimal(4400).ln()+Decimal('1e-120'))
        assert dec(logv)>Decimal(2200).ln()
        assert dec(logm)>Decimal(4400).ln()
    cases=[];type_count=0;probability_types=0;nonzero_mean_failure=0;nonzero_pair_failure=0
    for R in (64,128):
        m=R//2;data=[]
        for mixed in range(m+1):
            z=F(mixed,2*m);t=logv/m
            upper=z+t/2+up_root(t*t/4+t*z)
            U=min(F(1,4),upper);b=logm/(3*R)
            width=b+up_root(b*b+2*U*logm/R)
            for both_one in range(m-mixed+1):
                both_zero=m-mixed-both_one
                h=mixed+2*both_one
                coefficient=math.comb(m,mixed)*math.comb(m-mixed,both_one)*(1<<mixed)
                data.append((h,upper,width,coefficient))
        type_count+=len(data)
        for a in range(33):
            denominator=32**R;p=F(a,32);true_variance=p*(1-p)
            powers_a=[a**h for h in range(R+1)]
            powers_b=[(32-a)**h for h in range(R+1)]
            total=bad_pair=bad_mean=bad_joint=0
            for h,upper,width,coefficient in data:
                weight=coefficient*powers_a[h]*powers_b[R-h]
                total+=weight
                pair_failure=true_variance>upper
                interval_failure=abs(F(h,R)-p)>width
                bad_pair+=weight if pair_failure else 0
                bad_mean+=weight if interval_failure else 0
                bad_joint+=weight if pair_failure or interval_failure else 0
            assert total==denominator,(R,a,'probability normalisation')
            pair=F(bad_pair,denominator);mean=F(bad_mean,denominator);joint=F(bad_joint,denominator)
            assert pair<=av,(R,p,'pair coverage',pair)
            # Random-variance width is controlled jointly, not by the nominal
            # alpha_mean alone. Its actual finite-example error may be smaller.
            assert mean<=av+am,(R,p,'final interval coverage',mean)
            assert joint<=av+am,(R,p,'joint pair/mean coverage',joint)
            nonzero_pair_failure+=bool(pair);nonzero_mean_failure+=bool(mean)
            cases.append({'replicates':R,'probability_fraction':str(p),
                'enumerated_pair_types':len(data),'normalising_denominator':str(denominator),
                'pair_failure_probability_fraction':str(pair),
                'mean_interval_failure_probability_fraction':str(mean),
                'pair_or_mean_failure_probability_fraction':str(joint),
                'nominal_pair_failure_fraction':str(av),'nominal_joint_failure_fraction':str(av+am)})
            probability_types+=len(data)
    return {'status':'all exact Bernoulli pair-type coverage checks passed',
        'laws':len(cases),'distinct_pair_types':type_count,
        'law_pair_type_probability_terms':probability_types,
        'laws_with_nonzero_pair_failure':nonzero_pair_failure,
        'laws_with_nonzero_interval_failure':nonzero_mean_failure,
        'maximum_pair_failure_probability_fraction':str(max(F(r['pair_failure_probability_fraction']) for r in cases)),
        'maximum_joint_failure_probability_fraction':str(max(F(r['pair_or_mean_failure_probability_fraction']) for r in cases)),
        'probability_arithmetic':'Exact integer probability numerators; rational p=a/32, a=0,...,32; no simulated data',
        'confidence_log_source':'Actual retained 160-digit outward logarithms, independently checked at 300 digits' if log_variance_upper is not None else 'Separate 300-digit outward logarithms',
        'cases':cases,'scope':'Finite Bernoulli implementation examples, not a distribution-free proof'}


def retained_audit(root):
    root=Path(root);protocol_path=root/'results/rqmc_pricing_protocol.json';results_path=root/'results/rqmc_pricing_results.json'
    protocol=json.loads(protocol_path.read_text());results=json.loads(results_path.read_text())
    assert results['protocol_sha256']==sha(protocol_path)
    for source,digest in results['source_hashes'].items():
        assert (root/source).is_file(),('Missing retained source',source)
        assert sha(root/source)==digest,('Retained source hash mismatch',source)
        assert protocol['source_hashes'][source]==digest,('Protocol/result source hash mismatch',source)
    checked_hashes=dict(results['source_hashes'])
    selected={r['row']['case_id']:r for r in protocol['selected']}
    av=am=F(protocol['alpha_per_event_exact_fraction']);assert av==am==F(1,2200)
    checkpoints=protocol['checkpoint_replicates'];assert checkpoints==[1<<k for k in range(6,17)]
    assert len(checkpoints)*(av+am)==F(protocol['alpha_per_price_exact_fraction'])==F(1,100)
    records=[];count=0;empty=0
    tolerance=Decimal('1e-280')
    for price in results['configurations']:
        prior=selected[price['case_id']]
        rho=F(prior['finite_output_range_exact_fraction']);variance=F(prior['batch_variance_prior_exact_fraction']);mean=F(prior['batch_mean_prior_exact_fraction'])
        assert 0<=variance<=rho*rho/4 and 0<=mean<=rho<=1
        batch=list(map(int,price['batch_integer_sums']))
        assert len(batch)==price['replicates']
        last_fraction=F(0)
        for cindex,checkpoint in enumerate(price['checkpoints']):
            R=checkpoint['replicates'];assert R==checkpoints[cindex]
            statistics=checkpoint['exact_batch_statistics'];stored=checkpoint['confidence']
            D=int(statistics['batch_denominator']);active=batch[:R];m=R//2
            assert max(active,default=0)<=rho*D
            S=sum(active);S2=sum(x*x for x in active);Q=sum((active[2*i]-active[2*i+1])**2 for i in range(m))
            assert int(statistics['replicates'])==R and int(statistics['independent_pairs'])==m
            assert int(statistics['integer_batch_sum'])==S
            assert int(statistics['integer_batch_square_sum'])==S2
            assert int(statistics['integer_paired_difference_square_sum'])==Q
            assert F(statistics['finite_mean_exact_fraction'])==F(S,R*D)
            assert F(statistics['paired_variance_estimate_exact_fraction'])==F(Q,2*m*D*D)
            assert F(statistics['ordinary_replicate_sample_variance_exact_fraction'])==F(R*S2-S*S,R*(R-1)*D*D)
            ref=reference(rho,variance,mean,statistics,av,am)
            with localcontext() as ctx:
                ctx.prec=310
                assert dec(stored['log_variance_upper_fraction'])>ref['log_variance']
                assert dec(stored['log_mean_upper_fraction'])>ref['log_mean']
                assert dec(stored['pair_variance_upper_fraction'])+tolerance>=ref['pair']
                assert dec(stored['variance_upper_fraction'])+tolerance>=ref['variance']
                assert dec(stored['unclipped_mean_halfwidth_fraction'])+tolerance>=ref['width']
            assert F(stored['variance_upper_fraction'])==min(variance,F(stored['pair_variance_upper_fraction']))
            actual_width=F(stored['unclipped_mean_halfwidth_fraction']);sample=F(S,R*D)
            lo=max(F(0),sample-actual_width);hi=min(mean,sample+actual_width)
            assert F(stored['failure_budget_fraction'])==av+am
            assert F(stored['alpha_variance_fraction'])==av and F(stored['alpha_mean_fraction'])==am
            if lo>hi:
                empty+=1;assert stored['status']=='empty_intersection'
                assert stored['mean_lower_fraction'] is stored['mean_upper_fraction'] is None
                assert checkpoint['tolerance_met'] is False
            else:
                assert stored['status']=='nonempty'
                assert F(stored['mean_lower_fraction'])==lo and F(stored['mean_upper_fraction'])==hi
                assert F(stored['halfwidth_upper_fraction'])==(hi-lo)/2
                assert F(stored['estimate_fraction'])==(hi+lo)/2
                with localcontext() as ctx:
                    ctx.prec=310
                    assert dec(lo)<=ref['lower']+tolerance
                    assert dec(hi)+tolerance>=ref['upper']
                assert checkpoint['tolerance_met']==((hi-lo)/2<=F(prior['estimation_tolerance_exact_fraction']))
                last_fraction=(hi-lo)/2
            count+=1
        met=price['checkpoints'][-1]['tolerance_met']
        assert all(not c['tolerance_met'] for c in price['checkpoints'][:-1])
        assert met or price['replicates']==checkpoints[-1]
        assert (price['status']=='certified per-price tolerance met')==met
        records.append({'case_id':price['case_id'],'checkpoint_count':len(price['checkpoints']),
            'replicates':price['replicates'],'status':price['status'],
            'last_finite_mean_halfwidth_fraction':str(last_fraction)})
    return {'status':'all retained confidence checkpoint statistics and endpoints independently verified',
        'price_count':len(records),'checkpoint_count':count,'empty_intersections':empty,
        'reference_decimal_digits':300,'comparison_allowance_decimal':str(tolerance),
        'endpoint_replay':'Separate 300-digit logarithm and square-root evaluation; no rqmc_confidence imports',
        'protocol_sha256':sha(protocol_path),'pricing_results_sha256':sha(results_path),
        'checked_retained_source_hashes':checked_hashes,
        'records':records,
        'scope':'Verifies retained statistics and outward numerical formulae, not a claim that seeded pseudo-random shifts prove ideal coverage'}


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--root',default=str(Path(__file__).resolve().parents[1]));parser.add_argument('--out',default=str(Path(__file__).resolve().parents[1]/'results/rqmc_confidence_audit.json'))
    args=parser.parse_args();retained=retained_audit(args.root)
    first=json.loads((Path(args.root)/'results/rqmc_pricing_results.json').read_text())['configurations'][0]['checkpoints'][0]['confidence']
    bernoulli=exact_bernoulli_coverage(first['log_variance_upper_fraction'],first['log_mean_upper_fraction'])
    result={'status':'independent confidence audit passed','retained_checkpoint_audit':retained,
        'exact_bernoulli_coverage_audit':bernoulli,'audit_source_sha256':sha(__file__)}
    Path(args.out).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':result['status'],'prices':retained['price_count'],
        'checkpoints':retained['checkpoint_count'],'bernoulli_laws':bernoulli['laws'],
        'exact_probability_terms':bernoulli['law_pair_type_probability_terms'],
        'laws_with_nonzero_interval_failure':bernoulli['laws_with_nonzero_interval_failure']},indent=2))

if __name__=='__main__':main()
