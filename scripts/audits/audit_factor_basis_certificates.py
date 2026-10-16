"""Audit the complete registered rational-factor certificate pool.

No producer transformation or transfer routine is imported. Exact pair
invariants and derivatives use audit_factor_basis. The retained independently
audited one-dimensional Gaussian-cell library is replayed fresh for the
common raw transport masses; every candidate's D1/D2 and allowance is rebuilt.
"""
import argparse
import hashlib
import json
import sys
from fractions import Fraction as F
from pathlib import Path
import audit_factor_basis as independent


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def digest(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':'),default=str).encode()).hexdigest()


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--repo',type=Path,default=independent.DEFAULT_REPO)
    ap.add_argument('--certificate-dir',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args()
    sys.path.insert(0,str(args.repo))
    from flint import ctx
    import tensor_transfer as tt
    from boundary_grid_certificate import endpoints, serial, upper
    from boundary_derivative_certificate import payoff_norms
    ctx.prec=512
    pd=args.certificate_dir/'results'
    protocol=json.loads((pd/'factor_basis_protocol.json').read_text())
    result=json.loads((pd/'factor_basis_certificates.json').read_text())
    summary=json.loads((pd/'factor_basis_summary.json').read_text())
    original_by_id={m['model_id']:m for m in protocol['models']}
    assert len(original_by_id)==4
    expected_frozen=json.loads((args.repo/'results/four_survivor_quantum_protocol.json').read_text())['models']
    assert protocol['models']==expected_frozen
    assert protocol['financial_law']=={'L':'6','b_y':10,'q':28,'precision_bits':512}
    assert protocol['cap']=='1' and protocol['financial_price_answer_excluded'] is True
    assert protocol['registered_before_fresh_full_certificate_execution'] is True
    assert protocol['earlier_proxy_diagnostics_known'] is True
    assert result['financial_price_evaluations']==0
    assert result['protocol_sha256']==sha(pd/'factor_basis_protocol.json')
    assert result['source_sha256']==protocol['source_sha256']==sha(args.certificate_dir/'factor_basis_certificate.py')
    assert result['dependency_source_sha256']==protocol['dependency_source_sha256']
    assert result['input_sha256']==protocol['input_sha256']
    for name,hashvalue in protocol['dependency_source_sha256'].items():
        assert sha(args.repo/name)==hashvalue,name
    for name,hashvalue in protocol['input_sha256'].items():
        assert sha(args.repo/name)==hashvalue,name
    assert summary['certificate_sha256']==sha(pd/'factor_basis_certificates.json')
    # Fresh common cell transport, independently of the new producer.
    transport=tt.TransportMGF(6,10,28)
    b1=upper(transport.rho1(0))
    b2=upper(transport.rho2mass)
    fresh_transport_profile=transport.profile()
    retained_transport_profile=dict(result['transport_profile'])
    for profile in (fresh_transport_profile,retained_transport_profile):
        profile['node_certificate']=dict(profile['node_certificate'])
        profile['node_certificate'].pop('generation_seconds_diagnostic',None)
    assert retained_transport_profile==fresh_transport_profile
    assert result['actual_financial_node_integers']==[c['actual_node_integer'] for c in transport.cells]
    exact_first_dominating_mass=sum((max(abs(F(c['first_moment']['lower_exact_fraction'])),abs(F(c['first_moment']['upper_exact_fraction']))) for c in transport.cells),F(0))
    assert b1>=exact_first_dominating_mass
    expected_keys=[]
    for original in protocol['models']:
        p=original['outer_dimension']
        candidate_ids=['identity']
        for a in range(p):
            for b in range(a+1,p):
                candidate_ids.extend((f'plane_{a}_{b}_minus',f'plane_{a}_{b}_plus'))
        pool=protocol['candidate_pools'][original['model_id']]
        assert [c['id'] for c in pool]==candidate_ids
        assert len(pool)==1+p*(p-1)
        expected_keys.extend((original['model_id'],cid) for cid in candidate_ids)
    assert protocol['candidate_count']==len(expected_keys)==58
    actual_keys=[(r['original_model_id'],r['candidate_id']) for r in result['records']]
    assert actual_keys==expected_keys and len(set(actual_keys))==58
    audits=[]
    grouped={mid:[] for mid in original_by_id}
    norms={mid:payoff_norms(m['lambda'],F(1)) for mid,m in original_by_id.items()}
    coordinate_count=0
    for row in result['records']:
        mid,cid=row['original_model_id'],row['candidate_id']
        original=original_by_id[mid]
        candidate=next(c for c in protocol['candidate_pools'][mid] if c['id']==cid)
        assert candidate==row['candidate']
        p=original['outer_dimension']
        o=independent.identity(p)
        if candidate['plane'] is not None:
            a,b=candidate['plane'];c,s=F(candidate['cosine']),F(candidate['sine'])
            assert c==F(3,5) and s in (F(-4,5),F(4,5)) and c*c+s*s==1
            o[a][a]=o[b][b]=c;o[a][b]=-s;o[b][a]=s
        else:
            assert cid=='identity' and candidate['cosine']=='1' and candidate['sine']=='0'
        assert [[str(x) for x in rr] for rr in o]==row['invariants']['rotation_matrix']
        transformed=row['rotated_model']
        assert transformed['model_id']==mid if cid=='identity' else transformed['model_id']==mid+'__'+cid
        pair=independent.audit_pair(original,transformed,o)
        w,r,mu,lam,mg,rg,d,s=independent.parameters(transformed)
        invariant=row['invariants']
        gg=independent.gram(r)
        assert invariant['asset_outer_gram_sha256']==digest([[str(x) for x in rr] for rr in gg])
        assert invariant['full_covariance_sha256']==digest([[str(x) for x in rr] for rr in independent.covariance(r,lam)])
        assert invariant['centred_loading_gram_sha256']==digest([[str(x) for x in rr] for rr in independent.gram(s)])
        assert invariant['geometric_centred_cross_products']==[str(independent.dot(rg,rr)) for rr in s]
        assert F(invariant['geometric_outer_variance'])==independent.dot(rg,rg)
        assert F(invariant['effective_geometric_variance'])==lam*lam+independent.dot(rg,rg)
        assert invariant['unchanged_log_drifts']==mu or invariant['unchanged_log_drifts']==transformed['log_drifts']
        assert F(invariant['unchanged_geometric_mean_log_drift'])==mg
        assert invariant['finite_grid_invariance_claimed'] is False
        M1,M2=F(norms[mid]['M1_upper']),F(norms[mid]['M2_upper'])
        other=upper(transport.mass**(p-1))
        mass=serial(transport.mass**p)
        tail=upper(1-transport.mass**p)
        expected_rows=[];linear=quadratic=geom1=geom2=F(0)
        for j in range(p):
            A=max(abs(rr[j]) for rr in r)
            rj=rg[j]
            d1=M1*(A+abs(rj));d2=(M1+M2)*(A*A+rj*rj)
            ll=other*b1*d1;qq=other*b2*d2
            linear+=ll;quadratic+=qq;geom1+=A+abs(rj);geom2+=A*A+rj*rj
            expected_rows.append({'coordinate':j,'full_asset_loading_max':str(A),'geometric_loading':str(rj),
                'D1_upper':str(d1),'D2_upper':str(d2),'linear_error_upper':str(ll),'quadratic_error_upper':str(qq)})
            coordinate_count+=1
        actual=row['direct_global_transfer']
        assert actual['model_id']==transformed['model_id'] and actual['cap']=='1' and actual['payoff']=='capped' and actual['order']==1
        assert actual['coordinates']==expected_rows
        for key,value in {'linear_transport_upper':linear,'quadratic_transport_upper':quadratic,'tail_upper':tail,
                          'financial_transfer_upper':linear+quadratic+tail,'rho1_mass_upper':b1,'rho2_mass_upper':b2,
                          'other_coordinate_mass_upper':other,'geometric_B1':geom1,'geometric_B2':geom2}.items():
            assert F(actual[key])==value,(mid,cid,key)
        assert actual['mass_interval']==mass
        assert row['resource_selection_used_price_answer'] is False
        pair.update({'original_model_id':mid,'candidate_id':cid,'direct_transfer_coordinate_checks':p,
                     'direct_transfer_upper':str(linear+quadratic+tail)})
        audits.append(pair);grouped[mid].append(row)
    assert summary['analytic_candidates_retained']==58
    assert len(result['selected_models'])==4 and len(summary['records'])==4
    for mid,family in grouped.items():
        best=min(family,key=lambda x:(F(x['direct_global_transfer']['financial_transfer_upper']),x['candidate_id']))
        ident=next(row for row in family if row['candidate_id']=='identity')
        selected=next(row for row in result['selected_models'] if row['original_model_id']==mid)
        assert best==selected
        summary_row=next(row for row in summary['records'] if row['model_id']==mid)
        assert summary_row['analytic_candidates']==len(family) and summary_row['selected_candidate']==best['candidate_id']
        b,q=ident['direct_global_transfer'],best['direct_global_transfer']
        assert summary_row['identity_transfer']==b['financial_transfer_upper']
        assert summary_row['selected_transfer']==q['financial_transfer_upper']
        for key,field in (('transfer_reduction_fraction','financial_transfer_upper'),('B1_reduction_fraction','geometric_B1'),('B2_reduction_fraction','geometric_B2')):
            assert F(summary_row[key])==1-F(q[field])/F(b[field])
        assert F(summary_row['transfer_reduction_fraction'])>=0
        assert summary_row['screen_choice_is_global_resource_optimum'] is False
        assert summary_row['identity_retained_for_full_resource_comparison'] is True
    # Sanity lower bound independently of outward masses, in every basis.
    bound_info={'exact_first_dominating_mass':str(exact_first_dominating_mass),
                'replayed_first_mass_upper':str(b1),'replayed_second_mass_upper':str(b2)}
    output={'schema':'qbp-independent-complete-factor-basis-certificate-audit-v1','status':'passed',
        'registered_candidate_count':58,'direct_transfer_coordinate_checks':coordinate_count,
        'multinomial_exponent_checks':sum(a['multinomial_exponent_checks'] for a in audits),
        'pointwise_affine_checks':sum(a['pointwise_affine_checks'] for a in audits),
        'derivative_scalar_checks':sum(a['derivative_scalar_checks'] for a in audits),
        'selected_certificate_count':4,'pairs':audits,'raw_transport_check':bound_info,
        'finite_grid_counterexample':independent.finite_counterexample(),
        'source_sha256':sha(__file__),'mathematical_auditor_sha256':sha(Path(independent.__file__)),
        'protocol_sha256':sha(pd/'factor_basis_protocol.json'),
        'certificate_sha256':sha(pd/'factor_basis_certificates.json'),
        'summary_sha256':sha(pd/'factor_basis_summary.json'),
        'producer_source_sha256':sha(args.certificate_dir/'factor_basis_certificate.py'),
        'dependency_source_sha256':protocol['dependency_source_sha256'],'input_sha256':protocol['input_sha256'],
        'scope':'Exact model, covariance, derivative and complete direct-transfer reconstruction for every registered analytic basis. Shared previously audited cell primitives replayed fresh. Complete quantum resource reconstruction and originality are outside this record.'}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps({k:output[k] for k in ('status','registered_candidate_count','direct_transfer_coordinate_checks','multinomial_exponent_checks','derivative_scalar_checks')}))

if __name__=='__main__':
    main()
