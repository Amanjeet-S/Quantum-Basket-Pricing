"""Small bounded-support illustration; no financial or empirical benchmark."""
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal,localcontext,ROUND_HALF_EVEN
import argparse,hashlib,json,sys
parser=argparse.ArgumentParser()
parser.add_argument('--root',type=Path,required=True)
parser.add_argument('--output',type=Path,default=Path(__file__).with_name('weak_dispersion_range_examples.json'))
args=parser.parse_args();S=args.root.resolve()
sys.path.insert(0,str(S))
import iterative_estimator as tf
import full_batch_miqae as mi
out=[];endpoint_checks=0
zero_support_checks={}
for name,module in [('three_five',tf),('miqae',mi)]:
 p=module.plan(F(1,32),F(1,100),F(0))
 assert p.total_shots_upper==p.total_grover_calls_upper==0
 zero_support_checks[name]=True
for q in range(1,65):
 values={1,1<<q}
 for b in range(q+1):
  for a in ((1<<b)-1,1<<b,(1<<b)+1):
   if 1<=a<=1<<q:values.add(a)
 for a in values:
  c=(a-1).bit_length();R=F(1<<c,1<<q);A=F(a,1<<q)
  assert 0<=c<=q and A<=R<2*A
  assert c==0 if a==1 else c>0
  assert c==q if a==1<<q else c<=q
  endpoint_checks+=1
for n in range(1,21):
 eps=F(1,1<<n);q=2*n+3;Delta=F(1,1<<q)
 # Uniform input addresses yield V=-1,0,0,1 and centred h=(eps*V,-eps*V).
 # P(m)=m, m=0, H=L=1 are a mathematical test function, not a financial model.
 # For eps<=1/2, exp(eps)<2 follows by termwise comparison with sum (1/2)^k.
 # The proposition gives J<=eps**2. Nearest-grid arithmetic error is Delta/2.
 rho=eps*eps+Delta/2;scaled=rho*(1<<q);amax=scaled.numerator//scaled.denominator
 assert amax==8 and (amax-1).bit_length()==3
 R=F(8,1<<q);assert R==eps*eps
 t=eps*eps/32;tolerance=t/R;assert tolerance==F(1,32)
 with localcontext() as ctx:
  ctx.prec=220
  e=Decimal(eps.numerator)/Decimal(eps.denominator)
  gap=((e.exp()+(-e).exp())/2).ln()
  actual_a=int((gap*(1<<q)).to_integral_value(rounding=ROUND_HALF_EVEN))
  assert actual_a==4
 outputs=[actual_a,0,0,actual_a]
 mean=F(sum(outputs),4*(1<<q));var=sum((F(a,1<<q)-mean)**2 for a in outputs)/4
 assert mean==eps*eps/4 and var==eps**4/16
 original_flag_variance=mean*(1-mean)
 range_probability=mean/R;assert range_probability==F(1,4)
 caps={}
 for name,module in [('three_five',tf),('miqae',mi)]:
  p=module.plan(tolerance,F(1,100),F(1))
  caps[name]={'shots':p.total_shots_upper,'Q_calls':p.total_grover_calls_upper}
 out.append({'n':n,'epsilon':str(eps),'q':q,'rho':str(rho),'certificate_integer_maximum':amax,
  'attained_maximum_integer':actual_a,'range_R':str(R),'finite_mean_tolerance':str(t),'flag_tolerance':str(tolerance),
  'finite_mean':str(mean),'finite_variance':str(var),'canonical_flag_variance':str(original_flag_variance),
  'range_flag_probability':str(range_probability),'controller_caps':caps})
assert len({repr(x['controller_caps']) for x in out})==1
result={'status':'passed','date':'2026-10-06','edge_case_checks':endpoint_checks,'zero_support_controller_checks':zero_support_checks,'bounded_support_examples':out,
 'scope':'Exact rational range and controller planning. Gap-rounding values are checked at 220-digit Decimal precision as an illustration, not a new outward transcendental certificate or market-calibrated price. The analytic range bound and nearest-rounding error establish the stated conservative support certificate.',
 'proof_sha256':hashlib.sha256((S/'proof_sources/weak_dispersion_range_proof.tex').read_bytes()).hexdigest(),
 'controller_hashes':{n:hashlib.sha256((S/n).read_bytes()).hexdigest() for n in ('iterative_estimator.py','full_batch_miqae.py')},
 'check_source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 'limitations':'Mathematical test function P(m)=m only. No benchmark result, novelty assertion, uniform Gaussian range or total constant pricing cost is inferred.'}
args.output.parent.mkdir(parents=True,exist_ok=True)
args.output.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'status':result['status'],'edge_case_checks':endpoint_checks,'examples':len(out),'constant_controller_caps':out[0]['controller_caps'],'proof_sha256':result['proof_sha256']},indent=2))
