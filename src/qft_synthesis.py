"""Retained certified Clifford+T inverse QFT, with explicit global phases.

The production resource API uses only standard Python and retained certificates.
Wire 0 is the most significant computational-basis bit. Each gate tuple is in
chronological execution order. No controlled Grover operation is implemented here.
"""
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parent
@lru_cache(None)
def _data():
    p=ROOT/'results'/'qft_rotation_tapes.json'
    tapes=json.loads(p.read_text())
    cert=json.loads((ROOT/'results'/'qft_rotation_certificates.json').read_text())
    if hashlib.sha256(p.read_bytes()).hexdigest()!=cert['tapes_sha256']:
        raise ValueError('Rotation tape hash differs from its Acb certificate')
    bounds={r['key']:Fraction(r['operator_norm_upper_rational']) for r in cert['certificates']}
    if set(bounds)!=set(tapes['rotations']):raise ValueError('Incomplete certificates')
    if any(not r.get('passed') or r.get('precision_bits')!=384 for r in cert['certificates']):
        raise ValueError('Unaccepted rotation certificate')
    if any(not 0<=b<Fraction(1,2**17) for b in bounds.values()):
        raise ValueError('Invalid retained rotation error bound')
    return tapes,bounds

def rotation(k,sign=1):
    if sign not in (-1,1) or not 2<=k<=32:raise ValueError('Unsupported dyadic Rz')
    tapes,bounds=_data();key=f'{sign}:{k}';r=tapes['rotations'][key]
    return {'chronological_gates':tuple(r['chronological_gates']),
            'phase_pi':Fraction(r['phase_pi']),
            'norm_error_upper':bounds[key], 'k':k,'sign':sign}

def controlled_phase(distance,control,target,sign=-1):
    if not 1<=distance<=31 or control==target:raise ValueError('Invalid CP arguments')
    k=distance+1
    r=rotation(k,sign);s=rotation(k,-sign)
    tape=[]
    for g in r['chronological_gates']:tape.append((g,control))
    for g in r['chronological_gates']:tape.append((g,target))
    tape.append(('CX',control,target))
    for g in s['chronological_gates']:tape.append((g,target))
    tape.append(('CX',control,target))
    # When all three approximations are identities, the adjacent CNOTs cancel.
    if not r['chronological_gates'] and not s['chronological_gates']:tape=[]
    return {'gates':tuple(tape),
            'global_phase_correction_pi':Fraction(sign,2**(distance+2))+2*r['phase_pi']+s['phase_pi'],
            'norm_error_upper':2*r['norm_error_upper']+s['norm_error_upper'],
            'rotation_k':k,'identity_rotation_replacement':not r['chronological_gates']}

def inverse_qft_tape(phase_bits):
    m=int(phase_bits)
    if m!=phase_bits or not 1<=m<=32:raise ValueError('Supported inverse QFT widths are 1 to 32')
    tape=[];phase=Fraction(0);error=Fraction(0);cp=[]
    for j in range(m//2):
        k=m-1-j;tape.extend([('CX',j,k),('CX',k,j),('CX',j,k)])
    for j in range(m-1,-1,-1):
        for k in range(m-1,j,-1):
            record=controlled_phase(k-j,k,j,-1)
            tape.extend(record['gates']);phase+=record['global_phase_correction_pi'];error+=record['norm_error_upper']
            cp.append({'distance':k-j,'control':k,'target':j,'norm_error_upper':str(record['norm_error_upper']),'identity_rotation_replacement':record['identity_rotation_replacement']})
        tape.append(('H',j))
    return {'gates':tuple(tape),'global_phase_correction_pi':phase,'norm_error_upper':error,'controlled_phases':tuple(cp),'wire_convention':'wire 0 is the most significant basis bit; chronological gate order'}

@lru_cache(None)
def inverse_qft_profile(phase_bits,norm_budget=Fraction(1,64)):
    budget=Fraction(norm_budget);info=inverse_qft_tape(phase_bits)
    error=info['norm_error_upper']
    if budget<=0 or error>budget:raise ValueError('Retained QFT exceeds the requested operator norm budget')
    c=Counter(g[0]for g in info['gates']);nt=c['T']+c['Tdg'];depth=nt
    metrics={'T':nt,'measurements':0,'feedback_CZ':0,'reset_X':0,
             'explicit_T_depth_serial_upper':depth,
             'conservative_5layer_T_depth_serial_upper':depth,
             'expanded_primitive_operations':sum(c.values()),'uniform_H':0}
    metrics.update({f'gate_{g}':c[g]for g in ('H','S','Sdg','T','Tdg','CX')})
    return {'phase_register_bits':int(phase_bits),'norm_error_upper':error,
            'norm_error_upper_rational':str(error),'allocated_norm_budget':str(budget),
            'global_phase_correction_pi':str(info['global_phase_correction_pi']),
            'phase_relation':'exp(i*pi*global_phase_correction_pi) times the physical Clifford+T tape approximates the mathematical inverse QFT',
            'gate_counts':dict(c),'T_count':nt,'serial_T_depth_upper':depth,
            'metrics':metrics,'extra_ancillas':0,
            'controlled_phase_count':len(info['controlled_phases']),
            'identity_rotation_controlled_phases':sum(r['identity_rotation_replacement']for r in info['controlled_phases']),
            'proof':'Each retained Rz product has a fresh Acb operator-norm upper certificate. Unitary telescoping charges all three rotations per controlled phase. The inverse QFT is an uncontrolled final unitary, so its recorded aggregate scalar has no effect on outcome probabilities.'}

if __name__=='__main__':
    profiles=[inverse_qft_profile(m)for m in range(1,33)]
    (ROOT/'results'/'inverse_qft_profiles.json').write_text(json.dumps({'schema':'qbp-certified-inverse-qft-resources-v1','profiles':profiles},default=str,indent=2)+'\n')
    print('All 32 widths meet 1/64; m=32 T =',profiles[-1]['T_count'],'error =',float(profiles[-1]['norm_error_upper']))
