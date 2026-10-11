"""Generate retained actual tapes. The resource API itself needs no synth package."""
from pathlib import Path
from fractions import Fraction
import hashlib, json, sys, os
ROOT=Path(__file__).resolve().parent
VENDOR=Path(os.environ.get('QBP_SYNTHESIS_VENDOR',str(ROOT/'synthesis_vendor')))
sys.path.insert(0,str(VENDOR))
import mpmath as mp
from pygridsynth.gridsynth import gridsynth_gates
mp.mp.dps=100
out={"schema":"qbp-clifford-t-rz-tapes-v1","scope":"Chronological H/S/T and adjoint tapes. The exact rational multiple of pi in phase_pi multiplies the physical gate product and is part of the certified representation. These Rz rotations are used only as uncontrolled gates inside the inverse QFT.","generator":{"pygridsynth":"1.2.0","mpmath":"1.4.1","precision_decimal_digits":100,"seed":1729,"gridsynth_epsilon":"1/1048576","certified_required_upper":"1/131072","up_to_phase":False},"source_hashes":{},"rotations":{}}
for p in sorted((VENDOR/'pygridsynth').glob('*.py')):
    out['source_hashes'][p.name]=hashlib.sha256(p.read_bytes()).hexdigest()
for k in range(2,33):
    if k==2:
        chron=['T']; phase=Fraction(-1,8); provenance='Exact T identity'
    elif k>=19:
        chron=[]; phase=Fraction(0); provenance='Identity replacement, error charged explicitly'
    else:
        raw=gridsynth_gates(theta=mp.pi/(2**k),epsilon=mp.mpf(1)/(2**20),seed=1729,dps=100,up_to_phase=False)
        assert set(raw)<=set('HTSXW'),raw
        # The software displays a matrix product left to right. Reverse it for chronological execution.
        chron=list(reversed(raw.replace('W','')))
        # Replace X by H S S H, so the retained physical alphabet is H,S,T only.
        chron=[h for g in chron for h in (['H','S','S','H'] if g=='X' else [g])]
        phase=Fraction(raw.count('W'),4); provenance='pygridsynth matrix product converted to chronological execution'
    neg=[{'H':'H','S':'Sdg','T':'Tdg'}[g] for g in reversed(chron)]
    for sign,tape,p in [(1,chron,phase),(-1,neg,-phase)]:
        out['rotations'][f'{sign}:{k}']={'sign':sign,'k':k,'theta_pi':str(Fraction(sign,2**k)), 'phase_pi':str(p),'chronological_gates':tape,'gate_counts':{g:tape.count(g) for g in sorted(set(tape))},'T_count':sum(g in ('T','Tdg') for g in tape),'serial_T_depth':sum(g in ('T','Tdg') for g in tape),'provenance':provenance}
    print(k,len(chron),sum(g=='T' for g in chron),flush=True)
(ROOT/'results').mkdir(exist_ok=True)
(ROOT/'results'/'qft_rotation_tapes.json').write_text(json.dumps(out,indent=2)+'\n')
