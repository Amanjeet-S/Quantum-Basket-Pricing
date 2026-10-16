"""Additive literal-predicate depth and live-register resource check.

The primary independent graph/prior/T audit and all frozen outcomes are
preserved. This supplement replays only the already phase-checked literal
clean predicate at full 36/38-bit widths and derives depth/live caps from the
retained emitted blocks. It never invokes the producer cost frontend.
"""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
sys.dont_write_bytecode=True


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--step',type=Path,required=True)
    ap.add_argument('--repo',type=Path,required=True);ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args();q=args.step/'adversarial/quantum';source=q/'shared_threshold_oracle.py'
    datafile=q/'threshold_results/shared_threshold_resources.json';data=json.loads(datafile.read_text())
    primarypath=args.output.parent/'shared_threshold_financial_audit.json';primary=json.loads(primarypath.read_text())
    phasepath=args.output.parent/'shared_threshold_predicate_audit.json';phase=json.loads(phasepath.read_text())
    assert primary['status']==phase['status']=='passed'
    assert data['source_sha256']==primary['prototype_source_sha256']==phase['prototype_source_sha256']==sha(source)
    assert primary['resource_sha256']==sha(datafile)
    spec=importlib.util.spec_from_file_location('literal_checked_predicate',source)
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
    m.configure(args.repo,args.step/'resources/results/factor_basis_quantum_resources.json',args.output.parent/'configuration_only')
    facts=[]
    for row in data['rows']:
        for architecture,family in row['profiles'].items():
            profile=family['centred'];width=profile['integer_word_width']
            c=m.vo.ValueCircuit(store_gates=False,mode=architecture)
            lower,threshold,upper=(c.alloc(width) for _ in range(3));flag=c.alloc(1)[0]
            original=c.active_wires;m.emit_interval(c,lower,threshold,upper,flag)
            assert c.active_wires==original
            metrics=c.resources()
            extra=metrics['wires']-(3*width+1)
            assert profile['predicate_extra_qubits_upper']==extra
            assert profile['predicate_primitive_counts']==metrics['clifford_T_expanded_gate_counts']
            be=sum(b['explicit_T_depth_serial_upper'] for b in profile['blocks'])
            bf=sum(b['conservative_5layer_T_depth_serial_upper'] for b in profile['blocks'])
            assert profile['explicit_T_depth_serial_upper']==2*be+metrics['explicit_7T_sequence_T_depth']
            assert profile['conservative_5layer_T_depth_serial_upper']==2*bf+metrics['conservative_5layer_CCX_T_depth_upper']
            retained=profile['uniform_H']+1;peak=retained
            for block in profile['blocks']:
                peak=max(peak,retained+block['extra_qubits_upper']);retained+=block['output_bits']
            peak=max(peak,retained+extra)
            assert peak==profile['logical_qubits_upper']
            assert profile['reflection_core_qubits']==profile['active_qubits_after_clean_prep']==profile['uniform_H']+1
            for strategy in row['strategies']:
                if strategy['architecture']!=architecture:continue
                for event in strategy['result'].get('flags',[]):
                    mbits=event['plan']['phase_register_bits'];expectedpeak=max(peak,2*profile['reflection_core_qubits'])+mbits
                    assert event['complete_charge']['peak_logical_qubits']==expectedpeak
                    assert event['complete_charge']['metrics']['measurements']==event['plan']['repetitions']*expectedpeak
                    assert event['complete_charge']['metrics']['reset_X']==event['plan']['repetitions']*expectedpeak
            facts.append({'original_model_id':row['original_model_id'],'architecture':architecture,
                'predicate_width':width,'predicate_extra_qubits_upper':extra,
                'financial_live_peak':peak,'reflected_core':profile['reflection_core_qubits'],
                'explicit_predicate_T_depth':metrics['explicit_7T_sequence_T_depth'],
                'conservative_predicate_T_depth':metrics['conservative_5layer_CCX_T_depth_upper']})
    result={'schema':'qbp-independent-shared-threshold-depth-register-supplement-v1','status':'passed','profiles_checked':len(facts),
        'facts':facts,'source_sha256':sha(__file__),'prototype_source_sha256':sha(source),'resource_sha256':sha(datafile),
        'primary_financial_audit_sha256':sha(primarypath),'exact_phase_audit_sha256':sha(phasepath),
        'scope':'Additive replay of the phase-checked literal full-width predicate and independent serial-depth/live-register/reset composition. It preserves original outcomes and primary audits. Depths and live allocations are conservative caps, not parallel scheduling or physical hardware measurements.'}
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':'passed','profiles_checked':len(facts)}))

if __name__=='__main__':main()
