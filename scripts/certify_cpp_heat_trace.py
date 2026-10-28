"""Audit small native heat traces using exact rational discrete residuals.

This certifies an exported-coefficient discrete target, not the continuum price.
All cases are synthetic. Retained evidence must never be overwritten.
"""
from __future__ import annotations

import argparse
import bisect
from dataclasses import asdict
from fractions import Fraction
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import platform
import subprocess
import sys


MESHES = ((2, 1), (16, 3), (32, 16))


def load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def solve_exact(previous, lam, cap):
    """Small exact Gaussian elimination, independent of native Thomas arithmetic."""
    n = len(previous)
    matrix = [[1+2*lam if i == j else -lam if abs(i-j) == 1 else Fraction(0)
               for j in range(n)] for i in range(n)]
    rhs = list(previous)
    rhs[-1] += lam*cap
    for i in range(n):
        for j in range(i+1, n):
            multiplier = matrix[j][i]/matrix[i][i]
            if multiplier == 0:
                continue
            for k in range(i, n):
                matrix[j][k] -= multiplier*matrix[i][k]
            rhs[j] -= multiplier*rhs[i]
    result = [Fraction(0)]*n
    for i in reversed(range(n)):
        result[i] = (rhs[i]-sum(matrix[i][j]*result[j] for j in range(i+1,n)))/matrix[i][i]
    return result


def run(root: Path, executable: Path, output: Path):
    if output.exists():
        raise ValueError('Choose a new output directory; retained evidence cannot be overwritten.')
    if output.is_relative_to(root):
        raise ValueError('Generate fresh evidence outside the repository before review and installation.')
    engine = load(root/'src/capped_gbm.py', 'qbp_trace_contracts')
    cert = load(root/'src/cpp_residual_certificate.py', 'qbp_trace_certificate')
    sources = [Path(__file__).resolve(), root/'src/capped_gbm.py', root/'src/cpp_residual_certificate.py',
               root/'cpp/include/qbp/capped_gbm.hpp', root/'cpp/src/capped_gbm.cpp', root/'cpp/src/main.cpp']
    hashes = {str(path.relative_to(root)): digest(path) for path in sources}
    executable_hash = digest(executable)
    output.mkdir(parents=True)
    protocol = {
        'reporting_amendment': 'An initial run found four strict range overshoots. This version adds exact excess diagnostics and checks their residual coverage, without relaxing the strict range acceptance condition or removing failures.',
        'schema': 'qbp-native-trace-protocol-v1',
        'cases': [{'name': name, 'contract': asdict(c)} for name,c in engine.synthetic_contracts()],
        'meshes': [{'intervals': m,'time_steps': n} for m,n in MESHES],
        'target': 'Exact recurrence using binary64-parsed exported lambda, cap and initial vector; fixed Dirichlet boundaries 0/cap.',
        'scalar_target': 'Exact discrete final interpolation at the separately recorded binary64 index, weight and discount.',
        'independent_check': 'Exact rational Gaussian elimination on the one-interior-node, one-step meshes; compare both grid and scalar errors with exact residual bounds.',
        'acceptance': 'All dimensions/boundaries finite and valid; no observed range violations; all independent exact comparisons satisfy their rational bounds.',
        'limitations': 'No original-parameter, coordinate, payoff-construction, continuum, finite-domain, model or Greek certificate; largest meshes are trace audits, not exact target solves.',
        'arithmetic': 'Exact Fraction after parsing finite JSON values to binary64; upward-rounded finite binary64 summaries accompany exact fractions.',
    }
    (output/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
    summary, failures = [], []
    for name,c in engine.synthetic_contracts():
        lower, upper = engine.study_domain(c)
        for intervals, steps in MESHES:
            command = [str(executable),'pde','--spot',repr(c.spot),'--strike',repr(c.strike),
                       '--cap',repr(c.cap),'--maturity',repr(c.maturity),'--rate',repr(c.rate),
                       '--yield',repr(c.dividend),'--volatility',repr(c.volatility),
                       '--lower',repr(lower),'--upper',repr(upper),'--intervals',str(intervals),
                       '--time-steps',str(steps),'--dump-trace']
            payload = json.loads(subprocess.run(command,check=True,capture_output=True,text=True).stdout)
            nodes, trace = payload['grid']['nodes'], payload['heat_trace']
            if len(trace) != steps+1 or len(nodes) != intervals+1:
                raise ValueError('Native trace dimensions disagree with the requested protocol.')
            if any(not math.isfinite(x) for x in nodes) or any(b <= a for a,b in zip(nodes,nodes[1:])):
                raise ValueError('Published interpolation nodes must be finite and strictly increasing.')
            y = math.log(c.spot/c.strike)+c.log_drift*c.maturity
            if not nodes[0] <= y <= nodes[-1]:
                raise ValueError('The recorded scalar target lies outside the published domain.')
            left = min(bisect.bisect_right(nodes,y)-1, len(nodes)-2)
            weight = (y-nodes[left])/(nodes[left+1]-nodes[left])
            discount = math.exp(-c.rate*c.maturity)
            record = cert.certify_trace(trace,payload['grid']['lambda'],c.cap,
                                        price=payload['values']['price'],left_index=left,
                                        right_weight=weight,discount=discount)
            def exact(item):
                return Fraction(int(item['numerator']),int(item['denominator']))
            independent = None
            if intervals == 2:
                target = [Fraction.from_float(x) for x in trace[0][1:-1]]
                lam, cap = Fraction.from_float(payload['grid']['lambda']),Fraction.from_float(c.cap)
                for _ in range(steps):
                    target = solve_exact(target,lam,cap)
                actual_grid = max(abs(Fraction.from_float(x)-v) for x,v in zip(trace[-1][1:-1],target))
                full_target = [Fraction(0),*target,cap]
                theta, factor = Fraction.from_float(weight),Fraction.from_float(discount)
                target_scalar = factor*((1-theta)*full_target[left]+theta*full_target[left+1])
                actual_scalar = abs(Fraction.from_float(payload['values']['price'])-target_scalar)
                independent = {'grid_error':cert.record(actual_grid),'scalar_error':cert.record(actual_scalar),
                               'grid_bound_holds':actual_grid <= exact(record['final_grid_error_bound']),
                               'scalar_bound_holds':actual_scalar <= exact(record['scalar_price_error_bound'])}
                if not independent['grid_bound_holds'] or not independent['scalar_bound_holds']:
                    failures.append(f'{name}:{intervals}:{steps}:exact_check')
            if record['range_excess_within_residual_bound'] is not True:
                failures.append(f'{name}:{intervals}:{steps}:residual_range_coverage')
            if record['observed_range_violation']:
                failures.append(f'{name}:{intervals}:{steps}:range')
            stem = f'{name}_m{intervals}_n{steps}'
            # Timing is incidental and has no role in certificate acceptance.
            payload.pop('elapsed_seconds',None)
            document = {'contract':asdict(c),'native':payload,
                        'scalar_target_data':{'left_index':left,'right_weight':weight,'discount':discount,'coordinate':y},
                        'certificate':record,'independent_check':independent}
            (output/f'{stem}.json').write_text(json.dumps(document,indent=2,allow_nan=False)+'\n')
            summary.append({'case':name,'intervals':intervals,'steps':steps,
                            'grid_error_bound':record['final_grid_error_bound']['upper_binary64'],
                            'scalar_error_bound':record['scalar_price_error_bound']['upper_binary64'],
                            'observed_range_violation':record['observed_range_violation'],
                            'max_range_excess':record['max_observed_range_excess']['upper_binary64'],
                            'range_excess_within_residual_bound':record['range_excess_within_residual_bound']})
    for path in sources:
        if digest(path) != hashes[str(path.relative_to(root))]:
            raise RuntimeError('A source changed during the audit; retain diagnostics and rerun with fixed sources.')
    if digest(executable) != executable_hash:
        raise RuntimeError('The native executable changed during the audit.')
    manifest = {'schema':'qbp-native-trace-study-v1','status':'pass' if not failures else 'fail',
                'failures':failures,'records':summary,'trace_count':len(summary),
                'strict_range_failures':sum(row['observed_range_violation'] for row in summary),
                'residual_range_coverage_all':all(row['range_excess_within_residual_bound'] for row in summary),
                'max_observed_range_excess':max(row['max_range_excess'] for row in summary),
                'independent_exact_target_checks':sum(row['intervals']==2 for row in summary),
                'max_grid_error_bound':max(row['grid_error_bound'] for row in summary),
                'max_scalar_error_bound':max(row['scalar_error_bound'] for row in summary),
                'python':platform.python_version(),'source_sha256':hashes,'executable_sha256':executable_hash,
                'files':{p.name:{'sha256':digest(p),'bytes':p.stat().st_size} for p in sorted(output.iterdir())}}
    (output/'study_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    return manifest


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project-root',type=Path,default=Path(__file__).resolve().parents[1])
    parser.add_argument('--executable',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    report = run(args.project_root.resolve(),args.executable.resolve(),args.output.resolve())
    print(json.dumps({key:report[key] for key in ('status','trace_count','independent_exact_target_checks','max_grid_error_bound','max_scalar_error_bound')},indent=2))
    raise SystemExit(0 if report['status']=='pass' else 1)
