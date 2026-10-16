"""Exact block/count and constant-alias audit for the frozen financial family."""
import json
import hashlib
import sys
from decimal import Decimal, localcontext
from pathlib import Path

STAGE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(STAGE))
import encoding_resources as e
import generalisation as g
from reversible import Circuit


def active_internal_nodes(table, bits):
    def visit(lo, hi, depth):
        if not any(table[lo:hi]) or depth == 0:
            return 0
        mid = (lo + hi) // 2
        return 1 + visit(lo, mid, depth - 1) + visit(mid, hi, depth - 1)
    return visit(0, 1 << bits, bits)


def test_function_block(w, f, bits, dense):
    q = f + bits - 2  # span 4
    c = Circuit(store_gates=False)
    a = c.alloc(w)
    # Monotone rounded endpoints, either everywhere positive or tail-zero.
    endpoint = [j + (1 if dense else -((1 << bits)//2)) for j in range((1 << bits)+1)]
    endpoint = [max(0, x) for x in endpoint]
    table = [endpoint[j] | ((endpoint[j+1]-endpoint[j]) << w) for j in range(1 << bits)]
    node = e.Node('function', ('input',), 'general_cdf:0.2', bits)
    e.emit_node(c, node, [a], q, w, table)
    I = active_internal_nodes(table, bits)
    expected = 2*I + 2*f*(4*w+3*f) + 6*w
    actual = c.resources()['toffoli']
    assert actual == expected, (w, f, bits, dense, actual, expected)
    return dict(w=w, f=f, bits=bits, dense=dense, internal_nodes=I, actual=actual, formula=expected)


def alias_audit(models):
    original = e.Graph.const
    collisions = []
    context = {}
    def audited(self, value):
        value = Decimal(str(value))
        integer = e.rounded_scaled(value, self.q) % (1 << self.width)
        if integer in self.constants:
            name = self.constants[integer]
            previous = self.ideal_constants[name]
            if value != previous:
                with localcontext() as ctx:
                    ctx.prec = 220
                    difference = abs(value-previous)
                collisions.append({**context, 'q': self.q, 'previous': str(previous),
                                   'new': str(value), 'difference': str(difference)})
        return original(self, value)
    e.Graph.const = audited
    try:
        for model in models:
            for q in range(18,31):
                for cap in (1,10,25):
                    for method in ('direct','interval','full_global','full_cap'):
                        context.clear(); context.update(model=model['model_id'], cap=cap, method=method)
                        g.build_graph(model,method,q,8,8,q-1,4,4,4,cap)
    finally:
        e.Graph.const = original
    maximum = max((Decimal(row['difference']) for row in collisions), default=Decimal(0))
    return {'graphs_checked':len(models)*13*3*4, 'nonidentical_ideal_constant_collisions':len(collisions),
            'maximum_ideal_difference':str(maximum), 'examples':collisions[:8]}


def graph_count_audit(model, q=18, bits=6):
    tables, cache = {}, {}
    full = g.build_graph(model,'full_global',q,6,0,0,5,5,bits,10)
    interval = g.build_graph(model,'interval',q,6,0,q-1,5,5,bits,10)
    g.make_tables(full,tables); g.make_tables(interval,tables)
    fm = e.counted_composition(full,tables,cache)
    im = e.counted_composition(interval,tables,cache)
    strike = next(b['counts']['Toffoli'] for b in im['blocks'] if b['function']=='general_strike:10')
    ft = next(b['counts']['Toffoli'] for b in fm['blocks'] if b['function'] and b['function'].startswith('general_cap:'))
    it = next(b['counts']['Toffoli'] for b in im['blocks'] if b['function'] and b['function'].startswith('general_cdf:'))
    w = full.width
    actual = im['expanded_gate_counts']['Toffoli']-fm['expanded_gate_counts']['Toffoli']
    expected = 2*strike+24*w+4*(it-ft)
    assert actual == expected, (actual,expected)
    coarse = g.build_graph(model,'full_global',q,6,0,0,5,5,bits-1,10)
    g.make_tables(coarse,tables)
    cm = e.counted_composition(coarse,tables,cache)
    coarseft = next(b['counts']['Toffoli'] for b in cm['blocks'] if b['function'] and b['function'].startswith('general_cap:'))
    actual_coarsening = fm['expanded_gate_counts']['Toffoli']-cm['expanded_gate_counts']['Toffoli']
    assert actual_coarsening == 4*(ft-coarseft)
    return {'model':model['model_id'],'q':q,'probability_bits':bits,'width':w,
            'interval_minus_full_actual':actual,'interval_minus_full_formula':expected,
            'coarsening_actual':actual_coarsening,'coarsening_formula':4*(ft-coarseft)}


def main():
    rows = [test_function_block(w,f,b,dense) for w in (7,12,25,36)
            for f in (1,3,7) for b in (2,4) for dense in (False,True)
            if f+b-2 > 0 and f+b <= w and (1 << b) < (1 << w)]
    models=json.loads((STAGE/'results/generalisation_protocol.json').read_text())['models']
    result={'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), 'protocol_sha256':hashlib.sha256((STAGE/'results/generalisation_protocol.json').read_bytes()).hexdigest(), 'function_blocks_checked':len(rows),'function_blocks':rows,
            'constant_alias_audit':alias_audit(models),
            'complete_graph_cost_checks':[graph_count_audit(m) for m in models if m['assets']==2]}
    target=STAGE/'results/generalisation_compiler_audit.json'
    target.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k != 'function_blocks'},indent=2))


if __name__=='__main__':main()
