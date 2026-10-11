"""Finite quantum-value and dyadic-bucket contracts, without NumPy."""
from collections import Counter
from fractions import Fraction
from itertools import product
import unittest
import architecture as arch
import encoding_resources as er
import quantum_gates as qg
import value_oracle as vo


def fixture():
    graph=er.Graph(2);graph.width=3
    graph.output=graph.op('quantile',graph.inp('x',2),bits=2)
    return graph,{('quantile',2,2):[0,1,3,4]}


class ValueOracleTests(unittest.TestCase):
    def test_dirty_value_xor_clean_work_and_involution(self):
        graph,tables=fixture()
        for mode in arch.MODES:
            c,inputs,out=vo.compile_value(graph,tables,mode)
            for x,z in product(range(4),range(8)):
                state=c.run({inputs['x']:x,out:z})
                self.assertEqual(c.value(state,out),z^tables[('quantile',2,2)][x])
                self.assertTrue(all(bit==0 for w,bit in enumerate(state) if w not in inputs['x']+out))
            vo.value_into(c,graph,tables,inputs,out)
            for x,z in product(range(4),range(8)):
                state=c.run({inputs['x']:x,out:z})
                self.assertEqual(c.value(state,out),z)
                self.assertEqual(c.value(state,inputs['x']),x)

    def test_exhaustive_valid_bucket_predicates(self):
        for q in range(1,5):
            for b in range(q+1):
                c=vo.ValueCircuit(mode='exact7t');a=c.alloc(q+1)
                u=c.alloc(b+1 if b<q else 0);f=c.alloc(1)[0]
                vo.bucket_predicate_into(c,a,u,b,f)
                for value,U,dirty in product(range((1<<q)+1),range(1<<len(u)),(0,1)):
                    good=value==(1<<q) if b==q else (1<<b)<=value<(1<<(b+1)) and U<value
                    state=c.run({a:value,u:U,(f,):dirty})
                    self.assertEqual(state[f],dirty^bool(good))
                    self.assertTrue(all(bit==0 for w,bit in enumerate(state) if w not in a+u+(f,)))

    def test_exact_Born_count_reconstructs_finite_mean(self):
        graph,tables=fixture();mean=Fraction(0)
        for b in range(graph.q+1):
            c,inputs,u,f=vo.compile_bucket(graph,tables,b,'measured4t')
            good=sum(c.run({inputs['x']:x,u:U})[f] for x,U in product(range(4),range(1<<len(u))))
            probability=Fraction(good,4*(1<<len(u)))
            weight=Fraction(1) if b==graph.q else Fraction(1<<(b+1),1<<graph.q)
            mean+=weight*probability
        self.assertEqual(mean,Fraction(sum(tables[('quantile',2,2)]),4*(1<<graph.q)))

    def test_fused_and_blackbox_have_distinct_complete_costs(self):
        graph,tables=fixture()
        for mode in arch.MODES:
            cache={};value=vo.counted_value_composition(graph,tables,cache,mode)
            for b in range(graph.q+1):
                bucket=vo.counted_bucket_composition(graph,tables,cache,b,mode)
                predicate=bucket['predicate'];black=bucket['blackbox_value_oracle_alternative']
                self.assertEqual(bucket['logical_T_count'],value['logical_T_count']+predicate['logical_T_count'])
                self.assertEqual(black['logical_T_count'],2*value['logical_T_count']+predicate['logical_T_count'])
                expected=Counter({k:2*n for k,n in value['clifford_T_expanded_gate_counts'].items()})
                expected.update(predicate['clifford_T_expanded_gate_counts'])
                self.assertEqual(expected,Counter(black['clifford_T_expanded_gate_counts']))
                self.assertGreaterEqual(black['logical_qubits_serial_upper_bound'],bucket['logical_qubits_serial_upper_bound'])

    def test_literal_exact_CCX_has_correct_complex_phases(self):
        gates=qg.expand_toffoli([('CCX',0,1,2)])
        for x in range(8):
            state=[0j]*8;state[x]=1
            result=qg.replay(gates,state)
            expected=x^(4 if x&3==3 else 0)
            for j,z in enumerate(result):self.assertAlmostEqual(abs(z-(1 if j==expected else 0)),0,places=14)

    def test_retained_emitted_blocks_reproduce_fresh_flags(self):
        graph,tables=fixture()
        for mode in arch.MODES:
            old=arch.counted_composition(graph,tables,{},mode)
            body=vo.prepare_existing_body(old,graph.q)
            for b in range(graph.q+1):
                fast=vo.profile_from_prepared_body(body,b)
                fresh=vo.counted_bucket_composition(graph,tables,{},b,mode)
                self.assertEqual(fast['expanded_gate_counts'],fresh['expanded_gate_counts'])
                self.assertEqual(fast['grover_iterate_T_count'],fresh['grover_iterate_T_count'])
            for bits in range(graph.q+1):
                fast=vo.profile_from_prepared_range(body,bits)
                circuit,*_=vo.compile_range(graph,tables,bits,mode)
                self.assertEqual(fast['expanded_gate_counts'],circuit.resources()['primitive_gate_counts'])

    def test_clipped_flag_preserves_full_value_and_endpoint(self):
        graph,tables=fixture()
        for mode in arch.MODES:
            old=arch.counted_composition(graph,tables,{},mode)
            body=vo.prepare_existing_body(old,graph.q)
            for bits in range(graph.q+1):
                circuit,inputs,u,flag=vo.compile_clipped_range(graph,tables,bits,mode)
                fast=vo.profile_from_prepared_clipped_range(body,bits)
                self.assertEqual(fast['expanded_gate_counts'],old['expanded_gate_counts'])
                self.assertEqual(fast['expanded_gate_counts'],circuit.resources()['primitive_gate_counts'])
                successes=0
                for x,U in product(range(4),range(1<<bits)):
                    state=circuit.run({inputs['x']:x,u:U})
                    self.assertEqual(state[flag],int(U<tables[('quantile',2,2)][x]))
                    self.assertTrue(all(bit==0 for w,bit in enumerate(state) if w not in inputs['x']+u+(flag,)))
                    successes+=state[flag]
                self.assertEqual(successes,sum(min(a,1<<bits) for a in tables[('quantile',2,2)]))


if __name__=='__main__':unittest.main()
