"""Reject invalid certificate premises, rather than merely replaying formulas."""
import sys
from pathlib import Path
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from encoding_resources import build_graph
from resource_study import propagate


class CertificateContracts(unittest.TestCase):
    def test_reversed_cdf_operands_rejected(self):
        graph=build_graph('full',24,9,0,0,11,12,11,10)
        node=next(n for n in graph.nodes.values() if n.kind=='positive_sub')
        node.args=node.args[::-1]
        with self.assertRaisesRegex(AssertionError,'unproved ordered CDF subtraction'):
            propagate(graph,10)

    def test_unrelated_cdf_arguments_rejected(self):
        graph=build_graph('full',24,9,0,0,11,12,11,10)
        differences=[n for n in graph.nodes.values() if n.kind=='positive_sub']
        first=graph.nodes[differences[0].args[0]]
        log_input=next(name for name,node in graph.nodes.items()
                       if node.kind=='function' and node.meta=='log')
        # A CDF at an unrelated, already available argument retains a valid
        # address but invalidates the asserted ordering relation.
        first.args=(log_input,)
        with self.assertRaisesRegex(AssertionError,'unproved ordered CDF subtraction'):
            propagate(graph,10)

    def test_retained_formula_orders_established(self):
        for cap in (1,10,25):
            graph=build_graph('full',24,9,0,0,11,12,11,cap)
            cert=propagate(graph,cap)
            self.assertEqual(cert['ordered_cdf_subtractions_structurally_verified'],4)


if __name__=='__main__':unittest.main()
