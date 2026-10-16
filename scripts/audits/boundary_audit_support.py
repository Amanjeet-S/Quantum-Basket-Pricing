"""Retained exact integer reference interpreter, reused for boundary audits.

The node semantics are the unchanged checked unsigned and signed integer
compiler. No financial graph evaluate or Circuit.run method is called.
"""
from fractions import Fraction as F
import json

def exact_nearest_even(value):
    value=F(value);lower=value.numerator//value.denominator;rest=value-lower
    if rest>F(1,2) or rest==F(1,2) and lower%2:return lower+1
    return lower


def reference_scale(value,q):
    value=F(value)
    if value.denominator&(value.denominator-1)==0:
        return value.numerator,value.denominator.bit_length()-1
    return exact_nearest_even(value*(1<<q)),q


def reference_table_key(name,node):
    if node.kind=='quantile':return ('quantile',node.bits,str(node.meta))
    return ('scalar',json.dumps(node.meta.as_dict(),sort_keys=True),node.bits)


def interpret(graph,inputs,tables,*,domains=False):
    unit=1<<graph.q;modulus=1<<graph.width;values=dict(inputs);raw={};domain_records=[]
    for name,node in graph.nodes.items():
        arguments=[values[arg] for arg in node.args]
        if node.kind=='constant':value=exact_nearest_even(F(node.meta)*unit)
        elif node.kind=='add':value=arguments[0]+arguments[1]
        elif node.kind in ('sub','signed_sub'):value=arguments[0]-arguments[1]
        elif node.kind=='positive_sub':value=max(0,arguments[0]-arguments[1])
        elif node.kind=='mul':value=arguments[0]*arguments[1]//unit
        elif node.kind=='scale':
            integer,drop=reference_scale(node.meta,graph.q);value=arguments[0]*integer//(1<<drop)
        elif node.kind=='uniform':value=(2*arguments[0]+1)*(unit//(1<<(node.meta+1)))
        elif node.kind=='minimum':value=min(arguments)
        elif node.kind=='quantile':value=tables[reference_table_key(name,node)][arguments[0]]
        elif node.kind=='function':
            span=F(node.meta.span);exponent=span.numerator.bit_length()-span.denominator.bit_length()
            fractional=graph.q+exponent-node.bits
            index,remainder=divmod(arguments[0],1<<fractional)
            table=tables[reference_table_key(name,node)]
            assert 0<=index<len(table),('Function address outside retained table',name,index,len(table))
            packed=table[index];base=packed&(modulus-1);magnitude=(packed>>graph.width)&(modulus-1)
            decreasing=(packed>>(2*graph.width))&1
            adjustment=magnitude*remainder//(1<<fractional)
            value=base-adjustment if decreasing else base+adjustment
            assert value>=0,('Signed interpolant underflow',name,value)
            assert base>=magnitude if decreasing else base+magnitude<modulus
        elif node.kind in ('clamp','clamp_signed'):
            argument=arguments[0]
            if node.kind=='clamp_signed' and argument>=modulus//2:argument-=modulus
            value=max(0,min(argument,exact_nearest_even(F(node.meta)*unit)))
        else:raise ValueError(('Unknown independent node',node.kind))
        raw[name]=value
        if node.kind=='signed_sub':assert -modulus//2<=value<modulus//2,('Signed-domain wrap',name,value)
        else:assert 0<=value<1<<node.output_bits,('Unsigned-domain wrap',name,value,node.output_bits)
        values[name]=value% (1<<node.output_bits)
        if domains:domain_records.append((name,node.kind,value))
    return values[graph.output],values,raw,domain_records


