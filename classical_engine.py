"""Exact, vectorised classical evaluation of the retained finite pricing law.

All financial arithmetic is integer arithmetic. NumPy is an execution backend,
not a floating-point financial approximation. Timings and operation counts are
backend-specific diagnostics, not a comparison with physical quantum runtime.
"""
from __future__ import annotations

from collections import Counter
from fractions import Fraction
import hashlib
import json
import time

import numpy as np

U64 = np.uint64
M32 = U64((1 << 32) - 1)


def multiply_shift_u64(a, b, drop, width):
    """Return ((a*b) >> drop) modulo 2**width, without product overflow.

    Both inputs are unsigned 64-bit words. Four 32-bit products construct
    the exact 128-bit product. No Python object array is required.
    """
    if not 0 <= drop < 128 or not 1 <= width <= 64:
        raise ValueError('Invalid unsigned product window')
    a, b = np.asarray(a, dtype=U64), np.asarray(b, dtype=U64)
    if drop == 0:
        return (a*b) & U64((1 << width)-1)
    al, ah = a & M32, a >> U64(32)
    bl, bh = b & M32, b >> U64(32)
    ll, lh, hl, hh = al*bl, al*bh, ah*bl, ah*bh
    middle = (ll >> U64(32))+(lh & M32)+(hl & M32)
    low = (ll & M32) | ((middle & M32) << U64(32))
    high = hh+(lh >> U64(32))+(hl >> U64(32))+(middle >> U64(32))
    if drop < 64:
        out = (low >> U64(drop)) | (high << U64(64-drop))
    else:
        out = high >> U64(drop-64)
    return out & U64((1 << width)-1)


def fixed_multiply_u64(a, b, q, width):
    """Specialised exact q-bit truncation for this study's finite registers.

    Register widths are at most q+7 and q is at most 30. The largest result
    before masking has at most q+14 <= 44 bits, although the unshifted product
    can have 74 bits. The q-bit low limb products have at most 60 bits.
    """
    if not 1 <= q <= 30 or not 1 <= width <= q+7:
        raise ValueError('Outside the proved fixed-point product range')
    mask = U64((1 << q)-1)
    ah, al = a >> U64(q), a & mask
    bh, bl = b >> U64(q), b & mask
    out = ((ah*bh) << U64(q))+ah*bl+al*bh+((al*bl) >> U64(q))
    return out & U64((1 << width)-1)


class FiniteEvaluator:
    """Prepared exact interpreter for an encoding_resources.Graph.

    The graph and table contents determine the finite law. This constructor
    splits packed Python integer words once, records their hashes, and releases
    intermediate arrays at their final use. Inputs are one-dimensional equal
    length arrays in their declared unsigned address domains.
    """
    def __init__(self, graph, tables, require_normalised_output=True):
        started = time.perf_counter()
        if not 18 <= graph.q <= 30 or not graph.width <= graph.q+7:
            raise ValueError('Outside the retained study precision family')
        self.graph = graph
        self.require_normalised_output = require_normalised_output
        self.q, self.width = graph.q, graph.width
        self.mask = U64((1 << self.width)-1)
        self.prepared = {}
        self.table_manifest = []
        self.plan = []
        self.node_counts = Counter()
        self.ufunc_counts = Counter()
        self.used = Counter(v for node in graph.nodes.values() for v in node.args)
        self.used[graph.output] += 1
        for name, node in graph.nodes.items():
            self.node_counts[node.kind] += 1
            compiled = {'name':name, 'node':node}
            kind = node.kind
            if kind == 'scale':
                integer,drop,_ = node.meta
                if integer < 0 or integer >= 1 << self.width or drop not in (0,self.q):
                    raise ValueError('Constant multiplier outside the proved register family')
                if integer == 0:mode = 'zero'
                elif integer & (integer-1) == 0:
                    mode = 'power_two';compiled['constant_shift'] = integer.bit_length()-1-drop
                elif int(self.mask)*integer < 1 << 64:mode = 'native'
                elif drop == self.q and integer < 1 << self.q:mode = 'small_q_limb'
                elif drop == self.q:mode = 'full_q_limb'
                else:mode = 'native_modular'
                compiled['constant_product_mode'] = mode
            if kind in ('quantile', 'function'):
                key = ('quantile' if kind == 'quantile' else node.meta, node.bits, self.q)
                raw = tables[key]
                if len(raw) != 1 << node.bits:
                    raise ValueError(('Wrong table length',key))
                if key not in self.prepared:
                    if kind == 'quantile':
                        if not all(0 <= x <= int(self.mask) for x in raw):
                            raise ValueError('Quantile word exceeds its register')
                        data = np.asarray(raw,dtype=U64)
                    else:
                        if not all(0 <= x < 1 << (2*self.width) for x in raw):
                            raise ValueError('Function word exceeds its packed registers')
                        data = (np.fromiter((p & int(self.mask) for p in raw),U64,count=len(raw)),
                                np.fromiter((p >> self.width for p in raw),U64,count=len(raw)))
                    self.prepared[key] = data
                    self.table_manifest.append({'kernel':key[0],'address_bits':key[1],
                        'q':key[2],'entries':len(raw),
                        'sha256':hashlib.sha256(json.dumps(raw,separators=(',',':')).encode()).hexdigest(),
                        'numpy_storage_bytes':int(data.nbytes if kind == 'quantile' else sum(t.nbytes for t in data))})
                compiled['table'] = self.prepared[key]
                if kind == 'function':
                    from encoding_resources import function_span
                    span = function_span(node.meta)
                    if span < 1 or span & (span-1):
                        raise ValueError('Function span must be a positive binary power')
                    f = self.q+span.bit_length()-1-node.bits
                    if not 0 <= f < 64:
                        raise ValueError('Invalid interpolation fractional width')
                    compiled['fractional_bits'] = f
                    compiled['remainder_mask'] = U64((1 << f)-1)
                    max_delta = int(compiled['table'][1].max(initial=0))
                    compiled['native_product_safe'] = max_delta*((1 << f)-1) < 1 << 64
                    compiled['maximum_interpolation_product'] = str(max_delta*((1 << f)-1))
            self.plan.append(compiled)
            self._count_operations(compiled)
        self.preparation_seconds = time.perf_counter()-started

    def _count_operations(self, c):
        kind = c['node'].kind
        counts = Counter()
        if kind in ('add','sub','positive_sub'):counts[kind if kind == 'add' else 'sub'] += 1
        elif kind == 'mul':counts.update({'shift':4,'and':3,'mul':4,'add':3})
        elif kind == 'scale':
            integer,drop,_ = c['node'].meta
            mode = c['constant_product_mode']
            if mode == 'zero':pass
            elif mode == 'power_two':counts['shift'] += 1
            elif mode == 'small_q_limb':counts.update({'shift':2,'and':1,'mul':2,'add':1})
            elif mode == 'full_q_limb':counts.update({'shift':4,'and':3,'mul':4,'add':3})
            else:
                counts['mul'] += 1
                if drop:counts['shift'] += 1
        elif kind == 'shift':counts['shift'] += 1
        elif kind == 'quantile':counts['table_gather'] += 1
        elif kind == 'uniform':counts.update({'mul':2,'add':1})
        elif kind == 'function':
            counts.update({'shift':1,'and':1,'table_gather':2,'add':1})
            if c['native_product_safe']:counts.update({'shift':1,'mul':1})
            else:counts.update({'shift':9,'and':7,'mul':4,'add':5,'or':2})
        elif kind == 'clamp':counts.update({'less_than':1,'where':1,'minimum':1})
        elif kind != 'constant':raise ValueError(kind)
        if kind != 'constant':counts['and'] += 1
        self.ufunc_counts.update(counts)

    def operation_profile(self):
        """Per-sample arithmetic calls plus separately retained setup costs.

        NumPy calls operate on batches; these are scalar-lane operation counts,
        not CPU instructions. Validation, allocation, random generation and
        exact-moment reduction are outside this evaluator operation profile.
        """
        return {'q':self.q,'width':self.width,'inputs':dict(self.graph.inputs),
            'nodes_by_kind':dict(self.node_counts),'scalar_lane_operations_upper':dict(self.ufunc_counts),
            'scalar_lane_operation_upper_total':sum(self.ufunc_counts.values()),
            'table_manifest':self.table_manifest,'prepared_table_bytes':sum(x['numpy_storage_bytes'] for x in self.table_manifest),
            'preparation_seconds_diagnostic':self.preparation_seconds,
            'numpy_version':np.__version__,
            'interpolation_product_branches':[{'node':c['name'],'native_uint64_product_proved_safe':c['native_product_safe'],
                'maximum_interpolation_product':c['maximum_interpolation_product']}
                for c in self.plan if c['node'].kind == 'function'],
            'constant_product_branches':dict(Counter(c['constant_product_mode'] for c in self.plan if c['node'].kind == 'scale')),
            'qualification':'exact finite graph law; conservative scalar-lane upper counts also charge broadcast scalar operations to each lane; not an optimised classical lower bound or physical quantum runtime'}

    def evaluate(self, inputs):
        if set(inputs) != set(self.graph.inputs):
            raise ValueError('Input names do not match the graph')
        vals = {}
        size = None
        for name,bits in self.graph.inputs.items():
            raw = np.asarray(inputs[name])
            if raw.ndim != 1 or raw.dtype.kind not in 'iu' or (raw.dtype.kind == 'i' and np.any(raw < 0)):
                raise ValueError('Input addresses must be unsigned integer vectors')
            array = np.asarray(raw,dtype=U64)
            if np.any(array >= 1 << bits):raise ValueError(('Input out of address range',name))
            if size is None:size = len(array)
            elif len(array) != size:raise ValueError('Input lengths differ')
            vals[name] = array
        if size is None:raise ValueError('At least one sampled input is required')
        remaining = self.used.copy()
        for compiled in self.plan:
            name,node = compiled['name'],compiled['node']
            a = [vals[v] for v in node.args]
            kind = node.kind
            if kind == 'constant':v = U64(node.meta)
            elif kind == 'add':v = a[0]+a[1]
            elif kind in ('sub','positive_sub'):v = a[0]-a[1]
            elif kind == 'mul':v = fixed_multiply_u64(a[0],a[1],self.q,self.width)
            elif kind == 'scale':
                integer,drop,_ = node.meta
                mode = compiled['constant_product_mode']
                if mode == 'zero':v = np.zeros_like(a[0])
                elif mode == 'power_two':
                    shift = compiled['constant_shift']
                    v = a[0] << U64(shift) if shift >= 0 else a[0] >> U64(-shift)
                elif mode == 'small_q_limb':
                    high,low = a[0] >> U64(self.q),a[0] & U64((1 << self.q)-1)
                    v = high*U64(integer)+((low*U64(integer)) >> U64(self.q))
                elif mode == 'full_q_limb':v = fixed_multiply_u64(a[0],U64(integer),self.q,self.width)
                else:
                    v = a[0]*U64(integer)
                    if drop:v = v >> U64(drop)
            elif kind == 'shift':v = a[0] >> U64(node.meta)
            elif kind == 'quantile':v = compiled['table'][a[0]]
            elif kind == 'uniform':v = (U64(2)*a[0]+U64(1))*U64((1 << self.q) >> (node.bits+1))
            elif kind == 'function':
                f = compiled['fractional_bits']
                index = a[0] >> U64(f)
                base,delta = compiled['table']
                if np.any(index >= len(base)):
                    raise ValueError(('Function address out of domain',node.meta,name))
                rem = a[0] & compiled['remainder_mask']
                d = delta[index]
                product = (d*rem) >> U64(f) if compiled['native_product_safe'] else multiply_shift_u64(d,rem,f,self.width)
                v = base[index]+product
            elif kind == 'clamp':
                # A sign bit set represents a negative signed register. Its
                # unsigned value is removed before applying the upper cap.
                v = np.minimum(np.where(a[0] < (1 << (self.width-1)),a[0],U64(0)),U64(node.meta))
            else:raise ValueError(kind)
            vals[name] = v & self.mask
            for argument in node.args:
                remaining[argument] -= 1
                if remaining[argument] == 0:del vals[argument]
        out = vals[self.graph.output]
        if np.ndim(out) == 0:out = np.full(size,out,dtype=U64)
        if self.require_normalised_output and np.any(out > 1 << self.q):raise ValueError('The finite output is outside [0,1]')
        return out

    def sample(self, count, rng):
        """IID finite uniforms under the declared ideal random-bit model.

        NumPy's supplied generator is a reproducible pseudorandom diagnostic
        implementation. Its seed is not mathematical evidence of independence.
        """
        if not isinstance(count,int) or count < 0:raise ValueError('Invalid sample count')
        inputs = {name:rng.integers(0,1 << bits,size=count,dtype=U64) for name,bits in self.graph.inputs.items()}
        return self.evaluate(inputs)


class ExactMoments:
    """Mergeable exact sample count, integer sum and integer sum of squares."""
    def __init__(self,q):
        if not 1 <= q <= 30:raise ValueError('Invalid normalised output precision')
        self.q,self.count,self.total,self.total_squares = q,0,0,0

    def add(self, outputs):
        a = np.asarray(outputs,dtype=U64)
        if a.ndim != 1 or np.any(a > 1 << self.q):raise ValueError('Invalid normalised finite outputs')
        # At most 2**20 outputs per reduction. Low squares, high squares and
        # cross products all have uint64-safe sums. Assembly uses Python ints.
        for start in range(0,len(a),1 << 20):
            v = a[start:start+(1 << 20)]
            low,high = v & U64(65535),v >> U64(16)
            self.count += len(v)
            self.total += int(np.sum(v,dtype=U64))
            self.total_squares += (int(np.sum(low*low,dtype=U64))
                +(int(np.sum(low*high,dtype=U64)) << 17)
                +(int(np.sum(high*high,dtype=U64)) << 32))
        return self

    def merge(self, other):
        if self.q != other.q:raise ValueError('Cannot merge different precisions')
        self.count += other.count
        self.total += other.total
        self.total_squares += other.total_squares
        return self

    def as_dict(self):
        unit = 1 << self.q
        mean = Fraction(self.total,self.count*unit) if self.count else None
        numerator = self.count*self.total_squares-self.total*self.total
        variance = Fraction(numerator,self.count*(self.count-1)*unit*unit) if self.count > 1 else None
        return {'q':self.q,'count':self.count,'integer_sum':str(self.total),
            'integer_sum_of_squares':str(self.total_squares),'sample_variance_integer_numerator':str(numerator),
            'normalised_mean_exact_fraction':str(mean) if mean is not None else None,
            'normalised_sample_variance_exact_fraction':str(variance) if variance is not None else None,
            'normalised_mean_float_diagnostic':float(mean) if mean is not None else None,
            'normalised_sample_variance_float_diagnostic':float(variance) if variance is not None else None}
