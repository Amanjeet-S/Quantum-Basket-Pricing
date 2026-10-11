# Exact classical evaluation of the retained finite law

## 1. Scope and arithmetic contract

I use a prepared NumPy interpreter as an exact classical comparator for the finite financial graphs retained in the architecture study. Every input address is an unsigned integer in its declared domain. Every internal register has width $w\leq q+7$, with $18\leq q\leq30$. Addition and subtraction are reduced modulo $2^w$. Multiplication is unsigned multiplication followed by the graph's specified floor truncation and reduction modulo $2^w$. A signed clamp tests the register's highest bit, returns zero for a negative register, and otherwise caps its unsigned value at the graph's specified upper endpoint. In particular, the graph's operation named `positive_sub` is interpreted as its declared modular subtraction; I do not insert an additional clamp.

I split each packed interpolation word into its base and nonnegative slope once at preparation. The table content and its SHA-256 digest are retained in the operation profile. There is no floating-point exponential, logarithm, normal distribution function or interpolation calculation in the evaluator. The existing table generator remains a separate setup operation with its existing ordinate and quantile certificates. The fact that the classical evaluator and the quantum evaluator use the same finite table does not prove that this is the strongest classical method for the continuous model.

## 2. A product that cannot overflow after the required truncation

For $0\leq a,b<2^w$, write $a=2^q a_H+a_L$ and $b=2^q b_H+b_L$, where $0\leq a_L,b_L<2^q$. Exact distributivity and the fact that the first three terms are integers give

$$\left\lfloor\frac{ab}{2^q}\right\rfloor=2^q a_Hb_H+a_Hb_L+a_Lb_H+\left\lfloor\frac{a_Lb_L}{2^q}\right\rfloor.$$

The high limbs are smaller than $2^7$. The low-limb product is smaller than $2^{2q}\leq2^{60}$, so its multiplication does not overflow an unsigned 64-bit word. Each cross product is smaller than $2^{q+7}\leq2^{37}$. The shifted high product is smaller than $2^{q+14}\leq2^{44}$. More strongly, their sum is exactly $\lfloor ab/2^q\rfloor<2^{2w-q}\leq2^{q+14}\leq2^{44}$, so the additions are safe as well. The unshifted full product may have 74 bits; the implementation never calculates that full product in one unsigned word. Its final bit mask implements reduction modulo $2^w$ exactly.

For a known multiplier $k<2^q$, the formula reduces to $\lfloor ak/2^q\rfloor=a_Hk+\lfloor a_Lk/2^q\rfloor$, because $k_H=0$. If $k$ is a power of two, multiplication and the subsequent binary truncation reduce to a binary shift. A native product is used only when the preparation step has proved $(2^w-1)k<2^{64}$. For zero-drop multiplication, wrap modulo $2^{64}$ followed by the $w$-bit mask is exact even if the native product overflows, since $2^w$ divides $2^{64}$.

## 3. A general 128-bit product window

For other shifts I use 32-bit limbs. Put $B=2^{32}$, $a=a_0+Ba_1$ and $b=b_0+Bb_1$. Let $p_{ij}=a_i b_j$, define

$$m=\lfloor p_{00}/B\rfloor+(p_{01}\bmod B)+(p_{10}\bmod B),$$

and set

$$L=(p_{00}\bmod B)+B(m\bmod B),\qquad U=p_{11}+\lfloor p_{01}/B\rfloor+\lfloor p_{10}/B\rfloor+\lfloor m/B\rfloor.$$

Expanding the right-hand side gives $ab=L+2^{64}U$. Each limb product fits in 64 bits, $m<3B$ fits in 34 bits, and $U<2^{64}$ because $ab<2^{128}$. All partial additions used to form $U$ are nonnegative and at most $U$, so they also fit. For $0<d<64$, the low 64 bits of $\lfloor ab/2^d\rfloor$ are $(L\mathbin{\gg}d)\mathbin{\vert}(U\mathbin{\ll}(64-d))$, with the left shift reduced modulo $2^{64}$. The two bit ranges are disjoint, so bitwise OR equals their exact sum in the retained window. For $64\leq d<128$, the window is $U\mathbin{\gg}(d-64)$. The special case $d=0$ is the low native product word. A final mask selects the requested width. This establishes the helper's contract for every unsigned 64-bit pair, not only the financial values encountered in the numerical audit.

## 4. Interpolation, signed clipping and equality of laws

For an input $x$ and interpolation fractional width $f$, the stored table cell is $j=\lfloor x/2^f\rfloor$ and the remainder is $r=x\bmod2^f$. The evaluator returns $(b_j+\lfloor\delta_j r/2^f\rfloor)\bmod2^w$, which is the scalar graph rule and the emitted binary multiplication rule. A native product is used only if $\max_j\delta_j(2^f-1)<2^{64}$; otherwise the product-window helper in Section 3 is used. The address is checked before either table gather. A quantile node performs a gather of its retained integer word. A midpoint uniform node returns $(2i+1)2^{q-b-1}$, exactly as the graph specifies.

For the signed clamp, an unsigned word $a\geq2^{w-1}$ represents $a-2^w<0$, so the implementation's zero branch equals the signed lower clip. If $a<2^{w-1}$, the signed value is $a$ and its upper clip is simply the unsigned minimum with the cap. Every remaining graph operation agrees with its scalar integer rule. Induction in the graph's topological order therefore gives exactly the same output for every valid finite input tuple. Releasing an array after its final graph use changes only memory retention. It does not change a value consumed by any later node.

If the input addresses are independent discrete uniforms, the classical output has the same finite law as the coherent evaluator's uniformly prepared address registers. Its normalised output is $X=a/2^q$. The flag comparator gives conditional success probability $a/2^q$, so the unconditional quantum flag probability is $\mathbb E X$. An ordinary independent classical sample mean estimates that identical finite quantity without bias. The existing distribution and arithmetic certificates are required to relate it to the continuous financial target. NumPy's supplied seeded generator is a reproducible pseudorandom implementation of the sampling diagnostic. Its reproducibility does not prove mathematical independence or replace the ideal random-bit assumption in the confidence theorem.

## 5. Exact sample moments

The moment accumulator retains $N$, $S=\sum_i a_i$ and $Q=\sum_i a_i^2$ as Python integers. For the nonnegative outputs $a_i\leq2^q\leq2^{30}$, write $a_i=\ell_i+2^{16}h_i$ with $\ell_i<2^{16}$ and $h_i\leq2^{14}$. I reduce at most $2^{20}$ outputs in each NumPy chunk. The sums of low squares, high squares and cross products are respectively smaller than $2^{52}$, at most $2^{48}$ and smaller than $2^{50}$. The direct first-moment sum is at most $2^{50}$. They fit in unsigned 64-bit words, while their shifted assembly is done with Python integers. The identity

$$Q=\sum_i\ell_i^2+2^{17}\sum_i\ell_i h_i+2^{32}\sum_i h_i^2$$

therefore gives the exact squared sum even when $Q>2^{64}$. The exact normalised sample mean is $S/(N2^q)$, and the unbiased sample variance for $N>1$ is $(NQ-S^2)/(N(N-1)2^{2q})$. The output retains these fractions and their integer numerators. Floating renderings of them are labelled diagnostics. Sample variance is not substituted for a certified variance upper bound in a fixed confidence guarantee.

## 6. Cost and audit interpretation

The operation profile records financial node counts, conservative scalar-lane counts for this particular NumPy implementation, prepared table bytes, table hashes and separate preparation timing. Scalar broadcast operations are conservatively charged to every lane. Input validation, array allocation, random generation, moment reduction and table generation must be charged separately in any complete classical workflow. These counts are implementation descriptions. They are neither CPU instruction counts nor lower bounds for classical pricing, and a ratio to logical quantum T gates is not a physical speed-up.

The independent audit evaluates all distinct configurations represented by the retained selected architecture file against the scalar Python-integer graph. It also checks every retained function-table cell boundary and its adjacent register states, full unsigned product-window extremes, overflow-sensitive exact moments and the graph's signed negative branches. The formal equality of interpreters follows from the preceding identities; successful finite samples supplement that proof. I make no novelty claim for elementary limb arithmetic or integer moments.
