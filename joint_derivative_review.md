# 1. Checked theorem

I checked the proposed joint Hermite–Gaussian derivative certificate. Let ν be a positive finite Borel measure on real log strikes, with total mass H, λ>0, and define f(m)=∫Φ((m−x)/λ)ν(dx). For t>0 set g_m(t)=f(m+log t). The normalised payoff has ν mass one, so the H factor must be removed if f already denotes the normalised probability.

For an integer r≥1, the signed Stirling identity gives D_t^r=t^(−r)∑_{ℓ=1}^r s(r,ℓ)(tD_t)^ℓ. Since tD_t acts as differentiation in m+log t, and Φ^(ℓ)(z)=(-1)^(ℓ−1)He_(ℓ−1)(z)φ(z),

g_m^(r)(t)=t^(−r)∫A_(r,λ)((m+log t−x)/λ)φ((m+log t−x)/λ)ν(dx),

where A_(r,λ)(z)=∑_{ℓ=1}^r s(r,ℓ)λ^(−ℓ)(−1)^(ℓ−1)He_(ℓ−1)(z). Every derivative of the Gaussian distribution function is bounded, so differentiation under the finite measure is justified by dominated convergence. The finite sum does not introduce a convergence issue.

# 2. Global supremum certificate

For fixed r and λ, A is a nonzero polynomial of degree r−1. The leading coefficient is λ^(−r)(−1)^(r−1), since s(r,r)=1 and the probabilists' Hermite polynomial is monic. Therefore Q=A'−zA is a nonzero polynomial of degree r. The function Aφ tends to zero at both infinities and is nonzero somewhere, so its positive absolute maximum is attained at an interior point where Q=0. Roots of A alone can only give zero values. Thus D_(r,λ)=max_{z:Q(z)=0}|A(z)φ(z)|, with all real Q roots included.

For a rigorous numerical implementation I use exact rational polynomial coefficients whenever λ is rational. A complete Sturm isolation certificate must retain all real roots, including repeated roots after square-free reduction. Arb evaluation on each isolating interval then provides a global upper bound, and root refinement can reduce its width. A finite dense-grid maximum alone does not supply this certificate. At r=1 the root is exactly zero and D_(1,λ)=φ(0)/λ, providing a useful independent check.

# 3. Remainder and sharpness

Positivity gives |g_m^(r)(t)|≤H D_(r,λ)t^(−r). For t≥1 this is at most H D_(r,λ), so Taylor's finite Lagrange or integral remainder gives |g_m(T)−∑_{j=0}^k g_m^(j)(1)(T−1)^j/j!|≤H D_(k+1,λ)(T−1)^(k+1)/(k+1)! when T≥1. If T is random, E(T−1)^(k+1)<∞ justifies the expected remainder bound; the frozen lognormal finite mixture has moments of every finite order. No assertion of convergence as k tends to infinity follows.

The H D bound is sharp over the full class of positive strike measures of mass H and over m: choose a Dirac measure at a strike which places the Gaussian argument at a maximising root and take t=1. This does not establish sharpness for the cap's particular uniform measure over monetary strikes. Neither the extremum calculation nor this measure-class sharpness establishes theoretical priority or an improved best-published basket method without a source comparison.
