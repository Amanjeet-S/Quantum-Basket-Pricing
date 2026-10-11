# 1. Research question and access

I compare the proposed quantum Taylor-remainder importance state with published coherent reweighting, function loading and Gaussian-mixture preparation. The financial dynamics remain the established multivariate geometric Brownian motion model. The prospective new object is the pricing estimator and its preparation strategy. This is a source comparison, not a claim of theoretical priority.

Checked on 7 October 2026. The paper-index search failed with a DNS error. I used primary publisher and author-repository retrieval. An attempted Safari navigation did not complete, so I continued with retrieval. No purchase, subscription, data upgrade or message to a third party occurred.

# 2. Closest checked sources

## 2.1 Coherent rejection sampling

Ozols, M., Roetteler, M. and Roland, J. (2013). Quantum rejection sampling. *ACM Transactions on Computation Theory*, 5(3), article 11, 1–33. [DOI](https://doi.org/10.1145/2493252.2493256). The journal publication was verified through institutional metadata. I read the 19-page author-hosted conference/manuscript version at Microsoft, Section 4.1, equations (11) to (15) and the subsequent amplification discussion. I did not read the 33-page final journal body. [Consulted version](https://www.microsoft.com/en-us/research/wp-content/uploads/2016/02/Ozols20et20al201103.2774.pdf).

The method adds a controlled acceptance rotation to a prepared coherent sample and boosts its success by amplitude amplification. The post-selection norm determines the preparation overhead. This directly precedes converting a positive Gaussian proposal into a remainder-weighted proposal. The target density, Gaussian moment normaliser and financial remainder can distinguish the application, but coherent rejection and its quadratic improvement cannot be presented as new. The source assumes preparation-oracle access and does not pay for this project's Gaussian loader or reversible financial arithmetic.

## 2.2 Finite-series LCU preparation

Rosenkranz, M., Brunner, E., Marin-Sanchez, G., Fitzpatrick, N., Dilkes, S., Tang, Y., Kikuchi, Y. and Benedetti, M. (2025). Quantum state preparation for multivariate functions. *Quantum*, 9, 1703. [DOI](https://doi.org/10.22331/q-2025-04-09-1703). I read the final publisher PDF, Section 2.5, equations (25) to (32), Section 2.6, the Gaussian example in Section 4 and Appendix D. [Consulted final body](https://quantum-journal.org/papers/q-2025-04-09-1703/pdf/).

The paper prepares finite Fourier/Chebyshev representations by LCU, retaining complex coefficient phases and their one-norm normalisation. Its success probability is the squared target norm divided by the squared coefficient one-norm and the uniform-register size. Appendix D explicitly includes inverse-square-root success overhead. Gaussian state preparation and LCU coefficients are thus existing mechanisms. The prospective distinction is a closed Gaussian-tilt expansion tied to a bounded pricing remainder rather than a generic numerical function expansion. Their principal resource tables count CX gates and omit single-qubit gates, so they are not already fault-tolerant T budgets for this project.

## 2.3 Gaussian mixtures and financial applications

Zhuang, X.-N., Chen, Z.-Y., Xue, C., Xu, X.-F., Wang, C., Liu, H.-Y., Tan, M.-Y., Sun, T.-P., Wang, Y.-J., Zhang, J.-X., Wu, Y.-C. and Guo, G.-P. (2026). Statistics-informed parameterized quantum circuit: towards practical quantum state preparation and learning via maximum entropy principle. *npj Quantum Information*, 12, 45. [DOI](https://doi.org/10.1038/s41534-026-01191-5). I read the final typeset 13-page journal body at INSPIRE, pp. 6–8, Theorem 2, equations (16) to (18) and (21) to (25). Publication: 9 February 2026; version of record: 18 March 2026. [Consulted final body](https://inspirehep.net/files/2161536ac38147e4fec94d1be258cd2e).

Theorem 2 prepares a purified mixture retaining the parameter register. Its preparation cost includes a filling-rate factor. The paper demonstrates Gaussian mixtures and explicitly couples preparation to payoff evaluation and amplitude estimation for finance. Consequently, a purified Gaussian mixture or finance use alone is not new. A defensible difference would be an analytically chosen remainder proposal with a proved monetary error certificate and charged implementation. I did not inspect its supplementary theorem proofs, so the stated complexity has not been independently reproduced here.

## 2.4 Polynomial transformations and their filling ratio

Gonzalez-Conde, J., Watts, T. W., Rodriguez-Grasa, P. and Sanz, M. (2024). Efficient quantum amplitude encoding of polynomial functions. *Quantum*, 8, 1297. [DOI](https://doi.org/10.22331/q-2024-03-21-1297). I read the final publisher PDF, Section 2.1 on preprocessing, Sections 3.4.1–3.4.3, Lemma 2, equation (19) and the amplification discussion, and Section 3.5. [Consulted final body](https://quantum-journal.org/papers/q-2024-03-21-1297/pdf/).

QSVT transforms loaded amplitudes into a polynomial, followed by amplitude amplification. The filling ratio can make loading expensive; it is explicitly part of the resource analysis. Polynomial reweighting and charging its success overhead are established. The present proposal acts on a nonlinear exponential basket ratio and has an analytic finite Gaussian expansion, which is a different representation. Replacing one representation by the other needs a matched preparation-cost comparison, not an assertion that polynomial quantum loading is original.

## 2.5 LCU with a Taylor series

Berry, D. W., Childs, A. M., Cleve, R., Kothari, R. and Somma, R. D. (2015). Simulating Hamiltonian dynamics with a truncated Taylor series. *Physical Review Letters*, 114, 090502. [DOI](https://doi.org/10.1103/PhysRevLett.114.090502). I read the final publisher-version PDF at Macquarie University, the method summary, equations (3) to (6), and the LCU/amplification construction. [Consulted final body](https://research-management.mq.edu.au/ws/portalfiles/portal/62038617/Publisher%20version%20(open%20access).pdf).

The combination of truncated Taylor expansion and LCU is established. Its robust oblivious amplification concerns a truncated evolution operator close to unitary. The pricing filter is non-unitary and generally not close to unitary, so that near-deterministic operator conclusion does not transfer. Ordinary state-specific amplification is available, with preparation, adjoints and reflections charged. The proposed financial Taylor control and importance normalisation differ from Hamiltonian simulation, but that difference does not confer novelty on LCU.

## 2.6 Adaptive quantum importance sampling

Pyretzidis, K., Martínez de Lejarza, J. J. and Rodrigo, G. (2026). Unlocking multidimensional integration with quantum adaptive importance sampling. *Communications Physics*, 9, 269. [DOI](https://doi.org/10.1038/s42005-026-02684-7). Publisher metadata records publication on 22 May 2026 and the version of record on 14 August 2026. I read publisher-retrieved passages from Results, Algorithm summary, and the description of the finite-shot tiling correction. Direct full-page/PDF retrieval was blocked by the Nature authentication redirect, and an attempted Safari navigation did not complete. I did not inspect the complete final body or its proofs.

The checked passages describe trained PQCs used as proposals, with computational-basis measurements and classical statistical postprocessing. That differs from coherent amplitude estimation of a bounded analytically normalised remainder. It nevertheless establishes a published quantum importance-sampling antecedent. This limited inspection cannot exclude a closer result elsewhere in the article or literature.

# 3. Relevant supplementary gaps

McArdle, Gilyén and Berta's *Quantum state preparation without coherent arithmetic* now has a peer-reviewed successor: *Physical Review Letters*, 136(24), 240603 (2026), DOI [10.1103/ntvs-c48s](https://doi.org/10.1103/ntvs-c48s). I verified its publication through the authors' institutional record and PubMed. I have not read the final body during this task. It should replace a preprint-only reference if the Gaussian/QSVT loader is used.

Bauer, Deliyannis, Freytsis and Nachman's *Practical considerations for the preparation of multivariate Gaussian states on quantum computers*, arXiv:2109.10918, remains a preprint in its checked arXiv record. A search result associated DOI 10.1103/PhysRevResearch.6.023191 with it, but that DOI belongs to *Continuous optimization by quantum adaptive distribution search*. I rejected that false association. The preprint can supplement implementation comparisons, but I would not label it a journal publication.

Deiml and Peterseim's *Quantum sampling and moment estimation for transformed Gaussian random fields*, arXiv:2508.13879, is listed as a preprint by the institutional repository. It addresses generated random fields and bounded pointwise transformations. I retrieved its body but found no importance-sampling discussion in the inspected text. It is a wider related source, not evidence that the specific remainder proposal has been anticipated or that it is original.

# 4. Contribution boundary

I would call the candidate a quantum remainder-pricing construction whose proposal is determined by a proved analytic envelope. Its distinctive mathematical object is the coupling of an odd-order financial Taylor control, an even-power arithmetic-geometric ratio tilt, a closed Gaussian normaliser and a bounded correction observable. The generic loading, rejection, amplification and amplitude-estimation steps remain attributed to published methods. Neither a new acronym nor a new order of those steps establishes theoretical priority.

The decisive next test is a finite-register implementation with retained coefficients, probability normalisation, precision errors and resource counts, compared with both direct estimation and an equally informed classical sampler. A novel pricing construction can still be useful without a physical quantum advantage, but a continuous wavefunction calculation alone cannot show that a runnable quantum pricing method exists. The companion `positive_envelope_comparator.md` proves an immediate alternative to the signed-Gaussian LCU proposal and explains the relevant implementation charges.
