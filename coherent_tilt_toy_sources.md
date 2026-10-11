# 1. Portable source and historical evidence

I retain this experiment as evidence against adopting the signed Gaussian-tilt loader as the default signal preparation. It implements a finite classical terminal-model toy, its quantum signal circuit and amplification. It does not implement the financial remainder oracle, the complete pricing estimator or an end-to-end advantage comparison.

The initial evidence is preserved byte for byte in `history/tilt_toy_original_snapshot`, including its source, results, proof fragment, study and original manifest. The portable source revision retains the same mathematical target, gate tape, fixed synthesis seed and pinned generation versions. I regenerated every source-dependent result and certificate rather than changing historical hashes. `coherent_tilt_portability_manifest.json` records the original-to-current source mapping, regenerated result hashes and unchanged existing circuit dependency.

# 2. Replay and certification

Retained-gate replay needs NumPy and the repository's existing `quantum_gates.py`. It runs directly with `python audit_tilt_toy.py`. Independent rotation certification needs python-flint and runs with `python certify_tilt_toy.py`. These operations do not require mpmath, pygridsynth, a vendor directory, or a network connection. The primitive replay is a numerical diagnostic; the independent Acb384 matrix enclosures supply the rotation-error rigour.

For the revised dependency interface I repeated both operations in fresh processes with `QBP_SYNTHESIS_VENDOR` and `PYTHONPATH` unset. The actual signal and amplification counts remain 4,416 and 173,288 T gates, respectively. Their operator-norm upper bounds remain approximately 2.215×10⁻⁸ and 8.637×10⁻⁷. The corresponding directly prepared base-state signal uses 736 T gates. These are signal-only counts and are not full pricing costs.

# 3. Optional regeneration

Regeneration requires the versions in `requirements-tilt-synthesis.txt`: mpmath1.4.1 and pygridsynth1.2.0. I do not install these automatically. `QBP_SYNTHESIS_VENDOR` may point to an existing dependency target. Without that setting, the scripts consider `synthesis_vendor` beside the repository sources if it exists, then ordinary installed packages. No third-party libraries are copied into this package.

The generation sequence is `python coherent_tilt_toy.py`, `python synthesise_tilt_toy.py`, `python certify_tilt_toy.py`, then `python audit_tilt_toy.py`. The four mathematical test groups run with `python checks/test_coherent_tilt_toy.py` and use mpmath for the independently evaluated finite identities. Certification and retained-tape replay can still be performed when the optional generation packages are absent. Generation has NumPy available through the same numerical runtime used for the prototype, while certification uses the existing python-flint runtime.

# 4. Scope of the change

The original mathematical replay did not need a synthesis package. Its script imported the prototype module, which in turn unconditionally imported mpmath through a temporary absolute directory. I removed that accidental execution dependency by separating literal gate replay into `tilt_gate_replay.py` and reading retained model probabilities in the audit. Generation still needs numerical and synthesis libraries; their discovery is now explicit and repository-compatible. The independent certificate selects retained angle keys with exact decimal fractions, so it also avoids that temporary dependency.

The source changes improve reproducibility and portability. They do not alter the financial model, claimed contribution, statistical guarantee, or incomplete end-to-end scope of the original experiment.
