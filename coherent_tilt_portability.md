# 1 Retained evidence and dependency scope

I preserve the original toy sources, tapes and certificates byte for byte under history/tilt_toy_original_snapshot. The live files use explicit dependency discovery. The portability manifest records each original and regenerated source hash and verifies that the retained logical and Clifford and T tapes, rotation tapes and norm bounds have the same content.

The compiled-tape replay requires NumPy and the existing quantum_gates.py. The independent rotation certificate requires the existing Python-Flint runtime. Neither requires the optional synthesis vendor, mpmath or pygridsynth. The symbolic toy generator and its associated unit tests require the pinned optional mpmath version. Regenerating proposed rotation sequences also requires the pinned synthesis package. requirements-tilt-synthesis.txt records these optional versions. The runtime discovers an explicitly configured QBP_SYNTHESIS_VENDOR, a repository-local synthesis_vendor if present, or normal installed packages. No dependency library is copied into the repository.

# 2 Reproduction

The retained replay and independent certificate can be run in their respective existing runtimes:

```sh
python audit_tilt_toy.py
python certify_tilt_toy.py
```

With the pinned optional synthesis dependencies already available, the generator and its four test groups can be run:

```sh
python coherent_tilt_toy.py
python synthesise_tilt_toy.py
python -m unittest discover -s checks -p test_coherent_tilt_toy.py -v
```

# 3 Interpretation

This is a small, finite, forward-normalised classical toy model and a genuine compiled quantum signal. The adverse comparison concerns signal preparation weighted by the correction normalisation. The full financial observable and price estimator remain unimplemented, so the counts are not a complete pricing budget. The source and published-mechanism scope is recorded in coherent_tilt_toy_sources.md and quantum_method_sources.md.
