# First pass with inclusive top decays

The first 148,000 generated collisions gave **0.976 ± 0.355 fb** of tagged background in 123.3–126.7 GeV. The three electroweak samples had 16,000 collisions each. The tt̄, tt̄γ and tt̄γγ samples had 30,000, 30,000 and 40,000 respectively. A single tt̄γ event with a production photon and a photon radiated in a top decay contributed 0.333 fb. The result motivates the dedicated radiative top-decay samples in the [full campaign](../README.md).

The observed spectrum in this first simulation gives a median global 5σ exposure of **4.50 fb⁻¹ (5.40 design months)** and a 90% exposure of **7.60 fb⁻¹ (9.12 months)**. The matched GEM-background exposures are 4.23 and 7.08 fb⁻¹. The first-pass background error exceeds the 0.0408 fb target. [Spectrum and independent-event uncertainties](background.json.gz) · [Likelihood outputs](combined-mc.json)

To rerun this checkpoint, execute the following from the repository root.

```bash
code=discovery/mass_scan/gem_range/lepton_tag/background_mc
python3 "$code/refit.py" "$code/first-pass/background.json.gz" --output "$code/first-pass" --pilot --toys 20000
```
