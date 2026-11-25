# Diphoton spectra for the GEM search interval

The reconstructed diphoton mass spectra span 70–170 GeV in 1 GeV bins. The 80–160 GeV interval covers the GEM inclusive diphoton hypotheses and the extra 10 GeV on each side supplies data for fitting the continuum near the endpoints. The signal sample has a single generated Higgs mass of 125 GeV. The following files contain Monte Carlo predictions for one SSC experiment at 40 TeV.

| File | Contents |
| --- | --- |
| `rates.csv` | Signal and five background rates in fb in each 1 GeV bin, with the total background and its binwise Monte Carlo standard deviation |
| `pooled.npz` | Per-process rate vectors and full Monte Carlo covariance matrices for the pooled sample |
| `audit.json` | Exposures by process and generator phase space, file hashes, and comparisons with the published 105–160 GeV result |
| `baseline.npz` | Seven entries from the November 18 student release, including two generator domains for the gluon box |
| `fake_box.npz` | Individual results for 193 further dijet packets and 65 further gluon-box packets |
| `nlo_packets.npz` | Two results pooling the 489 further prompt diphoton and photon + jet packets, with the packet identities recorded in metadata |

Each compressed NumPy file can be read with `np.load(path, allow_pickle=False)`. The `sample` array labels the first axis of `rate_fb` and `covariance_fb2`. The `edges_gev` array gives the reconstructed mass boundaries. At integrated luminosity $L$ in inverse femtobarns, the expected bin counts are $L$ times the bin rates. To calculate the MC variance of a sum of bins, sum every element of the corresponding covariance submatrix and multiply by $L^2$. Collision counting fluctuations enter the discovery likelihood separately.

## Reproduce the pooled results

Run the following from the repository root with the discovery Python dependencies. It checks completeness and duplicates against the original discovery archive, validates the rate and covariance for each input over the common 55 bins and compares the combined results by process.

```bash
projection=discovery/mass_scan/gem_range/projection
python3 "$projection/pool.py" --output-dir /tmp/ssc40-gem-pool
```

## Rebuild the three inputs from event records on Timur

The raw projection uses NumPy, SciPy and Uproot. Run these commands from a checkout containing the scripts, with read access to the two frozen student source releases. The commands use a fresh temporary output directory.

```bash
projection=discovery/mass_scan/gem_range/projection
release=/home/shared/HiggsSSC/student-results/ssc40-discovery-20261123
archive="$release/analysis/discovery/data/spectra.tar.gz"
work=$(mktemp -d)
python3 "$projection/baseline.py" \
  --dataset /home/shared/HiggsSSC/diphoton/ssc40-20261118 \
  --tables "$release/analysis/hgammagamma/HOAnalysis/student_backgrounds/ssc40_v1" \
  --archive "$archive" --output "$work/baseline.npz"
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 "$projection/nlo_packets.py" \
  --release "$release" --workers 4 --output "$work/nlo_packets.npz"
python3 "$projection/fake_box.py" \
  --release "$release" --archive "$archive" --output "$work/fake_box.npz"
python3 "$projection/pool.py" --input-dir "$work" --output-dir "$work/pooled"
```

The projection requires two reconstructed photon candidates. All response hypotheses from a single generated collision are grouped before forming bin-to-bin covariances. Every signed NLO event uses its generator normalization. The signal uses the released iHixs production cross section and diphoton branching ratio. Dijet fakes use the nominal independent all-quark fake probabilities and the inclusive native cross section for each generator segment. Documented continuations enter once, with the original collision denominator. For the box contribution the broad sample supplies Born masses outside `[100,150)` GeV and the independent focused samples supply the complementary interval. Either sample can populate any reconstructed mass bin through detector response. The wider projection uses the frozen LHE and the validated list of discarded showers to recover the phase-space assignment for each saved detector event.
