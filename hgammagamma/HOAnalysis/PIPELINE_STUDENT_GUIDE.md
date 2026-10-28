# Full HO pipeline

Start with the [repository workflow](../../README.md): it installs private
dependencies and runs generation, showering, analysis and reporting. Learn the
LO route first, then select `--workflow gammagamma-ho`.

- [HO process, matching and normalization conventions](README.md).
- [What ihixs calculates](IHIXS_STUDENT_GUIDE.md).
- [Recovery](../../docs/TROUBLESHOOTING.md).
- [Recorded validation](../../VALIDATION.md).
- [Cut and XGBoost analyses](../../docs/ANALYSIS.md): their standard HTML reports
  automatically include diphoton mass spectra in 1 GeV bins, with PDF and CSV
  downloads. The cut report shows the spectrum after all other cuts; the
  classifier report also shows the independent test sample after its score cut.

The existing `hgammagamma/run_gammagamma_ho_pipeline.py` entry point remains
supported for existing installations. Keep its original `--work-dir` when
resuming those campaigns; the common command uses `inputs/pipeline/` for new
installations.
