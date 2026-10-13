# Full HO pipeline

Start with the [repository workflow](../../README.md): it installs private
dependencies and runs generation, showering, analysis and reporting. Learn the
LO route first, then select `--workflow gammagamma-ho`.

- [HO process, matching and normalization conventions](README.md).
- [What ihixs calculates](IHIXS_STUDENT_GUIDE.md).
- [Recovery](../../docs/TROUBLESHOOTING.md).
- [Recorded validation](../../VALIDATION.md).

The existing `hgammagamma/run_gammagamma_ho_pipeline.py` entry point remains
supported for existing installations. Keep its original `--work-dir` when
resuming those campaigns; the common command uses `inputs/pipeline/` for new
installations.
