#!/usr/bin/env bash
set -euo pipefail

# Compatibility entry point.  Campaign logic lives in the Python runner so
# the shell and Python interfaces cannot drift in sample definitions, PDF
# setup, response modes, or runtime-module handling.  Existing environment
# controls such as NEVENTS, RUN_TAG, RUN_SAMPLES, MG5_DIR, HERWIG_MODULE, and
# DRY_RUN are read directly by run_gammagamma_campaign.py.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [[ -n "${CAMPAIGN_PYTHON:-}" ]]; then
  PYTHON_EXE="${CAMPAIGN_PYTHON}"
elif command -v python3.11 >/dev/null 2>&1; then
  PYTHON_EXE="$(command -v python3.11)"
else
  PYTHON_EXE="$(command -v python3)"
fi

exec "${PYTHON_EXE}" "${SCRIPT_DIR}/run_gammagamma_campaign.py" "$@"
