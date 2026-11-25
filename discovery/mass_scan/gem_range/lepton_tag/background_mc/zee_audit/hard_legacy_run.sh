#!/usr/bin/env bash
set -eu
# Usage: MG5=/path/to/bin/mg5_aMC PYTHON=/path/to/python bash hard_legacy_run.sh OUTPUT
code=$(cd "$(dirname "$0")" && pwd)
out=${1:?Set an output directory}
mkdir -p "$out"
for process in stable_z z_decay; do
    seed=7342851
    if test "$process" = z_decay; then seed=7342852; fi
    "${PYTHON:-python3}" "$code/../ew_generate.py" --mg5 "${MG5:?Set MG5}" --process "$process" \
        --work "$out/$process" --card "$code/hard_$process.dat" --seed "$seed" \
        --events 60000 --cores 8 --set pta=10 --set etaa=-1 --set ptl=0 --set etal=-1 \
        --set dral=0 --set drll=0 --set mmll=0
done
