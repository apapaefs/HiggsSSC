#!/usr/bin/env bash
# Build the archived PYTHIA 5.6 / JETSET 7.3 control with the SSC detector functions.
set -euo pipefail
src="$(cd "$(dirname "$0")" && pwd)"
project="$(cd "$src/../../../../../.." && pwd)"
work="${1:?Usage: bash build_vintage.sh BUILD_DIRECTORY}"
mkdir -p "$work"
cd "$work"
for f in pythia56.f jetset73.f; do
    if [[ ! -f "$f" ]]; then
        curl -fLsS "https://pythia.org/download/jetset/$f" -o "$f"
    fi
done
python3 - "$src" <<'PY'
import hashlib, pathlib, sys
sys.path.insert(0,sys.argv[1])
from vintage_runs import SOURCES
for name, expected in SOURCES.items():
    if hashlib.sha256(pathlib.Path(name).read_bytes()).hexdigest() != expected:
        raise ValueError("Unexpected generator source: "+name)
PY
fc="${FC:-gfortran}"
cxx="${CXX:-c++}"
for f in pythia56 jetset73; do
    "$fc" -O2 -std=legacy -fallow-argument-mismatch -fno-automatic -c "$f.f" -o "$f.o"
done
"$fc" -O2 -c "$src/vintage.f90" -o main.o
"$cxx" -O2 -std=c++17 -I"$project" $(root-config --cflags) -c "$src/vintage_response.cc" -o response.o
"$cxx" main.o pythia56.o jetset73.o response.o $(root-config --libs) \
    -L"$(dirname "$("$fc" -print-file-name=libgfortran.a)")" -lgfortran -o vintage
