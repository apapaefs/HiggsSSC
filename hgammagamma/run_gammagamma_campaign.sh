#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

NEVENTS="${NEVENTS:-10000}"
EBEAM="${EBEAM:-20000}"
MG5_DIR="${MG5_DIR:-${REPO_ROOT}/MG5_aMC_v3_5_15}"
HERWIG="${HERWIG:-Herwig}"
HERWIG_PDF="${HERWIG_PDF:-NNPDF31_nnlo_as_0118}"
DEFAULT_HERWIG_ENV="${DEFAULT_HERWIG_ENV:-${HOME}/Projects/Herwig/Herwig-REAL-stable-gcc-full/bin/activate}"
HERWIG_ENV="${HERWIG_ENV:-}"
if [[ -z "${HERWIG_ENV}" && -f "${DEFAULT_HERWIG_ENV}" ]]; then
  HERWIG_ENV="${DEFAULT_HERWIG_ENV}"
fi
DEFAULT_COLLIER_DYLIB="${DEFAULT_COLLIER_DYLIB:-${HOME}/Projects/Herwig/Herwig-REAL-stable-gcc-full/opt/OpenLoops-2.1.4/lib/libcollier.2.dylib}"
COLLIER_DYLIB="${COLLIER_DYLIB:-}"
if [[ -z "${COLLIER_DYLIB}" && -f "${DEFAULT_COLLIER_DYLIB}" ]]; then
  COLLIER_DYLIB="${DEFAULT_COLLIER_DYLIB}"
fi
RUN_TAG="${RUN_TAG:-run_01}"
RUN_SAMPLES="${RUN_SAMPLES:-signal_gg_h_aa,bkg_prompt_aa}"
DRY_RUN="${DRY_RUN:-0}"
NB_CORE="${NB_CORE:-1}"
SEED_BASE="${SEED_BASE:-31122002}"

GEN_PHOTON_PT_MIN="${GEN_PHOTON_PT_MIN:-10.0}"
GEN_PHOTON_ETA_MAX="${GEN_PHOTON_ETA_MAX:-6.0}"
GEN_JET_PT_MIN="${GEN_JET_PT_MIN:-10.0}"
GEN_JET_ETA_MAX="${GEN_JET_ETA_MAX:-6.0}"
GEN_LEPTON_PT_MIN="${GEN_LEPTON_PT_MIN:-10.0}"
GEN_LEPTON_ETA_MAX="${GEN_LEPTON_ETA_MAX:-6.0}"
GEN_DRAA_MIN="${GEN_DRAA_MIN:-0.4}"
GEN_DRAJ_MIN="${GEN_DRAJ_MIN:-0.4}"
GEN_DRAL_MIN="${GEN_DRAL_MIN:-0.4}"
GEN_DRJJ_MIN="${GEN_DRJJ_MIN:-0.4}"
GEN_DRLL_MIN="${GEN_DRLL_MIN:-0.0}"

TEMPLATE_IN="${SCRIPT_DIR}/HW-template.in"
ANALYSIS_CODE_DIR="${SCRIPT_DIR}/LOAnalysis/Code"
ANALYSIS_EXE="${ANALYSIS_CODE_DIR}/HwSimPostAnalysis_gammagamma"

# name | category | model | process | madspin_decay | weight_scale
SAMPLES=(
  "signal_gg_h_aa|Signal|loop_sm_haa|g g > h [noborn=QCD]|h > a a|1.0"
  "bkg_prompt_aa|Backgrounds|sm|p p > a a||1.0"
  # Reducible fake-background starting points. Set the weight scale to the
  # desired total fake cross-section factor before enabling them.
  # "bkg_gamma_j|Backgrounds|sm|p p > a j||1.0"
  # "bkg_jj|Backgrounds|sm|p p > j j||1.0"
  # "bkg_ee|Backgrounds|sm|p p > e+ e-||1.0"
)

log() {
  printf '\n==> %s\n' "$*"
}

die() {
  printf 'ERROR: %s\n' "$*" >&2
  exit 1
}

is_dry_run() {
  [[ "${DRY_RUN}" == "1" || "${DRY_RUN}" == "true" || "${DRY_RUN}" == "TRUE" ]]
}

run_cmd() {
  printf '+'
  printf ' %q' "$@"
  printf '\n'
  if ! is_dry_run; then
    "$@"
  fi
}

run_in_dir() {
  local dir="$1"
  shift
  printf '+ cd %q &&' "${dir}"
  printf ' %q' "$@"
  printf '\n'
  if ! is_dry_run; then
    (cd "${dir}" && "$@")
  fi
}

ensure_dir() {
  if is_dry_run; then
    printf '+ mkdir -p %q\n' "$1"
  else
    mkdir -p "$1"
  fi
}

sample_requested() {
  local sample="$1"
  local requested=",${RUN_SAMPLES// /},"
  [[ "${RUN_SAMPLES}" == "all" || "${requested}" == *",${sample},"* ]]
}

require_inputs() {
  [[ -f "${TEMPLATE_IN}" ]] || die "missing Herwig template: ${TEMPLATE_IN}"
  [[ -x "${MG5_DIR}/bin/mg5_aMC" ]] || die "missing MG5 executable: ${MG5_DIR}/bin/mg5_aMC"
  if [[ -n "${HERWIG_ENV}" ]]; then
    [[ -f "${HERWIG_ENV}" ]] || die "missing Herwig activation script: ${HERWIG_ENV}"
  elif ! is_dry_run; then
    command -v "${HERWIG}" >/dev/null 2>&1 || die "Herwig not found on PATH; set HERWIG_ENV or HERWIG"
  fi
  if [[ -n "${COLLIER_DYLIB}" ]]; then
    [[ -f "${COLLIER_DYLIB}" ]] || die "missing COLLIER dynamic library: ${COLLIER_DYLIB}"
  fi
}

run_herwig_in_dir() {
  local dir="$1"
  shift

  if [[ -n "${HERWIG_ENV}" ]]; then
    local command_string
    local quoted_args=""
    local arg
    for arg in "$@"; do
      quoted_args+=" $(printf '%q' "${arg}")"
    done
    command_string="source $(printf '%q' "${HERWIG_ENV}"); $(printf '%q' "${HERWIG}")${quoted_args}"
    run_in_dir "${dir}" bash -lc "${command_string}"
  else
    run_in_dir "${dir}" "${HERWIG}" "$@"
  fi
}

prepare_collier_runtime() {
  local runtime_lib_dir="$1"
  local source_lib_dir

  [[ -n "${COLLIER_DYLIB}" ]] || return
  source_lib_dir="$(cd "$(dirname "${COLLIER_DYLIB}")" && pwd)"
  ensure_dir "${runtime_lib_dir}"

  if is_dry_run; then
    printf '+ cp -pL %q/*.dylib %q/\n' "${source_lib_dir}" "${runtime_lib_dir}"
    printf '+ patch OpenLoops dylib install names in %q\n' "${runtime_lib_dir}"
    printf '+ ensure %q exists\n' "${runtime_lib_dir}/libcollier.dylib"
    return
  fi

  local dylib dest tmp
  for dylib in "${source_lib_dir}"/*.dylib; do
    [[ -e "${dylib}" ]] || continue
    dest="${runtime_lib_dir}/$(basename "${dylib}")"
    tmp="${dest}.tmp.$$"
    cp -pL "${dylib}" "${tmp}"
    mv -f "${tmp}" "${dest}"
  done

  [[ -f "${runtime_lib_dir}/libcollier.dylib" ]] || \
    cp -pL "${COLLIER_DYLIB}" "${runtime_lib_dir}/libcollier.dylib"
  patch_openloops_dylibs "${runtime_lib_dir}"
}

patch_openloops_dylibs() {
  local runtime_lib_dir="$1"
  local dylib dep dep_base

  [[ "$(uname -s)" == "Darwin" ]] || return
  command -v install_name_tool >/dev/null 2>&1 || return
  command -v otool >/dev/null 2>&1 || return

  for dylib in "${runtime_lib_dir}"/*.dylib; do
    [[ -f "${dylib}" ]] || continue
    install_name_tool -id "@rpath/$(basename "${dylib}")" "${dylib}" 2>/dev/null || true
  done

  for dylib in "${runtime_lib_dir}"/*.dylib; do
    [[ -f "${dylib}" ]] || continue
    while IFS= read -r dep; do
      [[ "${dep}" == lib/*.dylib ]] || continue
      dep_base="$(basename "${dep}")"
      [[ -f "${runtime_lib_dir}/${dep_base}" ]] || continue
      install_name_tool -change "${dep}" "@loader_path/${dep_base}" "${dylib}" 2>/dev/null || true
    done < <(otool -L "${dylib}" | awk 'NR > 1 {print $1}')
  done
}

link_madloop_runtime_dirs() {
  local process_dir="$1"
  local runtime_lib_dir="$2"
  local subproc_dir

  [[ -n "${runtime_lib_dir}" ]] || return
  if is_dry_run; then
    printf '+ find %q/SubProcesses -maxdepth 1 -type d -name P\\* -exec ln -sfn %q {}/lib \\;\n' \
      "${process_dir}" "${runtime_lib_dir}"
    return
  fi

  while IFS= read -r subproc_dir; do
    ln -sfn "${runtime_lib_dir}" "${subproc_dir}/lib"
  done < <(find "${process_dir}/SubProcesses" -maxdepth 1 -type d -name 'P*' | sort)
}

restore_mg5_symmetry_factors() {
  local process_dir="$1"
  local subproc_dir

  if is_dry_run; then
    printf '+ restore missing symfact.dat from symfact_orig.dat under %q/SubProcesses/P*\n' "${process_dir}"
    return
  fi

  while IFS= read -r subproc_dir; do
    if [[ ! -f "${subproc_dir}/symfact.dat" && -f "${subproc_dir}/symfact_orig.dat" ]]; then
      cp -p "${subproc_dir}/symfact_orig.dat" "${subproc_dir}/symfact.dat"
    fi
  done < <(find "${process_dir}/SubProcesses" -maxdepth 1 -type d -name 'P*' | sort)
}

inject_mg5_runtime_rpath() {
  local process_dir="$1"
  local runtime_lib_dir="$2"
  local make_opts="${process_dir}/Source/make_opts"

  [[ -n "${runtime_lib_dir}" ]] || return
  if is_dry_run; then
    cat <<EOF
+ patch ${make_opts}
LDFLAGS += -Wl,-rpath,${runtime_lib_dir}
RPATH_LIBS += -Wl,-rpath,${runtime_lib_dir}
EOF
    return
  fi

  [[ -f "${make_opts}" ]] || die "missing MG5 make options: ${make_opts}"
  python3 - "${make_opts}" "${runtime_lib_dir}" <<'PY'
import sys
from pathlib import Path

make_opts = Path(sys.argv[1])
runtime_lib_dir = sys.argv[2]
marker = "# gamma-gamma campaign runtime library path"
block = (
    f"\n{marker}\n"
    f"LDFLAGS += -Wl,-rpath,{runtime_lib_dir}\n"
    f"RPATH_LIBS += -Wl,-rpath,{runtime_lib_dir}\n"
)

text = make_opts.read_text()
if marker not in text:
    make_opts.write_text(text.rstrip() + block)
PY
}

run_generate_events_in_dir() {
  local dir="$1"
  local runtime_lib_dir="$2"
  shift 2

  if [[ -n "${runtime_lib_dir}" ]]; then
    local quoted_args=""
    local arg
    for arg in "$@"; do
      quoted_args+=" $(printf '%q' "${arg}")"
    done
    local command_string
    command_string="export DYLD_LIBRARY_PATH=$(printf '%q' "${runtime_lib_dir}"):\${DYLD_LIBRARY_PATH:-}; "
    command_string+="export DYLD_FALLBACK_LIBRARY_PATH=$(printf '%q' "${runtime_lib_dir}"):\${DYLD_FALLBACK_LIBRARY_PATH:-}; "
    command_string+="./bin/generate_events${quoted_args}"
    run_in_dir "${dir}" bash -lc "${command_string}"
  else
    run_in_dir "${dir}" "./bin/generate_events" "$@"
  fi
}

write_mg5_process_card() {
  local card="$1"
  local model="$2"
  local process="$3"
  local output_dir="$4"

  if is_dry_run; then
    cat <<EOF
+ write MG5 process card ${card}
import model ${model}
define p = g u c d s u~ c~ d~ s~
define j = g u c d s u~ c~ d~ s~
define l+ = e+ mu+
define l- = e- mu-
define vl = ve vm vt
define vl~ = ve~ vm~ vt~
generate ${process}
output ${output_dir} -f
EOF
    return
  fi

  cat > "${card}" <<EOF
import model ${model}
define p = g u c d s u~ c~ d~ s~
define j = g u c d s u~ c~ d~ s~
define l+ = e+ mu+
define l- = e- mu-
define vl = ve vm vt
define vl~ = ve~ vm~ vt~
generate ${process}
output ${output_dir} -f
EOF
}

write_madspin_card() {
  local card="$1"
  local decay="$2"

  if is_dry_run; then
    cat <<EOF
+ write MadSpin card ${card}
set spinmode none
decay ${decay}
launch
EOF
    return
  fi

  cat > "${card}" <<EOF
set spinmode none
decay ${decay}
launch
EOF
}

patch_run_card() {
  local card="$1"
  local seed="$2"

  if is_dry_run; then
    cat <<EOF
+ patch run card ${card}
nevents=${NEVENTS}, ebeam1=${EBEAM}, ebeam2=${EBEAM}, iseed=${seed}
pta=${GEN_PHOTON_PT_MIN}, etaa=${GEN_PHOTON_ETA_MAX}, ptj=${GEN_JET_PT_MIN}, etaj=${GEN_JET_ETA_MAX}
ptl=${GEN_LEPTON_PT_MIN}, etal=${GEN_LEPTON_ETA_MAX}, draa=${GEN_DRAA_MIN}, draj=${GEN_DRAJ_MIN}, dral=${GEN_DRAL_MIN}
drjj=${GEN_DRJJ_MIN}, drll=${GEN_DRLL_MIN}
EOF
    return
  fi

  python3 - "$card" "$NEVENTS" "$EBEAM" "$seed" \
    "$GEN_PHOTON_PT_MIN" "$GEN_PHOTON_ETA_MAX" \
    "$GEN_JET_PT_MIN" "$GEN_JET_ETA_MAX" \
    "$GEN_LEPTON_PT_MIN" "$GEN_LEPTON_ETA_MAX" \
    "$GEN_DRAA_MIN" "$GEN_DRAJ_MIN" "$GEN_DRAL_MIN" "$GEN_DRJJ_MIN" "$GEN_DRLL_MIN" <<'PY'
import re
import sys
from pathlib import Path

card = Path(sys.argv[1])
updates = {
    "nevents": sys.argv[2],
    "ebeam1": sys.argv[3],
    "ebeam2": sys.argv[3],
    "iseed": sys.argv[4],
    "pta": sys.argv[5],
    "etaa": sys.argv[6],
    "ptj": sys.argv[7],
    "etaj": sys.argv[8],
    "ptl": sys.argv[9],
    "etal": sys.argv[10],
    "draa": sys.argv[11],
    "draj": sys.argv[12],
    "dral": sys.argv[13],
    "drjj": sys.argv[14],
    "drll": sys.argv[15],
}

pattern = re.compile(r"^(\s*)([^=]+?)(\s*=\s*)([A-Za-z0-9_]+)(.*)$")
lines = []
for line in card.read_text().splitlines(keepends=True):
    match = pattern.match(line)
    if not match:
        lines.append(line)
        continue
    indent, _value, sep, key, tail = match.groups()
    if key in updates:
        line = f"{indent}{updates[key]:>12}{sep}{key}{tail}\n"
    lines.append(line)

card.write_text("".join(lines))
PY
}

write_herwig_input() {
  local sample="$1"
  local lhe_file="$2"
  local seed="$3"
  local output="$4"

  if is_dry_run; then
    cat <<EOF
+ write Herwig input ${output}
replace LHEFILE.lhe.gz -> ${lhe_file}
set theGenerator:NumberOfEvents ${NEVENTS}
set theGenerator:RandomNumberGenerator:Seed ${seed}
set /Herwig/Partons/thePDFset:PDFName ${HERWIG_PDF}
saverun ${sample} theGenerator
EOF
    return
  fi

  python3 - "$TEMPLATE_IN" "$output" "$sample" "$lhe_file" "$NEVENTS" "$seed" "$HERWIG_PDF" <<'PY'
import re
import sys
from pathlib import Path

template, output, sample, lhe_file, nevents, seed, herwig_pdf = sys.argv[1:]
text = Path(template).read_text()
text = text.replace("LHEFILE.lhe.gz", lhe_file)
text = re.sub(r"set theGenerator:NumberOfEvents\s+\S+",
              f"set theGenerator:NumberOfEvents {nevents}", text)
text = re.sub(r"set theGenerator:RandomNumberGenerator:Seed\s+\S+",
              f"set theGenerator:RandomNumberGenerator:Seed {seed}", text)
text = re.sub(r"set /Herwig/Partons/thePDFset:PDFName\s+\S+",
              f"set /Herwig/Partons/thePDFset:PDFName {herwig_pdf}", text)
text = re.sub(r"set /Herwig/Analysis/HwSim:OutputLocation\s+\S+",
              "set /Herwig/Analysis/HwSim:OutputLocation events/", text)
text = re.sub(r"saverun\s+\S+\s+theGenerator",
              f"saverun {sample} theGenerator", text)
Path(output).write_text(text)
PY
}

find_lhe_file() {
  local process_dir="$1"
  local prefer_decayed="$2"
  local candidate
  local fallback=""

  while IFS= read -r candidate; do
    fallback="$candidate"
    if [[ "${prefer_decayed}" == "1" && "$candidate" == *decayed* ]]; then
      printf '%s\n' "$candidate"
      return 0
    fi
  done < <(find "${process_dir}/Events" -maxdepth 2 -type f \( -name 'unweighted_events.lhe.gz' -o -name 'unweighted_events.lhe' \) | sort)

  [[ -n "$fallback" ]] || return 1
  printf '%s\n' "$fallback"
}

write_root_input_list() {
  local herwig_events_dir="$1"
  local root_input="$2"
  local root_file
  local n_roots=0

  : > "${root_input}"
  while IFS= read -r root_file; do
    printf '%s\n' "${root_file}" >> "${root_input}"
    n_roots=$((n_roots + 1))
  done < <(find "${herwig_events_dir}" -type f -name '*.root' | sort)

  ((n_roots > 0)) || die "no HwSim ROOT files found in ${herwig_events_dir}"
}

run_sample() {
  local sample_index="$1"
  local name="$2"
  local category="$3"
  local model="$4"
  local process="$5"
  local madspin_decay="$6"
  local weight_scale="$7"
  local seed=$((SEED_BASE + sample_index))

  local sample_dir="${SCRIPT_DIR}/LOAnalysis/${category}/events/${name}"
  local card_dir="${sample_dir}/cards"
  local mg5_process_dir="${sample_dir}/mg5_process"
  local herwig_dir="${sample_dir}/herwig"
  local herwig_events_dir="${herwig_dir}/events"
  local runtime_lib_dir=""
  local root_input="${sample_dir}/${name}_hwsim_roots.input"
  local mg5_card="${card_dir}/${name}_proc_card.dat"
  local herwig_input="${herwig_dir}/${name}.in"
  local lhe_file

  log "Running sample ${name}"
  printf 'category=%s model=%s process=%s weight_scale=%s seed=%s\n' \
    "${category}" "${model}" "${process}" "${weight_scale}" "${seed}"

  ensure_dir "${card_dir}"
  ensure_dir "${herwig_events_dir}"
  if [[ -n "${COLLIER_DYLIB}" ]]; then
    runtime_lib_dir="${sample_dir}/lib"
    prepare_collier_runtime "${runtime_lib_dir}"
  fi

  write_mg5_process_card "${mg5_card}" "${model}" "${process}" "${mg5_process_dir}"
  run_cmd "${MG5_DIR}/bin/mg5_aMC" "${mg5_card}"
  inject_mg5_runtime_rpath "${mg5_process_dir}" "${runtime_lib_dir}"
  link_madloop_runtime_dirs "${mg5_process_dir}" "${runtime_lib_dir}"
  restore_mg5_symmetry_factors "${mg5_process_dir}"

  patch_run_card "${mg5_process_dir}/Cards/run_card.dat" "${seed}"
  if [[ -n "${madspin_decay}" ]]; then
    write_madspin_card "${mg5_process_dir}/Cards/madspin_card.dat" "${madspin_decay}"
  fi

  local generate_args=("${RUN_TAG}" "-f" "--nb_core=${NB_CORE}")
  run_generate_events_in_dir "${mg5_process_dir}" "${runtime_lib_dir}" "${generate_args[@]}"

  if is_dry_run; then
    if [[ -n "${madspin_decay}" ]]; then
      lhe_file="${mg5_process_dir}/Events/${RUN_TAG}_decayed_1/unweighted_events.lhe.gz"
    else
      lhe_file="${mg5_process_dir}/Events/${RUN_TAG}/unweighted_events.lhe.gz"
    fi
  else
    lhe_file="$(find_lhe_file "${mg5_process_dir}" "$([[ -n "${madspin_decay}" ]] && printf 1 || printf 0)")" \
      || die "could not locate MG5 LHE output for ${name}"
  fi

  write_herwig_input "${name}" "${lhe_file}" "${seed}" "${herwig_input}"
  run_herwig_in_dir "${herwig_dir}" read "${herwig_input}"
  run_herwig_in_dir "${herwig_dir}" run "${name}.run" "-N${NEVENTS}"

  if is_dry_run; then
    printf '+ find %q -type f -name *.root > %q\n' "${herwig_events_dir}" "${root_input}"
  else
    write_root_input_list "${herwig_events_dir}" "${root_input}"
  fi

  run_cmd "${ANALYSIS_EXE}" "${root_input}" -t "${RUN_TAG}" -w "${weight_scale}"
}

main() {
  require_inputs

  log "Building gamma-gamma post-analysis"
  run_cmd make -C "${ANALYSIS_CODE_DIR}" HwSimPostAnalysis_gammagamma

  local sample_index=0
  local entry name category model process madspin_decay weight_scale
  for entry in "${SAMPLES[@]}"; do
    IFS='|' read -r name category model process madspin_decay weight_scale <<< "${entry}"
    if ! sample_requested "${name}"; then
      continue
    fi
    run_sample "${sample_index}" "${name}" "${category}" "${model}" "${process}" "${madspin_decay}" "${weight_scale}"
    sample_index=$((sample_index + 1))
  done

  ((sample_index > 0)) || die "no samples selected by RUN_SAMPLES=${RUN_SAMPLES}"
}

main "$@"
