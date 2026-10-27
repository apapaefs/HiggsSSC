#!/usr/bin/env python3
"""Pinned VBF_H production, using the diphoton POWHEG stage/seed driver.

Setup and card preparation do not generate events. This CLI does, when run
explicitly by the student workflow in its activated Herwig/LHAPDF environment.
"""
import argparse
import json
import math
from pathlib import Path
import re
import subprocess
import sys
from types import SimpleNamespace

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from hgammagamma import run_gammagamma_ho_pipeline as shared
from hgammagamma import run_gammagamma_ho_campaign as ho
from hgammagamma import run_powheg_hjminnlo as driver
from hgammagamma.run_ihixs_normalization import write_same_or_new
from hfourlepton import mg5_guard

COMMIT = "b9afa1e52f8e90d20d645ea48226588e35642793"
PDF_ID = 335900
HIGGS_WIDTH = .004152
MASS_WINDOW_WIDTHS = 30.  # init_couplings.f in the pinned process
SCALE_GEV = 80.398       # set_fac_ren_scales in Born_phsp.f
LHE_NAME = "powheg-vbf-merged.lhe"
PROVENANCE = ("HiggsSSC POWHEG process VBF_H", "VBF_H commit " + COMMIT,
              "POWHEG core commit " + shared.POWHEG_COMMIT)
PHYSICS = dict(ih1=1, ih2=1, lhans1=PDF_ID, lhans2=PDF_ID,
               alphas_from_pdf=1,
               hwidth=HIGGS_WIDTH, higgsfixedwidth=1, complexpolescheme=0,
               hdecaymode=-1, renscfact=1, facscfact=1, bornonly=0,
               ptsupp=0, withdamp=1, withnegweights=1)


def version_include():
    # The older VBF makefile does not generate the include required by lhefwrite.f.
    return ''.join("      call pwhg_io_write(iun,'" + line.split(' ', 1)[0] + " '//\n"
                   "     & '" + line.split(' ', 1)[1] + "')\n" for line in PROVENANCE)


def build_command(config):
    # This makefile uses F77, not FC. No source edits or weakened physics tests.
    return ["make", "-j1", "pwhg_main", "ANALYSIS=default",
            f"F77={config['fc']} -fno-automatic -Jobj-gfortran -Iobj-gfortran",
            "FFLAGS=-fallow-argument-mismatch", f"CC={config['cc']}",
            f"CXX={config['cxx']}"]


def input_card(ebeam, mass=125., nevents=90):
    settings = dict(numevts=nevents, ebeam1=float(ebeam), ebeam2=float(ebeam), hmass=float(mass), **PHYSICS,
                    ncall1=20000, itmx1=2, ncall2=20000, itmx2=3,
                    foldcsi=1, foldy=1, foldphi=1, nubound=20000, xupbound=2,
                    fastbtlbound=1, storemintupb=1, storeinfo_rwgt=1, testplots=0,
                    rwl_file="'-'", rwl_format_rwgt=1, rwl_group_events=1000)
    text = '# VBF_H NLO, stable H; source defaults: muR=muF=mW=80.398 GeV.\n'
    text += '# Fixed-width Breit-Wigner over mH +/- 30 GammaH; no Born jet cuts.\n'
    text += ''.join(f'{key} {value}\n' for key, value in settings.items())
    text += "<initrwgt>\n<weightgroup name='scale_variation' combine='envelope'>\n"
    for index, (mur, muf) in enumerate(((1, 1), (.5, .5), (.5, 1), (1, .5), (1, 2), (2, 1), (2, 2)), 1001):
        text += f"<weight id='{index}'> renscfact={mur} facscfact={muf} </weight>\n"
    return text + '</weightgroup>\n</initrwgt>\n'


def parameter(prefix, key):
    matches = re.findall(rf'(?mi)^\s*{re.escape(key)}\s+(\S+)', prefix)
    if len(matches) != 1:
        raise ValueError(f'VBF_H requires one explicit POWHEG {key}')
    value = ho.number(matches[0])
    if not math.isfinite(value):
        raise ValueError(f'non-finite POWHEG VBF_H {key}')
    return value


def validate_header(prefix, args, non_nominal_higgs_mass=False):
    for line in PROVENANCE:
        if line not in prefix:
            raise ValueError('VBF LHE lacks pinned POWHEG VBF_H provenance: ' + line)
    expected = dict(PHYSICS, hmass=args.expected_higgs_mass,
                    ebeam1=args.ebeam, ebeam2=args.ebeam)
    for key, value in expected.items():
        if not math.isclose(parameter(prefix, key), value, rel_tol=1.e-10, abs_tol=1.e-12):
            raise ValueError(f'POWHEG VBF_H {key} differs from the campaign value {value}')


def mass_bounds(prefix, cfg):
    validate_header(prefix, SimpleNamespace(expected_higgs_mass=cfg.higgs_mass, ebeam=cfg.ebeam))
    offset = MASS_WINDOW_WIDTHS * HIGGS_WIDTH
    return max(0., cfg.higgs_mass - offset), min(2 * cfg.ebeam, cfg.higgs_mass + offset), HIGGS_WIDTH


def inspect_lhe(path, cfg, count):
    from hfourlepton import run_four_lepton_campaign as four
    info = ho.inspect_lhe(path, SimpleNamespace(signal=True, lhaid=PDF_ID),
        SimpleNamespace(ebeam=cfg.ebeam, expected_higgs_mass=cfg.higgs_mass, nevents=count),
        powheg_validator=validate_header)
    low, high, _ = mass_bounds(four.read_lhe_prefix(path), cfg)
    for is_event, event in four.lhe_blocks(path):
        if not is_event:
            continue
        header, particles = four.lhe_event_record(event)
        if len(particles) < 5:
            raise ValueError('VBF_H requires a complete H+2-parton event')
        scale = ho.number(header[3])
        if not math.isfinite(scale) or scale <= 0:
            raise ValueError('VBF_H requires a positive finite SCALUP')
        higgs = [row for row in particles if row[:2] == ['25', '1']]
        if len(higgs) != 1 or not low - 1.e-6 <= ho.number(higgs[0][10]) <= high + 1.e-6:
            raise ValueError('VBF_H event Higgs mass outside the pinned 30-width window')
    definitions = four.parse_lhe_optional_weight_definitions(path)
    if {d['id'] for d in definitions} != {str(i) for i in range(1001, 1008)}:
        raise ValueError('VBF_H requires the seven declared scale weights (1001..1007)')
    mapping = four.lhe_optional_weight_mapping(path, definitions)
    if mapping['absent_declared_ids']:
        raise ValueError('VBF_H scale weights are absent from the events')
    return info


def source_identity(process):
    def git(*args):
        return subprocess.check_output(['git', '-C', str(process), *args], text=True).strip()
    if git('rev-parse', 'HEAD') != COMMIT:
        raise ValueError('VBF_H source revision differs from the pinned version; run setup')
    if git('diff', 'HEAD', '--stat'):
        raise ValueError('VBF_H tracked sources changed; preserve this installation')
    core = subprocess.check_output(['git', '-C', str(process.parent), 'rev-parse', 'HEAD'], text=True).strip()
    if core != shared.POWHEG_COMMIT:
        raise ValueError('POWHEG core revision differs from the pinned version')
    changed = subprocess.check_output(['git', '-C', str(process.parent), 'diff', 'HEAD',
                                      '--ignore-submodules', '--stat'], text=True).strip()
    if changed:
        raise ValueError('POWHEG core tracked sources changed; preserve this installation')
    if (process / 'svn.version').read_text() != version_include():
        raise ValueError('VBF_H build provenance differs; run setup')
    binary = (process / 'pwhg_main').read_bytes()
    if any(line.encode() not in binary for line in PROVENANCE):
        raise ValueError('VBF_H executable lacks pinned build provenance; rebuild it through setup')
    return dict(process_commit=COMMIT, core_commit=core,
                executable_sha256=shared.digest(process / 'pwhg_main'),
                version_include_sha256=shared.digest(process / 'svn.version'))


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--powheg-dir', required=True, type=Path)
    p.add_argument('--run-dir', required=True, type=Path)
    p.add_argument('--nevents', type=driver.positive_int, required=True)
    p.add_argument('--jobs', type=driver.positive_int, default=4)
    p.add_argument('--ebeam', type=driver.positive_float, required=True)
    p.add_argument('--seed', type=driver.positive_int, required=True)
    args = p.parse_args(argv)
    directory, process = args.run_dir.resolve(), args.powheg_dir.resolve()
    card = input_card(args.ebeam, nevents=args.nevents)
    identity = dict(source=source_identity(process), card=card, nevents=args.nevents,
        jobs=args.jobs, seed=args.seed, driver_sha256=shared.digest(Path(driver.__file__)),
        adapter_sha256=shared.digest(Path(__file__)))
    stamp = directory / 'vbf_production.json'
    if stamp.exists():
        saved = json.loads(stamp.read_text())
        if saved['identity'] != identity:
            raise ValueError('VBF production configuration changed; use a new run tag')
        if saved.get('completed'):
            if shared.digest(directory / LHE_NAME) != saved['lhe_sha256']:
                raise ValueError('completed VBF LHE changed; preserve this run')
            return 0
    elif directory.exists() and any(directory.iterdir()):
        raise ValueError('unmanaged VBF production directory; use a new run tag')
    directory.mkdir(parents=True, exist_ok=True)
    shared.save(stamp, dict(identity=identity, completed=False))
    write_same_or_new(directory / 'vbf.input', card)
    nseeds = min(args.jobs, args.nevents)
    seeds = ''.join(str(1 + (args.seed - 1 + i * 104729) % 900000000) + '\n' for i in range(nseeds))
    write_same_or_new(directory / 'pwgseeds.dat', seeds)
    attempt = 1
    while (directory / f'driver-{attempt}.log').exists():
        attempt += 1
    command = [sys.executable, Path(driver.__file__), '--powheg-dir', process,
        '--run-dir', directory, '--card', directory / 'vbf.input', '--nevents', args.nevents,
        '--jobs', args.jobs, '--ebeam', args.ebeam, '--merged-lhe-name', LHE_NAME,
        '--skip-build', '--resume', '--no-herwig-env', '--no-herwig-module']
    mg5_guard.run_guarded(command, directory / f'driver-{attempt}.log',
                         directory / 'runtime.json', 'VBF_H', program='POWHEG')
    for stage in ('st1-xg1', 'st1-xg2', 'st2', 'st3', 'st4'):
        for seed in range(1, nseeds + 1):
            log = directory / f'run-{stage}-{seed}.log'
            if driver.log_has_failure(log) or not driver.log_is_complete(stage, seed, log, directory):
                raise ValueError(f'POWHEG VBF_H stage did not complete: {log}')
    info = inspect_lhe(directory / LHE_NAME, SimpleNamespace(ebeam=args.ebeam, higgs_mass=125.), args.nevents)
    if info['events'] != args.nevents:
        raise ValueError('VBF_H output count differs from the requested count')
    shared.save(stamp, dict(identity=identity, completed=True, lhe=info,
                           lhe_sha256=shared.digest(directory / LHE_NAME)))
    return 0


if __name__ == '__main__':
    sys.exit(main())
