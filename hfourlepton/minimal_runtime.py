"""Runtime controls and bounded postprocessing for the minimal HO preset.

Hard production groups run sequentially with the full core allocation. Once all
hard files exist, independent sample workers share that allocation. Only the
coordinator writes the main manifest; workers have private journals and logs.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from copy import deepcopy
from dataclasses import asdict, replace
import json
import os
from pathlib import Path
import re
import signal
import subprocess
import sys
import threading
import time

if str(Path(__file__).resolve().parents[1]) not in sys.path:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from hfourlepton import presets


def library_directory(mg5_dir, name):
    config = Path(mg5_dir) / "input/mg5_configuration.txt"
    values = re.findall(rf"(?m)^\s*{name}\s*=\s*([^\n#]*)", config.read_text()) if config.is_file() else []
    value = values[-1].strip() if values else "None"
    if not value or value.lower() == "none":
        return None
    path = Path(value).expanduser()
    return (path if path.is_absolute() else Path(mg5_dir) / path).resolve()


def collier_available(directory):
    # These are the archive and module locations accepted by MG5 3.5.15's
    # loop exporter; accepting an arbitrary .so/include layout fails at export.
    return directory is not None and (directory / "libcollier.a").is_file() and any(
        (parent / "collier.mod").is_file() for parent in (directory / "modules", directory.parent / "include"))


def validate_libraries(mg5_dir):
    ninja, collier = (library_directory(mg5_dir, name) for name in ("ninja", "collier"))
    if not collier_available(collier):
        raise ValueError("minimal requires COLLIER: run setup --workflow fourlepton-ho --preset minimal "
                         "(optionally --collier-dir /path/to/collier for an existing build)")
    if ninja is None or not all((ninja / name).is_file() for name in ("libninja.a", "libavh_olo.a")) or not (ninja.parent / "include/mninja.mod").is_file():
        raise ValueError("minimal requires Ninja/OneLOop: configure --ninja-dir during setup")
    return {"ninja": ninja, "collier": collier}


def production_settings(sample, cfg):
    if cfg.preset != "minimal":
        return None
    from hfourlepton import run_four_lepton_campaign as four
    settings = presets.description(cfg)
    libraries = {}
    for name in ("ninja", "collier"):
        directory = library_directory(cfg.mg5_dir, name)
        libraries[name] = {"directory": str(directory) if directory else None,
                          "libraries": {p.name: four.sha256_file(p) for p in sorted(directory.glob("lib*"))
                                        if p.is_file()} if directory else {}}
    loop = sample.perturbative_order == "LO_loop_induced"
    return dict(req_acc=settings["req_acc"] if sample.matching == "MC@NLO" else None,
                loop_survey_accuracy=settings["loop_survey_accuracy"] if loop else None,
                survey_splitting=cfg.nb_core if loop else None,
                nhel=1 if loop else None, output_dependencies="external", libraries=libraries,
                MLReductionLib="7|6|1" if sample.name.startswith("bkg_gg_") else "6|7|1",
                continuum_m4l_gev=None if sample.is_signal else settings["continuum_m4l_gev"],
                implementation_sha256=four.sha256_file(Path(__file__)))


def continuum_cut(nlo, lower=90., upper=180.):
    """MG5 3.5.15 dummy_cuts hook; QCD real radiation is excluded from m4l."""
    signature = "p,istatus,ipdg" if nlo else "p"
    declarations = ("      integer istatus(nexternal),ipdg(nexternal)\n" if nlo else "")
    selection = ("         if (istatus(i).eq.1.and.p(0,i).gt.0d0.and.\n"
                 "     &       (abs(ipdg(i)).eq.11.or.abs(ipdg(i)).eq.13))\n"
                 "     &       then\n" if nlo else "         if (i.gt.2) then\n")
    return f"""      logical function dummy_cuts({signature})
c     HiggsSSC minimal continuum window; all four charged leptons.
      implicit none
      include 'nexternal.inc'
{declarations}      double precision p(0:{4 if nlo else 3},nexternal),q(0:3),m2
      integer i,j,nlep
      q=0d0
      nlep=0
      do i=1,nexternal
{selection}            nlep=nlep+1
            do j=0,3
               q(j)=q(j)+p(j,i)
            enddo
         endif
      enddo
      m2=q(0)**2-q(1)**2-q(2)**2-q(3)**2
      dummy_cuts=nlep.eq.4.and.m2.ge.{lower**2:.8f}d0
     &           .and.m2.le.{upper**2:.8f}d0
      return
      end
"""


def validate_export_libraries(process_dir, sample):
    """Audit MG5 3.5.15's distinct NLO and loop-induced build layouts.

    The FKS exporter writes its main-program link flags in Source/make_opts
    and its virtual-amplitude flags/includes in SubProcesses/makefile_loop.
    The loop-induced exporter puts both in MadLoop_makefile_definitions.
    Read these files as text; evaluating a makefile could execute shell commands.
    """
    if sample.matching == "MC@NLO":
        definitions = (("Source/make_opts", "libcuttools"),
                       ("SubProcesses/makefile_loop", "LINKLIBS"))
    elif sample.matching == "LO":
        definitions = (("SubProcesses/MadLoop_makefile_definitions", "LINK_LOOP_LIBS"),)
    else:
        raise ValueError(f"no MG5 loop-export layout for {sample.matching}")
    for relative, variable in definitions:
        path = process_dir / relative
        if not path.is_file():
            raise ValueError(f"incomplete {sample.matching} MG5 export: missing {path}; inspect the export log")
        active = "\n".join(line.split("#", 1)[0] for line in path.read_text().splitlines())
        active = active.replace("\\\n", " ")
        flags = re.findall(rf"(?m)^[ \t]*{variable}[ \t]*[:+?]?=[ \t]*([^\n]*)", active)
        if "<Not_found_define_it_yourself>" in active:
            raise ValueError(f"MG5 could not find a loop-library include directory in {path}; inspect the export log")
        if not any(all(re.search(rf"(?:^|\s)-l{name}(?:\s|$)", value)
                       for name in ("collier", "ninja")) for value in flags):
            raise ValueError(f"MG5 did not enable both COLLIER and Ninja in {path} ({variable}); inspect the export log")


def configure_export(process_dir, sample, cfg):
    if cfg.preset != "minimal":
        return
    if cfg.dry_run:
        print(f"+ configure minimal loop libraries and continuum cuts: {process_dir}")
        return
    validate_export_libraries(process_dir, sample)
    card = process_dir / "Cards/MadLoopParams.dat"
    if not card.resolve().is_relative_to(process_dir.resolve()):
        raise ValueError(f"refusing to edit a shared MadLoop card: {card}")
    text = card.read_text()
    pattern = r"(?m)^#MLReductionLib\s*\n(?:[^#]*?)(?=^#|\Z)"
    match = re.search(pattern, text)
    if match is None:
        raise ValueError(f"cannot find MLReductionLib in {card}")
    # Comments explain the fallback order; change only the active parameter.
    lines = match.group().splitlines()
    comments = [line for line in lines[1:] if not line.strip() or line.lstrip().startswith("!")]
    order = "7|6|1" if sample.name.startswith("bkg_gg_") else "6|7|1"
    text = text[:match.start()] + "\n".join(["#MLReductionLib", order, *comments]) + "\n" + text[match.end():]
    card.write_text(text)
    if not sample.is_signal and not cfg.wide_continuum:
        hook = process_dir / "SubProcesses/dummy_fct.f"
        if not hook.resolve().is_relative_to(process_dir.resolve()):
            raise ValueError(f"refusing to edit a shared MG5 source: {hook}")
        original = hook.read_text()
        pattern = r"(?ims)^      logical\s+function\s+dummy_cuts\([^\n]*\).*?^      end\s*$"
        modified, count = re.subn(pattern, lambda _: continuum_cut(sample.matching == "MC@NLO"), original, count=1)
        if count != 1:
            raise ValueError(f"unrecognized MG5 dummy_cuts hook: {hook}")
        backup = hook.with_suffix(".before-minimal.f")
        if backup.exists():
            raise ValueError(f"minimal export already patched: {hook}")
        backup.write_text(original)
        hook.write_text(modified)


def _encode(value):
    if isinstance(value, Path):
        return {"__path__": str(value)}
    raise TypeError(type(value).__name__)


def _decode(value):
    return Path(value["__path__"]) if set(value) == {"__path__"} else value


def run_worker(job):
    from hfourlepton import run_four_lepton_campaign as four
    payload = json.loads(Path(job).read_text(), object_hook=_decode)
    cfg = four.Config(**payload["config"])
    manifest = payload["manifest"]
    samples = four.selected_samples(cfg)
    sample = next(s for s in samples if s.name == payload["sample"])
    entry = next(s for s in manifest["samples"] if s["name"] == sample.name)
    started = time.monotonic()
    try:
        four.run_sample(sample, samples, manifest, cfg)
    except BaseException as exc:
        entry.update(status="failed", error=str(exc))
        raise
    finally:
        entry["last_attempt_wall_seconds"] = time.monotonic() - started
        four.write_manifest(manifest, cfg)


def execute_campaign(samples, manifest, cfg):
    from hfourlepton import run_four_lepton_campaign as four
    groups = {}
    for sample in samples:
        groups.setdefault(sample.production_group_id, sample)
    timings = manifest.setdefault("production_wall_seconds", {})
    for name, sample in groups.items():
        four.log(f"Preparing hard production {name} ({cfg.nb_core} cores)")
        started = time.monotonic()
        four.obtain_production_lhe(sample, samples, cfg)
        timings[name] = time.monotonic() - started
        four.write_manifest(manifest, cfg)

    directory = cfg.run_dir / "workers" / cfg.detector_response
    directory.mkdir(parents=True, exist_ok=True)
    jobs = []
    for sample in samples:
        journal = directory / (sample.name + ".json")
        private = deepcopy(manifest)
        if journal.exists():
            previous = json.loads(journal.read_text())
            if previous.get("configuration_fingerprint") != manifest["configuration_fingerprint"]:
                raise ValueError(f"worker configuration differs: {journal}; use a new tag")
            saved = next(s for s in previous["samples"] if s["name"] == sample.name)
            private["samples"] = [saved if s["name"] == sample.name else s for s in private["samples"]]
        worker_cfg = replace(cfg, manifest_override=journal, postprocess_only=True)
        job = directory / (sample.name + ".job.json")
        job.write_text(json.dumps(dict(config=asdict(worker_cfg), sample=sample.name, manifest=private), default=_encode, indent=2) + "\n")
        jobs.append((sample.name, job, journal))

    active = {}
    lock = threading.Lock()
    stop = threading.Event()

    def launch(item):
        name, job, journal = item
        log = directory / (name + ".log")
        with log.open("a") as stream:
            with lock:
                if stop.is_set():
                    return
                process = subprocess.Popen([cfg.mg5_python or sys.executable, Path(__file__).resolve(),
                    "--worker", job], stdout=stream, stderr=subprocess.STDOUT, start_new_session=True,
                    env={**os.environ, "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1",
                         "PYTHONUNBUFFERED": "1"})
                active[name] = process
            try:
                code = process.wait()
            finally:
                with lock:
                    active.pop(name, None)
        if code:
            raise RuntimeError(f"sample {name} failed (exit {code}); inspect {log}")

    def merge(item):
        name, _, journal = item
        if journal.is_file():
            private = json.loads(journal.read_text())
            entry = next(s for s in private["samples"] if s["name"] == name)
            manifest["samples"] = [entry if s["name"] == name else s for s in manifest["samples"]]
            four.write_manifest(manifest, cfg)

    workers = presets.description(cfg)["postprocess_workers"]
    pool = ThreadPoolExecutor(max_workers=workers)
    futures = {pool.submit(launch, item): item for item in jobs}
    try:
        for future in as_completed(futures):
            future.result()
            merge(futures[future])
            print(f"Completed {futures[future][0]}; logs: {directory}", flush=True)
    except BaseException as error:
        stop.set()
        for future in futures:
            future.cancel()
        # A sample failure should not create partial outputs in healthy samples.
        # An explicit interrupt stops only process groups created by this run.
        if isinstance(error, (KeyboardInterrupt, SystemExit)):
            with lock:
                owned = list(active.values())
            for process in owned:
                try:
                    os.killpg(process.pid, signal.SIGTERM)
                except ProcessLookupError:
                    pass
            for process in owned:
                try:
                    process.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    try:
                        os.killpg(process.pid, signal.SIGKILL)
                    except ProcessLookupError:
                        pass
        else:
            print(f"{error}\nAllowing active sample workers to finish; queued workers are cancelled.", flush=True)
        raise
    finally:
        pool.shutdown(wait=True)
        for item in jobs:
            merge(item)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worker", required=True, type=Path)
    run_worker(parser.parse_args().worker)
