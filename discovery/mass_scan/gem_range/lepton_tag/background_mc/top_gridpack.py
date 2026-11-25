"""Generate new positive-weight events from a completed MG5 adaptive integration."""
from concurrent.futures import ThreadPoolExecutor
import fcntl
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict

from ew_generate import digest, lhe_info


def packet_events(path):
    with gzip.open(path, "rt") as stream:
        for line in stream:
            if line.strip() != "<event>":
                continue
            head = next(stream).split()
            particles = [next(stream) for _ in range(int(head[0]))]
            trailer = []
            for row in stream:
                trailer.append(row)
                if row.strip() == "</event>":
                    break
            else:
                raise ValueError("incomplete LHE event")
            yield head, particles, trailer


def kinematics(particles, observed):
    rows = [p.split() for p in particles]
    beam = [float(p[9]) for p in rows if p[1] == "-1"]
    photon = [[float(x) for x in p[6:10]] for p in rows if p[0] == "22" and p[1] == "1"]
    observed["sqrt_s_partonic"].append(2 * math.sqrt(beam[0]*beam[1]))
    observed["leading_photon_pt"].append(max(math.hypot(v[0], v[1]) for v in photon))
    if len(photon) == 2:
        px, py, pz, energy = map(sum, zip(*photon))
        observed["diphoton_mass"].append(math.sqrt(max(0., energy**2 - px**2 - py**2 - pz**2)))


def make_gridpack(base):
    target = base / "native-gridpack.tar.gz"
    with (base / "gridpack.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if not target.is_file():
            with tempfile.TemporaryDirectory(dir=base, prefix="gridpack-build-") as stage:
                scratch = Path(stage) / "madgraph"
                subprocess.run(["cp", "-a", "--reflink=auto", base / "madgraph", scratch], check=True)
                with (base / "create-gridpack.log").open("w") as log:
                    subprocess.run([sys.executable, "bin/madevent", "create_gridpack", "hard", "--nb_core=1"],
                                   cwd=scratch, stdout=log, stderr=subprocess.STDOUT, check=True)
                (scratch / "hard_gridpack.tar.gz").replace(target)
    return target


def run_lot(task):
    pack, out, n, seed, scratch_root = task
    out.mkdir(exist_ok=True, parents=True)
    native = out / "events.lhe.gz"
    if not native.is_file():
        with tempfile.TemporaryDirectory(dir=scratch_root or out, prefix=f"native-{seed}-") as stage:
            scratch = Path(stage)
            subprocess.run(["tar", "-xzf", pack, "-C", scratch], check=True)
            # One worker per native grid run. Multiple independent copies can run in parallel.
            with (scratch / "madevent/Cards/me5_configuration.txt").open("a") as card:
                card.write("\nnb_core = 1\nrun_mode = 0\n")
            with (out / "gridrun.log").open("w") as log:
                subprocess.run(["bash", "run.sh", str(n), str(seed)], cwd=scratch,
                               stdout=log, stderr=subprocess.STDOUT, check=True)
            candidate = scratch / "events.lhe.gz"
            info = lhe_info(candidate)
            if info["generated_events"] != n or not info["cross_section_pb"] > 0:
                raise ValueError(f"invalid native grid-run output in {out}: {info}")
            partial = out / "native-events.partial.gz"
            shutil.move(str(candidate), partial)
            partial.replace(native)
    info = lhe_info(native)
    if info["generated_events"] != n or list(map(int, info["lhe_init"][6:8])) != [335900, 335900]:
        raise ValueError(f"wrong event count or PDFs in {native}")
    return dict(lhe=str(native), seed=seed, lhe_sha256=digest(native), **info)


def generate(base, destination, n, seed, cores, scratch_root=None):
    base, destination = base.resolve(), destination.resolve()
    destination.mkdir(parents=True, exist_ok=True)
    source = json.loads((base / "hard.json").read_text())
    pack = make_gridpack(base)
    lots = min(cores, max(1, n // 3000))
    if scratch_root is not None:
        scratch_root.mkdir(mode=0o700, parents=True, exist_ok=True)
    # Preserve every native event weight in its original file. The composite LHE
    # uses a common positive unit of weight for all unweighted source events.
    tasks = [(pack, destination / "gridrun" / f"lot-{i:02}", n // lots + (i < n % lots), seed*100+i, scratch_root)
             for i in range(lots)]
    with ThreadPoolExecutor(max_workers=lots) as pool:
        native = list(pool.map(run_lot, tasks))
    sigma = source["cross_section_pb"]
    known = set()
    reference_kinematics, new_kinematics = defaultdict(list), defaultdict(list)
    for old in sorted(base.parent.glob("round-*/hard.json")):
        for _, particles, _ in packet_events(old.with_name("events.lhe.gz")):
            known.add(hashlib.sha256("".join(particles).encode()).digest())
            if old.parent == base:
                kinematics(particles, reference_kinematics)
    combined = destination / "events.partial.gz"
    with gzip.open(combined, "wt") as dst, gzip.open(base / "events.lhe.gz", "rt") as original:
        for line in original:
            if line.strip() == "<event>":
                break
            if "Number of Events" in line:
                line = re.sub(r"(Number of Events\s*:\s*)\d+", rf"\g<1>{n}", line)
            dst.write(line)
            if line.strip() == "<header>":
                dst.write(f"<!-- Native seeded gridpack continuation, {n} events. See hard.json for individual seeds. -->\n")
        for entry in native:
            reference_weight = None
            for head, particles, trailer in packet_events(Path(entry["lhe"])):
                weight = float(head[2])
                if reference_weight is None:
                    reference_weight = weight
                if weight <= 0. or not math.isclose(weight, reference_weight, rel_tol=1.e-7):
                    raise ValueError(f"Nonuniform weights in native unweighted lot {entry['seed']}")
                identity = hashlib.sha256("".join(particles).encode()).digest()
                if identity in known:
                    raise ValueError("repeated hard-event momenta in native gridpack generation")
                known.add(identity)
                kinematics(particles, new_kinematics)
                head[2] = f"{sigma:.10e}"
                dst.write("<event>\n" + " ".join(head) + "\n" + "".join(particles) + "".join(trailer))
        dst.write("</LesHouchesEvents>\n")
    info = lhe_info(combined)
    if info["generated_events"] != n or not math.isclose(info["cross_section_pb"], sigma, rel_tol=1.e-6):
        raise ValueError("merged LHE does not reproduce the parent normalization")
    avg = sum(l["cross_section_pb"] * l["generated_events"] / n for l in native)
    err = math.sqrt(sum((l["cross_section_error_pb"] * l["generated_events"] / n)**2 for l in native))
    pull = (avg-sigma)/math.hypot(err, source["cross_section_error_pb"])
    if abs(pull) > 5.:
        raise ValueError(f"native grid-run cross section differs from its parent by {pull:.2f} standard deviations")
    from scipy.stats import ks_2samp
    shape = {}
    for label, population in reference_kinematics.items():
        test = ks_2samp(population, new_kinematics[label], method="asymp")
        shape[label] = dict(ks=float(test.statistic), p=float(test.pvalue))
        if test.pvalue < 1.e-9:
            raise ValueError(f"hard-event distribution {label} differs from parent integration: {shape[label]}")
    final = destination / "events.lhe.gz"
    combined.replace(final)
    with gzip.open(final, "rt") as events, (destination / "banner.txt").open("w") as banner:
        for line in events:
            if line.strip() == "<event>":
                break
            banner.write(line)
    info.update(sample=source["sample"], hard_seed=seed, generated_events=n,
                cross_section_pb=sigma, cross_section_error_pb=source["cross_section_error_pb"],
                run_settings={**source["run_settings"], "nevents": n, "iseed": seed},
                integration_id=source["integration_id"], lhe_sha256=digest(final),
                mg5_version=source["mg5_version"], parent_hard_sha256=digest(base / "hard.json"),
                gridpack_sha256=digest(pack), native_weighted_mean_pb=avg, native_mean_error_pb=err,
                native_xsec_pull=pull, native_lots=native, kinematic_closure=shape,
                normalization="common parent integration and its shared cross-section uncertainty")
    for name in ("run_card.dat", "param_card.dat", "proc_card_mg5.dat"):
        # These are the immutable parent integration settings. Native lot
        # seeds and requested counts are recorded explicitly in hard.json.
        shutil.copy2(base / name, destination / name)
    info.update(run_card_sha256=digest(destination/"run_card.dat"),
                proc_card_sha256=digest(destination/"proc_card_mg5.dat"))
    temporary = destination / "hard.json.partial"
    temporary.write_text(json.dumps(info, indent=2) + "\n")
    os.replace(temporary, destination / "hard.json")
