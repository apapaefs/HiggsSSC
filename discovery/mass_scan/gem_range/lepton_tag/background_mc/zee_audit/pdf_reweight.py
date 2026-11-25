"""Replace the two incoming NNPDF densities by EHLQ1 on the frozen LO eeγ events."""
import argparse
from collections import defaultdict
import csv
import ctypes
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np

from vintage_runs import SOURCES

PDF_NAME = "NNPDF40_nlo_as_01180_qed"
WINDOW = (123.3, 126.7)
OBSERVABLES = ["hard_nnpdf_fb", "hard_ehlq_fb", "window_nnpdf_fb", "window_ehlq_fb"]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_events(path, count, beam_energy):
    # MadEvent's unweighted LO writer sets SCALUP to the common factorization scale.
    flavors = np.empty((count, 2), dtype=np.int32)
    energy = np.empty((count, 2))
    scale = np.empty(count)
    i = 0
    with gzip.open(path, "rt") as stream:
        for line in stream:
            if line.strip() != "<event>":
                continue
            header = next(stream).split()
            rows = [next(stream).split() for _ in range(int(header[0]))]
            incoming = sorted((r for r in rows if r[1] == "-1"), key=lambda r: -float(r[8]))
            if len(incoming) != 2 or any(abs(float(r[8])) != float(r[9]) for r in incoming):
                raise ValueError("Expected collinear, effectively massless incoming partons")
            flavors[i] = [int(r[0]) for r in incoming]
            energy[i] = [float(r[9]) for r in incoming]
            scale[i] = float(header[3])
            i += 1
    if i != count or np.any(flavors[:, 0] != -flavors[:, 1]) or not np.all(np.abs(flavors) <= 5):
        raise ValueError("Hard-event count or incoming quark assignment mismatch")
    x = energy / beam_energy
    if not np.all((x > 0.) & (x < 1.)) or not np.all(scale > 0.):
        raise ValueError("Unphysical x or factorization scale")
    return flavors, x, scale


def densities(library, pdf, flavors, x, scale):
    shape = flavors.shape
    ids = np.ascontiguousarray(flavors.ravel(), dtype=np.int32)
    xx = np.ascontiguousarray(x.ravel())
    qq = np.repeat(scale, 2)
    historical = np.empty_like(xx)
    double_array = np.ctypeslib.ndpointer(dtype="float64", ndim=1, flags="C_CONTIGUOUS")
    fn = library.ehlq1
    fn.argtypes = [ctypes.c_int, np.ctypeslib.ndpointer(dtype="int32", ndim=1, flags="C_CONTIGUOUS"),
                   double_array, double_array, double_array]
    fn.restype = None
    fn(len(ids), ids, xx, qq, historical)
    modern = np.fromiter((pdf.xfxQ(int(k), a, b) for k, a, b in zip(ids, xx, qq)),
                         dtype=float, count=len(ids))
    if np.any(modern <= 0.) or np.any(historical < 0.) or not np.all(np.isfinite(modern+historical)):
        raise ValueError("Unexpected nonpositive or nonfinite PDF density")
    return np.prod((historical / modern).reshape(shape), axis=1)


def recovered_indices(selected, flavors, x, beam_energy, receipt, replay_trees):
    import uproot
    ref = json.loads(receipt.read_text())
    start, stop = ref["source_first_event"], ref["source_stop_event"]
    old_tree = receipt.parent / "events/tagged.root"
    if old_tree.is_file():
        if digest(old_tree) != ref["root_sha256"]:
            raise ValueError("Original HwSim checksum mismatch")
        tree = old_tree
    else:
        record = json.loads((replay_trees / receipt.parent.name / "tree.json").read_text())
        tree = Path(record["root"])
    with uproot.open(tree) as file:
        branch = file["Data"]["incoming"].array(library="np")
    if len(branch) != ref["saved_events"]:
        raise ValueError("Mismatch in reconstructed collision count")
    pairs = {}
    for saved in selected:
        incoming = branch[saved]
        order = np.argsort(-incoming[3, :2])
        ids = incoming[4, order].astype(int)
        energies = incoming[0, order]
        match = (flavors[start:stop] == ids) & np.isclose(x[start:stop]*beam_energy, energies,
                                                           rtol=2.e-8, atol=1.e-7)
        matches = np.flatnonzero(match.all(axis=1))
        if len(matches) != 1:
            raise ValueError(f"No unique match for saved event {saved} in {receipt}")
        pairs[saved] = start + int(matches[0])
    if len(set(pairs.values())) != len(pairs):
        raise ValueError("Repeated hard collision")
    return pairs


def outputs(total, cov):
    result = {name:dict(rate_fb=value, se_fb=np.sqrt(cov[i,i])) for i,(name,value) in enumerate(zip(OBSERVABLES,total))}
    for name,i,j in (("hard_ehlq_over_nnpdf",1,0),("window_ehlq_over_nnpdf",3,2),
                     ("nnpdf_window_over_hard",2,0),("ehlq_window_over_hard",3,1)):
        ratio = total[i]/total[j]
        gradient = np.eye(4)[i]/total[j] - ratio*np.eye(4)[j]/total[j]
        result[name] = dict(value=ratio, paired_se=np.sqrt(max(0., gradient @ cov @ gradient)))
    difference = np.array([0., 0., -1., 1.])
    result["ehlq_minus_nnpdf_window_fb"] = dict(value=difference @ total, paired_se=np.sqrt(difference @ cov @ difference))
    acc_ratio = (total[3]/total[1]) / (total[2]/total[0])
    gradient = acc_ratio * np.array([1.,-1.,-1.,1.])/total
    result["ehlq_over_nnpdf_window_acceptance"] = dict(value=acc_ratio, paired_se=np.sqrt(gradient @ cov @ gradient))
    result["ehlq_plus_gem_gg_prescription_window_fb"] = 1.2*total[3]
    result["ehlq_plus_gem_gg_prescription_window_se_fb"] = 1.2*np.sqrt(cov[3,3])
    return result


def region(root, frozen, trees, library, pdf):
    with gzip.open(frozen, "rt") as stream:
        model = json.load(stream)
    shards = defaultdict(list)
    for shard in model["shards"]:
        if shard["sample"] == "z_ee_a":
            shards[shard["integration_id"]].append(shard)
    metadata = {json.loads(p.read_text())["integration_id"]:p for parent in (root, root/"ew")
                for p in parent.glob("z_ee_a*/hard.json")}
    frozen_zee = model["sources"]["z_ee_a"]
    nall = frozen_zee["attempted_events"]
    total, cov = np.zeros(4), np.zeros((4,4))
    nseen = selected_count = recovered = low_x = low_q = selected_low_x = 0
    event_rows = []
    for group, records in shards.items():
        meta_path = metadata[group]
        info = json.loads(meta_path.read_text())
        path = meta_path.with_name("events.lhe.gz")
        if digest(path) != records[0]["source_lhe_sha256"] or info["lhe_init"][6:8] != ["335900","335900"]:
            raise ValueError("Modern LHE hash or PDF ID mismatch")
        n = info["generated_events"]
        beam_energy = float(info["lhe_init"][2])
        flavors, x, scale = read_events(path, n, beam_energy)
        ratio = densities(library, pdf, flavors, x, scale)
        window = np.zeros(n)
        covered = np.zeros(n, dtype=bool)
        for row in records:
            receipt = root / row["receipt"]
            packet = receipt.with_name("packet.csv")
            if digest(receipt) != row["receipt_sha256"] or digest(packet) != row["packet_sha256"]:
                raise ValueError("Frozen packet hash mismatch")
            rec = json.loads(receipt.read_text())
            lo, hi = rec["source_first_event"], rec["source_stop_event"]
            if np.any(covered[lo:hi]):
                raise ValueError("Duplicate generated interval")
            covered[lo:hi] = True
            selected = defaultdict(float)
            with packet.open() as stream:
                for r in csv.DictReader(stream):
                    if WINDOW[0] <= float(r["mass_gev"]) < WINDOW[1]:
                        if float(r["evweight"]) != 1.:
                            raise ValueError("Expected unweighted positive collisions")
                        selected[int(r["source_event"])] += float(r["probability"])
            if not selected:
                continue
            if rec["saved_events"] == rec["attempted_events"]:
                backmap = {i:lo+i for i in selected}
            else:
                backmap = recovered_indices(selected, flavors, x, beam_energy, receipt, trees)
                recovered += len(selected)
            for saved, value in selected.items():
                idx = backmap[saved]
                window[idx] = value
                event_rows.append(dict(integration=group,event=idx,pdf_ratio=float(ratio[idx]),
                                       contribution_fb=value*1000.*info["cross_section_pb"]/nall))
        if not covered.all():
            raise ValueError("Incomplete shower coverage of hard source")
        nseen += n
        selected_count += np.count_nonzero(window)
        low_x += np.count_nonzero(np.min(x, axis=1)<1.e-4)
        selected_low_x += np.count_nonzero((np.min(x,axis=1)<1.e-4) & (window>0.))
        low_q += np.count_nonzero((scale**2<5.) | (scale**2>1.e8))
        observables = np.column_stack([np.ones(n), ratio, window, window*ratio])
        sigma = 1000.*info["cross_section_pb"]
        error = 1000.*info["cross_section_error_pb"]
        # Independent integrations of a single region are combined by their event counts.
        total += sigma*observables.sum(axis=0)/nall
        cov += (sigma/nall)**2*n*np.cov(observables, rowvar=False)
        mean_integral = observables.sum(axis=0)/nall
        cov += error**2 * np.outer(mean_integral,mean_integral)
    if nseen != nall or selected_count != frozen_zee["window_125"]["contributing_events"]:
        raise ValueError("Event counts disagree with frozen measurement")
    np.testing.assert_allclose([total[2], np.sqrt(cov[2,2])],
                               [frozen_zee["window_125"]["rate_fb"], frozen_zee["window_125"]["error_fb"]])
    return dict(spectrum_sha256=digest(frozen), attempted=nseen, window_events=selected_count,
                window_events_in_lossy_blocks_matched_by_incoming_flavour_and_momentum=recovered,
                generated_x_below_ehlq_fit=low_x, selected_x_below_ehlq_fit=selected_low_x,
                generated_q2_outside_ehlq_fit=low_q, **outputs(total,cov)), total, cov, event_rows


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("core_campaign","guard_campaign","core_spectrum","guard_spectrum","replay_trees","library","output"):
        p.add_argument(name, type=Path)
    p.add_argument("--weights", type=Path, help="Optional per-window PDF weights for paired diagnostics")
    a = p.parse_args()
    import lhapdf
    lhapdf.setVerbosity(0)
    pdf = lhapdf.mkPDF(PDF_NAME, 0)
    for source, expected in SOURCES.items():
        if digest(a.library.resolve().with_name(source)) != expected:
            raise ValueError(f"Historical source mismatch for {source}")
    lib = ctypes.CDLL(str(a.library.resolve()))
    report, weights = {}, {}
    totals, covariances = [], []
    for name,root,path in (("core",a.core_campaign,a.core_spectrum),("complement",a.guard_campaign,a.guard_spectrum)):
        report[name], values, covariance, weights[name] = region(root,path,a.replay_trees,lib,pdf)
        totals.append(values)
        covariances.append(covariance)
    cov = sum(covariances)
    report.update(wider=outputs(sum(totals), cov), observable_order=OBSERVABLES,
                  wider_covariance_fb2=cov.tolist(), modern_pdf=PDF_NAME,
                  historical_source_sha256=SOURCES,
                  historical_pdf="EHLQ set 1, internal PYTHIA 5.6 proton PYSTFU at the same factorization scale",
                  weight="Product of the two EHLQ1 incoming x*f values divided by the two NNPDF x*f values",
                  caveat="Changes the luminosity in the Born matrix element at fixed shower, detector response and coupling choices. Quoted errors describe MC and numerical integration.")
    a.output.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    if a.weights:
        a.weights.write_text(json.dumps(weights,separators=(",",":"),allow_nan=False)+"\n")
