"""Recover hard momenta for the reconstructed ee-gamma window in unambiguous shower blocks."""
import argparse
import csv
import gzip
import hashlib
import json
import sys
import tempfile
from collections import defaultdict
from pathlib import Path


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def match_saved(receipt, shower_dir, parents):
    import numpy as np
    import uproot
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from ew_shower import isolate
    from response_lhe import export
    tree = json.loads((shower_dir/receipt.parent.name/"tree.json").read_text())
    rec = json.loads(receipt.read_text())
    parent = parents[rec["integration_id"]].with_name("events.lhe.gz")
    with tempfile.TemporaryDirectory() as name:
        tmp = Path(name)
        digest = isolate(parent, tmp/"source.lhe.gz",rec["source_first_event"],rec["attempted_events"])
        if digest != rec["event_block_sha256"]:
            raise ValueError("Hard event block disagrees with frozen inputs")
        with uproot.open(tree["root"]) as f:
            incoming = np.asarray(f["Data"]["incoming"].array(library="np"))
        if len(incoming) != rec["saved_events"]:
            raise ValueError("Replayed event count differs from frozen input")
        export(tmp/"source.lhe.gz", incoming, tmp/"indices.csv")
        with (tmp/"indices.csv").open() as stream:
            return {int(s["saved_entry"]):int(s["lhe_entry"]) for s in csv.DictReader(stream)}


def region(root, frozen, shower_dir):
    with gzip.open(frozen, "rt") as stream:
        model = json.load(stream)
    records = [s for s in model["shards"] if s["sample"] == "z_ee_a"]
    parents = defaultdict(dict)
    attempted = mapped = omitted_candidates = 0
    event_rates = []
    identities = {}
    hard_paths = {json.loads(p.read_text())["integration_id"]:p for base in (root,root/"ew")
                  for p in base.glob("z_ee_a*/hard.json")}
    for shard in records:
        receipt = root / shard["receipt"]
        packet = receipt.with_name("packet.csv")
        if sha(receipt) != shard["receipt_sha256"] or sha(packet) != shard["packet_sha256"]:
            raise ValueError(f"Hash mismatch in {receipt}")
        r = json.loads(receipt.read_text())
        key = r["integration_id"]
        previous = identities.setdefault(key, shard["source_lhe_sha256"])
        if previous != r["source_lhe_sha256"] or previous != shard["source_lhe_sha256"]:
            raise ValueError("Hard parent identity changed")
        attempted += r["attempted_events"]
        lossless = r["saved_events"] == r["attempted_events"]
        if lossless or shower_dir is not None:
            mapped += r["attempted_events"]
        selected = defaultdict(float)
        with packet.open() as stream:
            for row in csv.DictReader(stream):
                if 123.3 <= float(row["mass_gev"]) < 126.7:
                    if float(row["evweight"]) != 1.:
                        raise ValueError("Expected unweighted positive hard collisions")
                    i = int(row["source_event"])
                    selected[i] += float(row["probability"]) * 1000 * r["cross_section_pb"]
        event_rates.extend(selected.values())
        if lossless or shower_dir is not None:
            if not lossless and selected:
                indices = match_saved(receipt,shower_dir,hard_paths)
                selected = {indices[i]: w for i,w in selected.items()}
            selected = {r["source_first_event"]+i: w for i,w in selected.items()}
            parent = parents[r["integration_id"]]
            if set(parent) & set(selected):
                raise ValueError("Duplicate source index")
            parent.update(selected)
        else:
            omitted_candidates += len(selected)
    s = model["sources"]["z_ee_a"]
    if abs(sum(event_rates)/attempted - s["window_125"]["rate_fb"]) > 1.e-12:
        raise ValueError("Window rate differs from frozen fit")
    events = []
    for key, indices in parents.items():
        if not indices:
            continue
        meta = hard_paths[key]
        if not sha(meta.with_name("events.lhe.gz")) == json.loads(meta.read_text())["lhe_sha256"] == identities[key]:
            raise ValueError("Parent LHE file changed")
        with gzip.open(meta.with_name("events.lhe.gz"), "rt") as stream:
            i = -1
            for line in stream:
                if line.strip() != "<event>":
                    continue
                i += 1
                header = next(stream).split()
                rows = [next(stream).split() for _ in range(int(header[0]))]
                if i not in indices:
                    continue
                incoming = sorted((r for r in rows if r[1] == "-1"), key=lambda r: -int(r[0]))
                outgoing = {int(r[0]):r for r in rows if r[1] == "1"}
                pid = int(incoming[0][0])
                if not (len(incoming) == 2 and 1 <= pid <= 4 and int(incoming[1][0]) == -pid
                        and set(outgoing) == {11,-11,22}):
                    raise ValueError("Unexpected external state")
                ordered = incoming + [outgoing[k] for k in (-11,11,22)]
                momenta = [[float(r[j]) for j in (9,6,7,8)] for r in ordered]
                events.append(dict(integration=key, event=i, hard_flavour=pid, p4=momenta,
                                   scale_gev=float(header[3]), contribution_fb=indices.pop(i)/mapped))
                if not indices:
                    break
        if indices:
            raise ValueError(f"Missing source event in {meta}")
    return dict(frozen_sha256=sha(frozen), attempted=attempted, mapped=mapped,
                unmapped_window_collisions=omitted_candidates, selected_fb=sum(event_rates)/attempted,
                integrations=s["integration_groups"], events=events)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("core_campaign","complement_campaign","core_spectrum","complement_spectrum","output"):
        p.add_argument(name, type=Path)
    p.add_argument("--showers",type=Path,help="Replay tree directory for exact matching of saved entries in partially written blocks")
    p.add_argument("--pdf-weights",type=Path,help="Embed independent full-event EHLQ/NNPDF ratios from pdf_reweight.py")
    a = p.parse_args()
    data = dict(core=region(a.core_campaign, a.core_spectrum,a.showers),
                complement=region(a.complement_campaign, a.complement_spectrum,a.showers))
    if a.pdf_weights is not None:
        pdf = json.loads(a.pdf_weights.read_text())
        for name, region in data.items():
            pairs = {(r["integration"],r["event"]):r for r in pdf[name]}
            for row in region["events"]:
                partner = pairs.pop((row["integration"],row["event"]))
                if abs(row["contribution_fb"] - partner["contribution_fb"]) > 1.e-13:
                    raise ValueError("PDF and matrix-element event weights differ")
                row["pdf_ratio"] = partner["pdf_ratio"]
            if pairs: raise ValueError("Unmapped PDF variation")
            region["pdf_weights_sha256"] = sha(a.pdf_weights)
    payload = (json.dumps(data,separators=(",", ":"),allow_nan=False)+"\n").encode()
    if a.output.suffix == ".gz":
        with a.output.open("wb") as dest, gzip.GzipFile(fileobj=dest,mode="wb",filename="",mtime=0) as out:
            out.write(payload)
    else:
        a.output.write_bytes(payload)
