"""Shower and reduce a disjoint event interval from an unweighted LO LHE integration."""
import argparse
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys

from ew_generate import digest


def isolate(parent, dst, start, count):
    # Hash uncompressed event records independently of the header, path and gzip timestamp.
    identity = hashlib.sha256()
    n = 0
    weights = []
    with gzip.open(parent, "rt") as src, dst.open("wb") as raw:
        with gzip.GzipFile(filename="", fileobj=raw, mode="wb", mtime=0) as data:
            in_header = True
            for line in src:
                if line.strip() == "<event>":
                    in_header = False
                    block = [line, next(src)]
                    for part in src:
                        block.append(part)
                        if part.strip() == "</event>":
                            break
                    else:
                        raise ValueError("truncated event")
                    if n >= start:
                        record = "".join(block).encode()
                        data.write(record)
                        identity.update(record)
                        weights.append(float(block[1].split()[2]))
                    n += 1
                    if n == start + count:
                        break
                elif in_header:
                    data.write(line.encode())
            if len(weights) != count:
                raise ValueError("LHE event interval out of bounds")
            data.write(b"</LesHouchesEvents>\n")
    if min(weights) <= 0 or not math.isclose(min(weights), max(weights), rel_tol=1e-7):
        raise ValueError("expected positive unweighted LO hard events")
    return identity.hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("hard","out","plugin","response","repo"):
        p.add_argument("--"+name, type=Path, required=True)
    p.add_argument("--sample", required=True)
    p.add_argument("--events", type=int, required=True)
    p.add_argument("--start", type=int, required=True)
    p.add_argument("--seed", type=int, required=True)
    p.add_argument("--keep-root", action="store_true")
    p.add_argument("--photon-matcher", type=Path)
    p.add_argument("--card", type=Path, help="Herwig input template with @LHE@, @SEED@, @PLUGIN@ and @EVENTS@")
    a=p.parse_args()
    out=a.out.resolve()
    meta=json.loads((a.hard/"hard.json").read_text())
    if meta["sample"] != a.sample or int(meta["run_settings"]["lhaid"]) != 335900:
        raise ValueError("unexpected source sample or PDF")
    if a.start < 0 or a.events < 1 or a.start+a.events > meta["generated_events"]:
        raise ValueError("invalid event block")
    out.mkdir(parents=True, exist_ok=True)
    if (out/"receipt.json").exists():
        raise FileExistsError(out/"receipt.json")
    (out/"events").mkdir(exist_ok=True)
    lhe=out/"source.lhe.gz"
    event_hash=isolate(a.hard/"events.lhe.gz",lhe,a.start,a.events)
    plugin=a.plugin.resolve()
    response=a.response.resolve()
    card=(a.card or Path(__file__).with_name("ew_cards").joinpath("lhe_shower.in")).read_text()
    for key,value in dict(LHE=lhe, SEED=a.seed, PLUGIN=plugin, EVENTS=a.events+2).items():
        card=card.replace(f"@{key}@",str(value))
    (out/"tagged.in").write_text(card)
    for name,command in (("read",["read","tagged.in"]),("run",["run","tagged.run",f"-N{a.events+2}"])):
        with (out/(name+".log")).open("w") as log:
            subprocess.run(["Herwig",*command,"-L",str(plugin.parent)],cwd=out,stdout=log,
                           stderr=subprocess.STDOUT,check=(name=="read"))
    sys.path.insert(0,str(a.repo.resolve()))
    from hgammagamma.ho_shower_completion import _parse_logs
    import uproot
    import numpy as np
    audit=_parse_logs(out,"tagged","source_exhausted")
    roots=list((out/"events").glob("*.root"))
    if len(roots)!=1 or audit["attempted_events"] != a.events:
        raise ValueError(f"incomplete Herwig output {audit}")
    with uproot.open(roots[0]) as data:
        signed=data["Data"]["evweight"].array(library="np")
    if len(signed) != audit["saved_events"] or not np.allclose(signed, 1., rtol=1e-6, atol=0):
        raise ValueError(f"unexpected shower weights {signed[:10]}")
    response_seed=a.seed+99901
    temporary=out/"packet.csv.tmp"
    cmd=[str(response),a.sample,str(roots[0]),str(temporary),str(response_seed)]
    extra={}
    if a.photon_matcher:
        sidecar=out/"hard-photons.csv"
        subprocess.run([sys.executable,str(a.photon_matcher.resolve()),str(lhe),str(roots[0]),str(sidecar)],check=True)
        cmd.append(str(sidecar))
        extra.update(photon_matcher_sha256=digest(a.photon_matcher),sidecar_sha256=digest(sidecar))
    subprocess.run(cmd,check=True)
    with temporary.open() as stream:
        if next(csv.reader(stream),None) is None:
            raise ValueError("missing packet header")
    temporary.replace(out/"packet.csv")
    info={k:meta[k] for k in ("sample","integration_id","cross_section_pb","cross_section_error_pb","hard_seed")}
    info.update(generation_region=meta.get("generation_region","nominal"),seed=a.seed, response_seed=response_seed, source_first_event=a.start,
                source_stop_event=a.start+a.events, source_lhe_sha256=meta["lhe_sha256"],
                lhe_sha256=digest(lhe), event_block_sha256=event_hash,
                signed_attempted=float(a.events),sumw2_attempted=float(a.events),
                signed_saved=float(signed.sum()),sumw2_saved=float((signed**2).sum()),
                attempted_events=a.events,saved_events=len(signed),herwig_completion=audit,
                packet_sha256=digest(out/"packet.csv"),plugin_sha256=digest(plugin),
                response_sha256=digest(response),root_sha256=digest(roots[0]),card_sha256=digest(out/"tagged.in"))
    info.update(extra)
    temp=out/"receipt.json.tmp"
    temp.write_text(json.dumps(info,indent=2)+"\n")
    temp.replace(out/"receipt.json")
    if not a.keep_root:
        roots[0].unlink()
        lhe.unlink()
        (out/"tagged.run").unlink()
    print(json.dumps({key:info[key] for key in ("sample","seed","attempted_events","saved_events")}),flush=True)


if __name__=="__main__":
    main()
