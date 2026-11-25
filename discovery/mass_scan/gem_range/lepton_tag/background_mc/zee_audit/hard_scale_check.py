"""Check that each hard incoming PDF scale equals SCALUP for this unmatched electroweak process."""
import argparse
import gzip
import json
from pathlib import Path
import re


def check(work):
    ref = json.loads((work/"hard.json").read_text())
    n = 0
    worst = 0.
    with gzip.open(work/"events.lhe.gz","rt") as stream:
        for line in stream:
            if line.strip() != "<event>":
                continue
            n += 1
            scale = float(next(stream).split()[3])
            q = {}
            for row in stream:
                m = re.match(r'\s*<pdfrwt beam="([12])">\s*(.*?)\s*</pdfrwt>',row)
                if m:
                    fields = m[2].split()
                    if int(fields[0]) != 1 or len(fields) != 4:
                        raise ValueError("Unexpected matched jet or extra PDF scale")
                    q[m[1]] = float(fields[-1])
                if row.strip() == "</event>": break
            if set(q) != {"1","2"}:
                raise ValueError("Missing per-beam PDF scales")
            worst = max(worst, *(abs(value/scale-1.) for value in q.values()))
    if n != ref["generated_events"] or worst > 2.e-6:
        raise ValueError("SCALUP differs from a factorization scale")
    return dict(events=n,max_relative_difference=worst,hard_lhe_sha256=ref["lhe_sha256"],
                settings=ref["run_settings"])


if __name__ == "__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("run",type=Path);p.add_argument("output",type=Path)
    a=p.parse_args()
    a.output.write_text(json.dumps(check(a.run), indent=2)+"\n")
