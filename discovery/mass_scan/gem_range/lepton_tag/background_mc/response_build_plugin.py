"""Add pre-recoil photon-copy momenta to a verified private HwSim ancestry build."""
import argparse
import hashlib
import json
import subprocess
from pathlib import Path

NAME = "HwSimPhotonHistory.so"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def replace_once(src, old, new):
    if src.count(old) != 1:
        raise ValueError(f"Unexpected HwSim source around {old}")
    return src.replace(old, new)


def build(base, out):
    proof = json.loads((base/"build.json").read_text())
    if digest(base/"HwSimAncestry.so") != proof["plugin_sha256"]:
        raise ValueError("Ancestry plugin differs from its build record")
    out.mkdir(parents=True, exist_ok=False)
    for name, sha in proof["source_sha256"].items():
        if digest(base/name) != sha:
            raise ValueError(f"Base source hash mismatch for {name}")
        text = (base/name).read_text()
        if name=="HwSim.h":
            text = replace_once(text,'return "HwSimAncestry.so"',f'return "{NAME}"')
            text = replace_once(text,"std::vector<int> truthPrevious, truthHard;",
                                "std::vector<int> truthPrevious, truthHard, truthPhotonIndex;\n"
                                "  std::vector<double> truthPhotonP4;")
        if name=="HwSim.cc":
            text = replace_once(text,"truthPrevious.clear(); truthHard.clear();",
                                "truthPrevious.clear(); truthHard.clear();\n"
                                "  truthPhotonIndex.clear(); truthPhotonP4.clear();")
            text = replace_once(text,"truthPdg.push_back(p->id());","""truthPdg.push_back(p->id());
    if(p->id()==22) {
      const auto q = p->momentum();
      truthPhotonIndex.push_back(i);
      truthPhotonP4.insert(truthPhotonP4.end(), {q.x()/GeV, q.y()/GeV, q.z()/GeV, q.e()/GeV});
    }""")
            text = replace_once(text,'Data->Branch("truthHard", &truthHard);','''Data->Branch("truthHard", &truthHard);
    Data->Branch("truthPhotonIndex", &truthPhotonIndex);
    Data->Branch("truthPhotonP4", &truthPhotonP4);''')
        (out/name).write_text(text)
    command = [NAME if word=="HwSimAncestry.so" else word for word in proof["command"]]
    with (out/"build.log").open("w") as log:
        subprocess.run(command, cwd=out, stdout=log, stderr=subprocess.STDOUT,check=True)
    receipt = dict(base_plugin_sha256=proof["plugin_sha256"], builder_sha256=digest(Path(__file__)),
                   plugin_sha256=digest(out/NAME), command=command,
                   source_sha256={name:digest(out/name) for name in proof["source_sha256"]})
    (out/"build.json").write_text(json.dumps(receipt,indent=2)+"\n")
    print(out/NAME)


if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for arg in ("base","output"):
        parser.add_argument(arg,type=Path)
    args = parser.parse_args()
    build(args.base.resolve(),args.output.resolve())
