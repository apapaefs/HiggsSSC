#!/usr/bin/env python3
"""Build a private, pinned HwSim ancestry plugin; never install over shared HwSim.

Run inside the intended Herwig/ROOT environment. The source repository may
have local changes: only the recorded Git revision is read.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shlex
import subprocess

BASE = "51ba33aea55ddbee0df8309d425b6fbdde1ceee3"
PATCH = Path(__file__).with_name("ancestry.patch")


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def flags(program, option):
    result = subprocess.run([program, option], text=True, capture_output=True)
    values = shlex.split(result.stdout)
    # ThePEG 2.3's shell helper returns the status of its final unrelated
    # `test`, even after successfully printing --cppflags.
    helper_quirk = (program == "thepeg-config" and option == "--cppflags"
                    and result.returncode == 1 and values
                    and all(v.startswith("-I") and Path(v[2:]).is_dir() for v in values))
    if result.returncode and not helper_quirk:
        raise RuntimeError(f"{program} {option} failed: {result.stderr}")
    return values


def build(source, output, cxx="c++"):
    identity = {"base_revision": BASE, "patch_sha256": digest(PATCH),
                "builder_sha256": digest(Path(__file__)), "cxx": cxx}
    proof = output / "build.json"
    if proof.exists():
        saved = json.loads(proof.read_text())
        if saved["identity"] != identity or digest(output / "HwSimAncestry.so") != saved["plugin_sha256"]:
            raise ValueError("private plugin build changed; use a new output directory")
        return saved
    if output.exists() and any(output.iterdir()):
        raise ValueError("nonempty build without proof; preserve it and use a new directory")
    output.mkdir(parents=True, exist_ok=True)
    names = subprocess.check_output(["git", "-C", str(source), "ls-tree", "--name-only", BASE], text=True).splitlines()
    names = [name for name in names if name.endswith(".h") or name == "HwSim.cc"]
    for name in names:
        (output / name).write_bytes(subprocess.check_output(
            ["git", "-C", str(source), "show", BASE + ":" + name]))
    subprocess.run(["git", "apply", "--unsafe-paths", "--directory=" + str(output),
                    str(PATCH.resolve())], cwd="/tmp", check=True)
    command = [cxx, "-O2", "-fPIC", "-shared", *flags("thepeg-config", "--cppflags"),
               *flags("root-config", "--cflags"), *flags("fastjet-config", "--cxxflags"),
               "HwSim.cc", "-o", "HwSimAncestry.so", *flags("fastjet-config", "--libs"),
               *flags("root-config", "--glibs"), "-lGenVector"]
    with (output / "build.log").open("w") as log:
        log.write(shlex.join(command) + "\n")
        log.flush()
        subprocess.run(command, cwd=output, stdout=log, stderr=subprocess.STDOUT, check=True)
    saved = {"identity": identity, "command": command, "plugin_sha256": digest(output / "HwSimAncestry.so"),
             "source_sha256": {p: digest(output / p) for p in names},
             "schema_version": 1}
    proof.write_text(json.dumps(saved, indent=2) + "\n")
    return saved


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-repo", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--cxx", default="c++")
    args = parser.parse_args()
    build(args.source_repo.resolve(), args.output_dir.resolve(), args.cxx)
    print(args.output_dir.resolve() / "HwSimAncestry.so")


if __name__ == "__main__":
    main()
