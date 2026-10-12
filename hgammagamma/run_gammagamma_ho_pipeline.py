#!/usr/bin/env python3
"""User-owned setup and complete SSC40 HO diphoton pipeline.

The existing generator/normalization/campaign wrappers remain responsible for
physics validation and resume checks. This entry point supplies their missing
installation, environment, ordering and first-run checks.
"""
from __future__ import annotations

import argparse
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import tarfile
import tempfile
from types import SimpleNamespace

try:
    from . import run_gammagamma_campaign as runtime
    from . import run_powheg_hjminnlo as powheg
    from . import ho_signal_normalization as norm
    from . import run_gammagamma_ho_campaign as ho_campaign
except ImportError:
    import run_gammagamma_campaign as runtime
    import run_powheg_hjminnlo as powheg
    import ho_signal_normalization as norm
    import run_gammagamma_ho_campaign as ho_campaign

REPO = Path(__file__).resolve().parent.parent
HO = REPO / "hgammagamma/HOAnalysis"
POWHEG_URL = "https://gitlab.com/POWHEG-BOX/V2/POWHEG-BOX-V2.git"
POWHEG_COMMIT = "d20df70ba6708da3978d5ab32d91b2a069457c1a"
HJ_COMMIT = "222b3d202bf0ef966a386e0e4d87c158f2a4681f"
MG5_COMMIT = "1d7559658bdc2b7eba476d4897aa594994fa82af"  # upstream v3.5.15
MG5_URL = f"https://codeload.github.com/mg5amcnlo/mg5amcnlo/tar.gz/{MG5_COMMIT}"
MG5_SHA256 = "e502f6ac9d2e35e00b9382230228284a7dbe105615bb57dfddb703793b6df03c"
CUBA_URL = "https://feynarts.de/cuba/Cuba-4.2.2.tar.gz"
CUBA_SHA256 = "8d9f532fd2b9561da2272c156ef7be5f3960953e4519c638759f1b52fe03ed52"
PDF_SOURCE = "https://lhapdfsets.web.cern.ch/current/"
PDF_SETS = ("NNPDF40_nlo_as_01180_qed", "NNPDF40_nnlo_as_01180_qed",
            "NNPDF40_an3lo_as_01180_qed_mhou", "PDF4LHC15_nnlo_100")
LO_PDF_SETS = ("NNPDF40_lo_as_01180",)
WORKFLOWS = ("gammagamma-lo", "gammagamma-ho", "fourlepton-lo")
CORE_PACKAGES = ("numpy", "matplotlib", "PyYAML", "six")
SAMPLES = ("signal_gg_h_aa", "bkg_prompt_aa", "bkg_gamma_j", "bkg_dy_ee")


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def save(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def frozen(path, data):
    if path.exists():
        if json.loads(path.read_text()) != data:
            raise ValueError(f"configuration differs from {path}; keep the original options or use a new directory/tag")
    else:
        save(path, data)


def repair_mg5_startup(path):
    """Repair the pinned 3.5.15 warning's missing early import on Python 3.9."""
    before = "import sys\nif sys.version_info < (3, 7):"
    after = "import logging\n" + before
    text = path.read_text()
    if after in text:
        return
    if text.count(before) != 1:
        raise ValueError(f"MG5 startup repair anchor changed: {path}")
    backup = path.with_name(path.name + ".before-pipeline")
    if not backup.exists():
        shutil.copy2(path, backup)
    path.write_text(text.replace(before, after))


def prepare_mg5_archive(directory):
    """Apply the release-tag archive's packaging step, keeping the original."""
    if not re.search(r"(?m)^version\s*=\s*3\.5\.15\s*$", (directory / "VERSION").read_text()):
        raise ValueError("MG5 archive is not the pinned 3.5.15 release")
    utility = directory / "bin/create_release.py"
    if utility.exists():
        # Upstream create_release.py removes this development-only utility.
        # Its presence enables an ANSI-coloured development banner in LHE
        # headers, which is invalid XML. Keep it outside bin for provenance.
        saved = directory / "source-packaging/create_release.py"
        saved.parent.mkdir(exist_ok=True)
        if saved.exists():
            raise ValueError(f"both active and archived MG5 release utilities exist: {directory}")
        utility.rename(saved)
    repair_mg5_startup(directory / "bin/mg5_aMC")


def pdf_files(directory, name):
    return [directory / name / (name + ".info"),
            *[directory / name / f"{name}_{member:04d}.dat" for member in range(101)]]


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("action", choices=("setup", "check", "run", "status"))
    p.add_argument("--work-dir", type=Path, default=HO / "inputs/pipeline",
                   help="private dependencies and normalization; default is inside this checkout")
    p.add_argument("--herwig-module", default=None, help="setup runtime; default on Linux: herwig/stable")
    p.add_argument("--herwig-env", type=Path, help="setup runtime: Herwig prefix or activation script")
    p.add_argument("--cc", default=None, help="setup C compiler; default /usr/bin/gcc")
    p.add_argument("--cxx", default=None, help="setup C++ compiler; default /usr/bin/g++")
    p.add_argument("--fc", default=None, help="setup Fortran compiler; default /usr/bin/gfortran")
    p.add_argument("--ninja-dir", type=Path, help="optional existing Ninja library directory for setup")
    p.add_argument("--jobs", type=int, default=4)
    p.add_argument("--run-tag", default="ho_pilot_01")
    p.add_argument("--nevents", type=int, default=30, help="30-event pilot by default; production requires a completed pilot")
    p.add_argument("--luminosity-fb", type=float, default=100.)
    p.add_argument("--signal-lhe", type=Path, help="optional existing signal; otherwise generate it with POWHEG")
    p.add_argument("--signal-normalization", type=Path,
                   help="optional validated record; otherwise calculate all 109 ihixs points locally")
    p.add_argument("--dry-run", action="store_true", help="print the plan; do not install or run programs")
    args = p.parse_args(argv)
    if args.jobs < 1 or args.nevents < 1 or not 0 < args.luminosity_fb < float("inf"):
        p.error("jobs, nevents and luminosity must be positive and finite")
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]*", args.run_tag):
        p.error("invalid run tag")
    if args.herwig_module and args.herwig_env:
        p.error("select either a module or an activation script")
    for key in ("work_dir", "signal_lhe", "signal_normalization", "herwig_env", "ninja_dir"):
        if getattr(args, key) is not None:
            setattr(args, key, getattr(args, key).expanduser().resolve())
    for path in (REPO, args.work_dir, args.signal_lhe, args.signal_normalization):
        if path and re.search(r"\s|[;#\"']", str(path)):
            p.error("generator paths must not contain whitespace or card metacharacters")
    return args


def configuration(args):
    activation = args.herwig_env
    if activation and activation.is_dir():
        activation /= "bin/activate"
    path = args.work_dir / "environment.json"
    if path.exists():
        config = json.loads(path.read_text())
        if config["repository"] != str(REPO) or config["work_dir"] != str(args.work_dir):
            raise ValueError("setup belongs to a different checkout; use a new --work-dir")
        for key in ("herwig_module", "herwig_env", "cc", "cxx", "fc", "ninja_dir"):
            value = activation if key == "herwig_env" else getattr(args, key)
            if value is not None and str(value) != str(config[key]):
                raise ValueError(f"{key} differs from saved setup; use a new --work-dir")
        return config
    return {"schema_version": 1, "repository": str(REPO), "work_dir": str(args.work_dir),
            "herwig_env": str(activation) if activation else None,
            "herwig_module": None if activation else (args.herwig_module or "herwig/stable"),
            "cc": args.cc or "/usr/bin/gcc", "cxx": args.cxx or "/usr/bin/g++",
            "fc": args.fc or "/usr/bin/gfortran",
            "ninja_dir": str(args.ninja_dir) if args.ninja_dir else None}


def activated_environment(config):
    from types import SimpleNamespace
    options = SimpleNamespace(herwig_env=config["herwig_env"], herwig_module=config["herwig_module"])
    commands = runtime.runtime_setup_commands(options)
    commands.append("python3 -c " + shlex.quote(
        "import json,os; print('HIGGSSSC_ENV='+json.dumps(dict(os.environ)))"))
    result = subprocess.run(["bash", "-c", "\n".join(commands)], env=runtime.clean_shell_env(),
                            text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if result.returncode or "HIGGSSSC_ENV=" not in result.stdout:
        raise ValueError("cannot activate Herwig runtime:\n" + result.stderr.strip())
    env = json.loads(result.stdout.rsplit("HIGGSSSC_ENV=", 1)[1])
    # Children inherit activated paths, not another wrapper's activation
    # defaults. Re-sourcing a site script could reset the private PDF path.
    for key in ("HERWIG_ENV", "HERWIG_MODULE", "DEFAULT_HERWIG_ENV", "HERWIG", "MG5_PYTHON"):
        env.pop(key, None)
    env.update(CC=config["cc"], CXX=config["cxx"], FC=config["fc"], F77=config["fc"])
    return env


class Pipeline:
    def __init__(self, args, config):
        self.args, self.config, self.work = args, config, args.work_dir
        self.env = activated_environment(config)
        self.prefix = Path(self.capture(["lhapdf-config", "--prefix"]))
        self.pdf_dir = self.work / "pdfs"
        self.env["LHAPDF_DATA_PATH"] = str(self.pdf_dir) + ":" + self.capture(["lhapdf-config", "--datadir"])
        self.python = self.work / "venv/bin/python"
        self.env["PATH"] = str(self.python.parent) + ":" + self.env["PATH"]
        self.env["PYTHONNOUSERSITE"] = "1"
        self.mg5 = self.work / "MG5_aMC_v3_5_15"
        self.powheg = self.work / "POWHEG-BOX-V2"
        self.cuba = self.work / "cuba"
        self.logs = self.work / "logs"
        self.boost = self.prefix / "include"
        if not (self.boost / "boost/version.hpp").is_file():
            self.boost = Path("/usr/include")

    def capture(self, argv):
        return subprocess.check_output([str(x) for x in argv], env=self.env, text=True).strip()

    def call(self, label, argv, cwd=None, extra_env=None, mg5=False):
        self.logs.mkdir(parents=True, exist_ok=True)
        log = self.logs / (label + ".log")
        print(f"[{label}] {shlex.join(str(x) for x in argv)}\n  Log: {log}", flush=True)
        # Append: failed attempts remain available when setup/run is resumed.
        with log.open("a+", errors="replace") as stream:
            stream.write("\nCOMMAND: " + shlex.join(str(x) for x in argv) + "\n")
            stream.flush()
            start = stream.tell()
            result = subprocess.run([str(x) for x in argv], cwd=cwd or REPO, env={**self.env, **(extra_env or {})},
                                    stdout=stream, stderr=subprocess.STDOUT)
            stream.seek(start)
            output = stream.read()
        if result.returncode or mg5 and re.search(r"interrupted with error:|Error detected in", output):
            tail = "\n".join(output.splitlines()[-50:])
            reason = f"exit status {result.returncode}" if result.returncode else "MG5 reported an error despite exit status 0"
            raise RuntimeError(f"{label} failed ({reason}); {log}\n{tail}")

    def download(self, url, filename, expected=None):
        target = self.work / "downloads" / filename
        target.parent.mkdir(parents=True, exist_ok=True)
        if not target.exists():
            temporary = target.with_suffix(target.suffix + ".part")
            self.call("download-" + filename,
                      ["curl", "--fail", "--location", "--retry", "2", url, "--output", temporary])
            temporary.replace(target)
        if expected and digest(target) != expected:
            raise ValueError(f"download checksum differs: {target}; preserve it for inspection")
        return target

    def unpack(self, archive, name, destination):
        if destination.exists():
            return
        self.work.mkdir(parents=True, exist_ok=True)
        staging = Path(tempfile.mkdtemp(prefix="unpack-", dir=self.work))
        with tarfile.open(archive) as tar:
            if not hasattr(tarfile, "data_filter"):
                raise ValueError("update Python to a security-patched release with tarfile.data_filter")
            tar.extractall(staging, filter="data")
        source = staging / name
        if not source.is_dir():
            raise ValueError(f"expected {name} in {archive}; extracted files retained at {staging}")
        source.rename(destination)
        staging.rmdir()

    def runtime_check(self):
        missing = [x for x in ("git", "curl", "make", "cmake", "Herwig", "lhapdf", "lhapdf-config",
                              "fastjet-config", "root-config", self.config["cc"], self.config["cxx"],
                              self.config["fc"]) if not shutil.which(x, path=self.env["PATH"])]
        if missing:
            raise ValueError("missing runtime tools: " + ", ".join(missing))
        self.call("check-pyroot", ["python3", "-c", "import ROOT; print(ROOT.gROOT.GetVersion())"])
        if not (self.boost / "boost/version.hpp").is_file():
            raise ValueError("Boost headers are missing from the Herwig prefix and /usr/include")
        probe = self.work / "checks"
        probe.mkdir(parents=True, exist_ok=True)
        card = probe / "hwsim.in"
        card.write_text("library HwSim.so\ncreate Herwig::HwSim /Herwig/Analysis/HiggsSSCPipelineProbe\n")
        self.call("check-hwsim", ["Herwig", "read", card], cwd=probe)

    @property
    def workflow(self):
        workflow = getattr(self.args, "workflow", "gammagamma-ho")
        if workflow not in WORKFLOWS:
            raise ValueError(f"unknown workflow: {workflow}")
        return workflow

    @property
    def pdf_sets(self):
        return PDF_SETS if self.workflow == "gammagamma-ho" else LO_PDF_SETS

    @property
    def with_ml(self):
        return getattr(self.args, "with_ml", False)

    def components(self):
        path = self.work / "components.json"
        return json.loads(path.read_text()) if path.is_file() else {}

    def record_component(self, name, metadata):
        components = self.components()
        components[name] = metadata
        save(self.work / "components.json", components)

    def setup_python(self):
        created = not self.python.exists()
        if created:
            self.call("python-venv", ["python3", "-m", "venv", "--system-site-packages", self.work / "venv"])
        components = self.components()
        # The original HO installer already recorded its working environment.
        # Do not let a later setup silently upgrade dependencies of saved runs.
        legacy = (self.work / "installed.json").is_file() and not created
        if created or "python_core" not in components and not legacy:
            requirements = REPO / "hgammagamma/HOAnalysis/requirements-pipeline.txt"
            self.call("python-packages", [self.python, "-m", "pip", "install", "-r", requirements])
        self.call("check-python", [self.python, "-c", "import ROOT,numpy,matplotlib,yaml,six"])
        if "python_core" not in components:
            self.record_component("python_core", {"packages": self.capture([self.python, "-m", "pip", "freeze"])})
        if self.with_ml:
            if created or "python_ml" not in components:
                # An optional classifier must not change the plotting/analysis
                # dependencies that produced an already validated pilot.
                versions = json.loads(self.capture([self.python, "-c",
                    "import importlib.metadata,json; print(json.dumps({name: importlib.metadata.version(name) "
                    f"for name in {CORE_PACKAGES!r}}}))"]))
                constraints = self.work / "ml-core-constraints.txt"
                constraints.write_text("".join(f"{name}=={versions[name]}\n" for name in CORE_PACKAGES))
                self.call("python-ml", [self.python, "-m", "pip", "install", "--constraint", constraints,
                          "--requirement", REPO / "requirements-ml.txt"])
            self.call("check-python-ml", [self.python, "-c", "import xgboost,sklearn,tqdm"])
            if "python_ml" not in components:
                self.record_component("python_ml", {"packages": self.capture([self.python, "-m", "pip", "freeze"])})

    def setup_cuba(self):
        if not (self.cuba / "lib/libcuba.a").is_file():
            archive = self.download(CUBA_URL, "Cuba-4.2.2.tar.gz", CUBA_SHA256)
            source = self.work / "Cuba-4.2.2"
            self.unpack(archive, "Cuba-4.2.2", source)
            self.call("cuba-configure", ["./configure", f"--prefix={self.cuba}",
                      f"--libdir={self.cuba}/lib", f"--includedir={self.cuba}/include", "--with-real=8"], source,
                      extra_env={"CFLAGS": "-O2 -fPIC"})
            self.call("cuba-build", ["make", "-j1", "lib"], source)
            self.call("cuba-install", ["make", "-j1", "install", "TOOLS_DEFAULT="], source)
        self.record_component("cuba", {"archive_sha256": CUBA_SHA256,
                                      "library_sha256": digest(self.cuba / "lib/libcuba.a")})

    def install_lo_model(self):
        source = REPO / "MG5_aMC_v3_5_15/models/loop_sm_haa"
        tracked = self.capture(["git", "-C", REPO, "ls-files", "--", source]).splitlines()
        files = [REPO / name for name in tracked if Path(name).suffix in (".py", ".dat")]
        if not files or not any(path.name == "parameters.py" for path in files):
            raise ValueError("tracked loop_sm_haa model source is missing")
        destination = self.mg5 / "models/loop_sm_haa"
        for path in files:
            target = destination / path.relative_to(source)
            if target.exists() and digest(target) != digest(path):
                raise ValueError(f"installed LO model differs from tracked source: {target}; use a new --work-dir")
        for path in files:
            target = destination / path.relative_to(source)
            if not target.exists():
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(path, target)
        self.record_component("lo_model", {str(path.relative_to(source)): digest(path) for path in files})

    def setup_mg5(self):
        if not (self.mg5 / "bin/mg5_aMC").is_file():
            archive = self.download(MG5_URL, "MG5-v3.5.15.tar.gz", MG5_SHA256)
            self.unpack(archive, f"mg5amcnlo-{MG5_COMMIT}", self.mg5)
        prepare_mg5_archive(self.mg5)
        candidates = ([Path(self.config["ninja_dir"])] if self.config["ninja_dir"] else
                      sorted(self.prefix.glob("opt/MG5*/HEPTools/ninja/lib")))
        ninja = next((p for p in candidates if (p / "libninja.a").is_file()
                      and (p / "libavh_olo.a").is_file() and (p.parent / "include/mninja.mod").is_file()), None)
        if self.config["ninja_dir"] and ninja is None:
            raise ValueError("--ninja-dir must contain Ninja/OneLOop libraries and ../include/mninja.mod")
        # These are this entry point's private MG5 files, never a shared MG5 installation.
        cfg = self.mg5 / "input/mg5_configuration.txt"
        # GitHub source archives omit the user configuration; MG5 supplies
        # defaults internally and writes this file on its first invocation.
        original = cfg.read_text() if cfg.exists() else ""
        text = original
        values = {"ninja": str(ninja) if ninja else "None", "collier": "None", "auto_update": "0",
                  "fortran_compiler": self.config["fc"], "cpp_compiler": self.config["cxx"],
                  "fastjet": self.capture(["which", "fastjet-config"]),
                  "lhapdf": self.capture(["which", "lhapdf-config"])}
        for key, value in values.items():
            text = re.sub(rf"(?m)^[ \t]*{key}[ \t]*=[^\n]*(?:\n|$)", "", text)
            text = text.rstrip() + f"\n{key} = {value}\n"
        if original != text:
            backup = cfg.with_name("mg5_configuration.before-pipeline.txt")
            if cfg.exists() and not backup.exists():
                shutil.copy2(cfg, backup)
            cfg.write_text(text)

        if self.workflow == "gammagamma-ho":
            # MG5 sometimes exits zero even when exporting a process failed. Test
            # that boundary before investing in integrations. Each attempt is kept.
            probe_marker = self.work / "mg5-export-checked.json"
            probe_identity = {"configuration_sha256": digest(cfg), "release_packaging": True,
                              "startup_sha256": digest(self.mg5 / "bin/mg5_aMC")}
            previous = json.loads(probe_marker.read_text()) if probe_marker.exists() else {}
            if previous.get("identity") != probe_identity:
                probe = Path(tempfile.mkdtemp(prefix="mg5-export-", dir=self.work / "checks"))
                card = probe / "process.mg5"
                card.write_text("set automatic_html_opening False\nimport model loop_sm-no_b_mass\n"
                                "generate p p > a a [QCD]\n" + f"output {probe / 'process'}\nquit\n")
                self.call("check-mg5-export", [self.python, self.mg5 / "bin/mg5_aMC", card], probe, mg5=True)
                if not (probe / "process/Cards/run_card.dat").is_file():
                    raise ValueError(f"MG5 did not export the probe process; inspect {probe}")
                save(probe_marker, {"process": str(probe / "process"), "identity": probe_identity})

        self.record_component("mg5", {"commit": MG5_COMMIT,
                              "configuration_sha256": digest(cfg),
                              "startup_sha256": digest(self.mg5 / "bin/mg5_aMC"),
                              "ninja": str(ninja) if ninja else None})
        if self.workflow == "gammagamma-lo":
            self.install_lo_model()

    def setup_powheg(self):
        if not self.powheg.exists():
            self.call("powheg-clone", ["git", "clone", "--filter=blob:none", "--sparse", "--no-checkout",
                                       POWHEG_URL, self.powheg])
        if not (self.work / "powheg-source-checked.json").exists():
            # Only an installation created under this managed work directory
            # reaches here. Resume an interrupted no-checkout clone, preserving
            # any conflicting changes through Git's ordinary checkout checks.
            self.call("powheg-sparse", ["git", "-C", self.powheg, "sparse-checkout", "set",
                                        "include", "svnversion", "MiNNLOStuff", "HJ"])
            self.call("powheg-pin", ["git", "-C", self.powheg, "checkout", POWHEG_COMMIT])
        if self.capture(["git", "-C", self.powheg, "rev-parse", "HEAD"]) != POWHEG_COMMIT:
            raise ValueError("existing private POWHEG checkout has a different commit; preserve it and use a new work directory")
        self.call("powheg-hj", ["git", "-C", self.powheg, "submodule", "update", "--init", "HJ"])
        if self.capture(["git", "-C", self.powheg / "HJ", "rev-parse", "HEAD"]) != HJ_COMMIT:
            raise ValueError("HJ submodule differs from the pinned version")
        save(self.work / "powheg-source-checked.json", {"powheg": POWHEG_COMMIT, "hj": HJ_COMMIT})
        process = self.powheg / "HJ/HJMiNNLO"
        (process / "obj-gfortran").mkdir(exist_ok=True)
        if not (process / "pwhg_main").is_file():
            self.call("powheg-build", ["make", "-j1", "pwhg_main", f"CXX={self.config['cxx']}",
                      f"CC={self.config['cc']}", f"FC={self.config['fc']}", "STDCLIB=-lstdc++"], process)
        self.record_component("powheg", {"commit": POWHEG_COMMIT, "hj_commit": HJ_COMMIT,
                                         "executable_sha256": digest(process / "pwhg_main")})

    def setup_analysis(self):
        if self.workflow == "fourlepton-lo":
            directory = REPO / "hfourlepton/LOAnalysis/Code"
            executable = "HwSimPostAnalysis_fourlepton"
            command = ["make", "-C", directory, executable, f"CXX={self.config['cxx']}"]
            component = "fourlepton_analysis"
        else:
            directory = REPO / "hgammagamma/LOAnalysis/Code"
            executable = "HwSimPostAnalysis_gammagamma_SSC"
            command = ["make", "-C", directory, executable, f"CPP={self.config['cxx']}",
                       f"CXX={self.config['cxx']}", "FASTJET_CPPFLAGS=", "FASTJET_LDFLAGS="]
            component = "gammagamma_analysis"
        self.call("analysis-build", command)
        self.record_component(component, {"executable_sha256": digest(directory / executable)})

    def setup(self):
        self.runtime_check()
        self.setup_python()
        self.pdf_dir.mkdir(parents=True, exist_ok=True)
        pdf_cmd = ["lhapdf", "--quiet", "--listdir", self.pdf_dir, "--pdfdir", self.pdf_dir, "--source", PDF_SOURCE]
        if not (self.pdf_dir / "pdfsets.index").is_file():
            self.call("pdf-index", [*pdf_cmd, "update"])
        self.install_pdfs(pdf_cmd, self.pdf_sets)
        for name in self.pdf_sets:
            self.record_component("pdf:" + name, {"members": 101,
                                  "info_sha256": digest(self.pdf_dir / name / (name + ".info"))})
        self.setup_mg5()
        if self.workflow == "gammagamma-ho":
            self.setup_cuba()
            self.setup_powheg()
            if not (REPO / "external/ihixs/CMakeLists.txt").is_file():
                self.call("ihixs-source", ["git", "submodule", "update", "--init", "external/ihixs"])
            self.record_component("ihixs", {"commit": self.capture(["git", "-C", REPO / "external/ihixs", "rev-parse", "HEAD"])})
        self.setup_analysis()
        self.check()
        # Capabilities grow in components.json. Keep the legacy provenance file
        # and every previously validated pilot byte-for-byte on extension.
        if not (self.work / "installed.json").exists():
            save(self.work / "installed.json", {"environment": self.config,
                 "downloads": {p.name: digest(p) for p in (self.work / "downloads").glob("*") if p.is_file()},
                 "components": self.components()})
        self.record_component(self.workflow, {"ready": True})
        print(f"Setup complete for {self.workflow}. Run its 30-event pilot next.")

    def install_pdfs(self, command, names=None):
        names = PDF_SETS if names is None else names
        missing = [name for name in names if any(not p.is_file() or not p.stat().st_size
                   for p in pdf_files(self.pdf_dir, name))]
        if missing:
            # LHAPDF otherwise skips a set found anywhere on its search path,
            # even when --pdfdir requests a different, user-owned destination.
            # Never reinstall a complete private set used by a saved run.
            self.call("pdf-install", [*command, "install", "--upgrade", *missing])
        for name in names:
            if any(not p.is_file() or not p.stat().st_size for p in pdf_files(self.pdf_dir, name)):
                raise ValueError(f"LHAPDF did not install all members of {name}; see pdf-install.log")

    def check(self):
        self.runtime_check()
        required = [self.python, self.mg5 / "bin/mg5_aMC", self.pdf_dir / "pdfsets.index"]
        if self.workflow == "gammagamma-ho":
            required.extend([self.cuba / "include/cuba.h", self.cuba / "lib/libcuba.a",
                             self.powheg / "HJ/HJMiNNLO/pwhg_main", REPO / "external/ihixs/CMakeLists.txt"])
        if self.workflow == "fourlepton-lo":
            required.append(REPO / "hfourlepton/LOAnalysis/Code/HwSimPostAnalysis_fourlepton")
        else:
            required.append(REPO / "hgammagamma/LOAnalysis/Code/HwSimPostAnalysis_gammagamma_SSC")
        if self.workflow == "gammagamma-lo":
            required.append(self.mg5 / "models/loop_sm_haa/parameters.py")
            for name, expected in self.components().get("lo_model", {}).items():
                path = self.mg5 / "models/loop_sm_haa" / name
                if not path.is_file() or digest(path) != expected:
                    raise ValueError(f"installed LO model changed: {path}; preserve existing runs and use a new --work-dir")
        for name in self.pdf_sets:
            required.extend(pdf_files(self.pdf_dir, name))
        missing = [str(p) for p in required if not p.is_file() or not p.stat().st_size]
        if missing:
            raise ValueError(f"setup is incomplete for {self.workflow}; run setup first. Missing:\n" + "\n".join(missing[:12]))
        imports = "import ROOT,numpy,matplotlib,yaml,six"
        if self.with_ml:
            imports += ",xgboost,sklearn,tqdm"
        self.call("check-python", [self.python, "-c", imports + "; print('Python dependencies ready')"])
        print(f"Dependency checks passed for {self.workflow}.")

    def preflight_inputs(self):
        """Reject supplied inputs before builds or an expensive ihixs calculation."""
        args = self.args
        record = args.signal_normalization
        validated = None
        if record:
            if not record.is_file():
                raise ValueError(f"supplied normalization is missing: {record}")
            validated = norm.load_record(record)
        if args.signal_lhe:
            if not args.signal_lhe.is_file():
                raise ValueError(f"supplied signal LHE is missing: {args.signal_lhe}")
            options = SimpleNamespace(nevents=args.nevents, ebeam=20000.,
                      expected_higgs_mass=validated["higgs_mass_gev"] if validated else 125.,
                      expected_native_parameters=validated.get("native_parameters", {}) if validated else {})
            ho_campaign.inspect_lhe(args.signal_lhe, ho_campaign.SAMPLES[0], options)
            ho_campaign.shower_completion.source_weight_summary(args.signal_lhe, args.nevents)

    def run(self):
        self.preflight_inputs()
        self.check()
        args = self.args
        if args.nevents > 30:
            verify_pilot(self.work)
        directory = HO / "runs" / args.run_tag
        marker = directory / "pipeline.json"
        if directory.exists() and not marker.exists() and any(directory.iterdir()):
            raise ValueError("this tag belongs to an existing campaign; use its original runner or a new pipeline tag")
        identity = {"environment": self.config, "nevents": args.nevents, "jobs": args.jobs,
                    "signal_lhe": str(args.signal_lhe) if args.signal_lhe else None,
                    "normalization": str(args.signal_normalization) if args.signal_normalization else None,
                    "powheg_card_sha256": digest(HO / "powheg-hjminnlo-ssc40-nnpdf40nnloqed.input")}
        frozen(marker, identity)
        record = args.signal_normalization or self.work / "normalization/ggf-ssc40-n3lo.json"
        if record.exists():
            norm.load_record(record)
        elif args.signal_normalization:
            raise ValueError(f"supplied normalization is missing: {record}")
        else:
            options = ["--work-dir", self.work / "normalization/ihixs", "--record", record,
                       "--lhapdf-dir", self.prefix, "--cuba-dir", self.cuba,
                       "--boost-dir", self.boost,
                       "--cc", self.config["cc"], "--cxx", self.config["cxx"], "--jobs", args.jobs, "--resume"]
            for stage in ("build", "benchmark", "calculate"):
                self.call("ihixs-" + stage, [self.python, REPO / "hgammagamma/run_ihixs_normalization.py",
                          "--stage", stage, *options, *(["--refine-failed"] if stage == "calculate" else [])])
            norm.load_record(record)
        lhe = args.signal_lhe or directory / "powheg/powheg-hjminnlo-merged.lhe"
        if not args.signal_lhe:
            card = directory / "powheg-pipeline.input"
            text = powheg.patched_input((HO / "powheg-hjminnlo-ssc40-nnpdf40nnloqed.input").read_text(),
                                       {"ncall2": 500000, "check_bad_st2": 1})
            if card.exists() and card.read_text() != text:
                raise ValueError("saved POWHEG card changed; preserve this run")
            card.write_text(text)
            self.call(args.run_tag + "-powheg", [self.python, REPO / "hgammagamma/run_powheg_hjminnlo.py",
                      "--nevents", args.nevents, "--jobs", args.jobs, "--powheg-dir", self.powheg / "HJ/HJMiNNLO",
                      "--run-dir", directory / "powheg", "--card", card,
                      "--no-herwig-env", "--no-herwig-module", "--skip-build", "--resume"])
        if not lhe.is_file():
            raise ValueError(f"signal LHE is missing: {lhe}")
        self.call(args.run_tag + "-campaign", [self.python, REPO / "hgammagamma/run_gammagamma_ho_campaign.py",
                  "--stage", "all", "--run-tag", args.run_tag, "--nevents", args.nevents,
                  "--nb-core", args.jobs, "--mg5-dir", self.mg5, "--mg5-python", self.python,
                  "--mg5-fortran", self.config["fc"], "--mg5-cxx", self.config["cxx"],
                  "--analysis-cxx", self.config["cxx"], "--no-herwig-module",
                  "--signal-normalization", record, "--signal-lhe", lhe, "--resume"])
        report = HO / "plots" / args.run_tag
        self.call(args.run_tag + "-report", [self.python, REPO / "hgammagamma/make_gammagamma_report.py",
                  "--analysis-root", directory, "--run-tag", args.run_tag, "--output-dir", report,
                  "--no-density", "--normalization", "event_xsec", "--luminosity-fb", args.luminosity_fb])
        artifacts = completed_artifacts(directory, report)
        if args.nevents <= 30:
            save(self.work / "pilot-complete.json", {"environment": self.config,
                 "artifacts": {name: value for name, value in artifacts.items()
                               if name != str(report / "index.html")}})
        print(f"Pipeline complete: {report / 'index.html'}")


def completed_artifacts(directory, report):
    paths = [report / "index.html"]
    for sample in SAMPLES:
        category = "Signal" if sample.startswith("signal_") else "Backgrounds"
        folder = directory / category / "events" / sample
        manifest = folder / "campaign.json"
        if not manifest.is_file() or not {"generate", "shower", "analyze"}.issubset(
                json.loads(manifest.read_text()).get("completed", [])):
            raise ValueError(f"analysis has not completed: {sample}")
        paths.append(manifest)
        paths.extend([folder / f"{sample}_hwsim_roots-{directory.name}.dat",
                      folder / f"{sample}_hwsim_roots-{directory.name}.top",
                      folder / f"{sample}_hwsim_roots-{directory.name}_var.root",
                      folder / f"normalization-{directory.name}.json"])
    return {str(p): digest(p) for p in paths}


def verify_pilot(work):
    path = work / "pilot-complete.json"
    if not path.is_file():
        raise ValueError("first complete a pilot: run --run-tag ho_pilot_01 --nevents 30")
    pilot = json.loads(path.read_text())
    if pilot["environment"] != json.loads((work / "environment.json").read_text()):
        raise ValueError("pilot used a different environment")
    for name, expected in pilot["artifacts"].items():
        if not Path(name).is_file() or digest(name) != expected:
            raise ValueError(f"pilot artifact changed: {name}; validate a new pilot before production")


def main(argv=None):
    args = parse_args(argv)
    if args.action == "status":
        directory = HO / "runs" / args.run_tag
        print(f"Campaign: {directory}")
        for path in sorted(directory.glob("*/events/*/campaign.json")):
            data = json.loads(path.read_text())
            print(data["sample"], ", ".join(data.get("completed", [])) or "prepared")
        print("Report:", (HO / "plots" / args.run_tag / "index.html").is_file())
        return 0
    config = configuration(args)
    if args.dry_run:
        stages = (["check shared runtime", "private Python/PDF/MG5/Cuba installation",
                   "NLO process export check", "pinned POWHEG build", "ihixs source", "analysis build"]
                  if args.action == "setup" else ["check installed runtime and dependencies"]
                  if args.action == "check" else [
                      "ihixs build/benchmark/109-point calculation (or validate record)",
                      "POWHEG signal (or validate supplied LHE)",
                      "MG5 backgrounds", "Herwig showers", "analysis", "report"])
        print(json.dumps({"action": args.action, "environment": config, "run_tag": args.run_tag,
                         "nevents": args.nevents, "pdf_sets": PDF_SETS,
                         "stages": stages}, indent=2))
        return 0
    if args.action != "setup" and not (args.work_dir / "environment.json").is_file():
        raise ValueError("run setup first")
    if args.action == "setup" and args.work_dir.exists() and not (args.work_dir / "environment.json").exists():
        if any(args.work_dir.iterdir()):
            raise ValueError("unmanaged non-empty work directory; select a new --work-dir")
    args.work_dir.mkdir(parents=True, exist_ok=True)
    with (args.work_dir / ".pipeline.lock").open("a") as lock:
        try:
            fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise ValueError("another pipeline command is using this work directory") from exc
        frozen(args.work_dir / "environment.json", config)
        pipeline = Pipeline(args, config)
        getattr(pipeline, args.action)()
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError, RuntimeError, subprocess.CalledProcessError) as exc:
        raise SystemExit(f"Error: {exc}")
