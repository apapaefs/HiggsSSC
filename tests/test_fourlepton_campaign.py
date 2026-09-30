import gzip
import io
import json
import math
import re
import subprocess
import sys
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

from hfourlepton import run_four_lepton_campaign as campaign


REPO_ROOT = Path(__file__).resolve().parents[1]


def config(*extra: str) -> campaign.Config:
    with patch.dict("os.environ", {}, clear=True):
        return campaign.parse_config(["--dry-run", "--no-herwig-module", *extra])


def lhe_text(cross_sections: list[float], weights: list[float] | None = None) -> str:
    subprocesses = "\n".join(
        f" {value:.12e} 0.0 1.0 {index + 1}"
        for index, value in enumerate(cross_sections)
    )
    events = ""
    for weight in weights or []:
        events += (
            "<event>\n"
            f" 2 1 {weight:.12e} 9.118800e+01 7.297e-03 1.18e-01\n"
            "</event>\n"
        )
    return f"""<LesHouchesEvents version="3.0">
<header>
<slha>
BLOCK MASS
  25 1.250900000e+02 # MH
DECAY 25 4.070000000e-03
</slha>
</header>
<init>
 2212 2212 2.000000e+04 2.000000e+04 0 0 331900 331900 3 {len(cross_sections)}
{subprocesses}
</init>
{events}</LesHouchesEvents>
"""


class FourLeptonSampleTests(unittest.TestCase):
    def test_irreducible_set_has_three_signal_and_six_continuum_samples(self) -> None:
        selected = campaign.selected_samples(config("--sample-set", "irreducible"))
        self.assertEqual(len(selected), 9)
        self.assertEqual(sum(sample.category == "signal" for sample in selected), 3)
        self.assertEqual(sum(sample.category == "irreducible" for sample in selected), 6)

    def test_loop_continuum_excludes_higgs_and_qq_uses_full_final_states(self) -> None:
        loop_samples = [
            sample for sample in campaign.SAMPLES if sample.name.startswith("bkg_gg_")
        ]
        qq_samples = [
            sample for sample in campaign.SAMPLES if sample.name.startswith("bkg_qq_")
        ]
        self.assertEqual(len(loop_samples), 3)
        self.assertTrue(all(" / h [noborn=QCD]" in sample.process for sample in loop_samples))
        self.assertEqual({sample.channel for sample in qq_samples}, {"4e", "4mu", "2e2mu"})
        card = campaign.process_card_text(qq_samples[0], Path("/tmp/output"))
        self.assertIn("define p = g u c d s b u~ c~ d~ s~ b~", card)

    def test_signal_channels_share_one_stable_production_group(self) -> None:
        signal = [sample for sample in campaign.SAMPLES if sample.is_signal]
        self.assertEqual({sample.production_group_id for sample in signal}, {"signal_gg_h_stable"})
        self.assertEqual({sample.process for sample in signal}, {"g g > h [noborn=QCD]"})
        self.assertTrue(all(sample.madspin_decay.startswith("h > ") for sample in signal))

    def test_ttz_strata_are_disjoint_and_share_stable_production(self) -> None:
        ttz = [sample for sample in campaign.SAMPLES if sample.name.startswith("bkg_ttz_")]
        self.assertEqual(len(ttz), 8)
        self.assertEqual(len({sample.prompt_decay_profile for sample in ttz}), 8)
        self.assertEqual({sample.production_group_id for sample in ttz}, {"bkg_ttz_stable"})
        self.assertEqual({sample.process for sample in ttz}, {"p p > t t~ z"})

    def test_semileptonic_mode_adds_biased_strata_but_retains_nominal_controls(self) -> None:
        cfg = config(
            "--sample-set",
            "reducible",
            "--heavy-flavour-bias",
            "semileptonic",
            "--nevents",
            "20000",
            "--bias-control-events",
            "250",
        )
        selected = campaign.selected_samples(cfg)
        nominal = [sample for sample in selected if not sample.is_hf_biased]
        biased = [sample for sample in selected if sample.is_hf_biased]
        self.assertEqual(len(nominal), len(biased))
        self.assertTrue(all(cfg.events_for(sample) == 250 for sample in nominal))
        self.assertTrue(all(cfg.events_for(sample) == 20000 for sample in biased))
        self.assertEqual({sample.bias_of for sample in biased}, {sample.name for sample in nominal})


class FourLeptonConfigurationTests(unittest.TestCase):
    def test_defaults_are_40_tev_lo_perfect(self) -> None:
        cfg = config()
        self.assertEqual(cfg.ebeam, 20000.0)
        self.assertEqual(cfg.detector_response, "perfect")
        self.assertEqual(cfg.source_backend, "mg5")
        self.assertEqual(campaign.LO_PDF_LHAID, 331900)
        self.assertEqual(cfg.analysis_target, "HwSimPostAnalysis_fourlepton")

    def test_corrected_signal_fractions_are_positive_and_normalized(self) -> None:
        cfg = config()
        self.assertEqual(cfg.signal_channel_fractions, campaign.DEFAULT_SIGNAL_CHANNEL_FRACTIONS)
        self.assertAlmostEqual(sum(cfg.signal_channel_fractions.values()), 1.0)
        self.assertAlmostEqual(cfg.signal_channel_fractions["4e"], 0.263127188)
        self.assertAlmostEqual(cfg.signal_channel_fractions["2e2mu"], 0.473745624)

    def test_invalid_signal_fraction_override_fails_closed(self) -> None:
        with patch("sys.stderr", new=io.StringIO()):
            with self.assertRaises(SystemExit):
                config("--signal-channel-fractions", "4e=.25,4mu=.25,2e2mu=.4")

    def test_per_sample_event_override_and_seed_are_stable(self) -> None:
        cfg = config(
            "--nevents",
            "100",
            "--sample-nevents",
            "bkg_qq_4e=17",
        )
        sample = campaign.SAMPLES_BY_NAME["bkg_qq_4e"]
        self.assertEqual(cfg.events_for(sample), 17)
        self.assertEqual(campaign.seed_for(sample.name, cfg), campaign.seed_for(sample.name, cfg))
        self.assertNotEqual(
            campaign.seed_for(sample.name, cfg),
            campaign.seed_for("bkg_qq_4mu", cfg),
        )

    def test_external_lhe_routes_by_sample_or_production_group(self) -> None:
        cfg = config(
            "--source-backend",
            "external-lhe",
            "--external-lhe",
            "signal_gg_h_stable=/tmp/stable.lhe.gz",
            "--run-samples",
            "signal_gg_h_4e",
        )
        sample = campaign.selected_samples(cfg)[0]
        self.assertEqual(
            campaign.resolve_external_lhe(sample, cfg),
            Path("/tmp/stable.lhe.gz"),
        )

    def test_runtime_setup_discards_only_unresolvable_compiler_exports(self) -> None:
        snippet = campaign.sanitize_toolchain_exports_snippet()
        probe = subprocess.run(
            [
                "bash",
                "-c",
                (
                    "FC='/definitely/missing/gfortran -O2'; "
                    "CXX='sh -x'; "
                    f"{snippet}; "
                    'test -z "${FC+x}" && test "$CXX" = "sh -x"'
                ),
            ],
            check=False,
        )
        self.assertEqual(probe.returncode, 0)

    def test_darwin_mg5_make_opts_uses_resolved_gcc_toolchain(self) -> None:
        original = """DEFAULT_CPP_COMPILER=clang
DEFAULT_F_COMPILER=gfortran
MACFLAG=-mmacosx-version-min=10.8
STDLIB=-lc++
STDLIB_FLAG=-stdlib=libc++
#end_of_make_opts_variables
"""
        rewritten = campaign.rewrite_mg5_make_opts(
            original,
            cxx_compiler="/opt/homebrew/bin/g++-16",
            fortran_compiler="/opt/homebrew/bin/gfortran-16",
            system="Darwin",
        )
        self.assertIn(
            "DEFAULT_CPP_COMPILER=/opt/homebrew/bin/g++-16",
            rewritten,
        )
        self.assertIn(
            "DEFAULT_F_COMPILER=/opt/homebrew/bin/gfortran-16",
            rewritten,
        )
        self.assertIn("STDLIB=-lstdc++", rewritten)
        self.assertIn("STDLIB_FLAG=\n", rewritten)
        self.assertIn("MACFLAG=-mmacosx-version-min=10.8", rewritten)
        self.assertEqual(
            campaign.rewrite_mg5_make_opts(
                rewritten,
                cxx_compiler="/opt/homebrew/bin/g++-16",
                fortran_compiler="/opt/homebrew/bin/gfortran-16",
                system="Darwin",
            ),
            rewritten,
        )

    def test_linux_mg5_make_opts_preserves_platform_library_flags(self) -> None:
        original = """DEFAULT_CPP_COMPILER=clang++
DEFAULT_F_COMPILER=gfortran
STDLIB=-lc++
STDLIB_FLAG=-stdlib=libc++
"""
        rewritten = campaign.rewrite_mg5_make_opts(
            original,
            cxx_compiler="/usr/bin/g++-13",
            fortran_compiler="/usr/bin/gfortran-13",
            system="Linux",
        )
        self.assertIn("DEFAULT_CPP_COMPILER=/usr/bin/g++-13", rewritten)
        self.assertIn("DEFAULT_F_COMPILER=/usr/bin/gfortran-13", rewritten)
        self.assertIn("STDLIB=-lc++", rewritten)
        self.assertIn("STDLIB_FLAG=-stdlib=libc++", rewritten)

    def test_mg5_make_opts_missing_required_assignment_fails_closed(self) -> None:
        with self.assertRaisesRegex(ValueError, "DEFAULT_F_COMPILER"):
            campaign.rewrite_mg5_make_opts(
                "DEFAULT_CPP_COMPILER=clang\n",
                cxx_compiler="/usr/bin/g++",
                fortran_compiler="/usr/bin/gfortran",
                system="Linux",
            )

    def test_matching_versioned_gfortran_is_preferred(self) -> None:
        cfg = config()
        cfg.dry_run = False
        with TemporaryDirectory() as tmpdir:
            directory = Path(tmpdir)
            cxx = directory / "g++-16"
            fortran = directory / "gfortran-16"
            cxx.touch(mode=0o755)
            fortran.touch(mode=0o755)
            with (
                patch.object(campaign, "capture_runtime_stdout", return_value=None),
                patch.object(campaign, "_compiler_works", return_value=True),
            ):
                resolved = campaign.mg5_fortran_compiler(cfg, cxx)
        self.assertEqual(resolved, fortran)

    def test_generated_make_opts_patch_resolves_both_compilers(self) -> None:
        cfg = config()
        cfg.dry_run = False
        with TemporaryDirectory() as tmpdir:
            process_dir = Path(tmpdir) / "process"
            source_dir = process_dir / "Source"
            source_dir.mkdir(parents=True)
            make_opts = source_dir / "make_opts"
            make_opts.write_text(
                "DEFAULT_CPP_COMPILER=clang\n"
                "DEFAULT_F_COMPILER=gfortran\n"
                "STDLIB=-lc++\n"
                "STDLIB_FLAG=-stdlib=libc++\n"
            )
            with (
                patch.object(
                    campaign,
                    "herwig_plugin_compiler",
                    return_value=Path("/toolchain/g++-16"),
                ),
                patch.object(
                    campaign,
                    "mg5_fortran_compiler",
                    return_value=Path("/toolchain/gfortran-16"),
                ),
                patch.object(campaign.platform, "system", return_value="Darwin"),
            ):
                campaign.patch_mg5_toolchain(process_dir, cfg)
            rewritten = make_opts.read_text()
        self.assertIn("DEFAULT_CPP_COMPILER=/toolchain/g++-16", rewritten)
        self.assertIn("DEFAULT_F_COMPILER=/toolchain/gfortran-16", rewritten)
        self.assertIn("STDLIB=-lstdc++", rewritten)
        self.assertIn("STDLIB_FLAG=\n", rewritten)


class FourLeptonCardTests(unittest.TestCase):
    def test_madspin_card_imports_stable_lhe_and_uses_direct_four_body_decay(self) -> None:
        sample = campaign.SAMPLES_BY_NAME["signal_gg_h_2e2mu"]
        text = campaign.madspin_card_text(Path("/tmp/stable.lhe.gz"), sample)
        self.assertIn("import /tmp/stable.lhe.gz", text)
        self.assertIn("set spinmode none", text)
        self.assertIn("decay h > e+ e- mu+ mu-", text)
        self.assertLess(text.index("set spinmode none"), text.index("import "))
        self.assertLess(text.index("import "), text.index("set seed "))

    def test_madspin_uses_the_dedicated_frontend(self) -> None:
        cfg = config()
        command = campaign.madspin_command(cfg, Path("/tmp/card.dat"))
        self.assertEqual(
            Path(command[0]),
            cfg.mg5_dir / "MadSpin" / "madspin",
        )

    def test_madspin_input_is_staged_uncompressed(self) -> None:
        with TemporaryDirectory() as tmpdir:
            directory = Path(tmpdir)
            source = directory / "stable.lhe.gz"
            destination = directory / "decay" / "stable_input.lhe"
            destination.parent.mkdir()
            with gzip.open(source, "wt") as stream:
                stream.write("<LesHouchesEvents>payload</LesHouchesEvents>\n")
            cfg = config()
            cfg.dry_run = False
            campaign.stage_madspin_input(source, destination, cfg)
            self.assertEqual(
                destination.read_text(),
                "<LesHouchesEvents>payload</LesHouchesEvents>\n",
            )
            self.assertFalse(destination.read_bytes().startswith(b"\x1f\x8b"))
            campaign.stage_madspin_input(source, destination, cfg)
            with gzip.open(source, "wt") as stream:
                stream.write("<LesHouchesEvents>changed</LesHouchesEvents>\n")
            with self.assertRaisesRegex(ValueError, "differs"):
                campaign.stage_madspin_input(source, destination, cfg)

    def test_banner_light_external_lhe_gets_minimal_madspin_metadata(
        self,
    ) -> None:
        text = """<LesHouchesEvents version="3.0">
<init>
 2212 2212 2.0e4 2.0e4 0 0 0 0 -4 1
 1.0 0.0 1.0 1
</init>
</LesHouchesEvents>
"""
        with TemporaryDirectory() as tmpdir:
            directory = Path(tmpdir)
            source = directory / "source.lhe"
            destination = directory / "prepared.lhe"
            source.write_text(text)
            cfg = config()
            cfg.dry_run = False
            plan = campaign.prepare_madspin_banner_input(
                source,
                destination,
                cfg,
            )
            prepared = destination.read_text()
            self.assertTrue(plan["created_header"])
            self.assertTrue(plan["added_mg5_process_card"])
            self.assertTrue(plan["added_mg_run_card"])
            self.assertTrue(plan["added_slha"])
            self.assertFalse(plan["replaced_slha"])
            self.assertIn("<MG5ProcCard>", prepared)
            self.assertIn("<MGRunCard>", prepared)
            self.assertIn("3.0 = lhe_version", prepared)
            self.assertIn("generate p p > h", prepared)
            self.assertIn("<slha>", prepared)
            self.assertRegex(prepared, r"(?m)^\s*25\s+1\.250900000e\+02")
            self.assertRegex(
                prepared,
                r"(?m)^\s*24\s+8\.0419[0-9]*e\+01",
            )
            self.assertRegex(
                prepared,
                r"(?m)^\s*DECAY\s+25\s+[1-9][0-9.]*e[-+][0-9]+",
            )
            self.assertLess(prepared.index("<header>"), prepared.index("<init>"))
            self.assertEqual(
                campaign.prepare_madspin_banner_input(
                    source,
                    destination,
                    cfg,
                ),
                plan,
            )

    def test_zero_width_existing_slha_is_replaced_without_replacing_mg_banner(
        self,
    ) -> None:
        process_card = """<MG5ProcCard>
<![CDATA[
import model loop_sm
generate g g > h [noborn=QCD]
]]>
</MG5ProcCard>
"""
        text = f"""<LesHouchesEvents version="3.0">
<header>
{process_card}<slha>
BLOCK MASS
 25 1.250900000e+02 # MH
DECAY 25 0.000000000e+00
</slha>
<custom>preserve-me</custom>
</header>
<init>
 2212 2212 2.0e4 2.0e4 0 0 0 0 -4 1
 1.0 0.0 1.0 1
</init>
</LesHouchesEvents>
"""
        with TemporaryDirectory() as tmpdir:
            directory = Path(tmpdir)
            source = directory / "source.lhe"
            destination = directory / "prepared.lhe"
            source.write_text(text)
            cfg = config()
            cfg.dry_run = False
            plan = campaign.prepare_madspin_banner_input(
                source,
                destination,
                cfg,
            )
            prepared = destination.read_text()

        self.assertFalse(plan["added_mg5_process_card"])
        self.assertFalse(plan["added_slha"])
        self.assertTrue(plan["replaced_slha"])
        self.assertIn(
            "nonpositive_decay_25",
            plan["slha_validation"]["reasons"],
        )
        self.assertEqual(prepared.count("<MG5ProcCard>"), 1)
        self.assertIn(process_card, prepared)
        self.assertIn("<custom>preserve-me</custom>", prepared)
        self.assertEqual(prepared.lower().count("<slha>"), 1)
        self.assertNotIn("DECAY 25 0.000000000e+00", prepared)
        self.assertRegex(
            prepared,
            r"(?m)^\s*DECAY\s+25\s+[1-9][0-9.]*e[-+][0-9]+",
        )

    def test_higgs_only_existing_slha_is_replaced_as_incomplete(self) -> None:
        text = """<LesHouchesEvents version="3.0">
<header>
<MG5ProcCard>
<![CDATA[
import model loop_sm
generate g g > h [noborn=QCD]
]]>
</MG5ProcCard>
<slha>
BLOCK MASS
 25 1.250900000e+02 # MH
DECAY 25 4.070000000e-03
</slha>
</header>
<init>
 2212 2212 2.0e4 2.0e4 0 0 0 0 -4 1
 1.0 0.0 1.0 1
</init>
</LesHouchesEvents>
"""
        with TemporaryDirectory() as tmpdir:
            directory = Path(tmpdir)
            source = directory / "source.lhe"
            destination = directory / "prepared.lhe"
            source.write_text(text)
            cfg = config()
            cfg.dry_run = False
            plan = campaign.prepare_madspin_banner_input(
                source,
                destination,
                cfg,
            )
            prepared = destination.read_text()

        self.assertFalse(plan["added_mg5_process_card"])
        self.assertFalse(plan["added_slha"])
        self.assertTrue(plan["replaced_slha"])
        self.assertFalse(plan["preserved_slha"])
        self.assertFalse(plan["slha_validation"]["valid"])
        self.assertIn(
            "missing_or_ambiguous_mass_24",
            plan["slha_validation"]["reasons"],
        )
        self.assertRegex(prepared, r"(?m)^\s*24\s+8\.0419[0-9]*e\+01")

    def test_standalone_complete_existing_mg5_slha_is_preserved(self) -> None:
        cfg = config()
        slha = campaign.canonical_madspin_slha(cfg)
        text = f"""<LesHouchesEvents version="3.0">
<header>
<MG5ProcCard>
<![CDATA[
import model loop_sm
generate g g > h [noborn=QCD]
]]>
</MG5ProcCard>
<MGRunCard>
<![CDATA[
 2 = nevents
 3.0 = lhe_version
]]>
</MGRunCard>
<slha>
{slha}</slha>
</header>
<init>
 2212 2212 2.0e4 2.0e4 0 0 0 0 -4 1
 1.0 0.0 1.0 1
</init>
</LesHouchesEvents>
"""
        with TemporaryDirectory() as tmpdir:
            directory = Path(tmpdir)
            source = directory / "source.lhe"
            destination = directory / "prepared.lhe"
            source.write_text(text)
            cfg.dry_run = False
            plan = campaign.prepare_madspin_banner_input(
                source,
                destination,
                cfg,
            )
            prepared = destination.read_text()

        self.assertFalse(plan["added_mg5_process_card"])
        self.assertFalse(plan["added_mg_run_card"])
        self.assertFalse(plan["added_slha"])
        self.assertFalse(plan["replaced_slha"])
        self.assertTrue(plan["preserved_slha"])
        self.assertTrue(plan["slha_validation"]["valid"])
        self.assertEqual(prepared, text)

    def test_wrong_higgs_mass_in_existing_slha_is_repaired(self) -> None:
        text = """<LesHouchesEvents version="3.0">
<header>
<MG5ProcCard>
<![CDATA[
import model loop_sm
generate g g > h [noborn=QCD]
]]>
</MG5ProcCard>
<slha>
BLOCK MASS
 25 1.250000000e+02 # stale MH
DECAY 25 4.070000000e-03
</slha>
</header>
<init>
 2212 2212 2.0e4 2.0e4 0 0 0 0 -4 1
 1.0 0.0 1.0 1
</init>
</LesHouchesEvents>
"""
        with TemporaryDirectory() as tmpdir:
            directory = Path(tmpdir)
            source = directory / "source.lhe"
            destination = directory / "prepared.lhe"
            source.write_text(text)
            cfg = config()
            cfg.dry_run = False
            plan = campaign.prepare_madspin_banner_input(
                source,
                destination,
                cfg,
            )
            prepared = destination.read_text()

        self.assertTrue(plan["replaced_slha"])
        self.assertIn("invalid_mass_25", plan["slha_validation"]["reasons"])
        self.assertNotIn("1.250000000e+02 # stale MH", prepared)
        self.assertRegex(prepared, r"(?m)^\s*25\s+1\.250900000e\+02")

    def test_patch_run_card_applies_lo_pdf_beams_and_loose_lepton_cuts(self) -> None:
        template = """ 100 = nevents
 6500 = ebeam1
 6500 = ebeam2
 1 = iseed
 nn23lo1 = pdlabel
 230000 = lhaid
 10 = ptl
 2.5 = etal
 0 = mmll
 0.4 = drll
"""
        with TemporaryDirectory() as tmpdir:
            card = Path(tmpdir) / "run_card.dat"
            card.write_text(template)
            cfg = config()
            cfg.dry_run = False
            campaign.patch_run_card(
                card,
                sample=campaign.SAMPLES_BY_NAME["bkg_qq_4e"],
                nevents=33,
                seed=77,
                cfg=cfg,
            )
            updated = card.read_text()
        self.assertIn("33 = nevents", updated)
        self.assertIn("20000.0 = ebeam1", updated)
        self.assertIn("331900 = lhaid", updated)
        self.assertIn("5.0 = ptl", updated)
        self.assertIn("3.0 = etal", updated)
        self.assertIn("4.0 = mmll", updated)

    def test_stable_process_run_card_does_not_require_lepton_cut_fields(self) -> None:
        template = """ 100 = nevents
 6500 = ebeam1
 6500 = ebeam2
 1 = iseed
 nn23lo1 = pdlabel
 230000 = lhaid
"""
        with TemporaryDirectory() as tmpdir:
            card = Path(tmpdir) / "run_card.dat"
            card.write_text(template)
            cfg = config()
            cfg.dry_run = False
            campaign.patch_run_card(
                card,
                sample=campaign.SAMPLES_BY_NAME["signal_gg_h_4e"],
                nevents=10,
                seed=11,
                cfg=cfg,
            )
            updated = card.read_text()
        self.assertNotIn("ptl", updated)
        self.assertIn("331900 = lhaid", updated)

    def test_param_card_higgs_mass_is_patched_to_12509(self) -> None:
        with TemporaryDirectory() as tmpdir:
            card = Path(tmpdir) / "param_card.dat"
            card.write_text(
                "BLOCK MASS\n  25 1.250000e+02 # MH\nDECAY 25 6.382e-03\n"
            )
            cfg = config()
            cfg.dry_run = False
            width = campaign.patch_param_card_higgs_mass(card, cfg)
            updated = card.read_text()
        self.assertIn("1.250900000e+02", updated)
        self.assertAlmostEqual(width, 6.382e-3)

    def test_hwsim_template_saves_every_required_record(self) -> None:
        text = campaign.Config.template_in.__get__(config(), campaign.Config).read_text()
        for setting in (
            "SaveObjects Yes",
            "SavePartons Yes",
            "SaveReconstructed Yes",
            "SaveOptionalWeights Yes",
        ):
            self.assertIn(setting, text)


class FourLeptonNormalizationTests(unittest.TestCase):
    def test_lhe_optional_weight_definitions_preserve_semantic_order(self) -> None:
        text = """<LesHouchesEvents version="3.0">
<header>
<initrwgt>
 <weightgroup name="scale variation" combine="envelope">
  <weight id="1002">  mur = 2.0   muf = 1.0 </weight>
  <weight id='1001'> mur = 0.5 muf = 1.0 </weight>
 </weightgroup>
</initrwgt>
</header>
<init>
 2212 2212 2.0e4 2.0e4 0 0 331900 331900 -4 1
 1.0 0.0 1.0 1
</init>
</LesHouchesEvents>
"""
        with TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "weights.lhe"
            path.write_text(text)
            definitions = campaign.parse_lhe_optional_weight_definitions(path)
        self.assertEqual(
            [definition["id"] for definition in definitions],
            ["1002", "1001"],
        )
        self.assertEqual(
            [definition["label"] for definition in definitions],
            [
                "lhe_id=1002:mur = 2.0 muf = 1.0",
                "lhe_id=1001:mur = 0.5 muf = 1.0",
            ],
        )

    def test_stable_higgs_validation_uses_event_records_without_banner_slha(
        self,
    ) -> None:
        event = """<LesHouchesEvents version="3.0">
<event>
 3 1 1.0 125.09 0.0 0.0
 21 -1 0 0 501 502 0 0 100 100 0 0 9
 21 -1 0 0 502 501 0 0 -20 20 0 0 9
 25 1 1 2 0 0 0 0 80 148.477567 125.09 0 9
</event>
</LesHouchesEvents>
"""
        cfg = config()
        with TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "stable.lhe"
            path.write_text(event)
            mass, width = campaign.validate_stable_higgs_lhe(path, cfg, 1)
            self.assertAlmostEqual(mass, cfg.higgs_mass)
            self.assertIsNone(width)
            path.write_text(event.replace(" 25 1 1 2", " 25 2 1 2"))
            with self.assertRaisesRegex(ValueError, "status-1 Higgs"):
                campaign.validate_stable_higgs_lhe(path, cfg, 1)

    def test_stable_higgs_validation_uses_events_when_banner_mass_is_stale(
        self,
    ) -> None:
        event = """<LesHouchesEvents version="3.0">
<header>
<slha>
BLOCK MASS
 25 1.250000000e+02
DECAY 25 4.070000000e-03
</slha>
</header>
<event>
 3 1 1.0 125.09 0.0 0.0
 21 -1 0 0 501 502 0 0 100 100 0 0 9
 21 -1 0 0 502 501 0 0 -20 20 0 0 9
 25 1 1 2 0 0 0 0 80 148.477567 125.09 0 9
</event>
</LesHouchesEvents>
"""
        cfg = config()
        with TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "stable.lhe"
            path.write_text(event)
            mass, width = campaign.validate_stable_higgs_lhe(path, cfg, 1)
        self.assertAlmostEqual(mass, cfg.higgs_mass)
        self.assertAlmostEqual(width, 4.07e-3)

    def test_lhe_init_cross_section_sums_subprocesses(self) -> None:
        with TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "events.lhe.gz"
            with gzip.open(path, "wt") as stream:
                stream.write(lhe_text([2.0, 3.5]))
            self.assertAlmostEqual(campaign.parse_lhe_init_cross_section_pb(path), 5.5)

    def test_lhe_init_parsers_ignore_nested_generator_metadata(self) -> None:
        with TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "events.lhe"
            text = lhe_text([2.0, 3.5]).replace(
                "</init>",
                "<generator name='MadGraph5_aMC@NLO' version='3.5.12'>"
                "MadGraph5_aMC@NLO</generator>\n</init>",
            )
            path.write_text(text)
            self.assertAlmostEqual(campaign.parse_lhe_init_cross_section_pb(path), 5.5)
            metadata = campaign.parse_lhe_init_metadata(path)
            self.assertEqual(metadata["n_subprocesses"], 2)
            idwtup, maxima = campaign.parse_lhe_weight_metadata(path)
            self.assertEqual(idwtup, 3)
            self.assertEqual(maxima, {1: 1.0, 2: 1.0})

    def test_signal_scale_corrects_madspin_instead_of_reapplying_br(self) -> None:
        sample = campaign.SAMPLES_BY_NAME["signal_gg_h_4e"]
        cfg = config()
        with TemporaryDirectory() as tmpdir:
            stable = Path(tmpdir) / "stable.lhe"
            decayed = Path(tmpdir) / "decayed.lhe"
            stable.write_text(lhe_text([100.0]))
            decayed.write_text(lhe_text([0.010]))
            normalization = campaign.normalization_from_lhe(
                sample,
                decayed,
                cfg,
                stable_signal_lhe=stable,
            )
        target_br = campaign.BR_H_TO_4L * cfg.signal_channel_fractions["4e"]
        self.assertAlmostEqual(normalization["target_cross_section_pb"], 100.0 * target_br)
        self.assertAlmostEqual(normalization["weight_scale"], 100.0 * target_br / 0.010)
        self.assertNotAlmostEqual(normalization["weight_scale"], target_br)
        self.assertIn("madspin_decayed_xsec", normalization["formula"])

    def test_external_cross_section_overrides_placeholder_lhe_values(self) -> None:
        sample = campaign.SAMPLES_BY_NAME["signal_gg_h_4e"]
        cfg = config(
            "--source-backend",
            "external-lhe",
            "--external-cross-section",
            "signal_gg_h_stable=125.0",
            "--external-cross-section",
            "signal_gg_h_4e=0.02",
        )
        with TemporaryDirectory() as tmpdir:
            stable = Path(tmpdir) / "stable.lhe"
            decayed = Path(tmpdir) / "decayed.lhe"
            stable.write_text(lhe_text([-1.0]))
            decayed.write_text(lhe_text([-1.0]))
            normalization = campaign.normalization_from_lhe(
                sample,
                decayed,
                cfg,
                stable_signal_lhe=stable,
            )
        target_br = campaign.BR_H_TO_4L * cfg.signal_channel_fractions["4e"]
        self.assertAlmostEqual(
            normalization["weight_scale"],
            125.0 * target_br / 0.02,
        )
        self.assertEqual(
            normalization["stable_production_cross_section_source"],
            "cli_external_cross_section:signal_gg_h_stable",
        )

    def test_predecay_denominator_retains_variable_lhe_weight_magnitudes(self) -> None:
        with TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "events.lhe"
            path.write_text(lhe_text([1.0], [2.5, -8.0, 0.0, 3.0]))
            summary = campaign.predecay_lhe_weight_summary(path, 4)
        self.assertEqual(
            summary["denominator_kind"], "predecay_lhe_handler_weight_sum"
        )
        self.assertEqual(summary["effective_xmaxup_source"], "full_file_max_abs_xwgtup")
        self.assertEqual(summary["denominator_sumw"], -0.3125)
        self.assertEqual(summary["denominator_sumabsw"], 1.6875)
        self.assertEqual(summary["denominator_sumw2"], 1.23828125)
        self.assertEqual(summary["denominator_npositive"], 2)
        self.assertEqual(summary["denominator_nnegative"], 1)


class FourLeptonDecayBiasTests(unittest.TestCase):
    def _b_decay_file(self, directory: Path) -> Path:
        path = directory / "HerwigBDecays.in"
        lines = []
        charge = {
            "B0": ("e+", "nu_e"),
            "B+": ("mu+", "nu_mu"),
            "Bbar0": ("e-", "nu_ebar"),
            "B-": ("mu-", "nu_mubar"),
            "B_s0": ("e+", "nu_e"),
            "B_sbar0": ("e-", "nu_ebar"),
        }
        for parent, (lepton, neutrino) in charge.items():
            lines.append(
                f"decaymode {parent}->D0,{lepton},{neutrino}; 0.1 1 /Herwig/Decays/HQET\n"
            )
            lines.append(
                f"decaymode {parent}->pi+,pi-; 0.9 1 /Herwig/Decays/Mambo\n"
            )
        path.write_text("".join(lines))
        return path

    def test_semileptonic_full_support_bias_installs_exact_reweighter(self) -> None:
        with TemporaryDirectory() as tmpdir:
            cfg = config(
                "--sample-set",
                "reducible",
                "--heavy-flavour-bias",
                "semileptonic",
            )
            cfg.herwig_b_decays = self._b_decay_file(Path(tmpdir))
            nominal = campaign.SAMPLES_BY_NAME["bkg_zbb_ee"]
            sample = campaign.biased_variant(nominal)
            block, metadata = campaign.decay_selection_block(sample, cfg)
        self.assertIn(
            "HiggsSSC::LHEHeavyFlavorBranchingRatioReweighter", block
        )
        self.assertIn("LHEBranchingRatioReweighter.so", block)
        self.assertIn("PostDecayHandlers", block)
        self.assertNotIn("PostHadronizationHandlers", block)
        self.assertNotIn("B0:SelectDecayModes", block)
        self.assertIn(
            "B0/B0->D0,nu_e,e+;:BranchingRatio 0.40000000000000002",
            block,
        )
        self.assertLess(
            block.index("LHEBranchingRatioReweighter"),
            block.index("B0/B0->D0,nu_e,e+;:BranchingRatio"),
        )
        self.assertTrue(metadata["decay_bias"]["enabled"])
        self.assertFalse(metadata["decay_bias"]["closure_validated"])
        self.assertTrue(metadata["decay_bias"]["full_tail_support"])
        self.assertIn("q_m =", metadata["decay_bias"]["proposal"])
        self.assertIn("b_m=4", metadata["decay_bias"]["bias_function"])
        proposal = metadata["decay_bias"]["parent_proposals"]["B0"]
        self.assertAlmostEqual(proposal["original_active_branching_sum"], 1.0)
        self.assertAlmostEqual(proposal["selected_branching_sum"], 0.1)
        self.assertAlmostEqual(proposal["proposal_normalization"], 1.3)

    def test_semileptonic_importance_ratio_closes_mode_probabilities(self) -> None:
        original = {"semileptonic": 0.1, "other": 0.9}
        bias = {"semileptonic": 4.0, "other": 1.0}
        normalizer = sum(original[mode] * bias[mode] for mode in original)
        weighted_probability = 0.0
        for mode in original:
            proposal = original[mode] * bias[mode] / normalizer
            importance_weight = normalizer / bias[mode]
            weighted_probability += proposal * importance_weight
        self.assertAlmostEqual(weighted_probability, 1.0)

    def test_runner_and_plugin_use_the_same_semileptonic_bias(self) -> None:
        plugin_source = (
            REPO_ROOT
            / "hfourlepton"
            / "Herwig"
            / "LHEBranchingRatioReweighter.cc"
        ).read_text()
        match = re.search(
            r"kHeavyFlavorSemileptonicBias\s*=\s*([-+0-9.eE]+)",
            plugin_source,
        )
        self.assertIsNotNone(match)
        self.assertAlmostEqual(
            float(match.group(1)),
            campaign.HEAVY_FLAVOUR_SEMILEPTONIC_BIAS_FACTOR,
        )

    def test_all_samples_complete_one_cli_dry_run(self) -> None:
        with TemporaryDirectory() as tmpdir:
            b_decays = self._b_decay_file(Path(tmpdir))
            result = subprocess.run(
                [
                    sys.executable,
                    str(
                        REPO_ROOT
                        / "hfourlepton"
                        / "run_four_lepton_campaign.py"
                    ),
                    "--dry-run",
                    "--sample-set",
                    "all",
                    "--heavy-flavour-bias",
                    "semileptonic",
                    "--herwig-b-decays",
                    str(b_decays),
                    "--nevents",
                    "1",
                    "--bias-control-events",
                    "1",
                    "--run-tag",
                    "unit_all_cards",
                    "--no-herwig-module",
                ],
                cwd=REPO_ROOT,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.PIPE,
                text=True,
            )
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_prompt_and_heavy_flavour_handlers_use_separate_stages(self) -> None:
        with TemporaryDirectory() as tmpdir:
            cfg = config(
                "--sample-set",
                "reducible",
                "--heavy-flavour-bias",
                "semileptonic",
            )
            cfg.herwig_b_decays = self._b_decay_file(Path(tmpdir))
            sample = campaign.biased_variant(
                campaign.SAMPLES_BY_NAME["bkg_ttbar"]
            )
            block, metadata = campaign.decay_selection_block(sample, cfg)
        self.assertIn("HiggsSSC::LHEBranchingRatioReweighter", block)
        self.assertIn(
            "HiggsSSC::LHEHeavyFlavorBranchingRatioReweighter", block
        )
        self.assertIn("PostHadronizationHandlers", block)
        self.assertIn("PostDecayHandlers", block)
        self.assertIn(
            "HiggsSSC::LHEBranchingRatioReweighter",
            metadata["branching_ratio_reweighter_implementations"],
        )
        self.assertIn(
            "HiggsSSC::LHEHeavyFlavorBranchingRatioReweighter",
            metadata["branching_ratio_reweighter_implementations"],
        )
        self.assertEqual(
            metadata["handlers"],
            [
                {
                    "class": "HiggsSSC::LHEBranchingRatioReweighter",
                    "stage": "PostHadronizationHandlers",
                    "scope": "prompt_abs_pdg_23_24",
                },
                {
                    "class": (
                        "HiggsSSC::"
                        "LHEHeavyFlavorBranchingRatioReweighter"
                    ),
                    "stage": "PostDecayHandlers",
                    "scope": (
                        "terminal_ground_state_B_abs_pdg_511_521_531"
                    ),
                },
            ],
        )

    def test_prompt_only_handler_runs_before_hadron_decays(self) -> None:
        sample = campaign.SAMPLES_BY_NAME["bkg_ttbar"]
        block, metadata = campaign.decay_selection_block(sample, config())
        self.assertIn("HiggsSSC::LHEBranchingRatioReweighter", block)
        self.assertNotIn("LHEHeavyFlavorBranchingRatioReweighter", block)
        self.assertIn("PostHadronizationHandlers", block)
        self.assertNotIn("PostDecayHandlers", block)
        self.assertFalse(metadata["decay_bias"]["enabled"])

    def test_nominal_zbb_has_no_forced_heavy_flavour_decay(self) -> None:
        sample = campaign.SAMPLES_BY_NAME["bkg_zbb_ee"]
        block, metadata = campaign.decay_selection_block(sample, config())
        self.assertNotIn("SelectDecayModes", block)
        self.assertFalse(metadata["branching_ratio_reweighter"])

    def test_prompt_decay_requires_undecayed_status_one_parents(self) -> None:
        ttz = campaign.SAMPLES_BY_NAME["bkg_ttz_zee_wpe_wme"]
        stable_event = """<LesHouchesEvents version="3.0">
<event>
 3 1 1.0 0.0 0.0 0.0
 6 1 0 0 0 0 0 0 0 172.5 172.5 0 9
 -6 1 0 0 0 0 0 0 0 172.5 172.5 0 9
 23 1 0 0 0 0 0 0 0 91.2 91.2 0 9
</event>
</LesHouchesEvents>
"""
        predecayed_event = stable_event.replace(
            " 23 1 0 0 0 0",
            " 23 2 0 0 0 0",
        )
        with TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "events.lhe"
            path.write_text(stable_event)
            metadata = campaign.validate_prompt_decay_parents(path, ttz, 1)
            self.assertEqual(metadata["required_status1_pdgs"], [6, -6, 23])
            path.write_text(predecayed_event)
            with self.assertRaisesRegex(ValueError, "pre-decayed"):
                campaign.validate_prompt_decay_parents(path, ttz, 1)

    def test_forced_decay_plugin_is_global_and_staged_for_read_and_run(self) -> None:
        sample = campaign.SAMPLES_BY_NAME["bkg_ttz_zee_wpe_wme"]
        cfg = config()
        cfg.dry_run = False
        with TemporaryDirectory() as tmpdir:
            temporary = Path(tmpdir)
            plugin = temporary / "build" / "LHEBranchingRatioReweighter.so"
            plugin.parent.mkdir()
            plugin.write_bytes(b"test plugin")
            work_dir = temporary / "sample"
            lhe = temporary / "events.lhe"
            lhe.write_text("placeholder")
            with (
                patch.object(
                    campaign.Config,
                    "decay_reweighter_plugin",
                    new=property(lambda _self: plugin),
                ),
                patch.object(campaign, "sample_work_dir", return_value=work_dir),
                patch.object(campaign, "run_runtime_command") as runtime,
                patch.object(campaign, "write_root_input_list", return_value=[]),
            ):
                campaign.run_herwig(sample, lhe, cfg)

            herwig_dir = work_dir / "herwig"
            staged = herwig_dir / plugin.name
            card = (herwig_dir / f"{sample.name}.in").read_text()
            self.assertEqual(staged.read_bytes(), b"test plugin")
            self.assertIn(f"globallibrary {staged}", card)
            commands = [call.args[0] for call in runtime.call_args_list]
            self.assertEqual(commands[0][:2], [cfg.herwig, "read"])
            self.assertEqual(
                commands[1][:3],
                [cfg.herwig, "run", f"{sample.name}.run"],
            )


class FourLeptonManifestAndCliTests(unittest.TestCase):
    def test_manifest_is_versioned_and_exposes_report_normalization_contract(self) -> None:
        cfg = config("--run-samples", "signal_gg_h_4e")
        manifest = campaign.build_manifest(cfg, campaign.selected_samples(cfg))
        self.assertEqual(manifest["schema_version"], 1)
        self.assertEqual(manifest["collider"]["sqrt_s_tev"], 40.0)
        entry = manifest["samples"][0]
        self.assertIn("native_cross_section_pb", entry["normalization"])
        self.assertIn("target_cross_section_pb", entry["normalization"])
        self.assertEqual(
            entry["normalization"]["weight_scale"],
            campaign.PENDING_SIGNAL_WEIGHT,
        )
        self.assertEqual(
            entry["normalization"]["denominator_kind"],
            "predecay_lhe_handler_weight_sum",
        )
        json.dumps(manifest)

    def test_manifest_persists_decay_handler_class_and_stage(self) -> None:
        cfg = config(
            "--sample-set",
            "reducible",
            "--run-samples",
            "bkg_ttbar",
        )
        manifest = campaign.build_manifest(cfg, campaign.selected_samples(cfg))
        handlers = manifest["samples"][0]["decay"]["handlers"]
        self.assertEqual(len(handlers), 1)
        self.assertEqual(
            handlers[0]["class"],
            "HiggsSSC::LHEBranchingRatioReweighter",
        )
        self.assertEqual(
            handlers[0]["stage"],
            "PostHadronizationHandlers",
        )

    def test_configuration_fingerprint_covers_code_and_software(self) -> None:
        cfg = config("--run-samples", "bkg_qq_4e")
        samples = campaign.selected_samples(cfg)
        tracked_a = {"runner": {"sha256": "a"}}
        tracked_b = {"runner": {"sha256": "b"}}
        versions_a = {"python": "3", "mg5": "1", "herwig": "7", "root": "6"}
        versions_b = {**versions_a, "herwig": "8"}
        baseline = campaign.configuration_fingerprint(
            cfg,
            samples,
            tracked_inputs=tracked_a,
            software_versions=versions_a,
        )
        self.assertNotEqual(
            baseline,
            campaign.configuration_fingerprint(
                cfg,
                samples,
                tracked_inputs=tracked_b,
                software_versions=versions_a,
            ),
        )
        self.assertNotEqual(
            baseline,
            campaign.configuration_fingerprint(
                cfg,
                samples,
                tracked_inputs=tracked_a,
                software_versions=versions_b,
            ),
        )

    def test_partial_herwig_stage_requires_matching_identity_and_roots(self) -> None:
        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir) / "events.root"
            marker = Path(tmpdir) / "stage.json"
            root.write_bytes(b"root")
            identity = {"sample": "test", "card_sha256": "abc"}
            marker.write_text(
                json.dumps(
                    {
                        "identity": identity,
                        "roots": [
                            {
                                "path": str(root.resolve()),
                                "size": root.stat().st_size,
                                "mtime_ns": root.stat().st_mtime_ns,
                            }
                        ],
                    }
                )
            )
            self.assertTrue(
                campaign.herwig_stage_matches(marker, identity, [root])
            )
            self.assertFalse(
                campaign.herwig_stage_matches(
                    marker,
                    {**identity, "card_sha256": "changed"},
                    [root],
                )
            )
            root.write_bytes(b"changed")
            self.assertFalse(
                campaign.herwig_stage_matches(marker, identity, [root])
            )

    def test_analyzer_command_uses_only_named_campaign_arguments(self) -> None:
        cfg = config("--run-tag", "unit", "--detector-response", "ssc")
        sample = campaign.SAMPLES_BY_NAME["bkg_qq_4e"]
        command = campaign.analysis_command(
            sample,
            Path("/tmp/roots.input"),
            123,
            1.0,
            cfg,
        )
        for option in (
            "--input-list",
            "--response-profile",
            "--seed",
            "--weight-scale",
            "--tag",
            "--output-dir",
            "--sample",
            "--category",
            "--channel",
        ):
            self.assertIn(option, command)
        self.assertEqual(
            campaign.expected_analysis_output(sample, cfg).name,
            "bkg_qq_4e_unit_ssc.root",
        )

    def test_manifest_records_versions_checksums_and_detector_constants(self) -> None:
        cfg = config("--run-samples", "bkg_qq_4e", "--detector-response", "ssc")
        manifest = campaign.build_manifest(cfg, campaign.selected_samples(cfg))
        software = manifest["software"]
        self.assertEqual(software["mg5_version"], "3.5.15")
        self.assertEqual(software["analyzer_schema_version"], 1)
        self.assertIsNotNone(manifest["tracked_inputs"]["runner"]["sha256"])
        detector = manifest["detector"]
        self.assertEqual(detector["electron"]["reconstruction_efficiency"], 0.90)
        self.assertEqual(detector["muon"]["table_identity"], "GEM_TDR_Fig4-20_approx_digitization")
        self.assertEqual(detector["trigger_efficiency"]["4e"], 0.98)

    def test_manifest_does_not_misattribute_fraction_override(self) -> None:
        cfg = config(
            "--run-samples",
            "signal_gg_h_4e",
            "--signal-channel-fractions",
            "4e=.25,4mu=.25,2e2mu=.5",
        )
        manifest = campaign.build_manifest(cfg, campaign.selected_samples(cfg))
        self.assertEqual(
            manifest["higgs"]["channel_fraction_provenance"],
            "cli_override",
        )
        self.assertNotIn(
            "Reference LO direct-1to4",
            manifest["higgs"]["channel_fraction_note"],
        )

    def test_no_pileup_noise_is_forwarded_to_analyzer(self) -> None:
        cfg = config("--no-pileup-noise")
        sample = campaign.SAMPLES_BY_NAME["bkg_qq_4e"]
        command = campaign.analysis_command(
            sample,
            Path("/tmp/roots.input"),
            123,
            1.0,
            cfg,
        )
        self.assertIn("--no-pileup-noise", command)
        self.assertFalse(campaign.build_manifest(cfg, [sample])["detector"]["pileup_noise_enabled"])


if __name__ == "__main__":
    unittest.main()
