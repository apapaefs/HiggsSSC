import array
import json
import math
import os
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

import ROOT


REPO_ROOT = Path(__file__).resolve().parents[1]
CODE_DIR = REPO_ROOT / "hfourlepton" / "LOAnalysis" / "Code"
ANALYZER = CODE_DIR / "HwSimPostAnalysis_fourlepton"
FIXTURE = CODE_DIR / "make_synthetic_hw_sim_fourlepton"


def tree_rows(path: Path, tree_name: str, branches: tuple[str, ...]) -> list[dict]:
    root_file = ROOT.TFile.Open(str(path), "READ")
    if root_file is None or root_file.IsZombie():
        raise OSError(f"could not open {path}")
    tree = root_file.Get(tree_name)
    if tree is None:
        root_file.Close()
        raise KeyError(f"{path} has no {tree_name} tree")
    rows = []
    for entry in tree:
        row = {}
        for branch in branches:
            value = getattr(entry, branch)
            if type(value).__name__ == "string":
                value = str(value)
            elif hasattr(value, "__len__") and not isinstance(value, str):
                value = [
                    str(item) if type(item).__name__ == "string" else item
                    for item in value
                ]
            row[branch] = value
        rows.append(row)
    root_file.Close()
    return rows


class FourLeptonAnalyzerIntegrationTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        if shutil.which("root-config") is None:
            raise unittest.SkipTest("ROOT development tools are unavailable")
        cls._temporary = tempfile.TemporaryDirectory()
        cls.work = Path(cls._temporary.name)
        build_environment = dict(os.environ)
        build_environment["CXX"] = "/definitely/missing/herwig-exported-cxx -std=c++14"
        subprocess.run(
            ["make", "-C", str(CODE_DIR), "clean", "all", "fixture"],
            check=True,
            capture_output=True,
            text=True,
            env=build_environment,
        )
        cls.full_input = cls.work / "synthetic-full.root"
        cls.minimal_input = cls.work / "synthetic-minimal.root"
        cls.unnamed_weights_input = cls.work / "synthetic-unnamed-weights.root"
        subprocess.run([str(FIXTURE), str(cls.full_input)], check=True)
        subprocess.run(
            [str(FIXTURE), str(cls.minimal_input), "minimal"], check=True
        )
        subprocess.run(
            [str(FIXTURE), str(cls.unnamed_weights_input), "unnamed-weights"],
            check=True,
        )

    @classmethod
    def tearDownClass(cls) -> None:
        subprocess.run(
            ["make", "-C", str(CODE_DIR), "clean"],
            check=False,
            capture_output=True,
            text=True,
        )
        cls._temporary.cleanup()

    def run_analyzer(
        self,
        input_path: Path,
        *,
        profile: str,
        sample: str,
        seed: int = 12345,
        weight_scale: float = 1.0,
        extra_args: tuple[str, ...] = (),
    ) -> Path:
        command = [
            str(ANALYZER),
            str(input_path),
            "--response-profile",
            profile,
            "--seed",
            str(seed),
            "--weight-scale",
            str(weight_scale),
            "--output-dir",
            str(self.work),
            "--sample",
            sample,
            *extra_args,
        ]
        subprocess.run(
            command,
            check=True,
            capture_output=True,
            text=True,
        )
        return self.work / f"{sample}_{profile}.root"

    def test_perfect_cut_masks_named_schema_and_signed_weights(self) -> None:
        output = self.run_analyzer(
            self.full_input,
            profile="perfect",
            sample="perfect",
            weight_scale=2.0,
        )
        hypotheses = tree_rows(
            output,
            "FourLepton",
            (
                "source_index",
                "cut_mask",
                "channel",
                "candidate_count",
                "event_weight",
                "m4l",
                "mZ1",
                "mZ2",
                "lepton_source_index",
                "lepton_dressed_photon_count",
                "n_jets",
                "n_bjets",
            ),
        )
        self.assertEqual(
            [int(row["cut_mask"]) for row in hypotheses],
            [1023, 1023, 1023, 1, 3, 7, 15, 63],
        )
        self.assertEqual([int(row["channel"]) for row in hypotheses[:3]], [3, 1, 2])
        self.assertEqual(int(hypotheses[1]["candidate_count"]), 5)
        self.assertEqual(
            list(hypotheses[1]["lepton_source_index"]),
            [0, 1, 2, 3],
        )
        self.assertEqual(list(hypotheses[0]["lepton_dressed_photon_count"]), [1, 0, 0, 0])
        self.assertAlmostEqual(float(hypotheses[0]["event_weight"]), 2.0)
        self.assertAlmostEqual(float(hypotheses[1]["event_weight"]), 4.0)
        self.assertAlmostEqual(float(hypotheses[2]["event_weight"]), -2.0)
        self.assertGreater(float(hypotheses[0]["m4l"]), 120.0)
        self.assertLess(float(hypotheses[0]["m4l"]), 140.0)
        self.assertGreater(float(hypotheses[0]["mZ1"]), 90.0)
        self.assertEqual(int(hypotheses[0]["n_jets"]), 1)
        self.assertEqual(int(hypotheses[1]["n_bjets"]), 1)

        sources = tree_rows(
            output,
            "SourceEvents",
            (
                "source_index",
                "generator_weight",
                "base_event_weight",
                "response_failure_probability",
                "response_probability_in_tree",
                "selected_event_weight",
                "response_closure_delta",
                "optional_weight_names",
                "optional_selected_event_weights",
            ),
        )
        self.assertEqual(len(sources), 8)
        self.assertEqual(
            [float(row["selected_event_weight"]) for row in sources[:3]],
            [2.0, 4.0, -2.0],
        )
        for row in sources:
            self.assertAlmostEqual(
                float(row["response_failure_probability"])
                + float(row["response_probability_in_tree"]),
                1.0,
                places=12,
            )
            self.assertAlmostEqual(float(row["response_closure_delta"]), 0.0)
        self.assertEqual(list(sources[0]["optional_weight_names"]), ["scale_down", "scale_up"])
        self.assertEqual(
            [round(float(value), 12) for value in sources[0]["optional_selected_event_weights"]],
            [1.8, 2.2],
        )

    def test_ssc_probabilities_missing_optional_branches_and_repeatability(self) -> None:
        first = self.run_analyzer(
            self.minimal_input, profile="ssc", sample="ssc_first", seed=9981
        )
        second = self.run_analyzer(
            self.minimal_input, profile="ssc", sample="ssc_second", seed=9981
        )
        branches = (
            "source_index",
            "cut_mask",
            "hypothesis_probability",
            "response_weight",
            "event_weight",
            "lepton_pt",
        )
        first_rows = tree_rows(first, "FourLepton", branches)
        second_rows = tree_rows(second, "FourLepton", branches)
        self.assertEqual(first_rows, second_rows)

        sources = tree_rows(
            first,
            "SourceEvents",
            (
                "response_failure_probability",
                "response_probability_in_tree",
                "selected_probability",
                "response_closure_delta",
                "optional_weights",
            ),
        )
        for row in sources:
            self.assertAlmostEqual(
                float(row["response_failure_probability"])
                + float(row["response_probability_in_tree"]),
                1.0,
                places=12,
            )
            self.assertAlmostEqual(float(row["response_closure_delta"]), 0.0)
            self.assertEqual(list(row["optional_weights"]), [])

        four_lepton_reco_probability = 0.90**2 * (0.85 * 0.95) ** 2
        self.assertAlmostEqual(
            float(sources[0]["response_probability_in_tree"]),
            four_lepton_reco_probability,
            places=12,
        )
        self.assertAlmostEqual(
            float(sources[0]["response_failure_probability"]),
            1.0 - four_lepton_reco_probability,
            places=12,
        )
        if int(first_rows[0]["cut_mask"]) == 1023:
            self.assertAlmostEqual(
                float(sources[0]["selected_probability"]),
                four_lepton_reco_probability * 0.99,
                places=12,
            )

    def test_invalid_profile_is_rejected(self) -> None:
        result = subprocess.run(
            [
                str(ANALYZER),
                str(self.minimal_input),
                "--response-profile",
                "not-a-profile",
            ],
            capture_output=True,
            text=True,
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("expected perfect or ssc", result.stderr)

    def test_empty_hw_sim_weight_names_get_stable_positional_labels(self) -> None:
        output = self.run_analyzer(
            self.unnamed_weights_input,
            profile="perfect",
            sample="unnamed_weights",
        )
        sources = tree_rows(
            output,
            "SourceEvents",
            ("optional_weight_names", "optional_weights"),
        )
        self.assertEqual(
            list(sources[0]["optional_weight_names"]),
            ["optional_0", "optional_1"],
        )
        self.assertEqual(len(sources[0]["optional_weights"]), 2)

    def test_optional_weight_names_file_recovers_semantics_and_metadata(self) -> None:
        names_file = self.work / "optional-weight-names.txt"
        names_file.write_text("# ordered HwSim weights\nscale_down\nscale_up\n")
        output = self.run_analyzer(
            self.unnamed_weights_input,
            profile="ssc",
            sample="semantic_weights",
            extra_args=("--optional-weight-names-file", str(names_file)),
        )
        sources = tree_rows(
            output,
            "SourceEvents",
            ("optional_weight_names", "optional_weights"),
        )
        self.assertEqual(
            list(sources[0]["optional_weight_names"]),
            ["scale_down", "scale_up"],
        )
        metadata = tree_rows(
            output,
            "AnalysisMetadata",
            (
                "optional_weight_names_source",
                "optional_weight_names_file",
                "detector_parameter_source",
                "cut_bit_zero_semantics",
                "electron_efficiency",
                "muon_efficiency",
                "trigger_4e_efficiency",
                "trigger_other_efficiency",
                "em_pileup_noise_enabled",
            ),
        )[0]
        self.assertEqual(metadata["optional_weight_names_source"], "configured_file")
        self.assertEqual(metadata["optional_weight_names_file"], str(names_file))
        self.assertIn("GEM_TDR", metadata["detector_parameter_source"])
        self.assertIn("not the baseline-accepted", metadata["cut_bit_zero_semantics"])
        self.assertAlmostEqual(float(metadata["electron_efficiency"]), 0.90)
        self.assertAlmostEqual(float(metadata["muon_efficiency"]), 0.85 * 0.95)
        self.assertAlmostEqual(float(metadata["trigger_4e_efficiency"]), 0.98)
        self.assertAlmostEqual(float(metadata["trigger_other_efficiency"]), 0.99)
        self.assertTrue(bool(metadata["em_pileup_noise_enabled"]))

        summary = json.loads(output.with_suffix(".summary.json").read_text())
        self.assertEqual(
            summary["diagnostics"]["optional_weight_names_from_file"], 8
        )
        self.assertEqual(
            summary["diagnostics"]["optional_weight_name_syntheses"], 0
        )

    def test_optional_weight_names_file_size_mismatch_is_rejected(self) -> None:
        names_file = self.work / "wrong-optional-weight-names.txt"
        names_file.write_text("only_one_name\n")
        result = subprocess.run(
            [
                str(ANALYZER),
                str(self.unnamed_weights_input),
                "--response-profile",
                "perfect",
                "--optional-weight-names-file",
                str(names_file),
                "--output-dir",
                str(self.work),
                "--sample",
                "wrong_semantic_weights",
            ],
            capture_output=True,
            text=True,
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("contains 1 names", result.stderr)
        self.assertIn("contains 2 weights", result.stderr)

    def test_selected_zero_weights_are_counted(self) -> None:
        output = self.run_analyzer(
            self.full_input,
            profile="perfect",
            sample="zero_selected_weights",
            weight_scale=0.0,
        )
        summary = json.loads(output.with_suffix(".summary.json").read_text())
        selected = summary["selected_scaled"]
        self.assertEqual(selected["entries"], 3)
        self.assertEqual(selected["zero_entries"], 3)
        self.assertEqual(selected["sumw"], 0)

    def test_more_than_ten_baseline_accepted_leptons_is_rejected(self) -> None:
        overflow_input = self.work / "accepted-lepton-overflow.root"
        root_file = ROOT.TFile(str(overflow_input), "RECREATE")
        tree = ROOT.TTree("Data", "accepted lepton overflow")
        numparticles = array.array("i", [11])
        objects = array.array("d", [0.0] * (8 * 10000))
        evweight = array.array("d", [1.0])
        tree.Branch("numparticles", numparticles, "numparticles/I")
        tree.Branch("objects", objects, "objects[8][10000]/D")
        tree.Branch("evweight", evweight, "evweight/D")
        coordinates = [
            (-2.0, 0.0),
            (-2.0, 1.5),
            (-2.0, 3.0),
            (-2.0, 4.5),
            (0.0, 0.0),
            (0.0, 1.5),
            (0.0, 3.0),
            (0.0, 4.5),
            (2.0, 0.0),
            (2.0, 1.5),
            (2.0, 3.0),
        ]
        for index, (eta, phi) in enumerate(coordinates):
            pt = 20.0
            px = pt * math.cos(phi)
            py = pt * math.sin(phi)
            pz = pt * math.sinh(eta)
            energy = math.sqrt(px * px + py * py + pz * pz)
            objects[index] = energy
            objects[10000 + index] = px
            objects[20000 + index] = py
            objects[30000 + index] = pz
            objects[40000 + index] = 11.0 if index % 2 == 0 else -11.0
        tree.Fill()
        tree.Write()
        root_file.Close()

        result = subprocess.run(
            [
                str(ANALYZER),
                str(overflow_input),
                "--response-profile",
                "perfect",
                "--output-dir",
                str(self.work),
                "--sample",
                "accepted_lepton_overflow",
            ],
            capture_output=True,
            text=True,
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("11 baseline-accepted leptons", result.stderr)
        self.assertIn("supports at most 10", result.stderr)

    def test_missing_required_branch_is_rejected(self) -> None:
        malformed = self.work / "missing-evweight.root"
        root_file = ROOT.TFile(str(malformed), "RECREATE")
        tree = ROOT.TTree("Data", "malformed")
        placeholder = array.array("d", [1.0])
        tree.Branch("placeholder", placeholder, "placeholder/D")
        tree.Fill()
        tree.Write()
        root_file.Close()
        result = subprocess.run(
            [str(ANALYZER), str(malformed)],
            capture_output=True,
            text=True,
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("missing required HwSim branch", result.stderr)


if __name__ == "__main__":
    unittest.main()
