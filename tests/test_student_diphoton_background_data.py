"""Physical integrity checks for the small event tables shipped to students."""
import hashlib
import json
from pathlib import Path
import unittest

import numpy as np

from read_root_varfiles import FEATURE_NAMES


class StudentDiphotonBackgroundDataTests(unittest.TestCase):
    def test_release_population_and_mass_window(self):
        directory = Path(__file__).resolve().parents[1] / "hgammagamma/HOAnalysis/student_backgrounds/ssc40_v1"
        manifest = json.loads((directory / "manifest.json").read_text())
        self.assertEqual(manifest["feature_names"], FEATURE_NAMES)
        self.assertEqual(manifest["energy_tev"], 40)
        reference = {"bkg_loop_gg_aa": (25000, 0.7390980753408),
                     "bkg_dijet_double_fake": (1320000, 3.01056023847953)}
        self.assertEqual(set(manifest["samples"]), set(reference))
        for name, sample in manifest["samples"].items():
            with self.subTest(name=name):
                path = directory / sample["file"]
                self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(), sample["sha256"])
                with np.load(path, allow_pickle=False) as data:
                    x, source, weights = (data[key] for key in ("variables", "sourceevent", "weight_pb"))
                    self.assertEqual(x.shape, (sample["rows"], len(FEATURE_NAMES)))
                    self.assertEqual(len(np.unique(source)), sample["represented_sources"])
                    self.assertEqual(sample["generated_sources"], reference[name][0])
                    self.assertTrue(np.all((source >= 0) & (source < sample["generated_sources"])))
                    self.assertTrue(np.all(np.isfinite(x)))
                    self.assertTrue(np.all(np.isfinite(weights) & (weights >= 0)))
                    self.assertAlmostEqual(weights.sum(), sample["sum_weight_pb"])
                    window = (x[:, 9] == 2) & (x[:, 0] >= 120) & (x[:, 0] <= 130)
                    self.assertAlmostEqual(weights[window].sum(), reference[name][1])
                    if name == "bkg_dijet_double_fake":
                        self.assertTrue(np.all(x[:, 9] == 2))
                        lower = data["lower_weight_pb"]
                        self.assertTrue(np.all((lower >= 0) & (lower <= weights)))
                        self.assertAlmostEqual(lower[window].sum(), 0.113334045152005)
                    else:
                        self.assertEqual(set(x[:, 9]), {0., 1., 2.})


if __name__ == "__main__":
    unittest.main()
