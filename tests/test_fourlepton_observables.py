from copy import deepcopy
import csv
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

import numpy as np

from hfourlepton import observables as obs


class ObservableTests(unittest.TestCase):
    spec = obs.Observable("m4l", "Mass [GeV]", "m4l", 0., 2., 2, "GeV")

    def row(self, source, mass, weight, channel=1):
        return dict(source_index=source, m4l=mass, event_weight=weight, channel=channel)

    def test_signed_weights_group_same_source_before_squaring(self):
        rows = [self.row(1, .2, .4), self.row(1, .4, .6), self.row(2, .3, -.25)]
        h = obs.histogram(rows, self.spec, 2.)
        np.testing.assert_allclose(h["cross_section_pb"], [0., 1.5, 0., 0.])
        np.testing.assert_allclose(h["variance_pb2"], [0., 4.25, 0., 0.])
        self.assertEqual(h["source_events"][1], 2)

    def test_shape_uncertainty_includes_cross_bin_correlation(self):
        rows = [self.row(1, .2, .25), self.row(1, 1.2, .75),
                self.row(2, .4, -.1), self.row(3, 1.5, .3)]
        h = obs.histogram(rows, self.spec, 1.)
        y, error = obs.density(h)
        source_bins = np.array([[.25, .75], [-.1, 0.], [0., .3]])
        w = source_bins.sum(axis=0)
        derivatives = (source_bins * w.sum() - w * source_bins.sum(axis=1)[:, None]) / w.sum()**2
        np.testing.assert_allclose(error**2, (derivatives**2).sum(axis=0))
        self.assertAlmostEqual(np.dot(y, np.diff(h["edges"])), 1.)
        one = obs.histogram(rows[:2], self.spec, 1.)
        np.testing.assert_allclose(obs.density(one)[1], 0., atol=1e-12)

    def test_flow_bins_and_signed_nonpositive_area(self):
        rows = [self.row(1, -1., 10.), self.row(2, 0., -2.),
                self.row(3, 1., 1.), self.row(4, 2., 20.)]
        h = obs.histogram(rows, self.spec, 1.)
        self.assertEqual(h["cross_section_pb"], [10., -2., 1., 20.])
        self.assertEqual(h["regular_total_pb"], -1.)
        self.assertIsNone(obs.density(h))

    def test_missing_data_are_unavailable_and_nonfinite_values_fail(self):
        self.assertEqual(obs.histogram([{}], self.spec, 1.)["status"], "branch_unavailable")
        self.assertEqual(obs.histogram([], self.spec, None)["status"], "normalization_unavailable")
        with self.assertRaises(ValueError):
            obs.histogram([self.row(1, float("nan"), 1.)], self.spec, 1.)

    def test_channel_correlations_and_signal_region_boundaries(self):
        rows = [self.row(1, 120., .1, 1), self.row(1, 125., .4, 1),
                self.row(1, 126., .6, 2), self.row(2, 130., .5, 3)]
        spec = obs.Observable("m4l", "Mass", "m4l", 100., 140., 1)
        with patch.object(obs, "OBSERVABLES", (spec,)):
            h = obs.sample_observables(rows, {"scale_pb_per_weight": 1.}, (120., 130.))
        inclusive = h["signal_region"]["inclusive"]["m4l"]
        self.assertEqual(inclusive["cross_section_pb"], [0., 1., 0.])
        self.assertEqual(inclusive["variance_pb2"], [0., 1., 0.])
        self.assertAlmostEqual(h["signal_region"]["4e"]["m4l"]["variance_pb2"][1], .16)
        self.assertAlmostEqual(h["signal_region"]["4mu"]["m4l"]["variance_pb2"][1], .36)

    def test_combines_independent_samples_and_exports_paired_plots(self):
        spec = obs.Observable("m4l", "Mass [GeV]", "m4l", 100., 140., 4, "GeV")
        rows = [self.row(1, 124., 1.), self.row(2, 125., -.2), self.row(3, 135., .1)]
        with patch.object(obs, "OBSERVABLES", (spec,)):
            values = obs.sample_observables(rows, {"scale_pb_per_weight": .001}, (120., 130.))
            samples = [dict(name=name, category=category, detector_profile="ssc", included_in_physics_totals=True,
                            observables=deepcopy(values)) for name, category in
                       (("signal_gg_h_4e", "signal"), ("bkg_qq_4e", "irreducible"))]
            combined = obs.combine([values["selected"]["inclusive"]["m4l"]] * 2)
            self.assertAlmostEqual(combined["regular_total_pb"], .0018)
            self.assertAlmostEqual(combined["regular_total_variance_pb2"], 2.1e-6)
            report = dict(metadata=dict(luminosity_fb=100., sqrt_s_tev=40., signal_region_gev=[120., 130.]),
                          samples=samples, audit={"errors": []})
            with TemporaryDirectory() as directory:
                result = obs.render(report, Path(directory))
                self.assertEqual(result["plots"], 2)
                for kind in ("png", "svg", "pdf"):
                    self.assertGreater((Path(directory) / "plots" / f"ssc_selected_m4l.{kind}").stat().st_size, 1000)
                with (Path(directory) / "data/m4l.csv").open() as stream:
                    records = list(csv.DictReader(stream))
                selected = [r for r in records if r["sample"] == "signal_gg_h_4e" and
                            r["region"] == "selected" and r["channel"] == "inclusive"]
                self.assertAlmostEqual(sum(float(r["events"]) for r in selected), 90.)
                self.assertIn("Unit-area shapes", (Path(directory) / "index.html").read_text())


if __name__ == "__main__":
    unittest.main()
