import math
from pathlib import Path
import tempfile
import unittest

import numpy as np
from scipy.stats import beta

from hgammagamma.combine_loop_gammagamma import (
    combine_statistics, embedded_run_card, herwig_rng_seed, normalized_herwig_card, parse_run_card, source_statistics,
)


class DisjointLoopCombinationTest(unittest.TestCase):
    def test_grouped_response_variance_and_full_source_denominator(self):
        # Two selected hypotheses from source0 must be squared together. The
        # fourth hard event has no saved response and still enters the exposure.
        row = source_statistics([True, True, False, True, False, False],
                                [.2, .3, .5, .4, .6, 1.], [0, 0, 0, 1, 1, 2],
                                attempted=4, cross_section=20.)
        self.assertAlmostEqual(row["selected_cross_section_pb"], 4.5)
        self.assertAlmostEqual(row["sum_squared_source_response"], .41)
        self.assertAlmostEqual(row["selected_cross_section_mc_error_pb"], 5*math.sqrt(4/3*(.41-.9**2/4)))
        self.assertEqual(row["selected_positive_sources"], 2)
        self.assertEqual(row["normalization_denominator"], 4)
        self.assertEqual(row["saved_sources"], 3)

    def test_zero_bound_and_independent_sum_have_simultaneous_coverage(self):
        a = source_statistics([False], [1.], [0], attempted=100, cross_section=20.)
        b = source_statistics([True], [1.], [0], attempted=50, cross_section=2.)
        self.assertAlmostEqual(a["positive_source_upper_bound"]["positive_response_probability_upper"], 1-.025**.01)
        combined = combine_statistics(a, b)
        expected = 20*(1-.025**.01) + 2*beta.ppf(.975, 2, 49)
        self.assertAlmostEqual(combined["selected_cross_section_upper_95_pb"], expected)
        self.assertAlmostEqual(combined["selected_cross_section_mc_error_pb"], b["selected_cross_section_mc_error_pb"])
        a["positive_source_upper_bound"]["confidence_level"] = .95
        with self.assertRaisesRegex(ValueError, "97.5%"):
            combine_statistics(a, b)

    def test_fractional_hypotheses_are_not_independent_binomial_trials(self):
        row = source_statistics([True, True, False], [.25, .25, .5], [0, 0, 0],
                                attempted=10, cross_section=100.)
        self.assertEqual(row["selected_positive_sources"], 1)
        self.assertEqual(row["selected_response_rows"], 2)
        self.assertAlmostEqual(row["positive_source_upper_bound"]["selected_cross_section_upper_pb"], 100*beta.ppf(.975, 2, 9))
        self.assertAlmostEqual(row["selected_cross_section_pb"], 5.)

    def test_nonunit_or_negative_response_fails(self):
        for weights in ([.7], [-1.]):
            with self.assertRaises(ValueError):
                source_statistics([True], weights, [0], attempted=10, cross_section=1.)

    def test_card_normalization_removes_only_nonphysical_controls(self):
        a = "set Reader:FileName /a\nset EventGenerator:NumberOfEvents 100\nset EventGenerator:RandomNumberGenerator:Seed 1\nset Shower:MaxPtIsMuF Yes\n"
        b = "set Reader:FileName /b\nset EventGenerator:NumberOfEvents 500\nset EventGenerator:RandomNumberGenerator:Seed 2\nset EventGenerator:DumpPeriod -1\nset Shower:MaxPtIsMuF Yes # comment\n"
        self.assertEqual(normalized_herwig_card(a), normalized_herwig_card(b))
        self.assertEqual(herwig_rng_seed(a), 1)
        self.assertEqual(herwig_rng_seed(b), 2)
        with self.assertRaisesRegex(ValueError, "exactly one"):
            herwig_rng_seed(a+a)
        self.assertNotEqual(normalized_herwig_card(a), normalized_herwig_card(b.replace("MaxPtIsMuF Yes", "MaxPtIsMuF No")))
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/"run_card.dat"
            path.write_text("100 = nevents ! count\n10.0 = pta # cut\n'lhapdf' = pdlabel\n")
            self.assertEqual(parse_run_card(path), {"nevents": "100", "pta": "10.0", "pdlabel": "'lhapdf'"})

    def test_generator_seed_is_proven_from_embedded_source_card(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/"events.lhe"
            path.write_text('<LesHouchesEvents>\n<header>\n<MGRunCard>\n181126033 = iseed ! actual input seed\n100 = nevents\n</MGRunCard>\n</header>\n')
            card, digest = embedded_run_card(path)
            self.assertEqual(card["iseed"], "181126033")
            self.assertEqual(len(digest), 64)
            path.write_text('<LesHouchesEvents>\n<event>\n')
            with self.assertRaisesRegex(ValueError, "embedded MG"):
                embedded_run_card(path)


if __name__ == "__main__":
    unittest.main()
