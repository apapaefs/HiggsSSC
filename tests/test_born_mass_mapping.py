import unittest
import numpy as np
from hgammagamma.born_mass_mapping import match_subsequence, region


class BornMassMappingTest(unittest.TestCase):
    def setUp(self):
        self.lhe = np.arange(80, dtype=float).reshape(10, 2, 4) + 1.

    def test_leading_internal_and_trailing_skips(self):
        indices = np.asarray([1, 2, 4, 5, 6, 8])
        raw = self.lhe[indices].copy()
        raw[2] = raw[2, ::-1]
        raw[3, 0, 0] += 5.e-8
        mapped, errors = match_subsequence(self.lhe, raw, 4)
        np.testing.assert_array_equal(mapped, indices)
        self.assertAlmostEqual(max(errors), 5.e-8)

    def test_ambiguous_duplicate_and_wrong_skip_count_fail(self):
        lhs = np.stack([self.lhe[0], self.lhe[0], self.lhe[1]])
        with self.assertRaisesRegex(ValueError, 'possible LHE matches'):
            match_subsequence(lhs, lhs[[0,2]], 1)
        with self.assertRaisesRegex(ValueError, 'population'):
            match_subsequence(self.lhe, self.lhe[:-1], 0)

    def test_reordered_or_recoiled_particles_fail(self):
        with self.assertRaisesRegex(ValueError, 'no row join is safe'):
            match_subsequence(self.lhe, self.lhe[::-1], 0)
        raw = self.lhe.copy();raw[4, 0, 1] += .01
        with self.assertRaisesRegex(ValueError, 'no row join is safe'):
            match_subsequence(self.lhe, raw, 0)

    def test_mass_partition_refuses_boundary_ambiguity(self):
        np.testing.assert_array_equal(region(np.asarray([99.,101.,149.,151.]),100.,150.,1.e-4),[False,True,True,False])
        with self.assertRaisesRegex(ValueError, 'serialization guard'):
            region(np.asarray([100.00001]),100.,150.,1.e-4)


if __name__ == '__main__':
    unittest.main()
