from __future__ import absolute_import
from ..PathogenStructure import Strain
import unittest


class TestStrain(unittest.TestCase):

    def test_to_sec(self):

        string_sec = {
                '22sec': 22.,
                '22.7sec': 22.7,
                '1min': 60.,
                '1.6min': 1.6*60,
                '1.6h': 1.6*60*60,
                '3d': 3*24*60*60,
                '3.1d': 3.1*24*60*60,
                }

        for s in string_sec:
            self.assertTrue(string_sec[s] == Strain.to_sec(s))


if __name__ == '__main__':
    unittest.main()
