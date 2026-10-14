import argparse
import io
from pathlib import Path
import tempfile
import unittest
from unittest import mock

import chira_utilities as utilities


class CigarTests(unittest.TestCase):
    def test_parse_all_sam_operations(self):
        self.assertEqual(utilities.parse_cigar('3H5S10M2I3D4N1P6=2X7S2H'),
                         [(3, 'H'), (5, 'S'), (10, 'M'), (2, 'I'), (3, 'D'),
                          (4, 'N'), (1, 'P'), (6, '='), (2, 'X'), (7, 'S'), (2, 'H')])

    def test_reject_malformed_cigars(self):
        for cigar in ('', '*', '0M', '01M', '10', 'M', '10Z', '10Mgarbage',
                      'junk10M', '10M-2I', '10M2S10M', '5S3H10M', '10M3H5S'):
            with self.subTest(cigar=cigar):
                with self.assertRaises(ValueError):
                    utilities.parse_cigar(cigar)

    def test_cigar_lengths_have_distinct_contracts(self):
        cigar = '3H5S10M2I3D4N1P6=2X7S2H'
        self.assertEqual(utilities.query_length(cigar, False), 37)
        self.assertEqual(utilities.query_length(cigar, True), 37)
        self.assertEqual(utilities.aligned_query_length(cigar), 20)
        self.assertEqual(utilities.reference_length(cigar), 25)
        self.assertEqual(utilities.alignment_length(cigar), 23)

    def test_match_positions(self):
        cases = [
            ('10M', (1, 10), (1, 10)),
            ('5S10M3S', (6, 15), (4, 13)),
            ('5S10=2X3S', (6, 17), (4, 15)),
            ('3H5S10M7S2H', (9, 18), (10, 19)),
            ('3H5S10M2I3D4N1P6=2X7S2H', (9, 28), (10, 29)),
            ('2I5M3I', (1, 10), (1, 10)),
            ('1M', (1, 1), (1, 1)),
        ]
        for cigar, forward, reverse in cases:
            with self.subTest(cigar=cigar):
                self.assertEqual(utilities.match_positions(cigar, False), forward)
                self.assertEqual(utilities.match_positions(cigar, True), reverse)

    def test_match_positions_require_reference_aligned_query_bases(self):
        for cigar in ('10S', '3H5S', '10I', '5D', '2S3I2S'):
            with self.subTest(cigar=cigar):
                with self.assertRaisesRegex(ValueError, 'no aligned query bases'):
                    utilities.match_positions(cigar, False)

    def test_reference_end_matches_bed_exclusive_end(self):
        # SAM start 101 -> BED start 100; reference span is 25 bases.
        for reverse in (False, True):
            self.assertEqual(utilities.alignment_end('101', '5S10M2I3D4N6=2X', reverse), 125)

    def test_chimera_overlap_counts_shared_endpoint(self):
        # Original-read arms are [1, 10] and [10, 20]: one shared base.
        self.assertFalse(utilities.is_chimeric('10M10S', '9S11M', False, False, 0))
        self.assertTrue(utilities.is_chimeric('10M10S', '9S11M', False, False, 1))
        self.assertTrue(utilities.is_chimeric('10M10S', '10S10M', False, False, 0))
        self.assertFalse(utilities.is_chimeric('10M10S', '11M9S', False, True, 0))


class UtilityTests(unittest.TestCase):
    def test_interval_conventions(self):
        self.assertEqual(utilities.overlap((1, 10), (10, 20)), 0)
        self.assertEqual(utilities.inclusive_overlap((1, 10), (10, 20)), 1)
        self.assertEqual(utilities.inclusive_overlap((5, 5), (5, 5)), 1)
        self.assertEqual(utilities.inclusive_overlap((1, 5), (6, 10)), 0)
        self.assertEqual(utilities.overlap((0, 10), (3, 7)), 4)
        self.assertEqual(utilities.inclusive_overlap((1, 10), (3, 7)), 5)

    def test_score_float(self):
        for value in ('0', '0.5', '1'):
            self.assertEqual(utilities.score_float(value), float(value))
        for value in ('-0.1', '1.1', 'nan', 'inf', '-inf'):
            with self.subTest(value=value):
                with self.assertRaises(argparse.ArgumentTypeError):
                    utilities.score_float(value)
        parser = argparse.ArgumentParser()
        parser.add_argument('--score', type=utilities.score_float)
        with mock.patch('sys.stderr', new=io.StringIO()):
            with self.assertRaises(SystemExit):
                parser.parse_args(['--score', 'not-a-number'])

    def test_median_sorted_input(self):
        self.assertEqual(utilities.median([1]), 1)
        self.assertEqual(utilities.median([1, 3, 9]), 3)
        self.assertEqual(utilities.median([1, 3, 5, 9]), 4)
        with self.assertRaisesRegex(ValueError, 'empty sequence'):
            utilities.median([])

    def test_extract_reference_lengths(self):
        with tempfile.TemporaryDirectory() as directory:
            fasta = Path(directory) / 'reference.fa'
            fasta.write_text('>tx1 description\nAC\nGT\n>tx2\nA\n')
            lengths = {'existing': 5}
            self.assertIsNone(utilities.extract_reflengths(fasta, lengths))
            self.assertEqual(lengths, {'existing': 5, 'tx1': 4, 'tx2': 1})


if __name__ == '__main__':
    unittest.main()
