import contextlib
import io
import math
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest import mock

import chira_quantify as quantify
import chira_merge as merging
import chira_utilities as utilities


def crl_row(segment, group, start, end, locus=None):
    locus = group if locus is None else locus
    return '\t'.join([segment, 'tx', str(locus), str(group), str(start), str(end), '+',
                      '{}M'.format(end - start), 'chr:1:{}:{}:+'.format(start, end),
                      'chr:1:{}:{}:+'.format(start, end), '1.00']) + '\n'


class GroupingTests(unittest.TestCase):
    def test_first_crl_is_considered_including_exact_similarity(self):
        for cutoff in (0, 0.7, 1):
            with self.subTest(cutoff=cutoff):
                groups = quantify.group_loci([{'a', 'b'}, {'a', 'b'}], True, cutoff, 1)
                self.assertEqual(len(groups), 1)
                self.assertEqual(set(groups[0]), {0, 1})

    def test_cutoff_one_keeps_nonidentical_read_sets_separate(self):
        groups = quantify.group_loci([{'a', 'b'}, {'a'}], True, 1, 1)
        self.assertEqual(len(groups), 2)

    def test_no_shared_reads_never_group_even_at_zero_cutoff(self):
        self.assertEqual(len(quantify.group_loci([{'a'}, {'b'}], True, 0, 1)), 2)

    def test_below_size_threshold_and_disabled_grouping(self):
        for enabled, minimum in ((True, 3), (False, 1)):
            self.assertEqual(len(quantify.group_loci([{'a', 'b'}, {'a', 'b'}], enabled, 0.7, minimum)), 2)

    def test_nonmatching_recent_crl_does_not_hide_earlier_match(self):
        groups = quantify.group_loci([set('abcdef'), set('xyzuv'), set('abcd')], True, 0.6, 1)
        self.assertEqual([set(group) for group in groups], [{0, 2}, {1}])

    def test_isolated_loci_have_one_membership_after_share_refinement(self):
        reads = [set('abcd'), set('abce'), set('abcf')]
        groups = quantify.group_loci(reads, True, 0.6, 1)
        self.assertEqual(set().union(*(set(group) for group in groups)), {0, 1, 2})
        for group in groups:
            union = set().union(*group.values())
            for locus_reads in group.values():
                self.assertGreaterEqual(len(locus_reads) / len(union), 0.6)

    def test_equal_best_memberships_are_retained(self):
        groups = quantify.group_loci([set('abc'), set('abd'), set('ab')], True, 0.6, 1)
        self.assertEqual([set(group) for group in groups], [{0, 2}, {1, 2}])

    def test_invalid_grouping_parameters(self):
        for cutoff in (-0.1, 1.1, math.nan):
            with self.assertRaises(ValueError):
                quantify.group_loci([{'a'}], True, cutoff, 1)
        with self.assertRaises(ValueError):
            quantify.group_loci([set()], True, 0.7, 1)


class NumericTests(unittest.TestCase):
    def quantify(self, rows, threshold=1e-10):
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / 'loci.txt'
            source.write_text(''.join(rows))
            return quantify.quantify_crls(str(source), threshold)

    def test_bed_lengths_determine_tpm(self):
        alpha, tpms = self.quantify([crl_row('read1|1', 0, 0, 1), crl_row('read2|1', 1, 0, 2)])
        self.assertEqual(alpha, {'read1|1': {'0': 1}, 'read2|1': {'1': 1}})
        self.assertAlmostEqual(tpms['0'], 2000000 / 3)
        self.assertAlmostEqual(tpms['1'], 1000000 / 3)
        self.assertAlmostEqual(sum(tpms.values()), 1000000)

    def test_symmetric_multimapping_stays_equal(self):
        alpha, tpms = self.quantify([crl_row('read1|1', 0, 0, 10), crl_row('read1|1', 1, 0, 10)])
        self.assertEqual(alpha['read1|1'], {'0': 0.5, '1': 0.5})
        self.assertEqual(tpms, {'0': 500000, '1': 500000})

    def test_unique_evidence_drives_multimap_probabilities(self):
        rows = [crl_row('unique1', 0, 0, 10), crl_row('unique2', 0, 0, 10),
                crl_row('unique3', 1, 0, 10), crl_row('multi', 0, 0, 10), crl_row('multi', 1, 0, 10)]
        alpha, tpms = self.quantify(rows)
        self.assertAlmostEqual(alpha['multi']['0'], 2 / 3, places=8)
        self.assertAlmostEqual(alpha['multi']['1'], 1 / 3, places=8)
        self.assertAlmostEqual(tpms['0'], 2000000 / 3, places=3)
        self.assertAlmostEqual(sum(sum(m.values()) for m in alpha.values()), 4)
        for memberships in alpha.values():
            self.assertAlmostEqual(sum(memberships.values()), 1)

    def test_duplicate_rows_do_not_increase_segment_counts(self):
        row = crl_row('read|UMI|100|1', 0, 0, 10)
        alpha, tpms = self.quantify([row, row])
        self.assertEqual(len(alpha), 1)
        self.assertEqual(tpms, {'0': 1000000})

    def test_tpm_uses_median_length_and_handles_zero_expression(self):
        self.assertEqual(quantify.tpm({'a': 1, 'b': 1}, {'a': {'x': 10, 'y': 30}, 'b': {'z': 20}}),
                         {'a': 500000, 'b': 500000})
        self.assertEqual(quantify.tpm({'a': 0}, {'a': {'x': 10}}), {'a': 0})
        self.assertEqual(quantify.tpm({}, {}), {})
        with self.assertRaises(ValueError):
            quantify.tpm({'a': 1}, {'a': {'x': 0}})

    def test_empty_input(self):
        self.assertEqual(self.quantify([]), ({}, {}))

    def test_invalid_locus_length(self):
        with self.assertRaisesRegex(ValueError, 'Invalid BED'):
            self.quantify([crl_row('read', 0, 5, 5)])

    def test_em_accepts_plain_dicts_and_reports_zero_abundance(self):
        alpha = {'read': {'a': 0.5, 'b': 0.5}}
        rho = {'a': 0.5, 'b': 0.5}
        self.assertIs(quantify.em(alpha, rho, 1, ['read'], 0), alpha)
        self.assertEqual(rho, {'a': 0.5, 'b': 0.5})
        with self.assertRaisesRegex(ValueError, 'Invalid total abundance'):
            quantify.em(alpha, {'a': 0, 'b': 0}, 1, ['read'], 1e-5)

    def test_nonconvergent_em_does_not_silently_return_partial_result(self):
        alpha = {'unique': {'a': 1}, 'multi': {'a': 0.5, 'b': 0.5}}
        with self.assertRaisesRegex(ValueError, 'did not converge within 1000'):
            quantify.em(alpha, {'a': 0.75, 'b': 0.25}, 2, ['multi'], 0)


class FileAndCliTests(unittest.TestCase):
    def test_merge_to_quantification_with_grouping(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / 'input.bed'
            source.write_text(''.join(utilities.bedentry(reference, '0', '10', read, '+', '10M') + '\n'
                                      for read in ('read2|AA|1', 'read10|TT|1') for reference in ('tx1', 'tx2')))
            merge_dir = root / 'merge'
            outdir = root / 'counts with spaces'
            with contextlib.redirect_stdout(io.StringIO()):
                merging.main(['-b', str(source), '-o', str(merge_dir)])
                status = quantify.main(['-b', str(merge_dir / 'segments.bed'), '-m', str(merge_dir / 'merged.bed'),
                                        '-o', str(outdir), '-crl', '-ls', '1', '-cs', '1'])
            self.assertEqual(status, 0)
            rows = [row.split('\t') for row in (outdir / 'loci.counts').read_text().splitlines()]
            self.assertEqual(len(rows), 4)
            self.assertEqual({row[3] for row in rows}, {'0'})
            self.assertEqual({tuple(row[-3:]) for row in rows}, {('1.00', '1.00', '1e+06')})
            self.assertEqual(rows[0][0], 'read2|AA|1|1')
            self.assertEqual(rows[-1][0], 'read10|TT|1|1')
            self.assertEqual({path.name for path in outdir.iterdir()}, {'loci.counts'})

    def test_empty_pipeline(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / 'empty'
            source.write_text('')
            output = root / 'new output'
            with contextlib.redirect_stdout(io.StringIO()):
                quantify.main(['-b', str(source), '-m', str(source), '-o', str(output)])
            self.assertEqual((output / 'loci.counts').read_text(), '')

    def test_sort_failure_retains_intermediate_files(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / 'empty'
            source.write_text('')
            output = root / 'out'
            stderr = io.StringIO()
            with mock.patch.object(quantify.subprocess, 'run', side_effect=subprocess.CalledProcessError(7, 'sort')):
                with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(stderr):
                    with self.assertRaises(SystemExit) as error:
                        quantify.main(['-b', str(source), '-m', str(source), '-o', str(output)])
            self.assertEqual(error.exception.code, 2)
            self.assertTrue((output / 'loci.txt').exists())
            self.assertTrue((output / 'loci.counts.temp').exists())

    def test_invalid_size_does_not_create_output(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / 'out'
            with contextlib.redirect_stderr(io.StringIO()):
                with self.assertRaises(SystemExit):
                    quantify.main(['-b', 'input', '-m', 'merged', '-o', str(output), '-ls', '0'])
            self.assertFalse(output.exists())

    def test_missing_segment_position_reports_error(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            bed = root / 'segments.bed'
            bed.write_text('')
            merged = root / 'merged.bed'
            merged.write_text('tx\t0\t10\t+\tread|1,tx,0,10,+,10M\n')
            with self.assertRaisesRegex(ValueError, 'No segment BED position'):
                quantify.build_crls(False, str(bed), str(merged), str(root / 'loci'), 0.7, 1)


if __name__ == '__main__':
    unittest.main()
