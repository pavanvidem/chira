import contextlib
import io
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest import mock

import chira_merge as merging
import chira_utilities as utilities


class LocusTests(unittest.TestCase):
    def merge(self, intervals, cutoff):
        descriptions = {'chr1': {interval: [str(interval)] for interval in intervals}}
        return merging.merge_overlapping_intervals(descriptions, 'chr1', cutoff)[0]

    def test_half_open_boundaries_and_fraction(self):
        self.assertEqual(self.merge([(0, 10), (10, 20)], 0), [(0, 10), (10, 20)])
        self.assertEqual(self.merge([(0, 10), (5, 15)], 0.5), [(0, 15)])
        self.assertEqual(self.merge([(0, 10), (5, 15)], 0.51), [(0, 10), (5, 15)])
        self.assertEqual(self.merge([(0, 10), (5, 6)], 1), [(0, 10)])

    def test_single_base_and_invalid_intervals(self):
        self.assertEqual(self.merge([(0, 1), (1, 2)], 0.7), [(0, 1), (1, 2)])
        for interval in ((1, 1), (2, 1), (-1, 2)):
            with self.assertRaisesRegex(ValueError, 'Invalid BED'):
                self.merge([interval], 0.7)

    def test_filter_does_not_mutate_input_and_deduplicates_references(self):
        alignments = {(1, 10, '+'): {'tx1,0,10,+': 'a', 'tx2,0,10,+': 'b', 'tx1,20,30,+': 'c'}}
        result = merging.filter_alignments(alignments, 0, [], [], False, 0.9)
        self.assertEqual(list(result[1, 10, '+'].values()), ['a', 'b'])
        self.assertEqual(len(alignments[1, 10, '+']), 3)

    def test_split_reference_pair_with_large_gap(self):
        alignments = {(1, 10, '+'): {'tx1,0,10,+': 'first'},
                      (20001, 20010, '+'): {'tx2,0,10,+': 'second'}}
        result = merging.filter_alignments(alignments, 0, ['tx1'], ['tx2'], True, 0.9)
        self.assertEqual(len(result), 2)

    def test_longest_alignment_ties_preserve_umi_read_id(self):
        alignments = ['read|AA|2|1,tx,0,5,+,5M', 'read|AA|2|2,tx,0,10,+,10M',
                      'read|AA|2|3,tx,20,30,+,10M', 'read|TT|2|1,tx,0,5,+,5M']
        self.assertEqual(merging.longest_alignments(alignments), alignments[1:])


class AnnotationTests(unittest.TestCase):
    def test_nested_gff_features_use_name_fallback(self):
        with tempfile.TemporaryDirectory() as directory:
            annotation = Path(directory) / 'annotation.gff'
            annotation.write_text('##gff-version 3\n'
                                  'chr1\ttest\tgene\t1\t20\t.\t+\t.\tID=g\n'
                                  'chr1\ttest\tmRNA\t1\t20\t.\t+\t.\tID=t;Parent=g\n'
                                  'chr1\ttest\texon\t1\t10\t.\t+\t.\tID=e;Parent=t;Name=tx\n')
            merging.parse_annotations(str(annotation), directory)
            self.assertEqual((Path(directory) / 'transcriptomic_exons.bed').read_text(),
                             'tx\t0\t10\ttx_e001\t1\t+\n')

    def test_interleaved_exons_are_deduplicated_and_sorted_in_transcript_orientation(self):
        with tempfile.TemporaryDirectory() as directory:
            annotation = Path(directory) / 'annotation.gtf'
            annotation.write_text(
                'chr1\ttest\texon\t201\t210\t.\t+\t.\ttranscript_id "plus";\n'
                'chr1\ttest\texon\t101\t110\t.\t-\t.\ttranscript_id "minus";\n'
                'chr1\ttest\texon\t101\t110\t.\t+\t.\ttranscript_id "plus";\n'
                'chr1\ttest\texon\t201\t210\t.\t-\t.\ttranscript_id "minus";\n'
                'chr1\ttest\texon\t101\t110\t.\t+\t.\ttranscript_id "plus";\n')
            merging.parse_annotations(str(annotation), directory)
            self.assertEqual((Path(directory) / 'genomic_exons.bed').read_text(),
                             'chr1\t200\t210\tminus_e001\t1\t-\n'
                             'chr1\t100\t110\tminus_e002\t1\t-\n'
                             'chr1\t100\t110\tplus_e001\t1\t+\n'
                             'chr1\t200\t210\tplus_e002\t1\t+\n')
            rows = (Path(directory) / 'transcriptomic_exons.bed').read_text().splitlines()
            self.assertEqual([row.split('\t')[1:3] for row in rows],
                             [['0', '10'], ['10', '20'], ['0', '10'], ['10', '20']])

    def test_missing_identifier_is_not_reused_from_previous_feature(self):
        with tempfile.TemporaryDirectory() as directory:
            annotation = Path(directory) / 'annotation.gtf'
            annotation.write_text('chr1\ttest\texon\t1\t10\t.\t+\t.\ttranscript_id "tx";\n'
                                  'chr1\ttest\texon\t21\t30\t.\t+\t.\tgene_id "gene";\n')
            with self.assertRaisesRegex(ValueError, 'Missing transcript_id/Name'):
                merging.parse_annotations(str(annotation), directory)

    def test_conversion_clips_fragments_and_writes_final_spliced_alignment(self):
        with tempfile.TemporaryDirectory() as directory:
            genomic_exons = Path(directory) / 'genomic.bed'
            genomic_exons.write_text('chr1\t100\t110\ttx_e001\t1\t+\n'
                                     'chr1\t200\t210\ttx_e002\t1\t+\n')
            target = Path(directory) / 'result.bed'
            rows = ('tx\t2\t5\tfirst\t1\t+\ttx\t0\t10\ttx_e001\t1\t+\n'
                    'tx\t8\t15\tlast\t1\t+\ttx\t0\t10\ttx_e001\t1\t+\n'
                    'tx\t8\t15\tlast\t1\t+\ttx\t10\t20\ttx_e002\t1\t+\n')
            def intersect(command, stdout=None):
                self.assertEqual(command[0], 'intersectBed')
                stdout.write(rows)
            with mock.patch.object(merging, 'run_command', side_effect=intersect):
                merging.transcript_to_genomic_pos(str(Path(directory) / 'input'), str(target),
                                                  str(genomic_exons), 'tx exons.bed')
            self.assertEqual(target.read_text(), 'chr1\t102\t105\tfirst\t1\t+\n'
                             'chr1\t200\t205\tlast\t1\t+\n')
            self.assertFalse((Path(directory) / 'input.overlap.txt').exists())

    def test_negative_exon_and_reverse_alignment_strands(self):
        with tempfile.TemporaryDirectory() as directory:
            genomic_exons = Path(directory) / 'genomic.bed'
            genomic_exons.write_text('chr1\t200\t210\ttx_e001\t1\t-\n')
            target = Path(directory) / 'result.bed'
            rows = ('tx\t2\t5\tforward\t1\t+\ttx\t0\t10\ttx_e001\t1\t-\n'
                    'tx\t2\t5\treverse\t1\t-\ttx\t0\t10\ttx_e001\t1\t-\n')
            with mock.patch.object(merging, 'run_command', side_effect=lambda command, stdout: stdout.write(rows)):
                merging.transcript_to_genomic_pos(str(Path(directory) / 'input.bed'), str(target),
                                                  str(genomic_exons), 'tx.bed')
            self.assertEqual(target.read_text(), 'chr1\t205\t208\tforward\t1\t-\n'
                             'chr1\t205\t208\treverse\t1\t+\n')

    def test_conversion_ties_are_deterministic_and_empty_intersection_is_valid(self):
        with tempfile.TemporaryDirectory() as directory:
            exons = Path(directory) / 'exons.bed'
            exons.write_text('chr1\t100\t110\te1\t1\t+\nchr1\t200\t210\te2\t1\t+\n')
            rows = ['tx\t8\t12\tlast\t1\t+\ttx\t0\t10\te1\t1\t+\n',
                    'tx\t8\t12\tlast\t1\t+\ttx\t10\t20\te2\t1\t+\n']
            target = Path(directory) / 'result.bed'
            for intersections, expected in ((rows, 'chr1\t108\t110\tlast\t1\t+\n'),
                                             (list(reversed(rows)), 'chr1\t108\t110\tlast\t1\t+\n'),
                                             ([], '')):
                with mock.patch.object(merging, 'run_command',
                                       side_effect=lambda command, stdout: stdout.writelines(intersections)):
                    merging.transcript_to_genomic_pos(str(Path(directory) / 'input.bed'), str(target),
                                                      str(exons), 'transcript.bed')
                self.assertEqual(target.read_text(), expected)


class MergeCliTests(unittest.TestCase):
    def test_full_overlap_merge_and_empty_input(self):
        for contents, count in (('', 0), (utilities.bedentry('tx', '0', '10', 'read|1', '+', '10M') + '\n', 1)):
            with self.subTest(count=count), tempfile.TemporaryDirectory() as directory:
                source = Path(directory) / 'input.bed'
                source.write_text(contents)
                output = Path(directory) / 'new output'
                with contextlib.redirect_stdout(io.StringIO()):
                    self.assertEqual(merging.main(['-b', str(source), '-o', str(output)]), 0)
                self.assertEqual(len((output / 'merged.bed').read_text().splitlines()), count)
                if count:
                    self.assertEqual((output / 'merged.bed').read_text(),
                                     'tx\t0\t10\t+\tread|1|1,tx,0,10,+,10M\n')

    def test_invalid_arguments_do_not_create_output(self):
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / 'input.bed'
            source.write_text('')
            output = Path(directory) / 'out'
            for options in (['-ls', '0'], ['-co', '-1'], ['-f2', 'ref.fa']):
                with contextlib.redirect_stderr(io.StringIO()):
                    with self.assertRaises(SystemExit) as error:
                        merging.main(['-b', str(source), '-o', str(output)] + options)
                self.assertEqual(error.exception.code, 2)
                self.assertFalse(output.exists())

    def test_blockbuster_failure_stops_before_parsing_output(self):
        with tempfile.TemporaryDirectory() as directory:
            with mock.patch.object(merging, 'run_command', side_effect=subprocess.CalledProcessError(7, 'sort')) as run:
                with self.assertRaises(subprocess.CalledProcessError):
                    merging.merge_loci_blockbuster(directory, 30, 10, 10, 0.1, 0.7)
                self.assertEqual(run.call_count, 1)


if __name__ == '__main__':
    unittest.main()
