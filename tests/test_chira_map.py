import contextlib
import io
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

import pysam

import chira_map as mapping


class CommandTests(unittest.TestCase):
    def test_bwa_honors_hit_limits_and_literal_paths(self):
        with mock.patch.object(mapping, 'run_pipeline') as pipeline:
            for align_type, hits in (('long', 7), ('short', 13)):
                mapping.align_with_bwa(align_type, 'index1', 'reads ; query.fa',
                                       'reference index', 'output dir', 12, 18, 1, 4, 6, 1, hits, 2)
                command, consumer, output = pipeline.call_args.args
                self.assertEqual(command[command.index('-h') + 1], str(hits))
                self.assertEqual(command[-2:], ['reference index', 'reads ; query.fa'])
                self.assertEqual(consumer, ['samtools', 'view', '-hb', '-'])
                self.assertEqual(output, 'output dir/index1.{}.bam'.format(align_type))

    def test_pipeline_streams_literal_arguments(self):
        with tempfile.TemporaryDirectory() as directory:
            output = str(Path(directory) / 'result ; output')
            value = 'literal ; $value with spaces'
            producer = [sys.executable, '-c', 'import sys; sys.stdout.write(sys.argv[1])', value]
            consumer = [sys.executable, '-c', 'import sys; sys.stdout.write(sys.stdin.read().upper())']
            with contextlib.redirect_stdout(io.StringIO()):
                mapping.run_pipeline(producer, consumer, output)
            self.assertEqual(Path(output).read_text(), value.upper())

    def test_pipeline_detects_either_command_failure(self):
        drain = [sys.executable, '-c', 'import sys; sys.stdin.read()']
        fail = [sys.executable, '-c', 'import sys; sys.exit(7)']
        for producer, consumer, failed_command in ((fail, drain, fail), (drain, fail, fail)):
            with self.subTest(producer=producer), tempfile.TemporaryDirectory() as directory:
                with contextlib.redirect_stdout(io.StringIO()):
                    with self.assertRaises(subprocess.CalledProcessError) as error:
                        mapping.run_pipeline(producer, consumer, str(Path(directory) / 'output'))
                self.assertEqual(error.exception.returncode, 7)
                self.assertEqual(error.exception.cmd, failed_command)

    def test_clan_failure_stops_before_output_conversion(self):
        with mock.patch.object(mapping, 'run_command', side_effect=subprocess.CalledProcessError(1, 'clan_search')) as run:
            with self.assertRaises(subprocess.CalledProcessError):
                mapping.align_with_clan('query', 'output', 'ref', 'idx', None, None, 2, 16, 'fw', 100, 1)
            self.assertEqual(run.call_count, 1)

    def test_sort_failure_retains_input(self):
        with tempfile.TemporaryDirectory() as directory:
            bed = Path(directory) / 'mapped.bed'
            bed.write_text('input\n')
            with mock.patch.object(mapping, 'run_command', side_effect=subprocess.CalledProcessError(1, 'sort')):
                with self.assertRaises(subprocess.CalledProcessError):
                    mapping.sort_bed(directory)
            self.assertEqual(bed.read_text(), 'input\n')


class BamConversionTests(unittest.TestCase):
    def convert(self, records, stranded='fw'):
        header = pysam.AlignmentHeader.from_dict({'HD': {'SO': 'queryname'},
                                                  'SQ': [{'SN': 'tx1', 'LN': 1000},
                                                         {'SN': 'tx2', 'LN': 1000}]})
        with tempfile.TemporaryDirectory() as directory:
            bam = Path(directory) / 'input.bam'
            bed = Path(directory) / 'output.bed'
            fasta = Path(directory) / 'unmapped.fa'
            with pysam.AlignmentFile(str(bam), 'wb', header=header) as handle:
                for record in records:
                    handle.write(pysam.AlignedSegment.fromstring(record, header))
            mapping.write_mapped_bed(str(bam), str(bed), str(fasta), stranded)
            return bed.read_text(), fasta.read_text()

    def test_primary_and_xa_use_same_query_span_with_deletions(self):
        # The deletion-containing alternative spans only 8 query bases.
        record = 'read1\t0\ttx1\t101\t60\t10M\t*\t0\t0\tAACCGGTTAA\t*\tXA:Z:tx2,+201,4M5D4M2S,0;tx2,+301,10M,0;'
        bed, fasta = self.convert([record])
        self.assertEqual([line.split('\t')[:3] for line in bed.splitlines()],
                         [['tx1', '100', '110'], ['tx2', '300', '310']])
        self.assertEqual(fasta, '')

    def test_reverse_xa_has_bed_coordinates_and_rescues_read(self):
        record = 'read1\t0\ttx1\t101\t60\t10M\t*\t0\t0\tAACCGGTTAA\t*\tXA:Z:tx2,-201,5M2D5M,0;'
        bed, fasta = self.convert([record], 'rc')
        self.assertEqual(bed, 'tx2\t200\t212\tread1,tx2,200,212,-,5M2D5M\t1\t-\n')
        self.assertEqual(fasta, '')

    def test_missing_secondary_sequence_and_last_unmapped_read(self):
        records = [
            'read1\t16\ttx1\t101\t60\t10M\t*\t0\t0\tAACCGGTTAA\t*',
            'read1\t272\ttx2\t201\t60\t10M\t*\t0\t0\t*\t*',
            'read2\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\t*',
        ]
        bed, fasta = self.convert(records)
        self.assertEqual(bed, '')
        self.assertEqual(fasta, '>read1\nTTAACCGGTT\n>read2\nACGT\n')

    def test_empty_bam(self):
        self.assertEqual(self.convert([]), ('', ''))

    def test_mapped_record_without_sequence_is_accepted(self):
        record = 'read1\t256\ttx1\t101\t60\t10M\t*\t0\t0\t*\t*'
        bed, fasta = self.convert([record])
        self.assertTrue(bed.startswith('tx1\t100\t110\t'))
        self.assertEqual(fasta, '')

    def test_unmapped_read_without_any_sequence_reports_id(self):
        record = 'read1\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*'
        with self.assertRaisesRegex(ValueError, 'read1'):
            self.convert([record])


class CliAndClanTests(unittest.TestCase):
    def test_bwa_orchestration_merges_sorts_and_cleans_generated_bams(self):
        for split_reference in (False, True):
            with self.subTest(split_reference=split_reference), tempfile.TemporaryDirectory() as directory:
                source = Path(directory) / 'query.fa'
                source.write_text('>read1\nACGT\n')
                output = Path(directory) / 'output with spaces'

                def fake_align(align_type, index_type, query, index, outdir, *parameters):
                    bam = Path(outdir) / '{}.{}.bam'.format(index_type, align_type)
                    reference = 'tx1' if index_type == 'index1' else 'tx2'
                    header = pysam.AlignmentHeader.from_dict({'SQ': [{'SN': reference, 'LN': 100}]})
                    with pysam.AlignmentFile(str(bam), 'wb', header=header) as handle:
                        handle.write(pysam.AlignedSegment.fromstring(
                            'read1\t0\t{}\t1\t60\t4M\t*\t0\t0\tACGT\t*'.format(reference), header))
                    return str(bam)

                options = ['-i', str(source), '-o', str(output), '-x1', 'index1']
                if split_reference:
                    options += ['-x2', 'index2']
                with mock.patch.object(mapping, 'align_with_bwa', side_effect=fake_align) as align:
                    with contextlib.redirect_stdout(io.StringIO()):
                        self.assertEqual(mapping.main(options), 0)
                self.assertEqual(align.call_count, 4 if split_reference else 2)
                rows = (output / 'sorted.bed').read_text().splitlines()
                self.assertEqual(len(rows), 2 if split_reference else 1)
                self.assertEqual({row.split('\t')[0] for row in rows},
                                 {'tx1', 'tx2'} if split_reference else {'tx1'})
                self.assertEqual((output / 'unmapped.fasta').read_text(), '')
                self.assertEqual({path.name for path in output.iterdir()},
                                 {'sorted.bam', 'sorted.bed', 'unmapped.fasta'})

    def test_cli_command_failure_stops_pipeline(self):
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / 'query.fa'
            source.write_text('>read1\nACGT\n')
            stderr = io.StringIO()
            with mock.patch.object(mapping, 'align_with_bwa', side_effect=subprocess.CalledProcessError(7, ['bwa', 'mem'])):
                with mock.patch.object(mapping.pysam, 'merge') as merge:
                    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(stderr):
                        with self.assertRaises(SystemExit) as error:
                            mapping.main(['-i', str(source), '-o', directory, '-x1', 'index'])
                    merge.assert_not_called()
            self.assertEqual(error.exception.code, 2)
            self.assertIn('exit status 7', stderr.getvalue())
            self.assertNotIn('Traceback', stderr.getvalue())

    def test_invalid_options_do_not_create_output(self):
        cases = [([], '--index1'), (['-b'], '--ref_fasta1'),
                 (['-b', '-x1', 'idx'], 'mutually exclusive'),
                 (['-x2', 'idx'], '--index1'),
                 (['-a', 'clan', '-x1', 'idx'], 'CLAN requires --ref_fasta1'),
                 (['-x1', 'idx', '-p', '0'], 'positive'),
                 (['-x1', 'idx', '-co', '-1'], 'non-negative')]
        for options, message in cases:
            with self.subTest(options=options), tempfile.TemporaryDirectory() as directory:
                source = Path(directory) / 'query.fa'
                source.write_text('>read\nACGT\n')
                output = Path(directory) / 'new output'
                stderr = io.StringIO()
                with contextlib.redirect_stderr(stderr):
                    with self.assertRaises(SystemExit) as error:
                        mapping.main(['-i', str(source), '-o', str(output)] + options)
                self.assertEqual(error.exception.code, 2)
                self.assertIn(message, stderr.getvalue())
                self.assertFalse(output.exists())

    def test_clan_second_index_requires_reference(self):
        args = mapping.create_parser().parse_args(['-i', 'query', '-o', 'out', '-a', 'clan',
                                                  '-f1', 'ref', '-x1', 'idx', '-x2', 'idx2'])
        with self.assertRaisesRegex(ValueError, 'together'):
            mapping.validate_args(args)

    def test_clan_conversion_preserves_clipping_and_reference_id(self):
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / 'out.map'
            source.write_text('header\nread1\t1\t3\t8\t10\ttx:1 description:101-106;\n')
            mapping.clan_to_bed(directory)
            self.assertEqual((Path(directory) / 'mapped.bed').read_text(),
                             'tx:1\t100\t106\tread1,tx:1,100,106,+,2S6M2S\t1\t+\n')

    def test_clan_empty_output(self):
        with tempfile.TemporaryDirectory() as directory:
            (Path(directory) / 'out.map').write_text('')
            mapping.clan_to_bed(directory)
            self.assertEqual((Path(directory) / 'mapped.bed').read_text(), '')


if __name__ == '__main__':
    unittest.main()
