import contextlib
import io
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest import mock
import subprocess

import chira_extract as extract


def count_row(read, segment, transcript, cigar, group, length=20, probability=1, share=1, tpm=10):
    start, end = 0, length
    return '\t'.join([read + '|' + str(segment), transcript, str(group), str(group), str(start), str(end), '+',
                      cigar, '{}:{}:{}:+'.format(transcript, start, end),
                      '{}:{}:{}:+'.format(transcript, start, end), str(share), str(probability), str(tpm)]) + '\n'


class ExtractionTests(unittest.TestCase):
    def test_partition_covers_every_read_once(self):
        for total in range(12):
            for workers in range(1, 15):
                ranges = extract.chunk_ranges(total, workers)
                self.assertEqual([i for start, end in ranges for i in range(start, end)], list(range(total)))

    def test_single_and_multi_worker_outputs_match(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fasta = root / 'ref.fa'
            fasta.write_text('>tx1\n' + 'A' * 30 + '\n>tx2\n' + 'T' * 30 + '\n')
            counts = root / 'loci.counts'
            counts.write_text(''.join(count_row('read{}|AA|1'.format(i), 1, 'tx1', '10M10S', 0, 10) +
                                      count_row('read{}|AA|1'.format(i), 2, 'tx2', '10S10M', 1, 10)
                                      for i in range(5)) + count_row('read5|TT|1', 1, 'tx1', '20M', 0))
            results = []
            for workers in (1, 2, 3, 10):
                output = root / str(workers)
                with contextlib.redirect_stdout(io.StringIO()):
                    self.assertEqual(extract.main(['-l', str(counts), '-o', str(output), '-f1', str(fasta),
                                                   '-p', str(workers), '-s']), 0)
                results.append(tuple((output / filename).read_text() for filename in ('chimeras', 'singletons', 'interactions')))
            self.assertTrue(all(result == results[0] for result in results))
            self.assertEqual(len(results[0][0].splitlines()), 6)
            self.assertEqual(len(results[0][1].splitlines()), 2)
            self.assertEqual(results[0][1].splitlines()[1].split('\t')[8], '30')

    def test_empty_input_produces_headers(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            counts = root / 'counts'
            counts.write_text('')
            fasta = root / 'ref'
            fasta.write_text('>tx\nACGT\n')
            output = root / 'out'
            with contextlib.redirect_stdout(io.StringIO()):
                extract.main(['-l', str(counts), '-o', str(output), '-f1', str(fasta), '-p', '4', '-s'])
            self.assertEqual(len((output / 'chimeras').read_text().splitlines()), 1)
            self.assertEqual(len((output / 'singletons').read_text().splitlines()), 1)
            self.assertEqual((output / 'interactions').read_text(), '')

    def test_tpm_cutoff_one_retains_maximum(self):
        with tempfile.TemporaryDirectory() as directory:
            counts = Path(directory) / 'counts'
            counts.write_text(count_row('read1', 1, 'tx', '20M', 0, tpm=5) +
                              count_row('read2', 1, 'tx', '20M', 1, tpm=10))
            self.assertEqual(extract.parse_counts_file(str(counts), 1), (2, 10))

    def test_annotation_lengths_are_deduplicated_and_state_is_independent(self):
        with tempfile.TemporaryDirectory() as directory:
            annotation = Path(directory) / 'annotation.gtf'
            annotation.write_text('chr1\ttest\texon\t1\t10\t.\t+\t.\ttranscript_id "tx"; gene_id "g";\n'
                                  'chr1\ttest\texon\t21\t30\t.\t+\t.\ttranscript_id "other";\n'
                                  'chr1\ttest\texon\t41\t50\t.\t+\t.\ttranscript_id "tx"; gene_id "g";\n'
                                  'chr1\ttest\texon\t1\t10\t.\t+\t.\ttranscript_id "tx"; gene_id "g";\n')
            first = extract.parse_annotations(str(annotation))
            second = extract.parse_annotations(str(annotation))
            self.assertEqual(first[1]['len']['tx'], 20)
            self.assertEqual(first, second)
            self.assertEqual(first[1]['gid']['other'], 'NA')

    def test_utr_without_cds_and_reference_colon(self):
        annotations = extract.empty_annotations()
        annotations[1]['UTR']['tx'] = [('chr:1', 0, 10, '+')]
        self.assertEqual(extract.guess_region('tx', 'chr:1:2:5:+', annotations), 'UTR')
        annotations[1]['CDS']['tx'] = [('chr:1', 10, 20, '+')]
        self.assertEqual(extract.guess_region('tx', 'chr:1:2:5:+', annotations), '5_prime_UTR')

    def test_reverse_hybrid_coordinates_and_end_positions(self):
        self.assertEqual(extract.hybridization_positions('..((..', '...))..'), (4, 5))
        self.assertEqual(extract.genomic_hybrid_position('chr:1:100:120:-', 2, 8), 'chr:1\t112\t118\t-')

    def test_worker_errors_are_propagated(self):
        with self.assertRaisesRegex(ValueError, 'Reference length unavailable'):
            extract.run_jobs(extract.reference_length, [('missing', {}, {})], 2)

    def test_split_reference_orientation_and_partial_length_fallback(self):
        annotations = extract.empty_annotations()
        annotations[1]['len']['tx2'] = 99
        chimeras, singletons = io.StringIO(), io.StringIO()
        rows = [count_row('read|1', 1, 'tx2', '10M10S', 1, 10),
                count_row('read|1', 2, 'tx1', '10S10M', 0, 10)]
        extract.extract_and_write('read|1', rows, set(), {'tx1': 30}, {'tx2': 40}, True,
                                  {}, 0, chimeras, singletons, False, annotations)
        fields = chimeras.getvalue().split('\t')
        self.assertEqual(fields[1:3], ['tx1', 'tx2'])
        self.assertEqual((fields[12], fields[16]), ('30', '99'))

    def test_reversed_summary_keeps_sequence_orientation_without_hybrid(self):
        with tempfile.TemporaryDirectory() as directory:
            chimeras, singletons = io.StringIO(), io.StringIO()
            rows = [count_row('read|1', 1, 'tx1', '10M10S', 0, 10),
                    count_row('read|1', 2, 'tx2', '10S10M', 1, 10)]
            extract.extract_and_write('read|1', rows, set(), {'tx1': 30, 'tx2': 30}, {}, None,
                                      {}, 0, chimeras, singletons, False, extract.empty_annotations())
            first = chimeras.getvalue().rstrip('\n').split('\t')
            first[29] = 'AAAA&CCCC'
            second = first[:]
            second[0] = 'read2|1'
            for left, right in ((1, 2), (7, 8), (20, 21), (24, 25), (26, 27)):
                second[left], second[right] = second[right], second[left]
            second[29] = 'CCCC&AAAA'
            (Path(directory) / 'chimeras').write_text('header\n' + '\t'.join(first) + '\n' + '\t'.join(second) + '\n')
            extract.write_interaction_summary(directory)
            result = (Path(directory) / 'interactions').read_text().split('\t')
            self.assertEqual(result[0], '2')
            self.assertEqual(result[9:11], ['AAAA', 'CCCC'])

    def test_checked_hybridization_failure(self):
        with mock.patch.object(extract.subprocess, 'run', side_effect=subprocess.CalledProcessError(7, 'IntaRNA')):
            with self.assertRaises(subprocess.CalledProcessError):
                extract.hybridize_with_intarna('AAAA', 'UUUU', [])

    def test_invalid_cli_arguments_do_not_create_output(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / 'out'
            for options in (['-p', '0'], ['-co', '-1'], ['-t', 'nan']):
                with contextlib.redirect_stderr(io.StringIO()):
                    with self.assertRaises(SystemExit):
                        extract.main(['-l', 'counts', '-f1', 'reference', '-o', str(output)] + options)
                self.assertFalse(output.exists())


@unittest.skipUnless(shutil.which('IntaRNA') and shutil.which('fastaFromBed'), 'IntaRNA and BEDTools required')
class HybridizationIntegrationTests(unittest.TestCase):
    def test_no_interaction_retains_sequences(self):
        params = ['--outMode', 'C', '--outCsvCols', 'id1,start1,start2,hybridDPfull,E', '--acc', 'N']
        self.assertEqual(extract.hybridize_with_intarna('A' * 20, 'A' * 20, params), ('NA', 'NA', 'NA'))
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / 'loci.fa.0').write_text('>tx1:0:20:+(+)\n' + 'A' * 20 + '\n>tx2:0:20:+(+)\n' + 'A' * 20 + '\n')
            chimeras = io.StringIO()
            rows = [count_row('read|1', 1, 'tx1', '20M20S', 0), count_row('read|1', 2, 'tx2', '20S20M', 1)]
            extract.extract_and_write('read|1', rows, set(), {'tx1': 20, 'tx2': 20}, {}, None,
                                      {}, 0, chimeras, io.StringIO(), True, extract.empty_annotations())
            (root / 'chimeras.0').write_text(chimeras.getvalue())
            extract.hybridize_and_write(directory, params, '0')
            fields = (root / 'chimeras-r.0').read_text().rstrip('\n').split('\t')
            self.assertEqual(fields[-4:], ['A' * 20 + '&' + 'A' * 20, 'NA', 'NA', 'NA'])

    def test_real_hybridization_and_summary(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fasta = root / 'reference with spaces.fa'
            fasta.write_text('>tx1\n' + 'G' * 20 + '\n>tx2\n' + 'C' * 20 + '\n')
            counts = root / 'counts'
            counts.write_text(count_row('read1|1', 1, 'tx1', '20M20S', 0) +
                              count_row('read1|1', 2, 'tx2', '20S20M', 1))
            output = root / 'output with spaces'
            with contextlib.redirect_stdout(io.StringIO()):
                extract.main(['-l', str(counts), '-o', str(output), '-f1', str(fasta), '-r', '-s'])
            row = (output / 'chimeras').read_text().splitlines()[1].split('\t')
            self.assertEqual(len(row), 33)
            self.assertNotEqual(row[30], 'NA')
            self.assertLess(float(row[32]), 0)
            self.assertEqual(len((output / 'interactions').read_text().splitlines()), 1)


if __name__ == '__main__':
    unittest.main()
