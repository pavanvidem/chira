import contextlib
import io
from pathlib import Path
import tempfile
import unittest

from chira_collapse import count_reads, main, write_fasta


def fastq(*sequences):
    return ''.join('@read{}\n{}\n+\n{}\n'.format(index, sequence, 'I' * len(sequence))
                   for index, sequence in enumerate(sequences, start=1))


class CollapseTests(unittest.TestCase):
    def collapse(self, contents, umi_len=0):
        output = io.StringIO()
        write_fasta(count_reads(io.StringIO(contents), umi_len), output)
        return output.getvalue()

    def test_counts_and_sorted_ids_without_umi(self):
        self.assertEqual(self.collapse(fastq('TT', 'AC', 'TT')),
                         '>1|1\nAC\n>2|2\nTT\n')

    def test_distinct_umis_are_kept_separate(self):
        self.assertEqual(self.collapse(fastq('TGAC', 'AAAC', 'TGAC', 'AATT'), 2),
                         '>1|AA|1\nAC\n>2|TG|2\nAC\n>3|AA|1\nTT\n')

    def test_output_is_independent_of_input_order(self):
        sequences = ['TGAC', 'AAAC', 'TGAC', 'AATT']
        self.assertEqual(self.collapse(fastq(*sequences), 2),
                         self.collapse(fastq(*reversed(sequences)), 2))

    def test_empty_input(self):
        self.assertEqual(self.collapse(''), '')

    def test_negative_umi_length(self):
        with self.assertRaisesRegex(ValueError, 'non-negative'):
            count_reads(io.StringIO(fastq('AC')), -1)

    def test_umi_must_be_shorter_than_read(self):
        for umi_len in (2, 3):
            with self.subTest(umi_len=umi_len):
                with self.assertRaisesRegex(ValueError, "Read 'read1' has length 2.*UMI length {}"
                                            .format(umi_len)):
                    self.collapse(fastq('AC'), umi_len)

    def test_empty_read_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "Read 'read1' has length 0"):
            self.collapse(fastq(''))

    def test_cli_writes_expected_output(self):
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / 'reads.fastq'
            output = Path(directory) / 'reads.fasta'
            source.write_text(fastq('TT', 'AC', 'TT'))
            with contextlib.redirect_stdout(io.StringIO()):
                status = main(['-i', str(source), '-o', str(output)])
            self.assertEqual(status, 0)
            self.assertEqual(output.read_text(), '>1|1\nAC\n>2|2\nTT\n')

    def test_cli_errors_preserve_existing_output(self):
        cases = [
            (fastq('AC'), ['-u', '-1'], 'non-negative'),
            (fastq('AACC', 'AA'), ['-u', '2'], "Read 'read2' has length 2"),
            (fastq('AACC', 'A'), ['-u', '2'], "Read 'read2' has length 1"),
            ('@broken\nAC\n+\nI\n', [], 'Lengths of sequence and quality values differ'),
            (None, [], 'No such file or directory'),
        ]
        for contents, options, message in cases:
            with self.subTest(message=message), tempfile.TemporaryDirectory() as directory:
                source = Path(directory) / 'reads.fastq'
                output = Path(directory) / 'reads.fasta'
                if contents is not None:
                    source.write_text(contents)
                output.write_text('existing output\n')
                stderr = io.StringIO()
                with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(stderr):
                    with self.assertRaises(SystemExit) as error:
                        main(['-i', str(source), '-o', str(output)] + options)
                self.assertEqual(error.exception.code, 2)
                self.assertIn(message, stderr.getvalue())
                self.assertNotIn('Traceback', stderr.getvalue())
                self.assertEqual(output.read_text(), 'existing output\n')

    def test_cli_reports_output_access_error(self):
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / 'reads.fastq'
            source.write_text(fastq('AC'))
            stderr = io.StringIO()
            with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(stderr):
                with self.assertRaises(SystemExit) as error:
                    main(['-i', str(source), '-o', directory])
            self.assertEqual(error.exception.code, 2)
            self.assertIn('Is a directory', stderr.getvalue())


if __name__ == '__main__':
    unittest.main()
