import contextlib
import io
from pathlib import Path
import shutil
import tempfile
import unittest

import chira_merge as merging
import chira_utilities as utilities


@unittest.skipUnless(shutil.which('intersectBed'), 'BEDTools is required')
class BedtoolsIntegrationTests(unittest.TestCase):
    def test_full_genomic_merge_on_both_strands_with_spliced_final_alignment(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            annotation = root / 'annotation with spaces.gtf'
            annotation.write_text(
                'chr1\ttest\texon\t101\t110\t.\t+\t.\ttranscript_id "plus";\n'
                'chr1\ttest\texon\t201\t210\t.\t+\t.\ttranscript_id "plus";\n'
                'chr2\ttest\texon\t301\t310\t.\t-\t.\ttranscript_id "minus";\n'
                'chr2\ttest\texon\t401\t410\t.\t-\t.\ttranscript_id "minus";\n')
            bed = root / 'input with spaces.bed'
            bed.write_text(
                utilities.bedentry('plus', '2', '5', 'read1|1', '+', '3M') + '\n' +
                utilities.bedentry('minus', '2', '5', 'read2|1', '+', '3M') + '\n' +
                utilities.bedentry('minus', '2', '5', 'read3|1', '-', '3M') + '\n' +
                utilities.bedentry('plus', '8', '15', 'read4|1', '+', '7M') + '\n')
            outdir = root / 'new output'
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(merging.main(['-b', str(bed), '-g', str(annotation), '-o', str(outdir)]), 0)
            rows = [row.split('\t') for row in (outdir / 'segments.bed').read_text().splitlines()]
            self.assertEqual([row[:3] + [row[5]] for row in rows],
                             [['chr1', '102', '105', '+'], ['chr2', '405', '408', '-'],
                              ['chr2', '405', '408', '+'], ['chr1', '200', '205', '+']])
            self.assertTrue(rows[-1][3].startswith('read4|1|1,'))
            self.assertEqual(len((outdir / 'merged.bed').read_text().splitlines()), 4)

    def test_empty_intersection(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / 'input.bed'
            source.write_text('')
            annotation = root / 'annotation.gtf'
            annotation.write_text('chr1\ttest\texon\t1\t10\t.\t+\t.\ttranscript_id "tx";\n')
            output = root / 'output'
            with contextlib.redirect_stdout(io.StringIO()):
                merging.main(['-b', str(source), '-g', str(annotation), '-o', str(output)])
            self.assertEqual((output / 'segments.bed').read_text(), '')
            self.assertEqual((output / 'merged.bed').read_text(), '')


@unittest.skipUnless(shutil.which('blockbuster.x'), 'Blockbuster is required')
class BlockbusterIntegrationTests(unittest.TestCase):
    def test_empty_input(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / 'input.bed'
            source.write_text('')
            output = root / 'output'
            with contextlib.redirect_stdout(io.StringIO()):
                merging.main(['-b', str(source), '-o', str(output), '-bb'])
            self.assertEqual((output / 'merged.bed').read_text(), '')

    def test_full_block_merge_preserves_alignment_descriptions(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / 'input with spaces.bed'
            source.write_text(''.join(utilities.bedentry('tx', '100', '120', 'read{}|1'.format(i), '+', '20M') + '\n'
                                      for i in range(1, 13)))
            output = root / 'new output'
            with contextlib.redirect_stdout(io.StringIO()):
                merging.main(['-b', str(source), '-o', str(output), '-bb', '-mc', '1', '-mb', '1'])
            rows = (output / 'merged.bed').read_text().splitlines()
            self.assertEqual(len(rows), 1)
            fields = rows[0].split('\t')
            self.assertEqual(fields[:4], ['tx', '100', '120', '+'])
            self.assertEqual(set(fields[4].split(';')),
                             {'read{}|1|1,tx,100,120,+,20M'.format(i) for i in range(1, 13)})


if __name__ == '__main__':
    unittest.main()
