from collections import defaultdict
import io
from pathlib import Path
import tempfile
import unittest

from chira_merge import filter_alignments, reads_to_segments, write_segments
import chira_utilities as utilities


def alignment(reference, cigar):
    return utilities.bedentry(reference, '0', str(utilities.reference_length(cigar)),
                              'read|1', '+', cigar) + '\n'


class ReadOverlapTests(unittest.TestCase):
    def test_filter_chimera_overlap_threshold(self):
        for cutoff, expected_count in ((0, 0), (1, 2)):
            with self.subTest(cutoff=cutoff):
                alignments = defaultdict(dict)
                alignments[1, 10, '+']['tx1,0,10,+'] = alignment('tx1', '10M10S')
                alignments[10, 20, '+']['tx2,0,11,+'] = alignment('tx2', '9S11M')
                filtered = filter_alignments(alignments, cutoff, [], [], True, 0.9)
                self.assertEqual(len(filtered), expected_count)

    def test_segment_overlap_fraction_includes_both_endpoints(self):
        # [1, 3] and [2, 4] share 2/3 of their bases, not 1/2.
        filtered = {(1, 3, '+'): {'tx1': alignment('tx1', '3M1S')},
                    (2, 4, '+'): {'tx2': alignment('tx2', '1S3M')}}
        output = io.StringIO()
        write_segments('read|1', filtered, 0.6, output)
        segment_ids = [line.split('\t')[3].split(',')[0]
                       for line in output.getvalue().splitlines()]
        self.assertEqual(segment_ids, ['read|1|1', 'read|1|1'])

    def test_single_base_segments_do_not_divide_by_zero(self):
        filtered = {(1, 1, '+'): {'tx1': alignment('tx1', '1M1S')},
                    (2, 2, '+'): {'tx2': alignment('tx2', '1S1M')}}
        output = io.StringIO()
        write_segments('read|1', filtered, 0.7, output)
        segment_ids = [line.split('\t')[3].split(',')[0]
                       for line in output.getvalue().splitlines()]
        self.assertEqual(segment_ids, ['read|1|1', 'read|1|2'])

    def test_equal_and_mismatch_cigars_flow_through_segment_conversion(self):
        with tempfile.TemporaryDirectory() as directory:
            bed = Path(directory) / 'input.bed'
            bed.write_text(alignment('tx1', '3H5S10=2X10S') + alignment('tx2', '3H17S10M'))
            reads_to_segments(str(bed), directory, 0.7, 0, [], [], True, 0.9)
            rows = (Path(directory) / 'segments.temp.bed').read_text().splitlines()
            self.assertEqual(len(rows), 2)
            self.assertEqual([row.split('\t')[3].split(',')[0] for row in rows],
                             ['read|1|1', 'read|1|2'])


if __name__ == '__main__':
    unittest.main()
