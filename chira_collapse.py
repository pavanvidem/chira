#!/usr/bin/env python
from collections import Counter
import argparse
from Bio import SeqIO


def count_reads(fastq_handle, umi_len=0):
    """Count reads by trimmed sequence and 5' UMI from a FASTQ handle."""
    if umi_len < 0:
        raise ValueError('UMI length must be non-negative')

    read_counts = Counter()
    for record in SeqIO.parse(fastq_handle, "fastq"):
        sequence = str(record.seq)
        if len(sequence) <= umi_len:
            raise ValueError("Read {!r} has length {}, which must be greater than UMI length {}"
                             .format(record.id, len(sequence), umi_len))
        umi = sequence[:umi_len]
        trimmed_sequence = sequence[umi_len:]
        read_counts[trimmed_sequence, umi] += 1
    return read_counts


def write_fasta(read_counts, fasta_handle):
    """Write counts in sequence/UMI order with deterministic, one-based IDs."""
    for sequence_id, (sequence, umi) in enumerate(sorted(read_counts), start=1):
        header = str(sequence_id)
        if umi:
            header += "|" + umi
        header += "|" + str(read_counts[sequence, umi])
        fasta_handle.write(">{}\n{}\n".format(header, sequence))


def non_negative_int(value):
    """Parse a non-negative integer command-line argument."""
    number = int(value)
    if number < 0:
        raise argparse.ArgumentTypeError('UMI length must be non-negative')
    return number


def create_parser():
    parser = argparse.ArgumentParser(description='Chimeric Read Annotator: collapse FASTQ reads to FASTA format',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-i', '--fastq', required=True, metavar='FASTQ',
                        help='Input FASTQ file')
    parser.add_argument('-o', '--fasta', required=True, metavar='FASTA',
                        help='Output FASTA file')
    parser.add_argument("-u", '--umi_len', type=non_negative_int, default=0,
                        help="Length of the UMI, if present. "
                        "It is trimmed from the 5' end of each read and appended to the tag id. "
                        "Must be non-negative and shorter than every read")
    parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.4.3')
    return parser


def main(argv=None):
    parser = create_parser()
    args = parser.parse_args(argv)
    print('Input FASTQ          : ' + args.fastq)
    print('Output FASTA         : ' + args.fasta)
    print('Length of the UMI    : ' + str(args.umi_len))

    try:
        with open(args.fastq) as fastq_handle:
            read_counts = count_reads(fastq_handle, args.umi_len)
        with open(args.fasta, "w") as fasta_handle:
            write_fasta(read_counts, fasta_handle)
    except (OSError, ValueError) as error:
        parser.exit(2, '{}: error: {}\n'.format(parser.prog, error))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
