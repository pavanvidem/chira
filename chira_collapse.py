#!/usr/bin/env python
from collections import defaultdict
import argparse
from Bio import SeqIO


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Chimeric Read Annotator: collapse FASTQ reads to FASTA format',
                                     usage='%(prog)s [-h] [-v,--version]',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-i', '--fastq', action='store', dest='fastq', required=True, metavar='',
                        help='Input fastq file')
    parser.add_argument('-o', '--fasta', action='store', dest='fasta', required=True, metavar='',
                        help='Output fasta file')
    parser.add_argument("-u", '--umi_len', action='store', type=int, default=0, help="Length of the UMI, if present."
                        "It is trimmed from the 5' end of each read and appended to the tag id")

    args = parser.parse_args()
    print('Input FASTQ          : ' + args.fastq)
    print('Output FASTA         : ' + args.fasta)
    print('Length of the UMI    : ' + str(args.umi_len))

    d_uniq_reads = defaultdict(lambda: defaultdict(int))
    with open(args.fastq) as fh_fastq:
        for record in SeqIO.parse(fh_fastq, "fastq"):
            umi = str(record.seq)[0:args.umi_len]
            sequence = str(record.seq)[args.umi_len:]
            d_uniq_reads[sequence][umi] += 1
    c = 1
    with open(args.fasta, "w") as fh_out:
        for sequence in sorted(d_uniq_reads.keys()):
            for umi in sorted(d_uniq_reads[sequence]):
                readcount = d_uniq_reads[sequence][umi]
                seqid = str(c)
                if umi:
                    seqid += "|" + umi
                seqid += "|" + str(readcount)
                fh_out.write(">" + seqid + "\n" + sequence + "\n")
                c += 1
