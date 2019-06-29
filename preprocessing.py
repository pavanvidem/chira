#!/usr/local/bin/python
import os
import sys
import time
from collections import defaultdict
import itertools
import gzip
import argparse
from Bio import SeqIO

p_cutadapt = '/usr/local/tools/cutadapt/1.8/iuc/package_cutadapt_1_8/980a47047f57/bin/cutadapt'


def fastq_to_fasta_umi(fastq, fasta):
    d_uniq_umi_reads = defaultdict(lambda: defaultdict(int))
    with gzip.open(fastq, 'rb') as f:
        sequences = itertools.islice(f, 1, None, 4)
        for sequence in sequences:
            d_uniq_umi_reads[sequence.rstrip("\n")[3:]][sequence[:3]] += 1
    f.close()

    fh_fasta = open(fasta, "w")
    c = 1
    for sequence in sorted(d_uniq_umi_reads.keys()):
        for umi in d_uniq_umi_reads[sequence].keys():
            readcount = d_uniq_umi_reads[sequence][umi]
            fh_fasta.write(">tag_" + str(c) + "_" + umi + "|" + str(readcount) + "\n")
            fh_fasta.write(sequence + "\n")
        c += 1
    fh_fasta.close()


def fastq_to_fasta(fastq, fasta):
    d_uniq_reads = defaultdict(int)
    with gzip.open(fastq, 'rb') as f:
        sequences = itertools.islice(f, 1, None, 4)
        for sequence in sequences:
            d_uniq_reads[sequence.rstrip("\n")] += 1
    f.close()

    fh_fasta = open(fasta, "w")
    c = 1
    for sequence in sorted(d_uniq_reads.keys()):
        readcount = d_uniq_reads[sequence]
        fh_fasta.write(">tag_" + str(c) + "|" + str(readcount) + "\n")
        fh_fasta.write(sequence + "\n")
        c += 1
    fh_fasta.close()


def run_cutadapt(fastq, adapter3p, adapter5p, umi, outputprefix):
    minlen=18
    if umi:
        minlen = 21
    cutadapt_call = \
        (p_cutadapt +
#            ' --mask-adapter'
#            ' --discard-untrimmed'
            ' -a ' + adapter3p +
            ' -g ' + adapter5p +
            ' -O 5'
            ' -m ' + str(minlen) +
            ' -q 28' +
            ' --output ' + outputprefix + '.trimmed.fastq.gz' +
            # TODO do not hard code file name
            ' --too-short-output ' + outputprefix + '.less18nt.fastq.gz' +
            ' ' + fastq +
            ' | tee ' + outputprefix + '.cutadapt.log')
    print(cutadapt_call)
    os.system(cutadapt_call)


def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return ''.join([complement[base] for base in seq[::-1]])


# TODO currently duplicates are recognized based on sequences, but UMI based duplicated should be considered
def seq_to_tags(matepair1_prefix, matepair2_prefix, umi_present):
    matepair2_dict = defaultdict(str)
    if matepair2_prefix:
        if os.path.exists(matepair2_prefix + ".trimmed.fastq.gz"):
            with gzip.open(matepair2_prefix + ".trimmed.fastq.gz", "rt") as matepair2_handle:
                for record in SeqIO.parse(matepair2_handle, "fastq"):
                    matepair2_dict[record.id] = str(record.seq)
            matepair2_handle.close()
    d_uniq_fragments = defaultdict(int)
    d_uniq_reads1 = defaultdict(int)
    d_uniq_reads2 = defaultdict(int)
    # capture all the possible fragments and reads from matepair1
    with gzip.open(matepair1_prefix + ".trimmed.fastq.gz", "rt") as matepair1_handle:
        for record in SeqIO.parse(matepair1_handle, "fastq"):
            fragment = str(record.seq)
            # TODO option to group UMIs together, not just splitting
            if umi_present:
                fragment = fragment[3:]
            if record.id in matepair2_dict:
                fragment = fragment + '#' + matepair2_dict[record.id]
                d_uniq_fragments[fragment] += 1
                del matepair2_dict[record.id]
            else:
                d_uniq_reads1[fragment] += 1
    matepair1_handle.close()
    # capture the remaining lonely reads from matepair2
    for seq in matepair2_dict.values():
        d_uniq_reads2[seq] += 1

    matepair1_fasta = open(matepair1_prefix + ".fasta", "w")
    matepair2_fasta = None
    if matepair2_prefix:
        matepair2_fasta = open(matepair2_prefix + ".fasta", "w")
    c = 1
    for sequence in sorted(d_uniq_fragments.keys()):
        readcount = d_uniq_fragments[sequence]
        seqs = sequence.split('#')
        matepair1_fasta.write(">tag_" + str(c) + "|" + str(readcount) + "\n")
        matepair1_fasta.write(seqs[0] + "\n")
        if matepair2_prefix:
            matepair2_fasta.write(">tag_" + str(c) + "|" + str(readcount) + "\n")
            matepair2_fasta.write(reverse_complement(seqs[1]) + "\n")
        c += 1
    for sequence in sorted(d_uniq_reads1.keys()):
        readcount = d_uniq_reads1[sequence]
        matepair1_fasta.write(">tag_" + str(c) + "|" + str(readcount) + "\n")
        matepair1_fasta.write(sequence + "\n")
        c += 1
    if matepair2_prefix:
        for sequence in sorted(d_uniq_reads2.keys()):
            readcount = d_uniq_reads2[sequence]
            matepair2_fasta.write(">tag_" + str(c) + "|" + str(readcount) + "\n")
            matepair2_fasta.write(reverse_complement(sequence) + "\n")
            c += 1

    matepair1_fasta.close()
    if matepair2_fasta:
        matepair2_fasta.close()


def main(argv):
    parser = argparse.ArgumentParser(description='Fastq to Fasta',
                                     usage='%(prog)s [-h] [-v,--version]',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-r1', '--matepair1', action='store', dest='matepair1', required=True, metavar='',
                        help='Input fastq file')
    parser.add_argument('-r2', '--matepair2', action='store', dest='matepair2', required=False, metavar='',
                        help='Input matepair fastq file')
    parser.add_argument('-a', '--adapter3p', action='store', dest='adapter3p', required=True, metavar='',
                        help='Output fasta file')
    parser.add_argument('-g', '--adapter5p', action='store', dest='adapter5p', required=True, metavar='',
                        help='Output fasta file')
    parser.add_argument("-u", '--umi', action='store_true', help="If umi present, it is trimmed from the read1")
    parser.add_argument('-o', '--outdir', action='store', dest='outdir', required=True, metavar='',
                        help='Output fasta file')

    args = parser.parse_args()
    print('Input matepair1  : ' + args.matepair1)
    if args.matepair2:
        print('Input matepair2  : ' + args.matepair2)
    print('3p               : ' + args.adapter3p)
    print('5p               : ' + args.adapter5p)
    print('UMI present?     : ' + str(args.umi))
    print('Output directory : ' + args.outdir)

    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    matepair1_name = os.path.splitext(os.path.basename(args.matepair1))[0].replace('.gz', '').replace('gzip', '')
    matepair1_prefix = os.path.join(args.outdir, matepair1_name)
    matepair2_prefix = None
    if args.matepair2:
        matepair2_name = os.path.splitext(os.path.basename(args.matepair2))[0].replace('.gz', '').replace('gzip', '')
        matepair2_prefix = os.path.join(args.outdir, matepair2_name)

    run_cutadapt(args.matepair1, args.adapter3p, args.adapter5p, args.umi, matepair1_prefix)
    # for matepair2 take reverse complement of the adapters and switch them
    if args.matepair2:
        adapter3p = reverse_complement(args.adapter5p)
        adapter5p = reverse_complement(args.adapter3p)
        # UMI option for read2 is always false considering that UMI is always present after 5' adapter.
        run_cutadapt(args.matepair2, adapter3p, adapter5p, False, matepair2_prefix)
    print(matepair1_prefix)
    # TODO: also consider cases where there is a UMI in single end reads
    seq_to_tags(matepair1_prefix, matepair2_prefix, args.umi)


if __name__ == "__main__":
    main(sys.argv[1:])
