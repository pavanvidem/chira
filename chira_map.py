#!/usr/bin/env python
import argparse
import os
import sys
from collections import defaultdict
import pysam
from multiprocessing import Process
import logging
import datetime
import math
import chira_utilities

def align_with_bwa(align_type, index_type, fasta, refindex, outdir, seed_length, align_score, processes):
    """
        Funtion that maps the reads to the transcriptome. Different parameters
        are used for long and short alignments.
        Parameters:
            align_type: String of alignmnet type (long or short) used to format the output file name.
            index_type: String of index1 type (index1 or index2) used to format the output file name.
            fasta: Path to query fasta file. For align_type short, long.unmapped.fa from the out_dir used.
            refindex: Path to the reference index file.
            outdir: Output directory path.
            seed_length: Minimum seed length
            align_score: Minimum alignment score
            processes: Number of processes to use
    """

    logging.info("| STRAT: Mapping " + align_type + "reads segments to " + index_type + " at " + refindex)

    query_fasta = mismatch_score = gap_o = gap_e = n_aligns = None
    bam = os.path.join(outdir, index_type + "." + align_type + ".bam")

    if align_type == "long":
        query_fasta = fasta
        mismatch_score = "4"
        n_aligns = "50"
        gap_o = "6"
        gap_e = "1"
    elif align_type == "short":
        query_fasta = os.path.join(outdir, "long.unmapped.fa")
        mismatch_score = "100"
        n_aligns = "100"
        gap_o = "100"
        gap_e = "100"

    bwa_params = ["-r 1",                       # look for internal seeds inside a seed longer than {-k} * {-r}
                  "-c 1000",
                  "-B " + mismatch_score,       # mismatch penalty
                  "-O " + gap_o,                # gap open penalty
                  "-E " + gap_e,                # gap extension penalty
                  "-L 0",                       # clipping penalty, we need soft clips
                  "-Y",                         # use soft clipping for supplementary alignments
                  "-h " + n_aligns,             # if there're -h hits with score >80% of the maxscore, output in XA
                  "-k " + str(seed_length),     # minimum seed length
                  "-T " + str(align_score),     # minimum alignment score
                  "-t " + str(processes),
                  refindex,
                  query_fasta
                  ]
    bwacall = ("bwa mem " + " ".join(bwa_params) + " | samtools view -hb - > " + bam)
    print(bwacall)
    os.system(bwacall)
    logging.info("| END: Mapping " + align_type + "reads segments to " + index_type + " at " + refindex)


def merge_bams(outdir, align_type, processes):
    logging.info("| START: merge and sort " + align_type + " alignments of both indices")
    # -f to force if file already exists
    pysam.merge("-f",  os.path.join(outdir, align_type + ".unsorted.bam"),
                os.path.join(outdir, "index1." + align_type + ".bam"),
                os.path.join(outdir, "index2." + align_type + ".bam"))
    pysam.sort("-m", "1G", "-@", str(processes),
               os.path.join(outdir, align_type + ".unsorted.bam"),
               "-T", os.path.join(outdir, align_type),
               "-o", os.path.join(outdir, align_type + ".bam"))
    pysam.index(os.path.join(outdir, align_type + ".bam"))
    os.path.join(outdir, align_type + ".unsorted.bam")
    if os.path.exists(os.path.join(outdir, "index1." + align_type + ".bam")):
        os.remove(os.path.join(outdir, "index1." + align_type + ".bam"))
    if os.path.exists(os.path.join(outdir, "index2." + align_type + ".bam")):
        os.remove(os.path.join(outdir, "index2." + align_type + ".bam"))

    logging.info("| END: merge and sort " + align_type + " alignments both indices")


def extract_singleton_reads(s, e, d_read_positions, co, d_mapped_reads, singletons_fasta, i):
    """
        For long alignments, this function to checks if there is at least an alignment per read
        with 2 non-overlapping arms. If there such an alignment, read qualifies as chimeric.
        Otherwise, as singleton and added to the given dictionary.
        Parameters:
            s: Start index of the chunk
            e: End index of the chunk
            d_read_positions: Dictionaries containing  the  [readid] -> [(read_sart, read_end; reference_id)]
            co: Maximum allowed overlap between the aligned read arms
            d_mapped_reads: dictionary containing all the mapped sequences
            singletons_fasta: File path to write singleton fasta sequences temporarily
            i: index to add to fasta file indicating the process count
    """
    d_singleton_reads = defaultdict(list)
    l_proper_chimeric_reads = []
    for readid in list(d_read_positions)[s:e]:
        if readid in l_proper_chimeric_reads:
            continue
        proper_chimera = False
        for (read_start, read_end, ref_id) in d_read_positions[readid]:
            # mark as unmapped if there is only one mapped segment that doesn't cover the entire read
            if readid in d_singleton_reads:
                for prev_pos in d_singleton_reads[readid]:
                    # same reference id
                    # TODO: accept as a proper chimera if its mapped parts are not overlapping
                    if ref_id == prev_pos[2]:
                        continue
                    # <= co indicates no overlap between read segments,
                    # a second non-overlapping soft clipped alignment found
                    if chira_utilities.overlap([prev_pos[0], prev_pos[1]], [read_start, read_end]) <= co:
                        proper_chimera = True
                        break
            if proper_chimera:
                l_proper_chimeric_reads.append(readid)
                if readid in d_singleton_reads:
                    del d_singleton_reads[readid]
            else:
                if [read_start, read_end, ref_id] not in d_singleton_reads[readid]:
                    d_singleton_reads[readid].append([read_start, read_end, ref_id])

    with open(singletons_fasta + "." + i, "w") as fh_singletons:
        for readid in d_singleton_reads.keys():
            fh_singletons.write(">" + readid + "\n" + d_mapped_reads[readid] + "\n")


def extract_unmapped(align_type, outdir, chimeric_overlap, stranded, processes):
    """
        Extracts the unmapped and singlton reads from the BAM alignments and writes them to a file
        Parameters:
            align_type: String of alignmnet type (long or short) used to get the BAM file path.
            outdir: Output directory path.
            chimeric_overlap: Maximum number of bases allowed between the chimeric segments of a read
            stranded: Strand specificity
            processes: Number of processes
    """
    logging.info("| START: extract unmapped reads of " + align_type + " alignments")
    fh_bam = pysam.Samfile(os.path.join(outdir, align_type + ".bam"), "rb")
    fh_mapped_bam = pysam.Samfile(os.path.join(outdir, align_type + ".mapped.bam"), "wb",
                                  template=fh_bam)
    fh_mapped_bed = open(os.path.join(outdir, align_type + ".mapped.bed"), "w")
    d_mapped_reads = defaultdict(str)
    d_unmapped_reads = defaultdict(str)
    d_read_positions = defaultdict(list)  # contains a tuple of matched base positions on the read and reference id
    for alignment in fh_bam.fetch(until_eof=True):
        readid = alignment.query_name
        readseq = None
        if alignment.is_reverse:
            readseq = chira_utilities.reverse_complement(alignment.query_sequence)
        else:
            readseq = alignment.query_sequence
        proper_stranded_align_found = False
        if alignment.is_unmapped or \
                (stranded == "fw" and alignment.is_reverse and not alignment.has_tag('XA')) or \
                (stranded == "rc" and not alignment.is_reverse and not alignment.has_tag('XA')):
            if readid not in d_mapped_reads:
                if readseq:
                    d_unmapped_reads[readid] = readseq
                else:
                    d_unmapped_reads[readid] = d_mapped_reads[readid]
        else:
            if readid in d_unmapped_reads:
                del d_unmapped_reads[readid]
            d_mapped_reads[readid] = readseq
            fh_mapped_bam.write(alignment)
            # in this case there might be alignments on wrong strand with XA tag
            # hence, write the alignment only if it mapped on desired strand
            optimal_match_len = 0
            if (stranded == "rc" and alignment.is_reverse) or (stranded == "fw" and not alignment.is_reverse):
                fh_mapped_bed.write(chira_utilities.bedentry(fh_bam.getrname(alignment.tid),
                                                             str(alignment.reference_start),
                                                             str(alignment.reference_end),
                                                             readid,
                                                             "-" if alignment.is_reverse else "+",
                                                             alignment.cigarstring) + "\n")
                d_read_positions[readid].append([alignment.query_alignment_start,
                                                alignment.query_alignment_end,
                                                fh_bam.getrname(alignment.tid)])
                proper_stranded_align_found = True
                optimal_match_len = chira_utilities.alignment_score(alignment.cigarstring, alignment.get_tag("NM"))
            # XA tag only present in bwa output
            if not alignment.has_tag('XA'):
                continue

            alt_alignments = alignment.get_tag('XA').rstrip(';').split(';')
            # either optimal_match_len already set or there must be at least a secondary alignment on desired strand
            if optimal_match_len == 0:
                for alt_alignment in alt_alignments:
                    f_alt_align = alt_alignment.split(',')
                    if f_alt_align[1].startswith('-'):
                        continue
                    alt_cigar = f_alt_align[2]
                    alt_nm = int(f_alt_align[3])
                    alt_match_len = chira_utilities.alignment_score(alt_cigar, alt_nm)
                    if alt_match_len > optimal_match_len:
                        optimal_match_len = alt_match_len

            # now examine XA tag for alternate alignments
            for alt_alignment in alt_alignments:
                f_alt_align = alt_alignment.split(',')
                # check if it mapped on desired strand
                if f_alt_align[1].startswith('-') and stranded == "fw" or \
                        f_alt_align[1].startswith('+') and stranded == "rc":
                    continue
                alt_referenceid = f_alt_align[0]
                alt_refstart = f_alt_align[1][1:]
                alt_refstrand = f_alt_align[1][0]
                alt_cigar = f_alt_align[2]
                alt_nm = int(f_alt_align[3])
                alt_refend = chira_utilities.alignment_end(alt_refstart, alt_cigar, f_alt_align[1].startswith("-"))

                alt_read_start, alt_read_end = chira_utilities.match_positions(alt_cigar,
                                                                               f_alt_align[1].startswith("-"))
                alt_match_len = chira_utilities.alignment_score(alt_cigar, alt_nm)
                # TODO an option to consider Suboptimal alignments by a threshold
                if alt_match_len < optimal_match_len:
                    continue
                fh_mapped_bed.write(chira_utilities.bedentry(alt_referenceid,
                                                             str(int(alt_refstart)-1),
                                                             str(alt_refend),
                                                             readid,
                                                             alt_refstrand,
                                                             alt_cigar) + "\n")
                d_read_positions[readid].append([alt_read_start, alt_read_end, alt_referenceid])
                proper_stranded_align_found = True
            # if there is not a single alignment on desired strand then consider the read as unmapped
            if not proper_stranded_align_found:
                if readid not in d_mapped_reads:
                    if readseq:
                        d_unmapped_reads[readid] = readseq
                    else:
                        d_unmapped_reads[readid] = d_mapped_reads[readid]
    fh_bam.close()
    fh_mapped_bam.close()
    fh_mapped_bed.close()

    if align_type == "short":
        with open(os.path.join(outdir, align_type + ".unmapped.fa"), "w") as fh_unmapped_fasta:
            for readid, readseq in d_unmapped_reads.items():
                fh_unmapped_fasta.write(">" + readid + "\n" + readseq + "\n")
        return
    # further process long mapped reads
    jobs = []
    print(str(datetime.datetime.now()), " START: multiprocessing")
    singletons_fasta = os.path.join(args.outdir, "siingleton_reads")
    for i in range(processes):
        s = i * math.ceil(len(d_read_positions)/processes)
        e = min(s + math.ceil(len(d_read_positions)/processes), len(d_read_positions))
        j = Process(target=extract_singleton_reads,
                    args=(s, e, d_read_positions, chimeric_overlap, d_mapped_reads, singletons_fasta, str(i)))
        jobs.append(j)
        j.start()
    for j in jobs:
        j.join()
    print(str(datetime.datetime.now()), " END: multiprocessing")
    # write all unmapped sequences to a fasta file
    with open(os.path.join(outdir, align_type + ".unmapped.fa"), "w") as fh_unmapped_fasta:
        for readid, readseq in d_unmapped_reads.items():
            fh_unmapped_fasta.write(">" + readid + "\n" + readseq + "\n")
        for n in range(processes):
            with open(singletons_fasta + "." + str(n)) as infile:
                for line in infile:
                    fh_unmapped_fasta.write(line)
            os.remove(singletons_fasta + "." + str(n))
    logging.info("| END: extract unmapped reads of " + align_type + " alignments")
    return


def merge_beds(long_mapped_bed, short_mapped_bed, mapped_bed):
    """
    Merge long and short mapped BED entries. There are chances that both long and
    short alignment settings lead to same alignments for some reads.

    Parameters:
       long_mapped_bed: Path to long mapped BED file
       short_mapped_bed: Path to short mapped BED file
       mapped_bed: Path to final output mapped BED file
    """

    lines1 = set(l for l in open(long_mapped_bed).read().splitlines())
    with open(mapped_bed, "w") as fh_out:
        for l in lines1:
            fh_out.write(l + "\n")
        with open(short_mapped_bed) as fh_short:
            for line in fh_short:
                if line.rstrip("\n") in lines1:
                    continue
                fh_out.write(line)


def align_with_clan(query_fasta, outdir, ref_fasta1, ref_index1, ref_fasta2, ref_index2,
                    chimeric_overlap, align_score, stranded, processes):
    logging.info("| START: Mapping reads using CLAN")

    ref_index = "-f " + ref_fasta1 + " -d " + ref_index1
    if ref_fasta2 and ref_index2:
        ref_index += " -F " + ref_fasta2 + " -D " + ref_index2
    n_aligns = 100
    clan_search_params = ["-r " + query_fasta,          # query fasta file
                          "-m " + str(n_aligns),        # number of maximum hits for each maximal fragment
                          "-l " + str(align_score),     # minimum length for each fragment
                          "-t " + str(processes),
                          "-v " + str(chimeric_overlap),
                          "-o " + os.path.join(outdir, "out.clan"),
                          ref_index]
    clan_search = ("clan_search " + " ".join(clan_search_params))
    print(clan_search)
    os.system(clan_search)
    logging.info("| END: Mapping reads using CLAN")

    ref_fasta = "-f " + ref_fasta1
    if ref_fasta2 and ref_index2:
        ref_fasta += " -F " + ref_fasta2

    clan_output_params = ["-r " + query_fasta,          # query fasta file
                          "-i " + os.path.join(outdir, "out.clan"),
                          "-o " + os.path.join(outdir, "out.map"),
                          ref_fasta]
    clan_output = ("clan_output " + " ".join(clan_output_params))
    print(clan_output)
    os.system(clan_output)


def clan_to_bed(outdir):
    with open(os.path.join(outdir, "out.map")) as fh_in, open(os.path.join(outdir, "mapped.bed"), "w") as fh_out:
        next(fh_in)
        for line in fh_in:
            [read_id,
             solution_id,
             read_mapped_begin,
             read_mapped_end,
             read_length,
             mapped_locations] = line.rstrip("\n").rstrip("\t").split("\t")
            lead_soft_clips = int(read_mapped_begin) - 1
            trail_soft_clips = int(read_length) - int(read_mapped_end)
            cigar = ""
            if lead_soft_clips != 0:
                cigar += str(lead_soft_clips) + "S"
            cigar += str(int(read_mapped_end) - int(read_mapped_begin) + 1) + "M"
            if trail_soft_clips != 0:
                cigar += str(trail_soft_clips) + "S"

            for mapped_location in mapped_locations.split(";"):
                d = mapped_location.split(":")
                # if header has spaces select the id only
                ref_id = ":".join(d[0:-1]).split(' ')[0]
                [ref_start, ref_end] = d[-1].split("-")
                # TODO: consider both strands
                fh_out.write("\t".join([ref_id, str(int(ref_start)-1), ref_end,
                                        ",".join([read_id, ref_id, str(int(ref_start)-1), ref_end, "+", cigar]),
                                        "1", "+"]) + "\n")


def score_float(x):
    x = float(x)
    if x < 0.0 or x > 1.0:
        raise argparse.ArgumentTypeError("%r not in range [0.0, 1.0]" % (x,))
    return x


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Chimeric Read Annotator: map reads to the reference',
                                     usage='%(prog)s [-h] [-v,--version]',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-a", '--aligner', type=str, choices=["bwa", "clan"], default='bwa', required=True,
                        dest='aligner', metavar='', help='Alignment program to use, bwa or clan')

    parser.add_argument('-i', '--query_fasta', action='store', dest='fasta', required=True,
                        metavar='', help='Path to query fasta file')

    parser.add_argument('-o', '--outdir', action='store', dest='outdir', required=True, metavar='',
                        help='Output directory path for the analysis')

    parser.add_argument('-x1', '--index1', action='store', dest='idx1', required=False,
                        metavar='', help='first prioroty index file')

    parser.add_argument('-x2', '--index2', action='store', dest='idx2', required=False,
                        metavar='', help='second priority index file')

    parser.add_argument('-f1', '--ref_fasta1', action='store', dest='ref_fasta1', required=False,
                        metavar='', help='First prioroty fasta file')

    parser.add_argument('-f2', '--ref_fasta2', action='store', dest='ref_fasta2', required=False,
                        metavar='', help='second priority fasta file')

    parser.add_argument("-b", '--build', action='store_true', dest='build_index',
                        help="Build indices from reference fasta files")

    parser.add_argument('-p', '--processes', action='store', type=int, default=1, metavar='',
                        dest='processes',
                        help='Number of processes to use')

    parser.add_argument("-s", '--stranded', type=str, choices=["fw", "rc", "both"], default='fw', metavar='',
                        dest='stranded',
                        help='''Strand-specificity of input samples.
                             fw = map to transcript strand;
                             rc = map to reverse compliment of transcript strand;
                             both = try to map on both strnads''')

    parser.add_argument("-l1", '--seed_length1', action='store', type=int, default=12, metavar='',
                        dest='seed_length1',
                        help='''Seed length for 1st mapping iteration.
                                bwa-mem parameter "-k"''')

    parser.add_argument("-l2", '--seed_length2', action='store', type=int, default=6, metavar='',
                        dest='seed_length2',
                        help='''Seed length for 2nd mapping iteration.
                                bwa-mem parameter "-k"''')

    parser.add_argument("-s1", '--align_score1', action='store', type=int, default=18, metavar='',
                        dest='align_score1',
                        help='''Minimum alignment score in 1st mapping iteration.
                                bwa-mem parameter "-T" and clan_search parameter "-l"''')

    parser.add_argument("-s2", '--align_score2', action='store', type=int, default=10, metavar='',
                        dest='align_score2',
                        help='''Minimum alignment score in 2nd mapping iteration.
                                It must be smaller than --align_score1 parameter.
                                bwa-mem parameter "-T" and clan_search parameter "-l"''')

    parser.add_argument('-co', '--chimeric_overlap', action='store', type=int, default=2, metavar='',
                        dest='chimeric_overlap',
                        help='Maximum number of bases allowed between the chimeric segments of a read')

    parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.0')

    args = parser.parse_args()
    print('Query fasta                          : ' + args.fasta)
    print('Output directory                     : ' + args.outdir)
    print('Aligner                              : ' + args.aligner)
    print('Build index?                         : ' + str(args.build_index))
    if args.idx1:
        print('1st priority BWA/CLAN index          : ' + args.idx1)
    if args.idx2:
        print('2nd priority BWA/CLAN  index         : ' + args.idx2)
    if args.ref_fasta1:
        print('1st priority reference fasta file    : ' + args.ref_fasta1)
    if args.ref_fasta2:
        print('2nd priority reference fasta file    : ' + args.ref_fasta2)
    print('Number of processes                  : ' + str(args.processes))
    print('Stranded                             : ' + args.stranded)
    print('Seed length                          : ' + str(args.seed_length1))
    if args.seed_length2:
        print('Seed length for 2nd iteration        : ' + str(args.seed_length2))
    print('Alignment score                      : ' + str(args.align_score1))
    if args.align_score2:
        print('Alignment score for 2nd iteration    : ' + str(args.align_score2))
    print('Chimeric overlap                     : ' + str(args.chimeric_overlap))
    print("===================================================================")

    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    if not args.idx1 and not args.ref_fasta1:
        sys.stderr.write("option -x1 or -f1 are required")
        sys.exit(1)

    if (args.idx1 or args.idx2) and args.build_index:
        sys.stderr.write("options -b and -x1 are mutually exclusive")
        sys.exit(1)
    if not args.idx1 and not args.idx2 and not args.build_index:
        sys.stderr.write("Either -b or -x1 is required")
        sys.exit(1)

    index1 = index2 = None
    if args.idx1:
        index1 = args.idx1
    if args.idx2:
        index2 = args.idx2
    if args.aligner == "clan":
        if args.build_index:
            index1 = os.path.join(args.outdir, "index1")
            os.system("clan_index -f " + args.ref_fasta1 + " -d " + index1)
            if args.ref_fasta2:
                index2 = os.path.join(args.outdir, "index2")
                os.system("clan_index -f " + args.ref_fasta2 + " -d " + index2)
        # use only align_score2 for mapping with CLAN
        align_with_clan(args.fasta, args.outdir, args.ref_fasta1, index1, args.ref_fasta2, index2,
                        args.chimeric_overlap, args.align_score2, args.stranded, args.processes)
        clan_to_bed(args.outdir)
    elif args.aligner == "bwa":
        # build indices
        if args.build_index:
            index1 = os.path.join(args.outdir, "index1")
            os.system("bwa index -p " + index1 + " " + args.ref_fasta1)
            if args.ref_fasta2:
                index2 = os.path.join(args.outdir, "index2")
                os.system("bwa index -p " + index2 + " " + args.ref_fasta2)
        # align with bwa
        align_with_bwa("long", "index1", args.fasta, index1,  args.outdir,
                       args.seed_length1, args.align_score1, args.processes)
        if index2:
            align_with_bwa("long", "index2", args.fasta, index2, args.outdir,
                           args.seed_length1, args.align_score1, args.processes)
            # if 2 indices given, merge bams from both indices and sort
            # TODO: handle cases where nothing maps to index1 or index2 and try to merge empty bam files
            merge_bams(args.outdir, "long", args.processes)
        else:
            logging.info("| START: sort long alignments")
            # if only  one index given just sort it
            pysam.sort("-m", "1G", "-@", str(args.processes),
                       os.path.join(args.outdir, "index1.long.bam"),
                       "-T", os.path.join(args.outdir, "index1.long"),
                       "-o", os.path.join(args.outdir, "long.bam"))
            pysam.index(os.path.join(args.outdir, "long.bam"))
            logging.info("| END: sort long alignments")
        print(str(datetime.datetime.now()), " START:extract unmapped long")
        extract_unmapped("long", args.outdir, args.chimeric_overlap, args.stranded, args.processes)
        if os.path.exists(os.path.join(args.outdir, "long.bam")):
            os.remove(os.path.join(args.outdir, "long.bam"))

        print(str(datetime.datetime.now()), " START:extract unmapped long")
        align_with_bwa("short", "index1", args.fasta, index1, args.outdir,
                       args.seed_length2, args.align_score2, args.processes)
        if index2:
            align_with_bwa("short", "index2", args.fasta, index2, args.outdir,
                           args.seed_length2, args.align_score2, args.processes)
            # if 2 indices given, merge bams from both indices and sort
            merge_bams(args.outdir, "short", args.processes)
        else:
            logging.info("| START: sort short alignments")
            # if only  one index given just sort it
            pysam.sort("-m", "1G", "-@", str(args.processes),
                       os.path.join(args.outdir, "index1.short.bam"),
                       "-T", os.path.join(args.outdir, "index1.short"),
                       "-o", os.path.join(args.outdir, "short.bam"))
            pysam.index(os.path.join(args.outdir, "short.bam"))
            logging.info("| END: sort short alignments")
        extract_unmapped("short", args.outdir, args.chimeric_overlap, args.stranded, args.processes)
        if os.path.exists(os.path.join(args.outdir, "short.bam")):
            os.remove(os.path.join(args.outdir, "short.bam"))

        logging.info("| START: merge both long and short alignments")
        merge_beds(os.path.join(args.outdir, "long.mapped.bed"),
                   os.path.join(args.outdir, "short.mapped.bed"),
                   os.path.join(args.outdir, "mapped.bed"))
        if os.path.exists(os.path.join(args.outdir, "long.mapped.bed")):
            os.remove(os.path.join(args.outdir, "long.mapped.bed"))
        if os.path.exists(os.path.join(args.outdir, "short.mapped.bed")):
            os.remove(os.path.join(args.outdir, "short.mapped.bed"))

        # -f to force if file already exists
        pysam.merge("-f", os.path.join(args.outdir, "mapped.bam"),
                    os.path.join(args.outdir, "long.mapped.bam"),
                    os.path.join(args.outdir, "short.mapped.bam"))
        if os.path.exists(os.path.join(args.outdir, "long.mapped.bam")):
            os.remove(os.path.join(args.outdir, "long.mapped.bam"))
        if os.path.exists(os.path.join(args.outdir, "short.mapped.bam")):
            os.remove(os.path.join(args.outdir, "short.mapped.bam"))

        logging.info("| END: merge both long and short alignments")
    else:
        sys.stderr.write("Unknown aligner!! Currently suppoted aligners: BWA-mem and CLAN")
        sys.exit(1)
