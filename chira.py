    #!/usr/local/bin/python
# coding=utf-8
import utilities
import quantifier
import multisample
import stats
import argparse
import os
import sys
from collections import defaultdict
import pysam
import multiprocessing
from multiprocessing import Process, Value, Lock, Pool, Manager
from functools import partial
import logging
import pkg_resources
import itertools
import datetime
import re
import pprint
sys.setrecursionlimit(10000)

d_gene_annotations = defaultdict(lambda: defaultdict(str))
d_transcript_annotations = defaultdict(lambda: defaultdict())


class ChiRA:

    def __init__(self, r1sample=None, r2sample=None, outputdir=None, index1=None, index2=None, f_gff=None,
                 f_bed=None, threads=None, f_goterms=None, f_mirfamilies=None, stranded=None, seed_length1=None,
                 seed_length2=None, align_score1=None, align_score2=None, min_looplen=None, tpm_cutoff=None,
                 score_cutoff=None, merge_overlap=None, crg_share_threshold=None, min_locus_size=None, em_threshold=None,
                 percent_segment_overlap=None, hybridize=None, current_matepair=None):
        self.r1sample = r1sample
        self.r2sample = r2sample
        self.outputdir = outputdir
        self.index1 = index1
        self.index2 = index2
        self.f_gff = f_gff
        self.f_bed = f_bed
        self.f_goterms = f_goterms
        self.f_mirfamilies = f_mirfamilies
        self.threads = threads
        self.stranded = stranded
        self.seed_length1 = seed_length1
        self.seed_length2 = seed_length2
        self.align_score1 = align_score1
        self.align_score2 = align_score2
        self.min_looplen = min_looplen
        self.tpm_cutoff = tpm_cutoff
        self.score_cutoff = score_cutoff
        self.merge_overlap = merge_overlap
        self.crg_share_threshold = crg_share_threshold
        self.min_locus_size = min_locus_size
        self.em_threshold = em_threshold
        self.percent_segment_overlap = percent_segment_overlap
        self.hybridize = hybridize
        self.current_matepair = current_matepair

    def map_reads_bwa(self, align_type, index):
        """
            Funtion that maps the reads to the miRNAs or transcriptome. Different parameters
            are used for miRNA and transcriptome mapping.
        """

        logging.info("Mapping " + align_type + "reads segments to " + index)

        query_fasta = refindex = None
        matepairdir = os.path.join(self.outputdir, self.current_matepair)
        bam = os.path.join(matepairdir, index + "." + align_type + ".bam")

        if align_type == "long":
            if self.current_matepair == "R1":
                query_fasta = self.r1sample
            elif self.current_matepair == "R2":
                query_fasta = self.r2sample
            seed_length = self.seed_length1
            align_score = self.align_score1
            mismatch_score = "3"
            n_aligns = "20"
            gap_o = "5"
            gap_e = "3"
        elif align_type == "short":
            query_fasta = os.path.join(matepairdir, "long.unmapped.fa")
            seed_length = self.seed_length2
            align_score = self.align_score2
            mismatch_score = "100"
            n_aligns = "50"
            gap_o = "100"
            gap_e = "100"

        if index == "index1":
            refindex = self.index1
        elif index == "index2":
            refindex = self.index2

        bwa_params = ["-r 1",                       # look for internal seeds inside a seed longer than {-k} * {-r}
                      "-A 1",                       # match score
                      "-B " + str(mismatch_score),  # mismatch penalty
                      "-O " + gap_o,                # gap open penalty
                      "-E " + gap_e,                # gap extension penalty
                      "-L 0",                       # clipping penalty, we need soft clips
                      "-Y",                         # use soft clipping for supplementary alignments
                      # "-h " + n_aligns,             # if there're -h hits with score >80% of the maxscore, output in XA
                      # "-a ",
                      "-k " + str(seed_length),     # minimum seed length
                      "-T " + str(align_score),     # minimum alignment score
                      "-t " + str(self.threads),
                      refindex,
                      query_fasta
                      ]
        bwacall = ("bwa mem " + " ".join(bwa_params) + " | samtools view -hb - > " + bam)
        print(bwacall)
        os.system(bwacall)

    def map_reads(self, align_type, index):
        """
            Funtion that maps the reads to the miRNAs or transcriptome. Different parameters
            are used for miRNA and transcriptome mapping.
        """

        logging.info("| STRAT: map " + align_type + " reads to " + index)

        query_fasta = refindex = None
        matepairdir = os.path.join(self.outputdir, self.current_matepair)
        bam = os.path.join(matepairdir, index + "." + align_type + ".bam")

        if align_type == "long":
            if self.current_matepair == "R1":
                query_fasta = self.r1sample
            elif self.current_matepair == "R2":
                query_fasta = self.r2sample
            seed_length = self.seed_length1
            align_score = self.align_score1
            mismatch_penalty = "2"
            n_aligns = "50"
            gap_o = "5"
            gap_e = "3"
            seed_mismatches = "1"
            seed_interval = str(int(seed_length * 0.5))
            # gap_o = "100"
            # gap_e = "100"
            # seed_mismatches = "0"
            # seed_interval = str(int(seed_length * 0.2))
        elif align_type == "short":
            query_fasta = os.path.join(matepairdir, "long.unmapped.fa")
            seed_length = self.seed_length2
            align_score = self.align_score2
            mismatch_penalty = "100"
            n_aligns = "100"
            gap_o = "100"
            gap_e = "100"
            seed_mismatches = "0"
            seed_interval = str(max(2, int(seed_length * 0.3)))
        if index == "index1":
            refindex = self.index1
        elif index == "index2":
            refindex = self.index2

        # bwa_params = ["-r 1",                       # look for internal seeds inside a seed longer than {-k} * {-r}
        #               "-A 1",                       # match score
        #               "-B " + str(mismatch_penalty),  # mismatch penalty
        #               "-O " + gap_o,                # gap open penalty
        #               "-E " + gap_e,                # gap extension penalty
        #               "-L 0",                       # clipping penalty, we need soft clips
        #               "-Y",                         # use soft clipping for supplementary alignments
        #               # "-h " + n_aligns,             # if there're -h hits with score >80% of the maxscore, output in XA
        #               "-a ",
        #               "-k " + str(seed_length),     # minimum seed length
        #               "-T " + str(align_score),     # minimum alignment score
        #               "-t " + str(self.threads),
        #               refindex,
        #               query_fasta
        #               ]
        # bwacall = ("bwa mem " + " ".join(bwa_params) + " | samtools view -hb - > " + bam)
        # print(bwacall)
        # os.system(bwacall)

        # minimap2_params = ["-A 1",
        #                    "-B " + str(mismatch_penalty),  # mismatch penalty
        #                    "-O " + gap_o,                # gap open penalty
        #                    "-E " + gap_e,                # gap extension penalty
        #                    "-a",
        #                    "-Y",
        #                    "--seed " + str(seed_length),     # minimum seed length
        #                    "-N " + str(n_aligns),
        #                    "-t " + str(self.threads),
        #                    "--sr",
        #                    "-s " + str(align_score),
        #                    # "-m " + str(align_score),
        #                    "--for-only",
        #                    "--frag yes",
        #                    refindex,
        #                    query_fasta
        #                    ]
        # minimap2call = ("minimap2 " + " ".join(minimap2_params) + " | samtools view -hb - > " + bam)
        # print(minimap2call)
        # os.system(minimap2call)

        if self.stranded == "fw":
            strand_option = "--norc"
        if self.stranded == "rc":
            strand_option = "--nofw"
        if self.stranded == "both":
            strand_option = ""
        bowtie2_params = ["--local",
                          "--ignore-quals",
                          strand_option,
                          "-k " + str(n_aligns),
                          "--gbar 1",
                          "-D 50",
                          "-R 10",
                          "-L " + str(seed_length),
                          "-N " + seed_mismatches,
                          "--ma 1",
                          "--np 0",
                          "--mp " + mismatch_penalty,
                          "--rdg " + str(gap_o) + "," + str(gap_e),
                          "--rfg " + str(gap_o) + "," + str(gap_e),
                          "-i L," + seed_interval + ",0",
                          "--score-min C," + str(align_score) + ",0",
                          "-p " + str(self.threads),
                          "-x " + refindex,
                          "-f -U " + query_fasta,
                          ]
        bowtie2call = ("bowtie2 " + " ".join(bowtie2_params) + " | samtools view -hb - > " + bam)
        print(bowtie2call)
        os.system(bowtie2call)

        logging.info("| END: map " + align_type + " reads to " + index)

    def merge_bams(self, align_type):
        logging.info("| START: merge and sort " + align_type + " alignments of both indices")
        matepairdir = os.path.join(self.outputdir, self.current_matepair)
        # -f to force if file already exists
        pysam.merge("-f",  os.path.join(matepairdir, align_type + ".unsorted.bam"),
                    os.path.join(matepairdir, "index1." + align_type + ".bam"),
                    os.path.join(matepairdir, "index2." + align_type + ".bam"))
        pysam.sort("-m", "1G", "-@", str(self.threads),
                   os.path.join(matepairdir, align_type + ".unsorted.bam"),
                   "-T", os.path.join(matepairdir, align_type),
                   "-o", os.path.join(matepairdir, align_type + ".bam"))
        pysam.index(os.path.join(matepairdir, align_type + ".bam"))
        # os.system("rm " + os.path.join(matepairdir, "index1." + align_type + ".bam")
        #           + " " + os.path.join(matepairdir, "index2." + align_type + ".bam"))

        logging.info("| END: merge and sort " + align_type + " alignments both indices")

    def extract_unmapped_bwa(self, align_type):
        matepairdir = os.path.join(self.outputdir, self.current_matepair)
        fh_bam = pysam.Samfile(os.path.join(matepairdir, align_type + ".bam"), "rb")
        fh_mapped_bam = pysam.Samfile(os.path.join(matepairdir, align_type + ".mapped.bam"), "wb",
                                      template=fh_bam)
        fh_mapped_bed = open(os.path.join(matepairdir, align_type + ".mapped.bed"), "w")
        # [read id -> sequence] which are unmapped or not mapped on desired strand
        d_unmapped_reads = defaultdict(str)
        d_short_singleton_reads = defaultdict(list)
        l_properly_mapped_reads = []
        # list of reads that are in d_unmapped_reads but have at least one alignment in desired strand
        d_mapped_reads = defaultdict(str)
        for alignment in fh_bam.fetch(until_eof=True):
            readid = alignment.query_name
            readseq = alignment.query_sequence
            read_start = alignment.query_alignment_start
            read_end = alignment.query_alignment_end

            # if it is unmapped or mapped to wrong strand, consider it as unmapped
            if alignment.is_unmapped or \
                    (self.stranded == "fw" and alignment.is_reverse) or \
                    (self.stranded == "rc" and not alignment.is_reverse):
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
                bedentry = '\t'.join([fh_bam.getrname(alignment.tid),
                                      str(alignment.reference_start),
                                      str(alignment.reference_end),
                                      ','.join([readid,
                                                fh_bam.getrname(alignment.tid),  # reference id
                                                str(alignment.reference_start),
                                                str(alignment.reference_end),
                                                "-" if alignment.is_reverse else "+",
                                                alignment.cigarstring]),
                                      "0",
                                      "-" if alignment.is_reverse else "+"])
                fh_mapped_bed.write(bedentry + "\n")

                # mark as unmapped if there is only one mapped segment that doesn't cover the entire read
                softclipped = False
                for (cigarop, cigarlength) in alignment.cigartuples:
                    if int(cigarop) == 4 and int(cigarlength) >= self.align_score2:
                        softclipped = True
                        break
                properly_mapped = False
                if softclipped and readid not in l_properly_mapped_reads:
                    for prev_pos in d_short_singleton_reads[readid]:
                        # <= 0 indicates no overlap between read segments,
                        # a second non-overlapping soft clipped alignment found
                        if max(0, min(prev_pos[1], read_end) - max(prev_pos[0], read_start)) == 0:
                                properly_mapped = True
                                break
                # it not at least one softclipped position, it is a long properly mapped singleton
                else:
                    properly_mapped = True

                if properly_mapped:
                    if readid not in l_properly_mapped_reads:
                        l_properly_mapped_reads.append(readid)
                    if readid in d_short_singleton_reads:
                        del d_short_singleton_reads[readid]
                else:
                    if [read_start, read_end] not in d_short_singleton_reads[readid]:
                        d_short_singleton_reads[readid].append([read_start, read_end])
        fh_bam.close()
        fh_mapped_bam.close()
        fh_mapped_bed.close()

        # write all unmapped sequences to a fasta file
        with open(os.path.join(matepairdir, align_type + ".unmapped.fa"), "w") as fh_unmapped_fasta:
            for readid, readseq in d_unmapped_reads.items():
                fh_unmapped_fasta.write(">" + readid + "\n" + readseq + "\n")
            for readid in d_short_singleton_reads.keys():
                fh_unmapped_fasta.write(">" + readid + "\n" + d_mapped_reads[readid] + "\n")
        logging.info("Mapping done")
        return

    def chunks(self, a, n):
        return (a[i::n] for i in range(n))

    def process_alignments(self, readids, d_mapped_read_aligns, d_short_singleton_reads):
        d_temp = defaultdict(list)
        l_properly_mapped_reads = []
        for readid in readids:
            for alignment in d_mapped_read_aligns[readid]:
                read_start = alignment.query_alignment_start
                read_end = alignment.query_alignment_end
                # mark as unmapped if there is only one mapped segment that doesn't cover the entire read
                softclipped = False
                for (cigarop, cigarlength) in alignment.cigartuples:
                    if int(cigarop) == 4 and int(cigarlength) >= self.align_score2:
                        softclipped = True
                        break
                properly_mapped = False
                if softclipped and readid not in l_properly_mapped_reads:
                    if readid in d_temp:
                        for prev_pos in d_temp[readid]:
                            # <= 0 indicates no overlap between read segments,
                            # a second non-overlapping soft clipped alignment found
                            if max(0, min(prev_pos[1], read_end) - max(prev_pos[0], read_start)) == 0:
                                properly_mapped = True
                                break
                # it not at least one softclipped position, it is a long properly mapped singleton
                else:
                    properly_mapped = True

                if properly_mapped:
                    if readid not in l_properly_mapped_reads:
                        l_properly_mapped_reads.append(readid)
                    if readid in d_temp:
                        del d_temp[readid]
                else:
                    if [read_start, read_end] not in d_temp[readid]:
                        d_temp[readid].append([read_start, read_end])
        d_short_singleton_reads.update(d_temp)

    def extract_unmapped(self, align_type):
        logging.info("| START: extract unmapped reads of " + align_type + " alignments")
        matepairdir = os.path.join(self.outputdir, self.current_matepair)
        fh_bam = pysam.Samfile(os.path.join(matepairdir, align_type + ".bam"), "rb")
        fh_mapped_bam = pysam.Samfile(os.path.join(matepairdir, align_type + ".mapped.bam"), "wb",
                                      template=fh_bam)
        fh_mapped_bed = open(os.path.join(matepairdir, align_type + ".mapped.bed"), "w")
        d_mapped_reads = defaultdict(str)
        d_unmapped_reads = defaultdict(str)
        d_mapped_read_aligns = defaultdict(list)
        for alignment in fh_bam.fetch(until_eof=True):
            readid = alignment.query_name
            readseq = alignment.query_sequence
            if alignment.is_unmapped or \
                    (self.stranded == "fw" and alignment.is_reverse) or \
                    (self.stranded == "rc" and not alignment.is_reverse):
                if readid not in d_mapped_reads:
                    if readseq:
                        d_unmapped_reads[readid] = readseq
                    else:
                        d_unmapped_reads[readid] = d_mapped_reads[readid]
                continue
            else:
                if readid in d_unmapped_reads:
                    del d_unmapped_reads[readid]
            d_mapped_reads[readid] = readseq
            d_mapped_read_aligns[readid].append(alignment)
            fh_mapped_bam.write(alignment)
            bedentry = '\t'.join([fh_bam.getrname(alignment.tid),
                                  str(alignment.reference_start),
                                  str(alignment.reference_end),
                                  ','.join([readid,
                                            fh_bam.getrname(alignment.tid),  # reference id
                                            str(alignment.reference_start),
                                            str(alignment.reference_end),
                                            "-" if alignment.is_reverse else "+",
                                            alignment.cigarstring]),
                                  "0",
                                  "-" if alignment.is_reverse else "+"])
            fh_mapped_bed.write(bedentry + "\n")
        fh_bam.close()
        fh_mapped_bam.close()
        fh_mapped_bed.close()

        m = Manager()
        # [read id -> sequence] which are unmapped or not mapped on desired strand
        d_short_singleton_reads = m.dict()

        processes = []
        print(str(datetime.datetime.now()), " START: multiprocessing")
        for readids in self.chunks(list(d_mapped_read_aligns.keys()), self.threads):
            p = Process(target=self.process_alignments,
                        args=(readids, d_mapped_read_aligns, d_short_singleton_reads))
            processes.append(p)
            p.start()
        for process in processes:
            process.join()
        print(str(datetime.datetime.now()), " END: multiprocessing")

        # write all unmapped sequences to a fasta file
        with open(os.path.join(matepairdir, align_type + ".unmapped.fa"), "w") as fh_unmapped_fasta:
            for readid, readseq in d_unmapped_reads.items():
                fh_unmapped_fasta.write(">" + readid + "\n" + readseq + "\n")
            for readid in d_short_singleton_reads.keys():
                fh_unmapped_fasta.write(">" + readid + "\n" + d_mapped_reads[readid] + "\n")
        logging.info("| END: extract unmapped reads of " + align_type + " alignments")
        return

    def reads_to_segments(self):
        matepairdir = os.path.join(self.outputdir, self.current_matepair)
        d_read_segments = defaultdict(lambda: defaultdict(set))
        with open(os.path.join(matepairdir, "mapped.bed")) as fh_in:
            with open(os.path.join(matepairdir, "mapped.segments.bed"), "w") as fh_out:
                for line in fh_in:
                    f = line.split("\t")
                    d = f[3].split(",")
                    readid = d[0]
                    current_cigar = d[-1]
                    current_match_start, current_match_end = utilities.match_positions(current_cigar)
                    similar_cigar_found = False
                    current_pos = current_match_start, current_match_end
                    matched_segment = ""
                    for segmentid in d_read_segments[readid]:
                        if current_pos in d_read_segments[readid][segmentid]:
                            similar_cigar_found = True
                            matched_segment = segmentid
                            break
                        l_read_pos = sorted(d_read_segments[readid][segmentid])
                        # check with the first position in this segment
                        previous_match_start = l_read_pos[0][0]
                        previous_match_end = l_read_pos[0][1]
                        overlap = max(0, min(previous_match_end, current_match_end)
                                      - max(previous_match_start, current_match_start))
                        overlap_percent1 = overlap / float(previous_match_end - previous_match_start)
                        overlap_percent2 = overlap / float(current_match_end - current_match_start)
                        if overlap_percent1 < self.percent_segment_overlap or \
                                overlap_percent2 < self.percent_segment_overlap:
                            continue

                        # check with the last position in this segment
                        previous_match_start = l_read_pos[-1][0]
                        previous_match_end = l_read_pos[-1][-1]
                        overlap = max(0, min(previous_match_end, current_match_end)
                                      - max(previous_match_start, current_match_start))
                        overlap_percent1 = overlap / float(previous_match_end - previous_match_start)
                        overlap_percent2 = overlap / float(current_match_end - current_match_start)
                        if overlap_percent1 >= self.percent_segment_overlap and \
                                overlap_percent2 >= self.percent_segment_overlap:
                            similar_cigar_found = True
                            matched_segment = segmentid
                            break
                    if not similar_cigar_found:
                        if d_read_segments[readid]:
                            matched_segment = max(d_read_segments[readid]) + 1
                        else:
                            matched_segment = 1
                    if current_pos not in d_read_segments[readid][matched_segment]:
                        d_read_segments[readid][matched_segment].add(current_pos)
                    fh_out.write("\t".join(f[0:3]) + "\t" +
                                 readid + "|" + str(matched_segment) + "," +
                                 ",".join(d[1:]) + "\t" +
                                 "\t".join(f[4:]))

    def merge_loci_and_quantify(self):
        matepairdir = os.path.join(self.outputdir, self.current_matepair)
        file_prefix = os.path.join(matepairdir, "mapped.segments")
        if self.f_gff:
            bed = file_prefix + ".genomic.bed"
        else:
            bed = file_prefix + ".bed"
        loci_groups_file = file_prefix + ".loci.groups.txt"
        loci_groups_counts_file = file_prefix + ".loci.groups.counts"

        # Merge the loci in the first part
        logging.info("| START: merge and create Common Read Loci (CRL)")
        quantifier.create_crgs(bed,
                               self.merge_overlap,
                               loci_groups_file,
                               self.crg_share_threshold,
                               self.min_locus_size)
        logging.info("| END: merge and create Common Read Loci (CRL)")
        logging.info("| START: quantify CRLs")
        d_read_group_fractions, d_group_tpms = quantifier.quantify_crgs(loci_groups_file, self.em_threshold)
        with open(loci_groups_file) as fh_in:
            with open(loci_groups_counts_file, "w") as fh_out:
                for line in fh_in:
                    f = line.rstrip("\n").split("\t")
                    readid = f[0]
                    groupid = f[3]
                    fh_out.write("\t".join([line.strip("\n"),
                                            "{:.6g}".format(d_read_group_fractions[readid][groupid]),
                                            "{:.6g}".format(d_group_tpms[groupid])]) + "\n")
        logging.info("| END: quantify CRLs")
        return

    def write_chimeras(self, readids, d_read_lines, filename):
        d_regions = defaultdict()
        fh_out = open(filename, "w")
        for readid in readids:
            loci = d_read_lines[readid]
            d_duplicate_pairs = {}
            line_pairs = itertools.combinations(loci, 2)
            for line1, line2 in line_pairs:
                [segmentid1, transcriptid1, locusid1, groupid1, tx_pos_start1, tx_pos_end1, tx_pos_strand1,
                 cigar1, genomic_pos1, locuspos1, locusshare1, prob1, tpm1] = line1.rstrip('\n').split('\t')
                [segmentid2, transcriptid2, locusid2, groupid2, tx_pos_start2, tx_pos_end2, tx_pos_strand2,
                 cigar2, genomic_pos2, locuspos2, locusshare2, prob2, tpm2] = line2.rstrip('\n').split('\t')
                # these are multimappings of the same segment
                if segmentid1 == segmentid2:
                    continue
                if utilities.cigars_overlapping(cigar1, cigar2):
                    continue
                if locuspos1 == locuspos2:
                    continue

                first_locus_score = float(prob1) * float(locusshare1)
                second_locus_score = float(prob2) * float(locusshare2)
                combined_score = first_locus_score * second_locus_score

                # check for duplicates
                if ",".join([readid,
                             genomic_pos1,
                             genomic_pos2,
                             transcriptid1,
                             transcriptid2,
                             str(combined_score)]) in d_duplicate_pairs or \
                        ",".join([readid,
                                  genomic_pos2,
                                  genomic_pos1,
                                  transcriptid2,
                                  transcriptid1,
                                  str(combined_score)]) in d_duplicate_pairs:
                    continue
                d_duplicate_pairs[",".join([readid,
                                            genomic_pos1,
                                            genomic_pos2,
                                            transcriptid1,
                                            transcriptid2,
                                            str(combined_score)])] = 1

                geneid1 = d_transcript_annotations['gid'][transcriptid1] if transcriptid1 in \
                                                                            d_transcript_annotations[
                                                                                'gid'] else 'NA'
                geneid2 = d_transcript_annotations['gid'][transcriptid2] if transcriptid2 in \
                                                                            d_transcript_annotations[
                                                                                'gid'] else 'NA'
                name1 = d_gene_annotations['name'][geneid1] if geneid1 in d_gene_annotations['name'] else 'NA'
                name2 = d_gene_annotations['name'][geneid2] if geneid2 in d_gene_annotations['name'] else 'NA'
                type1 = d_gene_annotations['type'][geneid1] if geneid1 in d_gene_annotations['type'] else 'NA'
                type2 = d_gene_annotations['type'][geneid2] if geneid2 in d_gene_annotations['type'] else 'NA'

                if type1 == 'miRNA' or type1 == 'tRNA':
                    if geneid1 in d_gene_annotations['family']:
                        name1 = d_gene_annotations['family'][geneid1]
                if type2 == 'miRNA' or type2 == 'tRNA':
                    if geneid2 in d_gene_annotations['family']:
                        name2 = d_gene_annotations['family'][geneid2]

                if transcriptid1 + '\t' + genomic_pos1 not in d_regions:
                    if self.f_gff:
                        region1 = utilities.guess_region(transcriptid1, genomic_pos1, d_transcript_annotations)
                    else:
                        region1 = "NA"
                    d_regions[transcriptid1 + '\t' + genomic_pos1] = region1
                if transcriptid2 + '\t' + genomic_pos2 not in d_regions:
                    if self.f_gff:
                        region2 = utilities.guess_region(transcriptid2, genomic_pos2, d_transcript_annotations)
                    else:
                        region2 = "NA"
                    d_regions[transcriptid2 + '\t' + genomic_pos2] = region2
                length1 = d_transcript_annotations['len'][transcriptid1] if transcriptid1 in d_transcript_annotations[
                    'len'] else 'NA'
                length2 = d_transcript_annotations['len'][transcriptid2] if transcriptid2 in d_transcript_annotations[
                    'len'] else 'NA'

                interaction = "\t".join([readid,
                                         transcriptid1,
                                         transcriptid2,
                                         geneid1,
                                         geneid2,
                                         name1,
                                         name2,
                                         type1,
                                         type2,
                                         d_regions[transcriptid1 + '\t' + genomic_pos1],
                                         d_regions[transcriptid2 + '\t' + genomic_pos2],
                                         str(tx_pos_start1),
                                         str(tx_pos_end1),
                                         tx_pos_strand1,
                                         str(length1),
                                         str(tx_pos_start2),
                                         str(tx_pos_end2),
                                         tx_pos_strand2,
                                         str(length2),
                                         genomic_pos1,
                                         genomic_pos2,
                                         cigar1,
                                         cigar2,
                                         locusid1,
                                         locusid2,
                                         locuspos1,
                                         locuspos2,
                                         groupid1,
                                         groupid2,
                                         str(float(tpm1) * float(locusshare1)),
                                         str(float(tpm2) * float(locusshare2)),
                                         str(first_locus_score),
                                         str(second_locus_score),
                                         str(combined_score)])
                fh_out.write(interaction + "\n")
        fh_out.close()

    def write_loci_pairs(self):
        matepairdir = os.path.join(self.outputdir, self.current_matepair)

        loci_groups_counts_file = os.path.join(matepairdir, "mapped.segments.loci.groups.counts")
        locipairs = os.path.join(matepairdir, "loci.pairs")

        d_read_lines = utilities.parse_counts_file(loci_groups_counts_file, self.tpm_cutoff, self.score_cutoff)

        processes = []
        print(str(datetime.datetime.now()), " START: multiprocessing")
        for i,  readids in enumerate(self.chunks(list(d_read_lines.keys()), self.threads)):
            p = Process(target=self.write_chimeras,
                        args=(readids, d_read_lines, locipairs + "." + str(i)))
            processes.append(p)
            p.start()
        for process in processes:
            process.join()

        # cleanup intermediate files
        with open(locipairs, 'w') as outfile:
            for i in range(self.threads):
                with open(locipairs + "." + str(i)) as infile:
                    for line in infile:
                        outfile.write(line)
                os.remove(locipairs + "." + str(i))

        return

    def hybridize_pairs(self, firstpart, secondpart):
        first_prefix, pairs_prefix, second_prefix = utilities.file_prefixes(os.path.join(self.outputdir, self.current_matepair), firstpart, secondpart)
        bed1 = pairs_prefix + "." + firstpart + ".bed"
        bed2 = pairs_prefix + "." + secondpart + ".bed"
        fasta1 = pairs_prefix + "." + firstpart + ".fa"
        fasta2 = pairs_prefix + "." + secondpart + ".fa"
        utilities.pairspos_to_bed(pairs_prefix, bed1, bed2)
        os.system(p_fastafrombed + " -s -name -fi " + f_genomicfasta + " -bed " + bed1 + " -fo " + fasta1)
        os.system(p_fastafrombed + " -s -name -fi " + f_genomicfasta + " -bed " + bed2 + " -fo " + fasta2)
        d_seq_records1 = utilities.fasta_to_dict(fasta1)
        d_seq_records2 = utilities.fasta_to_dict(fasta2)
        l_seqids = d_seq_records1.keys()
        chunk_size = len(l_seqids)/(self.threads-1)
        slices = utilities.chunks(l_seqids, chunk_size)
        print(len(l_seqids), str(chunk_size), len(slices))
        manager = multiprocessing.Manager()
        jobs = []
        for i, s in enumerate(slices):
            hybrid_file = pairs_prefix + ".hybrid."+str(i)
            print(len(l_seqids), len(s), hybrid_file)
            j = multiprocessing.Process(target=utilities.run_hybridmin,
                                        args=(d_seq_records1, d_seq_records2, s, hybrid_file, pairs_prefix, p_hybmin))
            jobs.append(j)
        for j in jobs:
            j.start()
        for j in jobs:
            j.join()
        os.system("cat " + pairs_prefix + ".hybrid.* > " + pairs_prefix + ".hybrids")
        os.system("rm " + pairs_prefix + ".hybrid.*")
        return

    def write_hybrids(self, firstpart, secondpart):
        first_prefix, pairs_prefix, second_prefix = utilities.file_prefixes(os.path.join(self.outputdir, self.current_matepair), firstpart, secondpart)
        locipairs = pairs_prefix + ".loci.pairs"
        hybrids = pairs_prefix + ".hybrids"
        interactome = pairs_prefix + ".interactome"

        d_hybrids = {}
        with open(hybrids) as fh_hybrids:
            for line in fh_hybrids:
                f = line.split("\t")
                d_hybrids[f[0]] = "\t".join(f[1:])

        fh_interactome = open(interactome, "w")
        with open(locipairs) as fh_locipairs:
            for line in fh_locipairs:
                fh_interactome.write(line.rstrip("\n") + "\t" + d_hybrids[line.split("\t")[0]])
        fh_interactome.close()

    def write_stats(self):
        d_mapped_readcounts_mirna = defaultdict()
        d_mapped_readcounts_other = defaultdict()
        d_unmapped_tagcounts = defaultdict(int)
        d_unmapped_readcounts = defaultdict(int)
        d_initial_tagcounts = defaultdict(int)

        matepairs = ['R1']
        # count whether a tag is seen each replicate
        d_initial_tagcounts = stats.tagcounts(self.r1sample, d_initial_tagcounts)
        if self.r2sample:
            matepairs.append('R2')
            d_initial_tagcounts = stats.tagcounts(self.r2sample, d_initial_tagcounts)

        for matepair in matepairs:
            self.current_matepair = matepair
            replicate_dir = os.path.join(self.outputdir, self.current_matepair)
            d_mapped_readcounts_mirna = stats.readcounts(
                os.path.join(replicate_dir, 'vs_miRNA.temp.bam'),
                d_mapped_readcounts_mirna)
            d_mapped_readcounts_other = stats.readcounts(
                os.path.join(replicate_dir, 'vs_transcriptome.temp.bam'),
                d_mapped_readcounts_other)
            d_unmapped_tagcounts = stats.unmappedcounts(
                os.path.join(replicate_dir, 'vs_transcriptome.unmapped.fa'),
                d_unmapped_tagcounts)

            d_hybrid_readcounts = stats.hybrid_readcounts(os.path.join(replicate_dir, "miRNA_vs_miRNA.loci.pairs"))
            d_hybrid_readcounts = stats.hybrid_readcounts(
                os.path.join(replicate_dir, "miRNA_vs_transcriptome.loci.pairs"))
            d_hybrid_readcounts = stats.hybrid_readcounts(
                os.path.join(replicate_dir, "transcriptom_vs_transcriptome.loci.pairs"))
            # d_hybrid_abundances = stats.hybrid_abundances(pairs_prefix + ".loci.pairs")

        n_mapped_readcount_mirna = sum(d_mapped_readcounts_mirna.values())
        n_mapped_readcount_other = sum(d_mapped_readcounts_other.values())
        diff = set(d_mapped_readcounts_mirna.keys()) - set(d_mapped_readcounts_other.keys())

        n_unmapped_readcount = 0
        n_initial_readcount = 0
        for tagid in d_initial_tagcounts:
            n_initial_readcount += int(tagid.split('|')[1])
        for tagid in d_unmapped_tagcounts:
            # if a tag is unmapped in all replicates which it initially present in
            # then it is not mapped at all
            if d_unmapped_tagcounts[tagid] == d_initial_tagcounts[tagid]:
                n_unmapped_readcount += int(tagid.split('|')[1])
        print(n_initial_readcount, n_mapped_readcount_mirna, n_mapped_readcount_other, len(diff), n_unmapped_readcount)

    def parse_arguments(self):
        def score_float(x):
            x = float(x)
            if x < 0.0 or x > 1.0:
                raise argparse.ArgumentTypeError("%r not in range [0.0, 1.0]" % (x,))
            return x

        parser = argparse.ArgumentParser(description='Chimeric read annotator for interactome data',
                                         usage='%(prog)s [-h] [-v,--version]',
                                         formatter_class=argparse.ArgumentDefaultsHelpFormatter)

        parser.add_argument('-r1', '--matepair1', action='store', dest='r1sample', required=True,
                            metavar='', help='Sample path for read1')

        parser.add_argument('-r2', '--matepair2', action='store', dest='r2sample', required=False,
                            metavar='', help='Sample path for read2')

        parser.add_argument('-o', '--out', action='store', dest='outputdir', required=True, metavar='',
                            help='Output directory path for the whole analysis')

        parser.add_argument('-x1', '--index1', action='store', dest='index1', required=True,
                            metavar='', help='first prioroty index file')

        parser.add_argument('-x2', '--index2', action='store', dest='index2', required=False,
                            metavar='', help='second priority index file')

        parser.add_argument('-g', '--gff', action='store', dest='f_gff', required=False,
                            metavar='', help='annotation file in GFF')

        parser.add_argument('-b', '--bed', action='store', dest='f_bed', required=False,
                            metavar='', help='annotation file in BED')

        parser.add_argument('-go', '--go', action='store', dest='f_goterms', required=False,
                            metavar='', help='File containing go terms')

        parser.add_argument('-mf', '--family', action='store', dest='f_mirfamilies', required=False,
                            metavar='', help='File containing miRNA family information')

        parser.add_argument('-p', '--threads', action='store', type=int, default=1, metavar='',
                            dest='threads',
                            help='Number of threads to use')

        parser.add_argument("-s", '--stranded', type=str, choices=["fw", "rc", "both"], default='fw', metavar='',
                            dest='stranded',
                            help='''Strand-specificity of input samples.
                                 fw = map to transcript strand;
                                 rc = map to reverse compliment of transcript strand;
                                 both = try to map on both strnads''')

        parser.add_argument("-sl1", '--seed_length1', action='store', type=int, default=12, metavar='',
                            dest='seed_length1',
                            help='''Seed length for 1st mapping iteration.
                                    bwa-mem parameter "-k"''')

        parser.add_argument("-sl2", '--seed_length2', action='store', type=int, default=6, metavar='',
                            dest='seed_length2',
                            help='''Seed length for 2nd mapping iteration.
                                    bwa-mem parameter "-k"''')

        parser.add_argument("-as1", '--align_score1', action='store', type=int, default=18, metavar='',
                            dest='align_score1',
                            help='''Minimum alignment score in 1st mapping iteration.
                                    It must be smaller than --align_score1 parameter.  
                                    bwa-mem parameter "-T"''')

        parser.add_argument("-as2", '--align_score2', action='store', type=int, default=12, metavar='',
                            dest='align_score2',
                            help='''Minimum alignment score in 2nd mapping iteration.
                                    It must be smaller than --align_score1 parameter.  
                                    bwa-mem parameter "-T"''')

        parser.add_argument('-m', '--min_looplen', action='store', type=int, default=3, metavar='',
                            dest='min_looplen',
                            help='Set minimum loop length in case of single RNA duplex')

        parser.add_argument('-tc', '--tpm_cutoff', action='store', type=score_float, default=0.1, metavar='',
                            dest='tpm_cutoff',
                            help='Transcripts with less than this percentile TPMs will be discarded in '
                                 'the final output. [0-1.0]')

        parser.add_argument('-sc', '--score_cutoff', action='store', type=score_float, default=0.0, metavar='',
                            dest='score_cutoff',
                            help='Hybrids with less than this score will be discarded in the final output. [0-1.0]')

        parser.add_argument('-mo', '--merge_overlap', action='store', type=score_float, default=0.8, metavar='',
                            dest='merge_overlap',
                            help='Minimum percentage overlap among BED entries inorder to merge. [0-1.0]')

        parser.add_argument('-cs', '--crg_share_threshold', action='store', type=score_float, default=0.7, metavar='',
                            dest='crg_share_threshold',
                            help='Minimum percentage overlap of a locus on a CRG inorder to merge it '
                                 'into a CRG. [0-1.0]')

        parser.add_argument('-ls', '--min_locus_size', action='store', type=int, default=5, metavar='',
                            dest='min_locus_size',
                            help='Minimum number of reads a locus should have in order to participate in CRG creation'
                                 'Always set this value relative to your sequencing depth. Setting this to lower leads'
                                 'CRGs of random multimappings Also consider setting the --crg_share_threshold option '
                                 'along with this')

        parser.add_argument('-e', '--em_threshold', action='store', type=score_float, default=0.05, metavar='',
                            dest='em_threshold',
                            help='The maximum difference of transcripts expression between two consecutive iterations '
                                 'of EM algorithm to converge.')

        parser.add_argument('-so', '--segment_overlap', action='store', type=score_float, default=0.8, metavar='',
                            dest='percent_segment_overlap',
                            help='Matching read positions with greater than this %% overlap are merged into a segment')

        parser.add_argument("-f", '--hybridize', action='store_true', help="Hybridize the predicted chimeras")

        parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.0')

        args = parser.parse_args()
        self.r1sample = args.r1sample
        self.r2sample = args.r2sample
        self.outputdir = args.outputdir
        self.index1 = args.index1
        self.index2 = args.index2
        self.f_gff = args.f_gff
        self.f_bed = args.f_bed
        self.f_goterms = args.f_goterms
        self.f_mirfamilies = args.f_mirfamilies
        self.threads = args.threads
        self.stranded = args.stranded
        self.seed_length1 = args.seed_length1
        self.seed_length2 = args.seed_length2
        self.align_score1 = args.align_score1
        self.align_score2 = args.align_score2
        self.min_looplen = args.min_looplen
        self.tpm_cutoff = args.tpm_cutoff  # we need percentile later
        self.score_cutoff = args.score_cutoff
        self.merge_overlap = args.merge_overlap
        self.crg_share_threshold = args.crg_share_threshold
        self.min_locus_size = args.min_locus_size
        self.em_threshold = args.em_threshold
        self.percent_segment_overlap = args.percent_segment_overlap
        self.hybridize = args.hybridize

        print('Sample R1                        : ' + args.r1sample)
        if args.r2sample:
            print('Sample R2                        : ' + args.r2sample)
        print('Output directory                 : ' + args.outputdir)
        if args.f_gff:
            print('Annotation GFF file              : ' + args.f_gff)
        if args.f_bed:
            print('Annotation BED file              : ' + args.f_bed)
        # print('BWA index                        : ' + args.index)
        if args.index1:
            print('1st priority BWA index           : ' + args.index1)
        if args.index2:
            print('2nd priority BWA  index          : ' + args.index2)
        if args.f_goterms:
            print('GO terms file                : ' + args.f_goterms)
        if args.f_mirfamilies:
            print('miRNA families file          : ' + args.f_mirfamilies)
        print('Number of threads                : ' + str(args.threads))
        print('Stranded                         : ' + args.stranded)
        print('Minimum loop length              : ' + str(args.min_looplen))
        print('TPM cutoff                       : ' + str(args.tpm_cutoff))
        print('Score cutoff                     : ' + str(args.score_cutoff))
        print('Merge overlap                    : ' + str(args.merge_overlap))
        print('CRG share                        : ' + str(args.crg_share_threshold))
        print('EM threshold                     : ' + str(args.em_threshold))
        print('Suboptimal hybrid threshold      : ' + str(args.percent_segment_overlap))


if __name__ == "__main__":

    chira = ChiRA()
    chira.parse_arguments()
    print("===================================================================")

    if not os.path.exists(chira.outputdir):
        os.makedirs(chira.outputdir)

    if chira.f_gff:
        # Parse the annotations and save them to dictionaries. Additionally write exons bed files to the outputdir
        print("Parsing the annotation file")
        utilities.parse_annotations(d_gene_annotations, d_transcript_annotations, chira.f_goterms,
                                    chira.f_mirfamilies, chira.f_gff, chira.f_bed, chira.outputdir)

    logging.basicConfig(level=logging.INFO,
                        filename=os.path.join(chira.outputdir, 'run.log'),  # log to this file
                        filemode='a',
                        format='%(asctime)s %(message)s')

    matepairs = ['R1']
    if chira.r2sample:
        matepairs.append('R2')

    for matepair in matepairs:
        # set this for each matepair and use the global class variable
        chira.current_matepair = matepair
        matepair_dir = os.path.join(chira.outputdir, matepair)
        if not os.path.exists(matepair_dir):
            os.makedirs(matepair_dir)
        chira.map_reads("long", "index1")
        # TODO: handle cases where nothing maps to index1 or index2 and try to merge empty bam files
        if chira.index2:
            chira.map_reads("long", "index2")
            chira.merge_bams("long")
        else:
            logging.info("| START: sort long alignments")
            pysam.sort("-m", "1G", "-@", str(chira.threads),
                       os.path.join(matepair_dir, "index1.long.bam"),
                       "-T", os.path.join(matepair_dir, "index1.long"),
                       "-o", os.path.join(matepair_dir, "long.bam"))
            pysam.index(os.path.join(matepair_dir, "long.bam"))
            logging.info("| END: sort long alignments")
        print(str(datetime.datetime.now()), " START:extract unmapped long")
        chira.extract_unmapped("long")
        print(str(datetime.datetime.now()), " START:extract unmapped long")
        chira.map_reads("short", "index1")
        if chira.index2:
            chira.map_reads("short", "index2")
            chira.merge_bams("short")
        else:
            logging.info("| START: sort short alignments")
            pysam.sort("-m", "1G", "-@", str(chira.threads),
                       os.path.join(matepair_dir, "index1.short.bam"),
                       "-T", os.path.join(matepair_dir, "index1.short"),
                       "-o", os.path.join(matepair_dir, "short.bam"))
            pysam.index(os.path.join(matepair_dir, "short.bam"))
            logging.info("| END: sort short alignments")
        chira.extract_unmapped("short")

        logging.info("| START: merge both long and short alignments")
        # # -f to force if file already exists
        pysam.merge("-f", os.path.join(matepair_dir, "mapped.bam"),
                    os.path.join(matepair_dir, "long.mapped.bam"),
                    os.path.join(matepair_dir, "short.mapped.bam"))
        os.system("cat " + os.path.join(matepair_dir, "long.mapped.bed")
                  + " " + os.path.join(matepair_dir, "short.mapped.bed")
                  + ">" + os.path.join(matepair_dir, "mapped.bed"))
        logging.info("| END: merge both long and short alignments")

        logging.info("| START: find and label exact read segments")
        chira.reads_to_segments()
        logging.info("| END: find and label exact read segments")

        if chira.f_gff:
            print("converting to genomic coordinates")
            logging.info("| START: convert to genomic coordinates")
            utilities.transcript_to_genomic_pos(os.path.join(matepair_dir, "mapped.segments.bed"),
                                                os.path.join(chira.outputdir, "genomic_exons.bed"),
                                                os.path.join(chira.outputdir, "transcriptomic_exons.bed"))
            logging.info("| END: convert to genomic coordinates")
        chira.merge_loci_and_quantify()
        logging.info("| START: write final interactions to a file")
        chira.write_loci_pairs()
        logging.info("| END: write final interactions to a file")
