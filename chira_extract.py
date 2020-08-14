#!/usr/bin/env python
import chira_utilities
import argparse
import os
import sys
from collections import defaultdict
from multiprocessing import Process
import itertools
import datetime
from BCBio import GFF
from Bio import SeqIO
import subprocess
import math


d_gene_annotations = defaultdict(lambda: defaultdict(str))
d_transcript_annotations = defaultdict(lambda: defaultdict())


def strandardize(strand):
    if strand == '-1':
        strand = '-'
    elif strand == '1':
        strand = '+'
    return strand


def guess_region(transcriptid, read_pos):
    [read_chr, read_start, read_end, read_strand] = read_pos.split(':')[-4:]
    region = 'NA'
    overlap_length = 0
    utr_start = utr_end = first_cds_start = last_cds_end = strand = None
    read_strand = strandardize(read_strand)
    if transcriptid in d_transcript_annotations['UTR']:
        for pos in d_transcript_annotations['UTR'][transcriptid]:
            chrom = pos[0]
            start = int(pos[1])
            end = int(pos[2])
            strand = strandardize(pos[3])
            if read_chr != chrom or read_strand != strand:
                continue

            if chira_utilities.overlap([start, end], [int(read_start), int(read_end)]) > overlap_length:
                overlap_length = chira_utilities.overlap([start, end], [int(read_start), int(read_end)])
                region = 'UTR'
                utr_start = start
                utr_end = end

    if transcriptid in d_transcript_annotations['CDS']:
        for pos in d_transcript_annotations['CDS'][transcriptid]:
            chrom = pos[0]
            start = int(pos[1])
            end = int(pos[2])
            strand = strandardize(pos[3])
            if read_chr != chrom or read_strand != strand:
                continue
            if chira_utilities.overlap([start, end], [int(read_start), int(read_end)]) > overlap_length:
                overlap_length = chira_utilities.overlap([start, end], [int(read_start), int(read_end)])
                region = 'CDS'
            if not first_cds_start or start < first_cds_start:
                first_cds_start = start
            if not last_cds_end or end > last_cds_end:
                last_cds_end = end

    # if region is still a UTR
    if region == 'UTR' and utr_end <= first_cds_start:
        region = '5_prime_UTR'
        if strand == '-':
            region = '3_prime_UTR'
    elif region == 'UTR' and utr_start >= last_cds_end:
        region = '3_prime_UTR'
        if strand == '-':
            region = '5_prime_UTR'

    # works for mirbase gff3
    if transcriptid in d_transcript_annotations['gid']:
        geneid = d_transcript_annotations['gid'][transcriptid]
        if geneid in d_transcript_annotations['mature']:
            for pos in d_transcript_annotations['mature'][geneid]:
                chrom = pos[0]
                start = int(pos[1])
                end = int(pos[2])
                strand = strandardize(pos[3])
                name = pos[4]
                if read_chr != chrom and read_strand != strand:
                    continue
                if chira_utilities.overlap([start, end], [int(read_start), int(read_end)]) > overlap_length:
                    overlap_length = chira_utilities.overlap([start, end], [int(read_start), int(read_end)])
                    if name.endswith("-3p"):
                        region = "3p_mature_mir"
                    elif name.endswith("-5p"):
                        region = "5p_mature_mir"
                    else:
                        region = "mature_mir"
    return region


def hybridize_with_intarna(seq1, seq2, intarna_params):
    # assuming first sequence query and second one is target
    output = os.popen(" ".join(["IntaRNA", intarna_params, "-q", seq1, "-t", seq2])).read()
    dotbracket = energy = pos = "NA"
    for alignment in output.split("\n"):
        if alignment.startswith('target'):
            dotbracket = alignment.split(";")[3]
            # exchange the dot bracket notion of query and target
            target_dotbracket = dotbracket.split("&")[0].replace("(", ")")
            query_dotbracket = dotbracket.split("&")[1].replace(")", "(")
            dotbracket = query_dotbracket + "&" + target_dotbracket
            energy = alignment.split(";")[4]
            # target is start1, query is start2
            pos = alignment.split(";")[2] + "&" + alignment.split(";")[1]
            break
    return dotbracket, pos, energy


def add_locus_to_set(locus, l_loci_bed):
    pos = locus.split(":")
    bed = "\t".join([":".join(pos[0:-3]), pos[-3], pos[-2], locus, "1", pos[-1]])
    if bed not in l_loci_bed:
        l_loci_bed.add(bed)


def update_best_hits(l_best_hits, hit_type):
    best_tpm = 0
    for i, hit in enumerate(l_best_hits):
        if hit_type == "chimera":
            tpm = float(hit[24]) + float(hit[25])
        elif hit_type == "singleton":
            tpm = float(hit[13])
        if tpm > best_tpm:
            best_tpm = tpm
    l_best_hits = [n for n in l_best_hits if n is not None]
    for i, hit in enumerate(l_best_hits):
        if hit_type == "chimera":
            tpm = float(hit[24]) + float(hit[25])
        elif hit_type == "singleton":
            tpm = float(hit[13])
        if tpm < best_tpm:
            l_best_hits[i] = None
    l_best_hits = [n for n in l_best_hits if n is not None]
    return l_best_hits


def extract_annotations(transcriptid, genomic_pos, d_regions, f_gtf):
    geneid = d_transcript_annotations['gid'][transcriptid] \
        if transcriptid in d_transcript_annotations['gid'] else 'NA'
    name = d_gene_annotations['name'][geneid] if geneid in d_gene_annotations['name'] else 'NA'
    biotype = d_gene_annotations['type'][geneid] if geneid in d_gene_annotations['type'] else 'NA'
    if biotype == 'miRNA' or biotype == 'tRNA':
        if geneid in d_gene_annotations['family']:
            name = d_gene_annotations['family'][geneid]

    if transcriptid + '\t' + genomic_pos in d_regions:
        region = d_regions[transcriptid + '\t' + genomic_pos]
    else:
        region = "NA"
        if f_gtf:
            region = guess_region(transcriptid, genomic_pos)
        if region == "NA":
            region = biotype
        d_regions[transcriptid + '\t' + genomic_pos] = region

    tx_length = d_transcript_annotations['len'][transcriptid] if transcriptid in d_transcript_annotations[
        'len'] else 'NA'
    return geneid, name, region, tx_length


def filter_alignments(lines, tpm_threshold, score_cutoff):
    l_read_aligns = []
    for line in lines:
        f = line.rstrip('\n').split('\t')
        crl_tpm = f[12]
        prob = f[10]
        locusshare = f[11]
        if float(crl_tpm) < tpm_threshold:
            continue
        locus_score = float(prob) * float(locusshare)
        if locus_score < score_cutoff:
            continue
        l_read_aligns.append(line)
    return l_read_aligns


def extract_and_write(readid, l_read_alignments, l_loci_bed, d_ref_lengths1, d_ref_lengths2, f_gtf,
                      d_regions, chimeric_overlap, fh_chimeras, fh_singletons, hybridize):
    chimera_found = False
    l_best_chimeras = []
    l_best_singletons = []
    alignment_pairs = list(itertools.combinations(l_read_alignments, 2))

    for alignment1, alignment2 in alignment_pairs:
        [segmentid1, transcriptid1, locusid1, crlid1, tx_pos_start1, tx_pos_end1, tx_pos_strand1,
         cigar1, genomic_pos1, locuspos1, locusshare1, prob1, tpm1] = alignment1.rstrip('\n').split('\t')
        [segmentid2, transcriptid2, locusid2, crlid2, tx_pos_start2, tx_pos_end2, tx_pos_strand2,
         cigar2, genomic_pos2, locuspos2, locusshare2, prob2, tpm2] = alignment2.rstrip('\n').split('\t')
        # these are multimappings of the same segment
        if segmentid1 == segmentid2 or locuspos1 == locuspos2 or crlid1 == crlid2 or\
                float(prob1) == 0 or float(prob2) == 0:
            continue
        # check these are chimeric arms
        if not chira_utilities.is_chimeric(cigar1, cigar2, tx_pos_strand1 == "-",
                                           tx_pos_strand2 == "-", chimeric_overlap):
            continue
        switch_alignments = False
        # if the second reference is provided then it is assumed to be a split reference
        if len(d_ref_lengths2) != 0:
            # check if both transcripts are from the same reference database
            if transcriptid1 in d_ref_lengths1 and transcriptid2 in d_ref_lengths1 or \
                    transcriptid1 in d_ref_lengths2 and transcriptid2 in d_ref_lengths2:
                continue
            # switch orientation of the alignments based on references
            if transcriptid1 in d_ref_lengths2 or transcriptid2 in d_ref_lengths1:
                switch_alignments = True
        # not a split reference
        else:
            # then switch alignments based on the reference id
            if transcriptid2 > transcriptid1:
                switch_alignments = True
            elif locuspos2 > locuspos1:
                switch_alignments = True
        if switch_alignments:
            [segmentid1, transcriptid1, locusid1, crlid1, tx_pos_start1, tx_pos_end1, tx_pos_strand1,
             cigar1, genomic_pos1, locuspos1, locusshare1, prob1, tpm1] = alignment2.rstrip('\n').split(
                '\t')
            [segmentid2, transcriptid2, locusid2, crlid2, tx_pos_start2, tx_pos_end2, tx_pos_strand2,
             cigar2, genomic_pos2, locuspos2, locusshare2, prob2, tpm2] = alignment1.rstrip('\n').split(
                '\t')
        first_locus_score = float("{:.2f}".format(float(prob1)))  # * float(locusshare1)))
        second_locus_score = float("{:.2f}".format(float(prob2)))  # * float(locusshare2)))
        combined_score = first_locus_score * second_locus_score
        tpm1 = "{:.2f}".format(float(tpm1))
        tpm2 = "{:.2f}".format(float(tpm2))

        geneid1, name1, region1, tx_len1 = extract_annotations(transcriptid1,
                                                               genomic_pos1,
                                                               d_regions,
                                                               f_gtf)
        geneid2, name2, region2, tx_len2 = extract_annotations(transcriptid2,
                                                               genomic_pos2,
                                                               d_regions,
                                                               f_gtf)

        arm1_start, arm1_end = chira_utilities.match_positions(cigar1, tx_pos_strand1 == "-")
        arm2_start, arm2_end = chira_utilities.match_positions(cigar2, tx_pos_strand2 == "-")
        read_length = chira_utilities.query_length(cigar1, tx_pos_strand1 == "-")

        if tx_len1 == "NA" or tx_len2 == "NA":
            if len(d_ref_lengths2) == 0:
                tx_len1 = d_reflen1[transcriptid1]
                tx_len2 = d_reflen1[transcriptid2]
            else:
                try:
                    tx_len1 = d_reflen1[transcriptid1]
                    tx_len2 = d_reflen2[transcriptid2]
                except KeyError:
                    tx_len1 = d_reflen2[transcriptid1]
                    tx_len2 = d_reflen1[transcriptid2]

        chimera = [readid, transcriptid1, transcriptid2,
                   geneid1, geneid2, name1, name2, region1, region2,
                   str(tx_pos_start1), str(tx_pos_end1), tx_pos_strand1, str(tx_len1),
                   str(tx_pos_start2), str(tx_pos_end2), tx_pos_strand2, str(tx_len2),
                   ",".join([str(arm1_start), str(arm1_end), str(arm2_start), str(arm2_end), str(read_length)]),
                   genomic_pos1, genomic_pos2,
                   locuspos1, locuspos2,
                   crlid1, crlid2,
                   tpm1, tpm2,
                   str(first_locus_score), str(second_locus_score), str(combined_score)]
        chimera_found = True
        l_best_chimeras.append(chimera)

    # if there are no pairs, then it is a singleton
    if not chimera_found:
        # singleton read
        for alignment in l_read_alignments:
            [segmentid, transcriptid, locusid, crlid, tx_pos_start, tx_pos_end, tx_pos_strand,
             cigar, genomic_pos, locuspos, locusshare, prob, tpm] = alignment.rstrip('\n').split('\t')

            geneid, name, region, tx_len = extract_annotations(transcriptid,
                                                               genomic_pos,
                                                               d_regions,
                                                               f_gtf)

            locus_score = "{:.2f}".format(float(prob) * float(locusshare))
            tpm = "{:.2f}".format(float(tpm))

            arm_start, arm_end = chira_utilities.match_positions(cigar, tx_pos_strand == "-")
            read_length = chira_utilities.query_length(cigar, tx_pos_strand == "-")

            singleton = [readid, transcriptid, geneid, name, region,
                         str(tx_pos_start), str(tx_pos_end), tx_pos_strand, str(tx_len),
                         ",".join([str(arm_start), str(arm_end), str(read_length)]),
                         genomic_pos,
                         locuspos,
                         crlid,
                         tpm,
                         locus_score]
            l_best_singletons.append(singleton)

    if len(l_best_chimeras) > 0:
        l_best_chimeras = update_best_hits(l_best_chimeras, "chimera")
        for a in l_best_chimeras:
            if hybridize:
                add_locus_to_set(a[20], l_loci_bed)
                add_locus_to_set(a[21], l_loci_bed)
            else:
                a.append("NA")
                a.append("NA")
                a.append("NA")
                a.append("NA")
            fh_chimeras.write("\t".join(a) + "\n")
    else:
        l_best_singletons = update_best_hits(l_best_singletons, "singleton")
        for b in l_best_singletons:
            fh_singletons.write("\t".join(b) + "\n")


def write_chimeras(chunk_start, chunk_end, total_read_count, d_ref_lengths1, d_ref_lengths2, hybridize,
                   chimeric_overlap, f_gtf, outdir, crl_file, tpm_threshold, score_cutoff, n):
    d_regions = defaultdict()
    l_loci_bed = set()
    file_chimeras = os.path.join(outdir, "chimeras." + n)
    file_singletons = os.path.join(outdir, "singletons." + n)
    # make bed entry for extracing locus sequence and hybridizing
    with open(file_chimeras, "w") as fh_chimeras, open(file_singletons, "w") as fh_singletons, open(crl_file) as fh_crl:
        prev_readid = None
        l_readlines = []
        read_count = 0
        for line in fh_crl:
            f = line.rstrip('\n').split('\t')
            # last field after | represents the segment id, rest of the string before is read id
            readid = '|'.join(f[0].split("|")[:-1])

            if prev_readid != readid:
                if chunk_start > read_count + 1:
                    read_count += 1
                    prev_readid = readid
                    l_readlines = []
                    continue
                if read_count > chunk_end:
                    break
                l_read_alignments = filter_alignments(l_readlines, tpm_threshold, score_cutoff)
                extract_and_write(prev_readid, l_read_alignments, l_loci_bed, d_ref_lengths1, d_ref_lengths2, f_gtf,
                                  d_regions, chimeric_overlap, fh_chimeras, fh_singletons, hybridize)

                read_count += 1
                prev_readid = readid
                l_readlines = []
            l_readlines.append(line)

        # write last read info
        if read_count >= total_read_count:
            l_read_alignments = filter_alignments(l_readlines, tpm_threshold, score_cutoff)
            extract_and_write(readid, l_read_alignments, l_loci_bed, d_ref_lengths1, d_ref_lengths2, f_gtf,
                              d_regions, chimeric_overlap, fh_chimeras, fh_singletons, hybridize)

    if hybridize:
        # loci sequences are neeeded to hybridize
        with open(os.path.join(outdir, "loci.bed.") + str(n), "w") as fh_bed:
            for bed_line in l_loci_bed:
                fh_bed.write(bed_line + "\n")


def hybridize_and_write(outdir, intarna_params, n):
    d_loci_seqs = defaultdict()
    fa_seq = SeqIO.parse(open(os.path.join(outdir, "loci.fa.") + n), 'fasta')
    for record in fa_seq:
        # ids of the form ENSMUST00000185852:1602:1618:+(+). Strip last 3 characters
        d_loci_seqs[record.id[:-3]] = str(record.seq).upper().replace('T', 'U')

    d_hybrids = defaultdict()
    file_chimeras = os.path.join(outdir, "chimeras." + n)
    with open(file_chimeras) as fh_chimeras, open(os.path.join(outdir, "chimeras-r." + n), "w") as fh_out:
        for line in fh_chimeras:
            seq1 = seq2 = dotbracket = pos = energy = "NA"
            a = line.rstrip("\n").split("\t")
            locuspos1 = a[20]
            locuspos2 = a[21]
            if locuspos1 in d_loci_seqs and locuspos2 in d_loci_seqs:
                seq1 = d_loci_seqs[locuspos1]
                seq2 = d_loci_seqs[locuspos2]
                if (locuspos1, locuspos2) in d_hybrids:
                    dotbracket, pos, energy = d_hybrids[locuspos1, locuspos2]
                else:
                    dotbracket, pos, energy = hybridize_with_intarna(seq1, seq2, intarna_params)
                    d_hybrids[locuspos1, locuspos2] = dotbracket, pos, energy
            a.append(seq1 + "&" + seq2)
            a.append(dotbracket)
            a.append(pos)
            a.append(energy)
            fh_out.write("\t".join(a) + "\n")
    os.remove(file_chimeras)


def parse_annotations(f_gtf):
    n_exon = 1
    exon_rel_start = 0
    exon_rel_end = 0
    exon_len = 0
    prev_transcript_id = None

    limit_info = dict(gff_type=["exon", "UTR", "CDS", "miRNA", "tRNA"])

    d_attributes = defaultdict(list)
    d_attributes['tid'] = ['transcript_id', 'Name']
    d_attributes['gid'] = ['gene_id', 'Alias']
    d_attributes['name'] = ['gene_name', 'Name']
    d_attributes['type'] = ['gene_biotype', 'Type']

    l_seen_exons = set()
    with open(f_gtf) as gtf_handle:
        # Chromosome seq level
        for rec in GFF.parse(gtf_handle, limit_info=limit_info, target_lines=1):
            # for each selected sub_feature
            for sub_feature in rec.features:

                # for some reason each qualifier is a list! take the first element
                for i in d_attributes['gid']:
                    if i in sub_feature.qualifiers:
                        gene_id = sub_feature.qualifiers[i][0]
                        break
                for i in d_attributes['tid']:
                    if i in sub_feature.qualifiers:
                        transcript_id = sub_feature.qualifiers[i][0]
                        break

                # NOTE: some mature miRs have multiple locations on genome
                # this is to select only the first location for mature miR

                if sub_feature.type == 'UTR':
                    if transcript_id not in d_transcript_annotations['UTR']:
                        d_transcript_annotations['UTR'][transcript_id] = []
                    d_transcript_annotations['UTR'][transcript_id].append([rec.id,
                                                                           str(sub_feature.location.start),
                                                                           str(sub_feature.location.end),
                                                                           str(sub_feature.location.strand)])
                elif sub_feature.type == 'CDS':
                    if transcript_id not in d_transcript_annotations['CDS']:
                        d_transcript_annotations['CDS'][transcript_id] = []
                    d_transcript_annotations['CDS'][transcript_id].append([rec.id,
                                                                           str(sub_feature.location.start),
                                                                           str(sub_feature.location.end),
                                                                           str(sub_feature.location.strand)])
                # remaining are exon, miRNA and tRNA lines
                else:
                    if transcript_id + "_e" + str(n_exon).zfill(3) in l_seen_exons:
                        continue

                    # reset the variables if it is a new transcript
                    if prev_transcript_id != transcript_id:
                        d_transcript_annotations['len'][transcript_id] = 0
                        n_exon = 1
                        exon_rel_start = 0
                        exon_rel_end = 0
                        exon_len = 0

                    d_transcript_annotations['gid'][transcript_id] = gene_id
                    # biotype
                    biotype = 'NA'
                    for i in d_attributes['type']:
                        if i in sub_feature.qualifiers:
                            biotype = sub_feature.qualifiers[i][0]
                    if sub_feature.type == "miRNA":
                        biotype = "miRNA"
                    if sub_feature.type == "tRNA":
                        biotype = "tRNA"
                    d_gene_annotations['type'][gene_id] = biotype

                    # name
                    try:
                        for i in d_attributes['name']:
                            if i in sub_feature.qualifiers:
                                d_gene_annotations['name'][gene_id] = sub_feature.qualifiers[i][0]
                    except KeyError:
                        d_gene_annotations['name'][gene_id] = "NA"

                    exon_len = sub_feature.location.end - sub_feature.location.start
                    exon_rel_end = exon_rel_start + exon_len
                    # TODO: check +1 to length ?
                    d_transcript_annotations['len'][transcript_id] += exon_len + 1
                    exon_rel_start = exon_rel_end
                    prev_transcript_id = transcript_id
                    l_seen_exons.add(transcript_id + "_e" + str(n_exon).zfill(3))
                    n_exon += 1


def parse_counts_file(crl_file, tpm_cutoff):
    d_crl_tpm = defaultdict(float)
    l_loci_bed = set()
    prev_readid = None
    read_count = 0
    with open(crl_file, "r") as fh_crl_file:
        for line in fh_crl_file:
            f = line.rstrip('\n').split('\t')
            readid = '|'.join(f[0].split("|")[:-1])
            if readid != prev_readid:
                read_count += 1
            crlid = f[3]
            crl_tpm = f[12]
            d_crl_tpm[crlid] = float(crl_tpm)
            b = f[9].split(":")
            locus_bed_entry = "\t".join([":".join(b[0:-3]), b[-3], b[-2], f[9], "1", b[-1]])
            if locus_bed_entry not in l_loci_bed:
                l_loci_bed.add(locus_bed_entry)
            prev_readid = readid

    uniq_tpms = sorted(list(set(d_crl_tpm.values())))
    tpm_threshold = uniq_tpms[int(tpm_cutoff * len(uniq_tpms))]

    return read_count, tpm_threshold


def merge_files(inprefix, outfile, header, r):
    temp_files = ""
    for i in range(r):
        temp_files += " " + inprefix + "." + str(i)
    os.system("cat " + temp_files +
              " | sort -u | sed '1s/^/" + header + "\\n/' > " +
              outfile)

    for i in range(r):
        os.remove(inprefix + "." + str(i))


def hybridization_positions(dotbracket1, dotbracket2):
    end1 = -1
    end2 = -1
    for i, c1 in enumerate(dotbracket1):
        if dotbracket1[i] == '(':
            end1 = i + 1
            dotbracket1[i] = "."
            for j, c2 in enumerate(reversed(dotbracket2)):
                if dotbracket2[j] == ')':
                    end2 = j + 1
                    dotbracket2[j] = '.'
                    break
    return end1, end2


def write_interaction_summary(outdir):
    d_interactions = defaultdict(lambda: defaultdict(list))
    with open(os.path.join(outdir, "chimeras")) as fh_in:
        next(fh_in)
        for line in fh_in:
            f = line.rstrip("\n").split("\t")
            (readid, ref1, ref2, region1, region2, locus1, locus2, tpm1, tpm2, score1, score2) = (f[0], f[1], f[2],
                                                                                                  f[7], f[8], f[20],
                                                                                                  f[21], f[24], f[25],
                                                                                                  f[26], f[27])
            tpm = str(float(tpm1) + float(tpm2))
            score = str(float(score1) + float(score2))
            sequence1 = sequence2 = "NA"
            if f[29] != "NA":
                [sequence1, sequence2] = f[29].split("&")
            dotbracket = f[30]
            hybrid_start_pos = f[31]
            mfe = f[32]
            interaction = "\t".join(locus1.split(":")) + "\t" + "\t".join(locus2.split(":"))
            hybridization_pos = interaction
            [refid1, ref_start1, re1, ref_strand1, refid2, ref_start2, re2, ref_strand2] = interaction.split("\t")
            hybridized_sequences = "NA\tNA"
            interaction_otherway = "\t".join(locus2.split(":")) + "\t" + "\t".join(locus1.split(":"))

            if interaction_otherway in d_interactions:
                interaction = interaction_otherway
                hybridization_pos = interaction
                (ref2, ref1, region2, region1, tpm2, tpm1, score2, score1) = (f[1], f[2], f[7], f[8], f[24], f[25],
                                                                              f[26], f[27])
            if dotbracket != "NA":
                hybrid_end1, hybrid_end2 = hybridization_positions(list(dotbracket.split("&")[0]),
                                                                   list(dotbracket.split("&")[1]))
                # decrease by 1 before adding to the reference start
                hybrid_start1 = int(hybrid_start_pos.split("&")[0]) - 1
                hybrid_start2 = int(hybrid_start_pos.split("&")[1]) - 1
                hybridization_pos1 = "\t".join([refid1, str(int(ref_start1) + hybrid_start1),
                                                str(int(ref_start1) + hybrid_end1), ref_strand1])
                hybridization_pos2 = "\t".join([refid2, str(int(ref_start2) + hybrid_start2),
                                               str(int(ref_start2) + hybrid_end2), ref_strand2])
                hybridized_sequence1 = sequence1[hybrid_start1:hybrid_end1]
                hybridized_sequence2 = sequence2[hybrid_start2:hybrid_end2]

                hybridization_pos = hybridization_pos1 + "\t" + hybridization_pos2
                hybridized_sequences = hybridized_sequence1 + "\t" + hybridized_sequence2

                if interaction_otherway in d_interactions:
                    hybrid_start_pos = f[31].split('&')[1] + '&' + f[31].split('&')[0]
                    hybridization_pos = hybridization_pos2 + "\t" + hybridization_pos1
                    hybridized_sequences = hybridized_sequence2 + "\t" + hybridized_sequence1
                    target_dotbracket = dotbracket.split("&")[0].replace("(", ")")
                    query_dotbracket = dotbracket.split("&")[1].replace(")", "(")
                    dotbracket = query_dotbracket + "&" + target_dotbracket
                    [sequence2, sequence1] = f[29].split("&")

            d_interactions[interaction]["readid"].append(readid)
            d_interactions[interaction]["ref1"].append(ref1)
            d_interactions[interaction]["ref2"].append(ref2)
            d_interactions[interaction]["region1"].append(region1)
            d_interactions[interaction]["region2"].append(region2)
            common_info = "\t".join([sequence1, sequence2, dotbracket, mfe, hybridized_sequences,
                                     hybrid_start_pos, hybridization_pos,
                                     tpm1, tpm2, tpm, score1, score2, score])
            d_interactions[interaction]["common"] = [common_info]

    with open(os.path.join(outdir, "interactions.temp"), "w") as fh_out:
        for interaction in d_interactions.keys():
            fh_out.write("\t".join([str(len(set(d_interactions[interaction]["readid"]))),
                                    interaction,
                                    d_interactions[interaction]["common"][0],
                                    ";".join(sorted(set(d_interactions[interaction]["region1"]))),
                                    ";".join(sorted(set(d_interactions[interaction]["region2"]))),
                                    ";".join(sorted(set(d_interactions[interaction]["ref1"]))),
                                    ";".join(sorted(set(d_interactions[interaction]["ref2"])))]) + "\n")

    os.system("sort -k 1nr,1 " + os.path.join(outdir, "interactions.temp") + " > " + os.path.join(outdir, "interactions"))
    os.remove(os.path.join(outdir, "interactions.temp"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Chimeric Read Annotator: extract chimeras',
                                     usage='%(prog)s [-h] [-v,--version]',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-l', '--loci', action='store', dest='crl_file', required=True,
                        metavar='', help='Input BED file with alignments')

    parser.add_argument('-o', '--out', action='store', dest='outdir', required=True, metavar='',
                        help='Path to output directory')

    parser.add_argument('-g', '--gtf', action='store', dest='f_gtf', required=False,
                        metavar='', help='Annotation GTF file')

    parser.add_argument('-p', '--processes', action='store', type=int, default=1, metavar='',
                        dest='processes',
                        help='Number of processes to use')

    parser.add_argument('-tc', '--tpm_cutoff', action='store', type=chira_utilities.score_float, default=0, metavar='',
                        dest='tpm_cutoff',
                        help='Transcripts with less than this percentile TPMs will be discarded in '
                             'the final output. [0-1.0]')

    parser.add_argument('-sc', '--score_cutoff', action='store', type=chira_utilities.score_float, default=0.0, metavar='',
                        dest='score_cutoff',
                        help='Hybrids with less than this score will be discarded in the final output. [0-1.0]')

    parser.add_argument('-co', '--chimeric_overlap', action='store', type=int, default=2, metavar='',
                        dest='chimeric_overlap',
                        help='Maximum number of bases allowed between the chimeric segments of a read')

    parser.add_argument("-r", '--hybridize', action='store_true', dest='hybridize',
                        help="Hybridize the predicted chimeras")

    parser.add_argument("-ns", '--no_seed', action='store_true', dest='no_seed',
                        help="Do not enforce seed interactions")

    parser.add_argument("-acc", '--accessibility', type=str, choices=["C", "N"], default='N', required=False,
                        dest='accessibility', metavar='', help='IntaRNA accessibility: C (compute) or N (not)')

    parser.add_argument("-m", '--intarna_mode', type=str, choices=["H", "M", "S"], default='H', required=False,
                        dest='intarna_mode', metavar='', help='IntaRNA mode: H (heuristic), M (exact), S (seed-only)')

    parser.add_argument('-t', '--temperature', action='store', type=float, default=37, metavar='',
                        dest='temperature',
                        help='IntaRNA temperature parameter in Celsius to setup the VRNA energy parameters')

    parser.add_argument('-sbp', '--seed_bp', action='store', type=int, default=5, metavar='',
                        dest='seed_bp', choices=range(2, 20),
                        help='IntaRNA --seedBP parameter: number of inter-molecular base pairs within the seed region')

    parser.add_argument('-smpu', '--seed_min_pu', action='store', type=chira_utilities.score_float, default=0,
                        metavar='', dest='seed_min_pu',
                        help='IntaRNA --seedMinPu parameter: minimal unpaired probability '
                             '(per sequence) a seed region may have')

    parser.add_argument('-accw', '--acc_width', action='store', type=int, default=150, metavar='',
                        dest='acc_width', choices=range(0, 99999),
                        help='IntaRNA --accW parameter:  sliding window size for accessibility computation')

    parser.add_argument('-f1', '--ref_fasta1', action='store', dest='ref_fasta1', required=True,
                        metavar='', help='First prioroty fasta file')

    parser.add_argument('-f2', '--ref_fasta2', action='store', dest='ref_fasta2', required=False,
                        metavar='', help='second priority fasta file')

    parser.add_argument('-f', '--ref', action='store', dest='f_ref', required=False,
                        metavar='', help='Reference fasta file')

    parser.add_argument("-s", '--summerize', action='store_true', dest='summerize',
                        help="Summerize interactions at loci level")

    parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.4.1')

    args = parser.parse_args()

    print('CRL file                             : ' + args.crl_file)
    print('Outpur direcoty                      : ' + args.outdir)
    if args.f_gtf:
        print('Annotation file                      : ' + args.f_gtf)
    print('Number of processes                  : ' + str(args.processes))
    print('TPM cutoff                           : ' + str(args.tpm_cutoff))
    print('Score cutoff                         : ' + str(args.score_cutoff))
    print('Chimeric overlap                     : ' + str(args.chimeric_overlap))
    print('Hybridize chimeric loci?             : ' + str(args.hybridize))
    print('Do not enforce seed interaction      : ' + str(args.no_seed))
    print('1st priority reference fasta file    : ' + args.ref_fasta1)
    if args.ref_fasta2:
        print('2nd priority reference fasta file    : ' + args.ref_fasta2)
    if args.f_ref:
        print('Reference genomic fasta file         : ' + args.f_ref)
    print('Summerize interactions at loci level : ' + str(args.summerize))
    print("===================================================================")

    if args.hybridize and args.f_gtf and not args.f_ref:
        sys.stderr.write("Need the reference fasta file to hybridize. Make sure to provide the genomic fasta file"
                         " in case you already provided a GTF file.\n")
        sys.exit(1)

    if args.temperature < 0 or args.temperature > 100:
        sys.stderr.write("IntaRNA tempertature must be between 0 and 100!\n")
        sys.exit(1)
    if args.f_gtf:
        # Parse the annotations and save them to dictionaries. Additionally write exons bed files to the outputdir
        print("Parsing the annotation file")
        parse_annotations(args.f_gtf)

    d_reflen1 = defaultdict()
    d_reflen2 = defaultdict()
    # extract lengths of references from the reference fasta files, keys are reference ids used to find chimeric reads
    chira_utilities.extract_reflengths(args.ref_fasta1, d_reflen1)
    if args.ref_fasta2:
        chira_utilities.extract_reflengths(args.ref_fasta2, d_reflen2)

    print("Parsing CRLs file")
    no_of_reads, tpm_cutoff_value = parse_counts_file(args.crl_file, args.tpm_cutoff)
    print("Done")

    print(str(datetime.datetime.now()), " START: multiprocessing")

    jobs = []
    # write chimeric reads
    for k in range(args.processes):
        s = k * math.ceil(no_of_reads / args.processes)
        e = min(s + math.floor(no_of_reads / args.processes), no_of_reads)
        print(k, s, e, no_of_reads)
        j = Process(target=write_chimeras, args=(s, e, no_of_reads, d_reflen1, d_reflen2, args.hybridize,
                                                 args.chimeric_overlap, args.f_gtf, args.outdir,
                                                 args.crl_file, tpm_cutoff_value, args.score_cutoff, str(k)))
        jobs.append(j)

    for j in jobs:
        j.start()
    for j in jobs:
        j.join()

    # file name prefixes
    chimeras_prefix = chimeras_file = os.path.join(args.outdir, "chimeras")
    singletons_prefix = os.path.join(args.outdir, "singletons")

    # if hybridization required
    if args.hybridize:
        chimeras_prefix = os.path.join(args.outdir, "chimeras-r")
        f_reference = "NA"
        if args.f_ref:
            f_reference = args.f_ref
        else:
            f_reference = os.path.join(args.outdir, 'merged_reference.fa')
            l_ref = [args.ref_fasta1]
            if args.ref_fasta2:
                l_ref.append(args.ref_fasta2)
            with open(f_reference, 'w') as fh_out_ref:
                for fname in l_ref:
                    with open(fname) as infile:
                        for lin in infile:
                            fh_out_ref.write(lin)

        for k in range(args.processes):
            # before starting hybridization, extract fasta sequences
            # do not run them in parallel as they need a lot memory
            bed = os.path.join(args.outdir, "loci.bed.") + str(k)
            fa = os.path.join(args.outdir, "loci.fa.") + str(k)
            # loci.bed should exist already, now extract loci sequences
            process = subprocess.Popen(['fastaFromBed', '-s', '-nameOnly',
                                        '-fi', f_reference,
                                        '-bed', bed,
                                        '-fo', fa],
                                       stdout=subprocess.PIPE,
                                       stderr=subprocess.STDOUT,
                                       universal_newlines=True)
            for l in sorted(set(process.stdout.readlines())):
                print(l, end="")

        noseed_param = ""
        if args.no_seed:
            noseed_param = "--noSeed"
        common_intarna_params = " ".join(["--outMode C", "--outCsvCols id1,start1,start2,hybridDPfull,E",
                                          noseed_param, "-m", args.intarna_mode, "--acc", args.accessibility,
                                          "--temperature", str(args.temperature), "--seedBP", str(args.seed_bp),
                                          "--seedMinPu", str(args.seed_min_pu), "--accW", str(args.acc_width)])
        jobs = []
        # hybridize chimeric reads
        for k in range(args.processes):
            j = Process(target=hybridize_and_write, args=(args.outdir, common_intarna_params, str(k)))
            jobs.append(j)

        for j in jobs:
            j.start()
        for j in jobs:
            j.join()

        for k in range(args.processes):
            os.remove(os.path.join(args.outdir, "loci.fa.") + str(k))
            os.remove(os.path.join(args.outdir, "loci.bed.") + str(k))

        # remove the temporary reference file
        if not args.f_ref:
            if os.path.exists(f_reference):
                os.remove(f_reference)
            if os.path.exists(f_reference + ".fai"):
                os.remove(f_reference + ".fai")

    # cleanup intermediate files
    # header fields
    header_chimeras = "\t".join(["tagid",
                                 "txid1",
                                 "txid2",
                                 "geneid1",
                                 "geneid2",
                                 "symbol1",
                                 "symbol2",
                                 "region1",
                                 "region2",
                                 "tx_pos_start1",
                                 "tx_pos_end1",
                                 "tx_pos_strand1",
                                 "length1",
                                 "tx_pos_start2",
                                 "tx_pos_end2",
                                 "tx_pos_strand2",
                                 "length2",
                                 "read_info",
                                 "genomic_pos1",
                                 "genomic_pos2",
                                 "locus1",
                                 "locus2",
                                 "groupid1",
                                 "groupid2",
                                 "tpm1",
                                 "tpm2",
                                 "score1",
                                 "score2",
                                 "score",
                                 "sequences",
                                 "hybrid",
                                 "hybrid_pos",
                                 "mfe"])
    merge_files(chimeras_prefix, chimeras_file, header_chimeras, args.processes)
    header_singletons = "\t".join(["tagid",
                                   "txid",
                                   "geneid",
                                   "symbol",
                                   "region",
                                   "tx_pos_start",
                                   "tx_pos_end",
                                   "tx_pos_strand",
                                   "length",
                                   "read_info",
                                   "genomic_pos",
                                   "locus",
                                   "groupid",
                                   "tpm",
                                   "score"])
    merge_files(singletons_prefix, singletons_prefix, header_singletons, args.processes)
    print(str(datetime.datetime.now()), " END: multiprocessing")
    if args.summerize:
        write_interaction_summary(args.outdir)

