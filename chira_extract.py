#!/usr/bin/env python
import chira_utilities
import argparse
import os
import sys
from collections import defaultdict
from multiprocessing import Process
import logging
import itertools
import datetime
import math
from BCBio import GFF
import sqlite3
import csv

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


def hybridize_with_intarna(seq1, seq2):
    output = os.popen("IntaRNA -t " + seq1 + " -q " + seq2 + " --outMode C").read()
    alignment = output.split("\n")[1]
    dotbracket = energy = pos = "NA"
    if alignment:
        dotbracket = alignment.split(";")[7]
        energy = alignment.split(";")[8]
        pos = alignment.split(";")[1] + "&" + alignment.split(";")[4]
    return dotbracket, pos, energy


def hybridize_with_rnahybrid(id1, id2, seq1, seq2, outdir, n):
    dotbracket = energy = pos = "NA"
    hybrid_prefix = os.path.join(outdir, ".".join([n, id1, id2]))
    file1 = os.path.join(outdir, n + "." + id1)
    file2 = os.path.join(outdir, n + "." + id2)
    with open(file1, "w") as f1:
        f1.write(seq1)
    with open(file2, "w") as f2:
        f2.write(seq2)
    os.system("hybrid-min -o " + hybrid_prefix + " " + file1 + " " + file2 + " > /dev/null")
    with open(hybrid_prefix + ".dG") as fh_dg:
        next(fh_dg)
        for line in fh_dg:
            energy = line.split("\t")[1]
    l_i = []
    l_j = []
    with open(hybrid_prefix + ".37.plot") as fh_plot:
        next(fh_plot)
        for line in fh_plot:
            f = line.split("\t")
            l_i.append(int(f[0]))
            l_j.append(int(f[1]))
    for index in range(1, len(seq1)+1):
        if index in l_i:
            dotbracket += "("
        else:
            dotbracket += "."
    dotbracket += "&"
    for index in range(1, len(seq2)+1):
        if index in l_j:
            dotbracket += ")"
        else:
            dotbracket += "."
    if not l_i:
        l_i.append(-1)
        l_j.append(-1)
    pos = str(l_i[0]) + "&" + str(l_j[-1])
    return dotbracket, pos, energy


def update_best_hits(l_best_hits):
    longest = 0
    for hit in l_best_hits:
        hit_length = int(hit[10]) - int(hit[9])
        if hit[13] != "NA" and hit[14] != "NA":
            hit_length += int(hit[14]) - int(hit[13])
        if hit_length > longest:
            longest = hit_length
    best_tpm = 0
    for i, hit in enumerate(l_best_hits):
        length = int(hit[10]) - int(hit[9])
        tpm = float(hit[25])
        if hit[13] != "NA" and hit[14] != "NA":
            length += int(hit[14]) - int(hit[13])
            tpm += float(hit[26])
        if length < 0.9 * longest:
            l_best_hits[i] = None
        else:
            # if tpm < best_tpm:
            #     l_best_hits[i] = None
            # else:
            #     best_tpm = tpm
            if tpm > best_tpm:
                best_tpm = tpm
    l_best_hits = [n for n in l_best_hits if n is not None]
    for i, hit in enumerate(l_best_hits):
        tpm = float(hit[25])
        if hit[26] != "NA":
            tpm += float(hit[26])
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


def write_chimeras(l_readids, l_loci, d_loci_seqs, refids1, refids2, chimeric_overlap, hyb_program,
                   f_gtf, file_chimeras, file_singletons, n, outdir):
    d_hybrids = defaultdict()
    d_regions = defaultdict()
    with open(file_chimeras + "." + n, "w") as fh_chimeras, open(file_singletons + "." + n, "w") as fh_singletons:
        for i, readid in enumerate(l_readids):
            loci = l_loci[i]
            chimera_found = False
            l_best_chimeras = []
            l_best_singletons = []
            alignment_pairs = list(itertools.combinations(loci, 2))
            # if there are no pairs, then it is a long singleton
            for alignment1, alignment2 in alignment_pairs:
                [segmentid1, transcriptid1, locusid1, groupid1, tx_pos_start1, tx_pos_end1, tx_pos_strand1,
                 cigar1, genomic_pos1, locuspos1, locusshare1, prob1, tpm1] = alignment1.rstrip('\n').split('\t')
                [segmentid2, transcriptid2, locusid2, groupid2, tx_pos_start2, tx_pos_end2, tx_pos_strand2,
                 cigar2, genomic_pos2, locuspos2, locusshare2, prob2, tpm2] = alignment2.rstrip('\n').split('\t')
                # these are multimappings of the same segment
                if segmentid1 == segmentid2 or locuspos1 == locuspos2 or groupid1 == groupid2:
                    continue
                switch_order = chira_utilities.switch_alignments(cigar1, cigar2, tx_pos_strand1 == "-",
                                                                 tx_pos_strand2 == "-", chimeric_overlap)
                # switch_alignments returns None if the cigars are overlapping. Hence not a chimeric pair
                if not switch_order:
                    continue
                # switch_alignments returns "y" if the cigars are not in linear, then switch the alignments
                if switch_order == "y":
                    [segmentid1, transcriptid1, locusid1, groupid1, tx_pos_start1, tx_pos_end1, tx_pos_strand1,
                     cigar1, genomic_pos1, locuspos1, locusshare1, prob1, tpm1] = alignment2.rstrip('\n').split('\t')
                    [segmentid2, transcriptid2, locusid2, groupid2, tx_pos_start2, tx_pos_end2, tx_pos_strand2,
                     cigar2, genomic_pos2, locuspos2, locusshare2, prob2, tpm2] = alignment1.rstrip('\n').split('\t')

                # if it is a split reference
                if len(refids2) != 0:
                    # check if both transcripts are from the same reference database
                    if transcriptid1 in refids1 and transcriptid2 in refids1 or \
                            transcriptid1 in refids2 and transcriptid2 in refids2:
                        continue

                first_locus_score = float(prob1) * float(locusshare1)
                second_locus_score = float(prob2) * float(locusshare2)
                combined_score = first_locus_score + second_locus_score

                geneid1, name1, region1, tx_len1 = extract_annotations(transcriptid1,
                                                                       genomic_pos1,
                                                                       d_regions,
                                                                       f_gtf)
                geneid2, name2, region2, tx_len2 = extract_annotations(transcriptid2,
                                                                       genomic_pos2,
                                                                       d_regions,
                                                                       f_gtf)

                chimera = [readid,
                           transcriptid1, transcriptid2,
                           geneid1, geneid2,
                           name1, name2,
                           region1, region2,
                           str(tx_pos_start1), str(tx_pos_end1), tx_pos_strand1, str(tx_len1),
                           str(tx_pos_start2), str(tx_pos_end2), tx_pos_strand2, str(tx_len2),
                           genomic_pos1, genomic_pos2,
                           cigar1, cigar2,
                           locuspos1, locuspos2,
                           groupid1, groupid2,
                           str(float(tpm1) * float(locusshare1)), str(float(tpm2) * float(locusshare2)),
                           str(first_locus_score), str(second_locus_score), str(combined_score)]
                chimera_found = True
                l_best_chimeras.append(chimera)

            if not chimera_found:
                # singleton read
                for alignment in loci:
                    [segmentid, transcriptid, locusid, groupid, tx_pos_start, tx_pos_end, tx_pos_strand,
                     cigar, genomic_pos, locuspos, locusshare, prob, tpm] = alignment.rstrip('\n').split('\t')

                    geneid, name, region, tx_len = extract_annotations(transcriptid,
                                                                       genomic_pos,
                                                                       d_regions,
                                                                       f_gtf)

                    first_locus_score = combined_score = float(prob) * float(locusshare)
                    seq = "NA"
                    if locuspos in d_loci_seqs:
                        seq = d_loci_seqs[locuspos]

                    singleton = [readid,
                                 transcriptid, "NA",
                                 geneid, "NA",
                                 name, "NA",
                                 region, "NA",
                                 str(tx_pos_start), str(tx_pos_end), tx_pos_strand, str(tx_len),
                                 "NA", "NA", "NA", "NA",
                                 genomic_pos, "NA",
                                 cigar, "NA",
                                 locuspos, "NA",
                                 groupid, "NA",
                                 str(float(tpm) * float(locusshare)), "NA",
                                 str(first_locus_score), "NA",
                                 str(combined_score),
                                 seq, "NA", "NA", "NA"]
                    l_best_singletons.append(singleton)

            if len(l_best_chimeras) > 0:
                l_best_chimeras = update_best_hits(l_best_chimeras)
                for a in l_best_chimeras:
                    seq1 = seq2 = dotbracket = pos = energy = "NA"
                    locuspos1 = a[21]
                    locuspos2 = a[22]
                    if (locuspos1, locuspos2) not in d_hybrids:
                        if locuspos1 in d_loci_seqs and locuspos2 in d_loci_seqs:
                            seq1 = d_loci_seqs[locuspos1]
                            seq2 = d_loci_seqs[locuspos2]
                            if hyb_program == "intarna":
                                dotbracket, pos, energy = hybridize_with_intarna(seq1, seq2)
                            elif hyb_program == "hybrid-min":
                                dotbracket, pos, energy = hybridize_with_rnahybrid(locuspos1, locuspos2,
                                                                                   seq1, seq2, outdir, n)
                    a.append(seq1 + "&" + seq2)
                    a.append(dotbracket)
                    a.append(pos)
                    a.append(energy)
                    fh_chimeras.write("\t".join(a) + "\n")
            else:
                l_best_singletons = update_best_hits(l_best_singletons)
                for b in l_best_singletons:
                    fh_singletons.write("\t".join(b) + "\n")


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
                    # reset the variables if it is a new transcript
                    if prev_transcript_id != transcript_id:
                        d_transcript_annotations['len'][transcript_id] = 0
                        n_exon = 1
                        exon_rel_start = 0
                        exon_rel_end = 0
                        exon_len = 0

                    if transcript_id + "_e" + str(n_exon).zfill(3) in l_seen_exons:
                        continue

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


def parse_loci_fasta(f_loci, d_loci_seq):
    with open(f_loci) as fh_loci_fa:
        seqid = seq = ""
        for line in fh_loci_fa:
            if line.startswith('>'):
                if seq:
                    # convert DNA to RNA
                    d_loci_seq[seqid] = seq.upper().replace('T', 'U')
                seq = ""
                # cut out extra strand from the header
                # 17:84770948:84770962:-(-)
                seqid = line.lstrip('>')[:-4]
            else:
                seq += line.rstrip('\n')
    return


def extract_refids(s_refids, ref_fasta):
    with open(ref_fasta) as fh_ref_fasta:
        for line in fh_ref_fasta:
            if line.startswith('>'):
                s_refids.add(line[1:].rstrip("\n"))
    return


def parse_counts_file(crl_file, tpm_cutoff, score_cutoff, hybridize, f_ref, outdir):
    d_group_tpm = defaultdict(float)

    l_loci_bed = set()
    with open(crl_file, "r") as fh_crl_file:
        for line in fh_crl_file:
            f = line.rstrip('\n').split('\t')
            groupid = f[3]
            group_tpm = f[12]
            d_group_tpm[groupid] = float(group_tpm)
            b = f[9].split(":")
            if hybridize:
                locus_bed_entry = "\t".join([b[0], b[1], b[2], f[9], "0", b[3]])
                if locus_bed_entry not in l_loci_bed:
                    l_loci_bed.add(locus_bed_entry)

    if hybridize:
        with open(os.path.join(outdir, "loci.bed"), "w") as fh_bed:
            for bed in l_loci_bed:
                fh_bed.write(bed + "\n")
        os.system("fastaFromBed -s -nameOnly -fi " + f_ref + " -bed " +
                  os.path.join(outdir, "loci.bed") + " -fo " + os.path.join(outdir, "loci.fa"))

    uniq_tpms = sorted(list(set(d_group_tpm.values())))
    tpm_threshold = uniq_tpms[int(tpm_cutoff * len(uniq_tpms))]

    d_read_aligns = defaultdict(list)
    with open(crl_file, "r") as fh_crl_file:
        for line in fh_crl_file:
            f = line.rstrip('\n').split('\t')
            # last field after | represents the segment id, rest of the string before is read id
            readid = '|'.join(f[0].split("|")[:-1])
            group_tpm = f[12]
            prob = f[10]
            locusshare = f[11]
            if float(group_tpm) < tpm_threshold:
                continue
            locus_score = float(prob) * float(locusshare)
            if locus_score < score_cutoff:
                continue
            d_read_aligns[readid].append(line)
    return d_read_aligns


def score_float(x):
    x = float(x)
    if x < 0.0 or x > 1.0:
        raise argparse.ArgumentTypeError("%r not in range [0.0, 1.0]" % (x,))
    return x


def tsv_to_sqlite(in_file, out_file):
    conn = sqlite3.connect(out_file)
    cursor = conn.cursor()

    header = True
    with open(in_file) as fh_in:
        reader = csv.reader(fh_in, delimiter='\t')
        columns = next(reader)
        columns = [h.strip() for h in columns]
        if header:
            sql_command = 'CREATE TABLE IF NOT EXISTS Chimeras(%s)' % ', '.join(['%s' % column for column in columns])
            cursor.execute(sql_command)
            header = False

        query = 'insert into Chimeras({0}) values ({1})'
        query = query.format(','.join(columns), ','.join('?' * len(columns)))
        cursor = conn.cursor()
        for row in reader:
            cursor.execute(query, row)
        conn.commit()
        conn.close()


def merge_files(outfile, r):
    # header fields
    header = "\t".join(["tagid",
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
                        "genomic_pos1",
                        "genomic_pos2",
                        "cigar1",
                        "cigar2",
                        "locus1",
                        "locus2",
                        "groupid1",
                        "groupid2",
                        "tpm1",
                        "tpm2",
                        "score1",
                        "score2",
                        "score",
                        "sequence1",
                        "sequence2",
                        "hybrid"])
    with open(outfile, 'w') as fh_out:
        fh_out.write(header + "\n")
        for i in range(r):
            with open(outfile + "." + str(i)) as infile:
                for line in infile:
                    fh_out.write(line)
            os.remove(outfile + "." + str(i))


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

    parser.add_argument('-tc', '--tpm_cutoff', action='store', type=score_float, default=0, metavar='',
                        dest='tpm_cutoff',
                        help='Transcripts with less than this percentile TPMs will be discarded in '
                             'the final output. [0-1.0]')

    parser.add_argument('-sc', '--score_cutoff', action='store', type=score_float, default=0.0, metavar='',
                        dest='score_cutoff',
                        help='Hybrids with less than this score will be discarded in the final output. [0-1.0]')

    parser.add_argument('-co', '--chimeric_overlap', action='store', type=int, default=2, metavar='',
                        dest='chimeric_overlap',
                        help='Maximum number of bases allowed between the chimeric segments of a read')

    parser.add_argument("-q", '--create-sqlite', action='store_true', dest='create_sqlitedb',
                        help="Hybridize the predicted chimeras")

    parser.add_argument("-r", '--hybridize', action='store_true', dest='hybridize',
                        help="Hybridize the predicted chimeras")

    parser.add_argument("-t", '--hyb_program', type=str, choices=["intarna", "hybrid-min"], default='bwa',
                        required=False, dest='hyb_program', metavar='',
                        help='Program to hybridize, intarna or hybrid-min')

    parser.add_argument('-f1', '--ref_fasta1', action='store', dest='ref_fasta1', required=True,
                        metavar='', help='First prioroty fasta file')

    parser.add_argument('-f2', '--ref_fasta2', action='store', dest='ref_fasta2', required=False,
                        metavar='', help='second priority fasta file')

    parser.add_argument('-f', '--ref', action='store', dest='f_ref', required=False,
                        metavar='', help='Reference fasta file')

    parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.0')

    args = parser.parse_args()

    print('CRL file                             : ' + args.crl_file)
    print('Outpur direcoty                      : ' + args.outdir)
    print('Create sqlite database?              : ' + str(args.create_sqlitedb))
    if args.f_gtf:
        print('Annotation file                      : ' + args.f_gtf)
    print('Number of processes                  : ' + str(args.processes))
    print('TPM cutoff                           : ' + str(args.tpm_cutoff))
    print('Score cutoff                         : ' + str(args.score_cutoff))
    print('Chimeric overlap                     : ' + str(args.chimeric_overlap))
    print('Hybridize chimeric loci?             : ' + str(args.hybridize))
    print('1st priority reference fasta file    : ' + args.ref_fasta1)
    if args.ref_fasta2:
        print('2nd priority reference fasta file    : ' + args.ref_fasta2)
    if args.f_ref:
        print('Reference genomic fasta file         : ' + args.f_ref)
    print("===================================================================")

    if args.hybridize and args.f_gtf and not args.f_ref:
        sys.stderr.write("Need the reference fasta file to hybridize. Make sure to provide the genomic fasta file"
                         " in case you already provided a GTF file. \n")
        sys.exit(1)

    if args.f_gtf:
        # Parse the annotations and save them to dictionaries. Additionally write exons bed files to the outputdir
        print("Parsing the annotation file")
        parse_annotations(args.f_gtf)

    s_refids1 = set()
    s_refids2 = set()
    extract_refids(s_refids1, args.ref_fasta1)
    if args.ref_fasta2:
        extract_refids(s_refids2, args.ref_fasta2)

    d_read_alignments = parse_counts_file(args.crl_file, args.tpm_cutoff, args.score_cutoff,
                                          args.hybridize, args.f_ref, args.outdir)

    d_loci_sequences = defaultdict()
    if args.hybridize:
        parse_loci_fasta(os.path.join(args.outdir, "loci.fa"), d_loci_sequences)

    print(str(datetime.datetime.now()), " START: multiprocessing")
    chimeras_file = os.path.join(args.outdir, "chimeras")
    singletons_file = os.path.join(args.outdir, "singletons")
    jobs = []
    for k in range(args.processes):
        s = k * math.ceil(len(d_read_alignments) / args.processes)
        e = min(s + math.ceil(len(d_read_alignments) / args.processes), len(d_read_alignments))
        readids = []
        l_locus = []
        for read in list(d_read_alignments)[s:e]:
            l_locus.append(d_read_alignments[read])
            readids.append(read)
        j = Process(target=write_chimeras, args=(readids, l_locus, d_loci_sequences, s_refids1, s_refids2,
                                                 args.chimeric_overlap, args.hyb_program, args.f_gtf,
                                                 chimeras_file,  singletons_file, str(k), args.outdir))
        jobs.append(j)
    d_read_alignments.clear()
    for j in jobs:
        j.start()
    for j in jobs:
        j.join()
    # cleanup intermediate files
    merge_files(chimeras_file, args.processes)
    merge_files(singletons_file, args.processes)
    print(str(datetime.datetime.now()), " END: multiprocessing")
    # os.system("gzip -f " + chimeras_file)
    # os.system("gzip -f " + singletons_file)
    if args.create_sqlitedb:
        if os.path.exists(chimeras_file + '.chira.sqlite'):
            os.system("rm " + chimeras_file + '.chira.sqlite')
        tsv_to_sqlite(chimeras_file, chimeras_file + '.chira.sqlite')
    logging.info("| END: write final interactions to a file")
