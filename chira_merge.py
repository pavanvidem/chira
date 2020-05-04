#!/usr/bin/env python
import chira_utilities
import argparse
import os
import sys
from collections import defaultdict
import datetime
import itertools
import re
from BCBio import GFF


def reads_to_segments(bed, outdir, segment_overlap_fraction, chimeric_overlap, refids1, refids2, chimeric_only):
    prev_readid = None
    prev_read_alignments = defaultdict(lambda: defaultdict())
    with open(bed) as fh_in, open(os.path.join(outdir, "segments.temp.bed"), "w") as fh_out:
        for line in fh_in:
            desc = line.split("\t")[3].split(",")
            readid = desc[0]

            # next read encountered , process the previous read alignments
            alignments1 = defaultdict(lambda: defaultdict())
            alignments2 = defaultdict(lambda: defaultdict())

            if prev_readid != readid:
                alignments_pairs = []
                if len(refids1) > 0 and len(refids2) > 0:
                    for readpos in prev_read_alignments:
                        for refid in prev_read_alignments[readpos]:
                            if refid in refids1:
                                alignments1[readpos][refid] = prev_read_alignments[readpos][refid]
                            elif refid in refids2:
                                alignments2[readpos][refid] = prev_read_alignments[readpos][refid]
                    alignments_pairs = list(itertools.product(alignments1.keys(), alignments2.keys()))
                else:
                    alignments_pairs = list(itertools.combinations(prev_read_alignments, 2))

                filtered_alignments = defaultdict(lambda: defaultdict())
                longest_chimeric = 0
                longest_singleton = 0
                chimera_found = False
                for readpos1, readpos2 in alignments_pairs:
                    if chira_utilities.overlap(readpos1, readpos2) <= chimeric_overlap:
                        chimera_found = True
                        chimeric_length = readpos1[1] - readpos1[0] + readpos2[1] - readpos2[0] + 2
                        if chimeric_length > longest_chimeric:
                            longest_chimeric = chimeric_length

                if not chimera_found and not chimeric_only:
                    for readpos in prev_read_alignments:
                        singleton_length = readpos[1] - readpos[0] + 1
                        if singleton_length > longest_singleton:
                            longest_singleton = singleton_length
                    for readpos in prev_read_alignments:
                        singleton_length = readpos[1] - readpos[0] + 1
                        # TODO: make this an option
                        if singleton_length < longest_singleton * 0.8:
                            continue
                        for refid, bedline in prev_read_alignments[readpos].items():
                            filtered_alignments[readpos][refid] = bedline
                else:
                    for readpos1, readpos2 in alignments_pairs:
                        if chira_utilities.overlap(readpos1, readpos2) <= chimeric_overlap:
                            chimeric_length = readpos1[1] - readpos1[0] + readpos2[1] - readpos2[0] + 2
                            if chimeric_length < longest_chimeric * 0.8:
                                continue
                            if len(refids1) > 0 and len(refids2) > 0:
                                if readpos1 not in filtered_alignments:
                                    for refid, bedline in alignments1[readpos1].items():
                                        filtered_alignments[readpos1][refid] = bedline
                                if readpos2 not in filtered_alignments:
                                    for refid, bedline in alignments2[readpos2].items():
                                        filtered_alignments[readpos2][refid] = bedline
                            else:
                                if readpos1 not in filtered_alignments:
                                    for refid, bedline in prev_read_alignments[readpos1].items():
                                        filtered_alignments[readpos1][refid] = bedline
                                if readpos2 not in filtered_alignments:
                                    for refid, bedline in prev_read_alignments[readpos2].items():
                                        filtered_alignments[readpos2][refid] = bedline
                d_read_segments = defaultdict(list)
                for (current_match_start, current_match_end, current_strand) in filtered_alignments:
                    similar_cigar_found = False
                    matched_segment = ""
                    for segmentid in d_read_segments:
                        if (current_match_start, current_match_end, current_strand) in d_read_segments[segmentid]:
                            similar_cigar_found = True
                            matched_segment = segmentid
                            break
                        l_read_pos = sorted(d_read_segments[segmentid])
                        # check with the first position in this segment
                        first_match_start = l_read_pos[0][0]
                        first_match_end = l_read_pos[0][1]
                        overlap = max(0, min(first_match_end, current_match_end)
                                      - max(first_match_start, current_match_start))
                        overlap_percent1 = overlap / float(first_match_end - first_match_start)
                        overlap_percent2 = overlap / float(current_match_end - current_match_start)
                        # reciprocal overlap is needed, continue if one fails
                        if overlap_percent1 < segment_overlap_fraction or \
                                overlap_percent2 < segment_overlap_fraction:
                            continue

                        # check with the last position in this segment
                        last_match_start = l_read_pos[-1][0]
                        last_match_end = l_read_pos[-1][1]
                        overlap = max(0, min(last_match_end, current_match_end)
                                      - max(last_match_start, current_match_start))
                        overlap_percent1 = overlap / float(last_match_end - last_match_start)
                        overlap_percent2 = overlap / float(current_match_end - current_match_start)
                        if overlap_percent1 >= segment_overlap_fraction and \
                                overlap_percent2 >= segment_overlap_fraction:
                            similar_cigar_found = True
                            matched_segment = segmentid
                            break
                    if not similar_cigar_found:
                        if d_read_segments:
                            matched_segment = max(d_read_segments) + 1
                        else:
                            matched_segment = 1
                    if (current_match_start, current_match_end, current_strand) not in d_read_segments[matched_segment]:
                        d_read_segments[matched_segment].append((current_match_start, current_match_end, current_strand))

                    for refid, alignment in filtered_alignments[current_match_start, current_match_end, current_strand].items():
                        f = alignment.split('\t')
                        d = f[3].split(",")
                        fh_out.write("\t".join(f[0:3]) + "\t" +
                                     prev_readid + "|" + str(matched_segment) + "," +
                                     ",".join(d[1:]) + "\t" +
                                     "\t".join(f[4:]))
                prev_read_alignments.clear()

            cigar = desc[-1]
            strand = desc[-2]
            referenceid = desc[-5]
            match_start, match_end = chira_utilities.match_positions(cigar, strand == "-")
            prev_read_alignments[match_start, match_end, strand][referenceid] = line
            prev_readid = readid


def merge_loci_overlap(outdir, alignment_overlap_fraction):
    f_bed = os.path.join(outdir, "segments.bed")
    # Merge the aligned loci
    d_desc = defaultdict(lambda: defaultdict(list))
    print(str(datetime.datetime.now()), "Reading segments BED file")
    with open(f_bed) as fh_bed:
        for line in fh_bed:
            f = line.rstrip('\n').split('\t')
            d_desc[f[0] + "\t" + f[5]][(int(f[1]), int(f[2]))].append(f[3])
    print(str(datetime.datetime.now()), "Read segments BED file")

    print(str(datetime.datetime.now()), "merging overlapping alignments")
    merged_bed = os.path.join(outdir, "merged.bed")
    with open(merged_bed, "w") as fh_out:
        # TODO: consider cases where a read belongs to different overlapping loci
        for chrom in sorted(d_desc.keys()):
            startpos_sorted = sorted(d_desc[chrom], key=lambda tup: tup[0])
            t_merged = []
            d_mergeddesc = defaultdict(list)
            for currentpos in startpos_sorted:
                if not t_merged:
                    t_merged.append(currentpos)
                    d_mergeddesc[currentpos[0]].extend(d_desc[chrom][currentpos])
                else:
                    prevpos = t_merged[-1]
                    if currentpos[0] <= prevpos[1]:
                        overlap_len = min(currentpos[1], prevpos[1]) - currentpos[0] + 1
                        if overlap_len / float(currentpos[1] - currentpos[0] + 1) >= alignment_overlap_fraction \
                                or overlap_len / float(prevpos[1] - prevpos[0] + 1) >= alignment_overlap_fraction:
                            end = max(prevpos[1], currentpos[1])
                            t_merged[-1] = (prevpos[0], end)  # replace by merged interval
                            d_mergeddesc[(prevpos[0])].extend(d_desc[chrom][currentpos])
                        else:
                            t_merged.append(currentpos)
                            d_mergeddesc[currentpos[0]].extend(d_desc[chrom][currentpos])
                            continue
                    else:
                        t_merged.append(currentpos)
                        d_mergeddesc[currentpos[0]].extend(d_desc[chrom][currentpos])
            del d_desc[chrom]
            for pos in t_merged:
                d_longest_alignments = defaultdict(int)
                # choose per locus per transcript only one longest alignment
                for alignment in d_mergeddesc[pos[0]]:
                    [segmentid, transcriptid, start, end, tx_strand, cigar] = alignment.split(',')
                    # cutting out the segment id gives the read id
                    readid = '|'.join(segmentid.split('|')[:-1])
                    if int(end) - int(start) >= d_longest_alignments[readid + "\t" + transcriptid]:
                        d_longest_alignments[readid + "\t" + transcriptid] = int(end) - int(start)
                l_alignments = []
                for alignment in d_mergeddesc[pos[0]]:
                    [segmentid, transcriptid, start, end, tx_strand, cigar] = alignment.split(',')
                    # cutting out the segment id gives the read id
                    readid = '|'.join(segmentid.split('|')[:-1])
                    if int(end) - int(start) >= d_longest_alignments[readid + "\t" + transcriptid]:
                        l_alignments.append(alignment)
                [chromid, strand] = chrom.split("\t")
                fh_out.write("\t".join([chromid, str(pos[0]), str(pos[1]), strand, ";".join(sorted(l_alignments))]) + "\n")
    d_desc.clear()
    return


def merge_loci_blockbuster(outdir, distance, min_cluster_height, min_block_height, scale):
    os.system("sort -k 1,1 -k 6,6 -k 2n,2 -k 3n,3 " + os.path.join(outdir, "segments.bed") + " > " +
              os.path.join(outdir, "segments.sorted.bed"))
    os.system("blockbuster.x " +
              " -distance " + str(distance) +
              " -minClusterHeight " + str(min_cluster_height) +
              " -minBlockHeight " + str(min_block_height) +
              " -scale " + str(scale) +
              " -print 2 " +
              os.path.join(outdir, "segments.sorted.bed") + " > " +
              os.path.join(outdir, "segments.blockbuster"))

    bed_entry = None
    merged_bed = os.path.join(outdir, "merged.bed")

    with open(merged_bed, "w") as fh_out:
        with open(os.path.join(outdir, "segments.blockbuster")) as fh_blockbuster:
            for line in fh_blockbuster:
                if line.startswith('>'):
                    if bed_entry:
                        fh_out.write(bed_entry.rstrip(";") + "\n")
                    f = line.rstrip('\n').split('\t')
                    bed_entry = '\t'.join([f[1], f[2], f[3], f[4]]) + "\t"
                else:
                    bed_entry += line.split("\t")[3] + ";"
        fh_out.write(bed_entry.rstrip(";") + "\n")
    return


def parse_annotations(gtf, outdir):
    n_exon = 1
    exon_rel_start = 0
    exon_rel_end = 0
    exon_len = 0
    prev_transcript_id = None

    limit_info = dict(gff_type=["exon", "UTR", "CDS", "miRNA", "tRNA"])

    d_attributes = defaultdict(list)
    d_attributes['tid'] = ['transcript_id', 'Name']
    d_attributes['gid'] = ['gene_id', 'Alias']

    genomic_exons = os.path.join(outdir, "genomic_exons.bed")
    transcriptomic_exons = os.path.join(outdir, "transcriptomic_exons.bed")

    l_seen_exons = set()
    with open(gtf) as gff_handle, \
            open(genomic_exons, "w") as fh_genomic_exons, \
            open(transcriptomic_exons, "w") as fh_transcriptomic_exons:
        # Chromosome seq level
        for rec in GFF.parse(gff_handle, limit_info=limit_info, target_lines=1):
            # for each selected sub_feature
            for sub_feature in rec.features:
                # each qualifier is a list! take the first element
                for i in d_attributes['tid']:
                    if i in sub_feature.qualifiers:
                        transcript_id = sub_feature.qualifiers[i][0]
                        break
                # remaining are exon, miRNA and tRNA lines
                if sub_feature.type == 'exon' or sub_feature.type == 'miRNA' or sub_feature.type == 'tRNA':
                    # reset the variables if it is a new transcript
                    if prev_transcript_id != transcript_id:
                        n_exon = 1
                        exon_rel_start = 0
                        exon_rel_end = 0
                        exon_len = 0
                    if transcript_id + "_e" + str(n_exon).zfill(3) in l_seen_exons:
                        continue

                    gene_exon_bed_entry = '\t'.join([rec.id,
                                                     str(sub_feature.location.start),
                                                     str(sub_feature.location.end),
                                                     transcript_id + "_e" + str(n_exon).zfill(3),
                                                     "1",
                                                     "-" if str(sub_feature.location.strand) == "-1" else "+"
                                                     ])
                    exon_len = sub_feature.location.end - sub_feature.location.start
                    exon_rel_end = exon_rel_start + exon_len
                    tx_exon_bed_entry = '\t'.join([transcript_id,
                                                   str(exon_rel_start),
                                                   str(exon_rel_end),
                                                   transcript_id + "_e" + str(n_exon).zfill(3),
                                                   "1",
                                                   "-" if str(sub_feature.location.strand) == "-1" else "+"
                                                   ])
                    fh_genomic_exons.write(gene_exon_bed_entry + "\n")
                    fh_transcriptomic_exons.write(tx_exon_bed_entry + "\n")

                    exon_rel_start = exon_rel_end
                    prev_transcript_id = transcript_id
                    l_seen_exons.add(transcript_id + "_e" + str(n_exon).zfill(3))
                    n_exon += 1
    return


def transcript_to_genomic_pos(transcriptomic_bed, genomic_bed, f_geneexonbed, f_txexonbed):
    id2chr = {}
    exid2start = {}
    exid2end = {}
    print("Read in genomic exon features .. ")
    fh_genomic_exon_bed = open(f_geneexonbed, "r")
    for line in fh_genomic_exon_bed:
        # chrName, s, e, exid, pol = line.rstrip('\n').split('\t')[0,1,2,3,5]
        f = line.rstrip('\n').split('\t')
        chrname = f[0]
        s = int(f[1])
        e = int(f[2])
        exid = f[3]
        match = re.match("(.+?)_e", exid)
        transcriptid = match.group(1)
        id2chr[transcriptid] = chrname
        exid2start[exid] = s
        exid2end[exid] = e
    fh_genomic_exon_bed.close()

    print("done")
    print("Transcripts with exons in annotation:   " + str(len(id2chr)))

    overlapout = transcriptomic_bed.replace(".bed", ".overlap.txt")
    print("Calculating bed overlap .. ")
    intersect_command = ("intersectBed -a " + transcriptomic_bed + " -b " + f_txexonbed + " -wb > " + overlapout)
    print(intersect_command)
    os.system(intersect_command)

    # ENSMUST00000023614      5301    5319    2175|tag_1004650|1,ENSMUST00000023614,ENSMUSG00000022897,5301,5319 \
    #       0       +       ENSMUST00000023614      1817    5754    ENSMUST00000023614_e012 0       +
    prev_readid = None
    l_withjunctions = []
    with open(overlapout) as fh_overlapout, open(genomic_bed, "w") as fh_wojunctions:
        for line in fh_overlapout:
            # transcriptid, s, e, readid, exonstart, exonend, exonid, pol
            f = line.rstrip('\n').split('\t')
            transcriptid = f[0]
            txstart = int(f[1])
            txend = int(f[2])
            readid = f[3]
            exonstart = int(f[7])
            exonid = f[9]
            pol = f[11]
            if prev_readid and prev_readid != readid:
                if len(set(l_withjunctions)) == 1:
                    fh_wojunctions.write(set(l_withjunctions).pop())
                l_withjunctions.clear()

            # Calculate genomic hit positions.
            # Plus strand case.
            if pol == "+":
                genomicstart = exid2start[exonid] + (txstart - exonstart)
                genomicend = exid2start[exonid] + (txend - exonstart)
            # Minus strand case.
            elif pol == "-":
                genomicstart = exid2end[exonid] - (txend - exonstart)
                genomicend = exid2end[exonid] - (txstart - exonstart)
            else:
                sys.exit("ERROR: Check polarity " + pol + "in line: " + line + "\n")
            b = "\t".join([id2chr[transcriptid], str(genomicstart), str(genomicend), readid, "1", pol]) + "\n"
            l_withjunctions.append(b)
            prev_readid = readid
        # write the last segment
        if len(set(l_withjunctions)) == 1:
            fh_wojunctions.write(set(l_withjunctions).pop())

    print("done")
    if os.path.exists(overlapout):
        os.remove(overlapout)

    return


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Chimeric Read Annotator: merge alignments and convert coordinates',
                                     usage='%(prog)s [-h] [-v,--version]',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-bb', '--block_based', action='store_true', dest='block_based')

    parser.add_argument('-b', '--bed', action='store', dest='bed', required=True,
                        metavar='', help='Input BED file with alignments')

    parser.add_argument('-o', '--outdir', action='store', dest='outdir', required=True, metavar='',
                        help='Output directory path for the whole analysis')

    parser.add_argument('-g', '--gtf', action='store', dest='gtf', required=False,
                        metavar='', help='Annotation GTF file')

    parser.add_argument('-ao', '--alignment_overlap', action='store', type=chira_utilities.score_float, default=0.7, metavar='',
                        dest='alignment_overlap_fraction',
                        help='Minimum percentage overlap among BED entries inorder to merge. [0-1.0]')

    parser.add_argument('-so', '--segment_overlap', action='store', type=chira_utilities.score_float, default=0.7, metavar='',
                        dest='segment_overlap_fraction',
                        help='Matching read positions with greater than this %% overlap are merged into a segment')

    parser.add_argument('-d', '--distance', action='store', type=int, default=40, metavar='',
                        dest='distance',
                        help='Blockbuster parameter distance')

    parser.add_argument('-mc', '--min_cluster_height', action='store', type=float, default=50, metavar='',
                        dest='min_cluster_height',
                        help='Blockbuster parameter minClusterHeight')

    parser.add_argument('-mb', '--min_block_height', action='store', type=float, default=1, metavar='',
                        dest='min_block_height',
                        help='Blockbuster parameter minBlockHeight')

    parser.add_argument('-sc', '--scale', action='store', type=chira_utilities.score_float, default=0.5, metavar='',
                        dest='scale',
                        help='Blockbuster parameter scale')

    parser.add_argument('-co', '--chimeric_overlap', action='store', type=int, default=2, metavar='',
                        dest='chimeric_overlap',
                        help='Maximum number of bases allowed between the chimeric segments of a read')

    parser.add_argument('-f1', '--ref_fasta1', action='store', dest='ref_fasta1', required=False,
                        metavar='', help='First prioroty fasta file')

    parser.add_argument('-f2', '--ref_fasta2', action='store', dest='ref_fasta2', required=False,
                        metavar='', help='second priority fasta file')

    parser.add_argument("-c", '--chimeric_only', action='store_true', dest='chimeric_only',
                        help="Consider chimeric reads only for merging")

    parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.2.0')

    args = parser.parse_args()

    print('Input BED file                       : ' + args.bed)
    print('Output directory                     : ' + args.outdir)
    print('Segment overlap fraction             : ' + str(args.segment_overlap_fraction))
    if args.gtf:
        print('Annotation file                      : ' + args.gtf)
    if args.block_based:
        print('Merge method                         : blockbuser based')
        print('Blockbuster distance                 : ' + str(args.distance))
        print('Blockbuster minClusterHeight         : ' + str(args.min_cluster_height))
        print('Blockbuster minBlockHeight           : ' + str(args.min_block_height))
        print('Blockbuster scale                    : ' + str(args.scale))
    else:
        print('Merge method                         : overlap based')
        print('Alignment overlap fraction           : ' + str(args.alignment_overlap_fraction))
    print('Chimeric overlap                     : ' + str(args.chimeric_overlap))
    if args.ref_fasta1:
        print('1st priority reference fasta file    : ' + args.ref_fasta1)
    if args.ref_fasta2:
        print('2nd priority reference fasta file    : ' + args.ref_fasta2)
    print("===================================================================")

    d_reflen1 = defaultdict()
    d_reflen2 = defaultdict()
    if args.ref_fasta1 and args.ref_fasta2:
        chira_utilities.extract_reflengths(args.ref_fasta1, d_reflen1)
        if args.ref_fasta2:
            chira_utilities.extract_reflengths(args.ref_fasta2, d_reflen2)

    chira_utilities.print_w_time("START: Reads to segments")
    reads_to_segments(args.bed, args.outdir, args.segment_overlap_fraction, args.chimeric_overlap,
                      d_reflen1.keys(), d_reflen2.keys(), args.chimeric_only)
    chira_utilities.print_w_time("END: Reads to segments")
    if args.gtf:
        # Parse the annotations and save them to dictionaries. Additionally write exons bed files to the outdir
        chira_utilities.print_w_time("START: Parse the annotation file")
        parse_annotations(args.gtf, args.outdir)
        chira_utilities.print_w_time("END: Parse the annotation file")
        chira_utilities.print_w_time("START: Convert to genomic coordinates")
        transcript_to_genomic_pos(os.path.join(args.outdir, "segments.temp.bed"),
                                  os.path.join(args.outdir, "segments.bed"),
                                  os.path.join(args.outdir, "genomic_exons.bed"),
                                  os.path.join(args.outdir, "transcriptomic_exons.bed"))
        chira_utilities.print_w_time("END: Convert to genomic coordinates")
        os.remove(os.path.join(args.outdir, "segments.temp.bed"))
    else:
        os.system(" ".join(["mv", os.path.join(args.outdir, "segments.temp.bed"),
                            os.path.join(args.outdir, "segments.bed")]))
    if args.block_based:
        chira_utilities.print_w_time("START: blockcuster based merging")
        merge_loci_blockbuster(args.outdir, args.distance, args.min_cluster_height, args.min_block_height, args.scale)
        chira_utilities.print_w_time("END: blockcuster based merging")
    else:
        chira_utilities.print_w_time("START: overlap based merging")
        merge_loci_overlap(args.outdir, args.alignment_overlap_fraction)
        chira_utilities.print_w_time("END: overlap based merging")

