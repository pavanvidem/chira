#!/usr/bin/env python
import chira_utilities
import argparse
import os
import sys
from collections import defaultdict
import datetime
import re
from BCBio import GFF


def reads_to_segments(bed, outdir, segment_overlap_fraction):
    d_read_segments = defaultdict(lambda: defaultdict(set))
    with open(bed) as fh_in, open(os.path.join(outdir, "segments.temp.bed"), "w") as fh_out:
        for line in fh_in:
            f = line.split("\t")
            d = f[3].split(",")
            readid = d[0]
            current_cigar = d[-1]
            current_strand = d[-2]
            current_match_start, current_match_end = chira_utilities.match_positions(current_cigar,
                                                                                     current_strand == "-")
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
                # reciprocal overlap is needed, continue if one fails
                if overlap_percent1 < segment_overlap_fraction or \
                        overlap_percent2 < segment_overlap_fraction:
                    continue

                # check with the last position in this segment
                previous_match_start = l_read_pos[-1][0]
                previous_match_end = l_read_pos[-1][-1]
                overlap = max(0, min(previous_match_end, current_match_end)
                              - max(previous_match_start, current_match_start))
                overlap_percent1 = overlap / float(previous_match_end - previous_match_start)
                overlap_percent2 = overlap / float(current_match_end - current_match_start)
                if overlap_percent1 >= segment_overlap_fraction and \
                        overlap_percent2 >= segment_overlap_fraction:
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


def merge_loci_and_quantify(outdir, alignment_overlap_fraction):
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
                                                     '0',
                                                     '-' if str(sub_feature.location.strand) == '-1' else '+'
                                                     ])
                    exon_len = sub_feature.location.end - sub_feature.location.start
                    exon_rel_end = exon_rel_start + exon_len
                    tx_exon_bed_entry = '\t'.join([transcript_id,
                                                   str(exon_rel_start),
                                                   str(exon_rel_end),
                                                   transcript_id + "_e" + str(n_exon).zfill(3),
                                                   '0',
                                                   '-' if str(sub_feature.location.strand) == '-1' else '+'
                                                   ])
                    fh_genomic_exons.write(gene_exon_bed_entry + "\n")
                    fh_transcriptomic_exons.write(tx_exon_bed_entry + "\n")

                    exon_rel_start = exon_rel_end
                    prev_transcript_id = transcript_id
                    l_seen_exons.add(transcript_id + "_e" + str(n_exon).zfill(3))
                    n_exon += 1
    return


def transcript_to_genomic_pos(transcriptomic_bed, genomic_bed, f_geneexonbed, f_txexonbed):
    id2pol = {}
    id2chr = {}
    exid2start = {}
    exid2end = {}
    print("Read in genomic exon features .. ")
    fh_genomic_exonsne_exon_bed = open(f_geneexonbed, "r")
    for line in fh_genomic_exonsne_exon_bed:
        # chrName, s, e, exid, pol = line.rstrip('\n').split('\t')[0,1,2,3,5]
        f = line.rstrip('\n').split('\t')
        chrname = f[0]
        s = int(f[1])
        e = int(f[2])
        exid = f[3]
        pol = f[5]
        match = re.match("(.+?)_e", exid)
        transcriptid = match.group(1)
        id2pol[transcriptid] = pol
        id2chr[transcriptid] = chrname
        exid2start[exid] = s
        exid2end[exid] = e
    fh_genomic_exonsne_exon_bed.close()
    print("done")
    print("Transcripts with exons:   " + str(len(id2pol)))

    overlapout = transcriptomic_bed.replace(".bed", ".overlap.txt")
    print("Calculating bed overlap .. ")
    intersect_command = ("intersectBed -a " + transcriptomic_bed + " -b " + f_txexonbed + " -wb > " + overlapout)
    print(intersect_command)
    os.system(intersect_command)

    fh_overlapout = open(overlapout, "r")

    d_idsseen = {}
    d_idsdouble = {}
    d_ids_plus = {}
    d_ids_minus = {}
    d_withjunctions = defaultdict(list)

    # ENSMUST00000023614      5301    5319    2175|tag_1004650|1,ENSMUST00000023614,ENSMUSG00000022897,5301,5319 \
    #       0       +       ENSMUST00000023614      1817    5754    ENSMUST00000023614_e012 0       +
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
        # Calculate genomic hit positions.
        # Plus strand case.
        if pol == "+":
            genomicstart = exid2start[exonid] + (txstart - exonstart)
            genomicend = exid2start[exonid] + (txend - exonstart)
            d_ids_plus[readid] = d_ids_plus.get(readid, 0) + 1
        # Minus strand case.
        elif pol == "-":
            genomicstart = exid2end[exonid] - (txend - exonstart)
            genomicend = exid2end[exonid] - (txstart - exonstart)
            d_ids_minus[readid] = d_ids_minus.get(readid, 0) + 1
        else:
            sys.exit("ERROR: Check polarity " + pol + "in line: " + line + "\n")

        # Store ID.
        if readid in d_idsseen:
            d_idsdouble[readid] = d_idsdouble.get(readid, 0) + 1
        else:
            d_idsseen[readid] = d_idsseen.get(readid, 0) + 1

        d_withjunctions[readid].append(id2chr[transcriptid] + "\t" +
                                       str(genomicstart) + "\t" +
                                       str(genomicend) + "\t" +
                                       readid +
                                       "\t0\t" +
                                       pol + '\n')
        # there are cases that a single mature miRNA presnt on different strands at different genomic positions
        # these are usually derived from different precusrsors
        # if str(id2pol[transcriptid]) != pol:
        #     print(transcriptid, readid, str(id2pol[transcriptid]), pol)
        # sys.exit("Different strands at " + line)
    fh_overlapout.close()

    print("Total alignments:        " + str(len(d_ids_plus)+len(d_ids_minus)))
    print("Hits on plus strand:     " + str(len(d_ids_plus)))
    print("Hits on minus strand:    " + str(len(d_ids_minus)))

    # Write alignments spanning exon borders to a file
    fh_wojunctions = open(genomic_bed, "w")
    n_final = 0
    print("Filtering out exon border reads .. ")
    for readid, lines in d_withjunctions.items():
        # Alignments over exon borders have more than 2 lines in hash
        # This way we can identify them and filter them out.
        # TODO: Considering them as 2 seperate alignments at genomic level doesn't harm. Just need to make sure that
        # they should not be treated as multi-mapped read. Give them a special ID while quantification
        if len(lines) > 1:
            continue
        n_final += 1

        for line in lines:
            fh_wojunctions.write(line)
    fh_wojunctions.close()

    print("done")
    if os.path.exists(overlapout):
        os.remove(overlapout)

    return


def score_float(x):
    x = float(x)
    if x < 0.0 or x > 1.0:
        raise argparse.ArgumentTypeError("%r not in range [0.0, 1.0]" % (x,))
    return x


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Chimeric Read Annotator: merge alignments and convert coordinates',
                                     usage='%(prog)s [-h] [-v,--version]',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-b', '--bed', action='store', dest='bed', required=True,
                        metavar='', help='Input BED file with alignments')

    parser.add_argument('-o', '--outdir', action='store', dest='outdir', required=True, metavar='',
                        help='Output directory path for the whole analysis')

    parser.add_argument('-g', '--gtf', action='store', dest='gtf', required=False,
                        metavar='', help='Annotation GTF file')

    parser.add_argument('-ao', '--alignment_overlap', action='store', type=score_float, default=0.7, metavar='',
                        dest='alignment_overlap_fraction',
                        help='Minimum percentage overlap among BED entries inorder to merge. [0-1.0]')

    parser.add_argument('-so', '--segment_overlap', action='store', type=score_float, default=0.7, metavar='',
                        dest='segment_overlap_fraction',
                        help='Matching read positions with greater than this %% overlap are merged into a segment')

    parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.0')

    args = parser.parse_args()

    print('Input BED file                       : ' + args.bed)
    print('Output directory                     : ' + args.outdir)
    if args.gtf:
        print('Annotation file                      : ' + args.gtf)
    print('Alignment overlap fraction           : ' + str(args.alignment_overlap_fraction))
    print('Segment overlap fraction             : ' + str(args.segment_overlap_fraction))
    print("===================================================================")
    reads_to_segments(args.bed, args.outdir, args.segment_overlap_fraction)
    if args.gtf:
        # Parse the annotations and save them to dictionaries. Additionally write exons bed files to the outdir
        print("Parsing the annotation file")
        parse_annotations(args.gtf, args.outdir)
        print("converting to genomic coordinates")
        transcript_to_genomic_pos(os.path.join(args.outdir, "segments.temp.bed"),
                                  os.path.join(args.outdir, "segments.bed"),
                                  os.path.join(args.outdir, "genomic_exons.bed"),
                                  os.path.join(args.outdir, "transcriptomic_exons.bed"))
    else:
        os.system(" ".join(["mv", os.path.join(args.outdir, "segments.temp.bed"),
                            os.path.join(args.outdir, "segments.bed")]))
    merge_loci_and_quantify(args.outdir, args.alignment_overlap_fraction)
