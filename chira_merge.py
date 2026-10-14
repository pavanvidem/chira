#!/usr/bin/env python
import chira_utilities
import argparse
import os
import shlex
import subprocess
from collections import defaultdict
import datetime
import itertools
from BCBio import GFF


def run_command(command, stdout=None):
    """Run literal arguments and stop on an external command failure."""
    print(shlex.join(command))
    subprocess.run(command, stdout=stdout, check=True)


def filter_alignments(prev_read_alignments, chimeric_overlap, refids1, refids2, chimeric_only, lt):
    """Select longest qualifying singletons or closest chimeric pairs."""
    refids1, refids2 = set(refids1), set(refids2)
    split_reference = bool(refids1 and refids2)

    alignments1 = defaultdict(list)
    alignments2 = defaultdict(list)
    alignments_pairs = []
    # choose first alignment if there are multiple alignmnets on same reference postion
    # if split reference
    if split_reference:
        # split the alignments based on their references
        l_pos1 = []
        l_pos2 = []
        for readpos in prev_read_alignments:
            for refpos in prev_read_alignments[readpos]:
                if refpos.split(',')[0] in refids1:
                    alignments1[refpos].append(readpos)
                    if readpos not in l_pos1:
                        l_pos1.append(readpos)
                elif refpos.split(',')[0] in refids2:
                    alignments2[refpos].append(readpos)
                    if readpos not in l_pos2:
                        l_pos2.append(readpos)
        alignments_pairs = list(itertools.product(l_pos1, l_pos2))
    else:
        alignments_pairs = list(itertools.combinations(prev_read_alignments, 2))

    filtered_alignments = defaultdict(lambda: defaultdict())
    longest_chimeric = 0
    longest_singleton = 0
    d_closest_chimeric = defaultdict(lambda: float('inf'))
    # find longest and closest chimeric
    for readpos1, readpos2 in alignments_pairs:
        if chira_utilities.inclusive_overlap(readpos1, readpos2) <= chimeric_overlap:
            # total length of chimeric segments
            chimeric_length = readpos1[1] - readpos1[0] + readpos2[1] - readpos2[0] + 2
            if chimeric_length > longest_chimeric:
                longest_chimeric = chimeric_length
            # distance between chimeric segments
            if readpos1[0] < readpos2[0]:
                chimeric_distance = readpos2[0] - readpos1[1]
            else:
                chimeric_distance = readpos1[0] - readpos2[1]
            for refpos1 in prev_read_alignments[readpos1]:
                for refpos2 in prev_read_alignments[readpos2]:
                    if split_reference:
                        if refpos1.split(',')[0] not in refids1 or refpos2.split(',')[0] not in refids2:
                            continue
                        if readpos1 not in alignments1[refpos1] or readpos2 not in alignments2[refpos2]:
                            continue
                    if chimeric_distance < d_closest_chimeric[refpos1, refpos2]:
                        d_closest_chimeric[refpos1, refpos2] = chimeric_distance
    # find longest singleton
    for readpos in prev_read_alignments:
        singleton_length = readpos[1] - readpos[0] + 1
        if singleton_length > longest_singleton:
            longest_singleton = singleton_length

    # singleton is longer than chimeric, a possible singleton
    if longest_singleton > longest_chimeric:
        if not chimeric_only:  # not chimera_found and
            for readpos in prev_read_alignments:
                singleton_length = readpos[1] - readpos[0] + 1
                if singleton_length < longest_singleton * lt:
                    continue
                seen_references = set()
                for refpos, bedline in prev_read_alignments[readpos].items():
                    # for each read position and reference consider 1st alignmnet only
                    reference = refpos.split(',')[0]
                    if reference not in seen_references:
                        filtered_alignments[readpos][refpos] = bedline
                        seen_references.add(reference)
    # chimeric read
    else:
        for readpos1, readpos2 in alignments_pairs:
            if chira_utilities.inclusive_overlap(readpos1, readpos2) > chimeric_overlap:
                continue
            chimeric_length = readpos1[1] - readpos1[0] + readpos2[1] - readpos2[0] + 2
            if chimeric_length < longest_chimeric * lt:
                continue
            # distance between chimeric segments
            if readpos1[0] < readpos2[0]:
                chimeric_distance = readpos2[0] - readpos1[1]
            else:
                chimeric_distance = readpos1[0] - readpos2[1]
            for refpos1 in prev_read_alignments[readpos1]:
                for refpos2 in prev_read_alignments[readpos2]:
                    qualified = False
                    if split_reference:
                        if refpos1.split(',')[0] in refids1 and refpos2.split(',')[0] in refids2:
                            if readpos1 in alignments1[refpos1] and readpos2 in alignments2[refpos2]:
                                qualified = True
                    else:
                        qualified = True
                    # for each reference position choose the chimeric segments with minimal distance
                    if qualified and chimeric_distance == d_closest_chimeric[refpos1, refpos2]:
                        filtered_alignments[readpos1][refpos1] = prev_read_alignments[readpos1][refpos1]
                        filtered_alignments[readpos2][refpos2] = prev_read_alignments[readpos2][refpos2]
    return filtered_alignments


def reciprocal_read_overlap(first, second, cutoff):
    """Check inclusive overlap as a fraction of both read spans."""
    shared = chira_utilities.inclusive_overlap(first, second)
    return (shared / float(first[1] - first[0] + 1) >= cutoff and
            shared / float(second[1] - second[0] + 1) >= cutoff)


def write_segments(prev_readid, filtered_alignments, segment_overlap_fraction, fh_out):
    d_read_segments = defaultdict(list)
    for (current_match_start, current_match_end, current_strand) in sorted(filtered_alignments):
        similar_cigar_found = False
        matched_segment = ""
        for segmentid in d_read_segments:
            if (current_match_start, current_match_end, current_strand) in d_read_segments[segmentid]:
                similar_cigar_found = True
                matched_segment = segmentid
                break
            l_read_pos = sorted(d_read_segments[segmentid])
            # check with the first position in this segment
            current = (current_match_start, current_match_end)
            # reciprocal overlap is needed, continue if one fails
            if not reciprocal_read_overlap(l_read_pos[0], current, segment_overlap_fraction):
                continue

            # check with the last position in this segment
            if reciprocal_read_overlap(l_read_pos[-1], current, segment_overlap_fraction):
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


def update_cigar(bedline, merged_readpos, merged_refpos):
    [readid, refid, start, end, strand, cigar] = bedline.split('\t')[3].split(',')
    read_len = chira_utilities.query_length(cigar, strand == "-")
    new_cigar = ""
    # TODO check if the cigar is correct for reverse strandede alignment
    if merged_readpos[0] > 0:
        new_cigar += str(merged_readpos[0] - 1) + "S"
    new_cigar += str(merged_readpos[1] - merged_readpos[0] + 1) + "M"
    if read_len - merged_readpos[1] > 0:
        new_cigar += str(read_len - merged_readpos[1]) + "S"
    #TODO check if nw cigar has correct length
    new_bedline = "\t".join([refid, str(merged_refpos[0]), str(merged_refpos[1]),
                            ",".join([readid, refid, str(merged_readpos[0]), str(merged_readpos[1]), strand, new_cigar]),
                            "1", strand]) + "\n"
    return new_bedline


def stitch_alignments(alignments_by_ref, read_alignments, prev_readid):
    for refid in alignments_by_ref.keys():
        l_merged_refpos = []
        l_merged_readpos = []
        if len(alignments_by_ref[refid]) < 2:
            continue
        if len(alignments_by_ref[refid]) != len(set(tuple(l) for l in alignments_by_ref[refid].values())):
            continue
        l_delete = []
        l_add = defaultdict()
        for curr_ref_pos in sorted(alignments_by_ref[refid].keys()):
            if not l_merged_refpos:
                l_merged_refpos.append(curr_ref_pos)
            else:
                prev_ref_pos = l_merged_refpos[-1]
                if curr_ref_pos[0] <= prev_ref_pos[1] + 1:
                    # same reference position different read positions
                    # ENSMUST00000137264	502	516	tag_118603|1,ENSMUST00000137264,502,516,+,25S14M20S	1	+
                    # ENSMUST00000137264	502	516	tag_118603|1,ENSMUST00000137264,502,516,+,45S14M	1	+
                    if len(alignments_by_ref[refid][curr_ref_pos]) > 1:
                        continue
                    curr_read_pos = alignments_by_ref[refid][curr_ref_pos][0]
                    if not l_merged_readpos:
                        l_merged_readpos.append(curr_read_pos)
                    else:
                        prev_read_pos = l_merged_readpos[-1]
                        # both the reference and reads are mergable
                        if curr_read_pos[0] <= prev_read_pos[1] + 1:
                            refpos_end = max(prev_ref_pos[1], curr_ref_pos[1])
                            readpos_end = max(prev_read_pos[1], curr_read_pos[1])
                            l_merged_refpos[-1] = (prev_ref_pos[0], refpos_end)
                            l_merged_readpos[-1] = (prev_read_pos[0], readpos_end)

                            prev_reference = ",".join([refid[0], str(prev_ref_pos[0]), str(prev_ref_pos[1]), refid[1]])
                            curr_reference = ",".join([refid[0], str(curr_ref_pos[0]), str(curr_ref_pos[1]), refid[1]])
                            if alignments_by_ref[refid][prev_ref_pos]:
                                l_delete.append((alignments_by_ref[refid][prev_ref_pos][0], prev_reference))
                            l_delete.append((curr_read_pos, curr_reference))

                            new_bedline = update_cigar(read_alignments[curr_read_pos][curr_reference], l_merged_readpos[-1], l_merged_refpos[-1])
                            new_reference = ",".join([refid[0], str(l_merged_refpos[-1][0]), str(l_merged_refpos[-1][1]), refid[1]])
                            new_read = (l_merged_readpos[-1][0], l_merged_readpos[-1][1], refid[1])
                            l_add[new_read, new_reference] = new_bedline

        for read, ref in l_delete:
            del read_alignments[read][ref]
        for read, ref in l_add:
            read_alignments[read][ref] = l_add[(read, ref)]


def reads_to_segments(bed, outdir, segment_overlap_fraction, chimeric_overlap, refids1, refids2, chimeric_only, lt):
    """Convert contiguous read groups from name-sorted BED into segments."""
    with open(bed) as fh_in, open(os.path.join(outdir, "segments.temp.bed"), "w") as fh_out:
        for readid, lines in itertools.groupby(fh_in, key=lambda line: line.split('\t')[3].split(',')[0]):
            read_alignments = defaultdict(dict)
            for line in lines:
                desc = line.split('\t')[3].split(',')
                if len(desc) != 6:
                    raise ValueError('Invalid alignment description for {!r}'.format(readid))
                strand, cigar = desc[-2:]
                match_start, match_end = chira_utilities.match_positions(cigar, strand == '-')
                referencepos = ','.join(desc[1:5])
                read_alignments[match_start, match_end, strand][referencepos] = line.rstrip('\n') + '\n'
            filtered = filter_alignments(read_alignments, chimeric_overlap, refids1, refids2, chimeric_only, lt)
            write_segments(readid, filtered, segment_overlap_fraction, fh_out)


def merge_overlapping_intervals(d_desc, chrom, alignment_overlap_fraction):
    """Greedily merge half-open intervals with sufficient overlap on either side."""
    startpos_sorted = sorted(d_desc[chrom])
    t_merged = []
    d_mergeddesc = defaultdict(list)
    for currentpos in startpos_sorted:
        if currentpos[0] < 0 or currentpos[1] <= currentpos[0]:
            raise ValueError('Invalid BED interval: {!r}'.format(currentpos))
        if not t_merged:
            t_merged.append(currentpos)
            d_mergeddesc[currentpos[0]].extend(d_desc[chrom][currentpos])
        else:
            prevpos = t_merged[-1]
            if currentpos[0] < prevpos[1]:
                overlap_len = chira_utilities.overlap(currentpos, prevpos)
                if overlap_len / float(currentpos[1] - currentpos[0]) >= alignment_overlap_fraction \
                        or overlap_len / float(prevpos[1] - prevpos[0]) >= alignment_overlap_fraction:
                    end = max(prevpos[1], currentpos[1])
                    t_merged[-1] = (prevpos[0], end)  # replace by merged interval
                    d_mergeddesc[prevpos[0]].extend(d_desc[chrom][currentpos])
                else:
                    t_merged.append(currentpos)
                    d_mergeddesc[currentpos[0]].extend(d_desc[chrom][currentpos])
                    continue
            else:
                t_merged.append(currentpos)
                d_mergeddesc[currentpos[0]].extend(d_desc[chrom][currentpos])
    return t_merged, d_mergeddesc


def longest_alignments(alignments):
    """Keep all ties for the longest reference span per read/transcript pair."""
    parsed = []
    longest = defaultdict(int)
    for alignment in alignments:
        segment, transcript, start, end, strand, cigar = alignment.split(',')
        key = (segment.rsplit('|', 1)[0], transcript)
        length = int(end) - int(start)
        longest[key] = max(longest[key], length)
        parsed.append((alignment, key, length))
    return [alignment for alignment, key, length in parsed if length == longest[key]]


def merge_loci_overlap(outdir, alignment_overlap_fraction, min_locus_size):
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
        for chrom in sorted(d_desc.keys()):
            t_merged, d_mergeddesc = merge_overlapping_intervals(d_desc, chrom, alignment_overlap_fraction)
            del d_desc[chrom]
            for pos in t_merged:
                l_alignments = longest_alignments(d_mergeddesc[pos[0]])
                [chromid, strand] = chrom.split("\t")
                # ignore locus if has not enough alignments
                if len(l_alignments) < min_locus_size:
                    continue
                fh_out.write("\t".join([chromid, str(pos[0]), str(pos[1]), strand, ";".join(sorted(l_alignments))]) + "\n")
    d_desc.clear()
    return


def write_merged_pos(d_mergeddesc, prev_chrom_strand, fh_out):
    for merged_pos in sorted(d_mergeddesc):
        l_alignments = d_mergeddesc[merged_pos]
        [chromid, strand] = prev_chrom_strand.split("\t")
        fh_out.write(
            "\t".join([chromid, str(merged_pos[0]), str(merged_pos[1]), strand, ";".join(sorted(l_alignments))]) + "\n")
    d_mergeddesc.clear()


def merge_loci_blockbuster(outdir, distance, min_cluster_height, min_block_height, scale, alignment_overlap_fraction):
    with open(os.path.join(outdir, 'segments.sorted.bed'), 'w') as output:
        run_command(['sort', '-k', '1,1', '-k', '6,6', '-k', '2n,2', '-k', '3n,3',
                     os.path.join(outdir, 'segments.bed')], stdout=output)
    with open(os.path.join(outdir, 'segments.blockbuster'), 'w') as output:
        run_command(['blockbuster.x', '-distance', str(distance),
                     '-minClusterHeight', str(min_cluster_height), '-minBlockHeight', str(min_block_height),
                     '-scale', str(scale), '-print', '1', os.path.join(outdir, 'segments.sorted.bed')], stdout=output)

    d_desc = defaultdict(lambda: defaultdict(list))
    with open(os.path.join(outdir, "segments.blockbuster")) as fh_blockbuster:
        for line in fh_blockbuster:
            f = line.rstrip('\n').split('\t')
            if line.startswith('>'):
                continue
            d_desc[f[1] + "\t" + f[4]][(int(f[2]), int(f[3]))].append(f[0])

    d_merged = defaultdict()
    for chrom in sorted(d_desc.keys()):
        t_merged, d_mergeddesc = merge_overlapping_intervals(d_desc, chrom, alignment_overlap_fraction)
        d_merged[chrom] = t_merged

    merged_bed = os.path.join(outdir, "merged.bed")
    d_mergeddesc = defaultdict(list)
    prev_chrom_strand = None
    with open(merged_bed, "w") as fh_out, open(os.path.join(outdir, "segments.sorted.bed")) as fh_in:
        for line in fh_in:
            f = line.rstrip('\n').split('\t')
            chrom_strand = f[0] + "\t" + f[5]
            start = f[1]
            end = f[2]
            if prev_chrom_strand != chrom_strand:
                write_merged_pos(d_mergeddesc, prev_chrom_strand, fh_out)
            if chrom_strand in d_merged:
                for merged_pos in d_merged[chrom_strand]:
                    # within the merged psotion range
                    if int(start) >= merged_pos[0] and int(end) <= merged_pos[1]:
                        d_mergeddesc[merged_pos].append(f[3])
                        break
            prev_chrom_strand = chrom_strand
        # write last chromosome positions
        write_merged_pos(d_mergeddesc, prev_chrom_strand, fh_out)
    os.remove(os.path.join(outdir, "segments.blockbuster"))

    return


def parse_annotations(gtf, outdir):
    """Build paired exon BEDs independently of annotation feature ordering."""
    transcripts = defaultdict(set)

    def visit(features, chromosome):
        for feature in features:
            if feature.type in ('exon', 'miRNA', 'tRNA'):
                ids = feature.qualifiers.get('transcript_id') or feature.qualifiers.get('Name')
                if not ids:
                    raise ValueError('Missing transcript_id/Name for {} on {}'.format(feature.type, chromosome))
                strand = feature.location.strand
                if strand not in (-1, 1):
                    raise ValueError('Missing or invalid annotation strand for {!r}'.format(ids[0]))
                interval = (chromosome, int(feature.location.start), int(feature.location.end), strand)
                if interval[2] <= interval[1]:
                    raise ValueError('Empty annotation interval for {!r}'.format(ids[0]))
                for transcript_id in ids:
                    transcripts[transcript_id].add(interval)
            visit(feature.sub_features, chromosome)

    with open(gtf) as handle:
        for record in GFF.parse(handle, limit_info={'gff_type': ['exon', 'miRNA', 'tRNA']}, target_lines=1):
            visit(record.features, record.id)

    with open(os.path.join(outdir, 'genomic_exons.bed'), 'w') as genomic, \
            open(os.path.join(outdir, 'transcriptomic_exons.bed'), 'w') as transcriptomic:
        for transcript_id, exons in sorted(transcripts.items()):
            if len({(exon[0], exon[3]) for exon in exons}) != 1:
                raise ValueError('Conflicting chromosome/strand for {!r}'.format(transcript_id))
            strand = '-' if next(iter(exons))[3] == -1 else '+'
            offset = 0
            for number, (chromosome, start, end, _) in enumerate(
                    sorted(exons, key=lambda exon: (exon[1], exon[2]), reverse=strand == '-'), start=1):
                exon_id = '{}_e{:03d}'.format(transcript_id, number)
                genomic.write('\t'.join([chromosome, str(start), str(end), exon_id, '1', strand]) + '\n')
                transcriptomic.write('\t'.join([transcript_id, str(offset), str(offset + end - start),
                                                exon_id, '1', strand]) + '\n')
                offset += end - start


def transcript_to_genomic_pos(transcriptomic_bed, genomic_bed, f_geneexonbed, f_txexonbed):
    """Convert exon intersections, keeping the longest fragment per alignment."""
    exons = {}
    with open(f_geneexonbed) as handle:
        for line in handle:
            chromosome, start, end, exon_id, score, strand = line.rstrip('\n').split('\t')
            exons[exon_id] = chromosome, int(start), int(end)
    overlapout = transcriptomic_bed + '.overlap.txt'
    with open(overlapout, 'w') as output:
        run_command(['intersectBed', '-a', transcriptomic_bed, '-b', f_txexonbed, '-wb'], stdout=output)
    with open(overlapout) as intersections, open(genomic_bed, 'w') as output:
        for description, lines in itertools.groupby(intersections, key=lambda line: line.split('\t')[3]):
            fragments = set()
            for line in lines:
                fields = line.rstrip('\n').split('\t')
                txstart, txend = int(fields[1]), int(fields[2])
                exonstart, exonend = int(fields[7]), int(fields[8])
                txstart, txend = max(txstart, exonstart), min(txend, exonend)
                if txend <= txstart:
                    continue
                chromosome, start, end = exons[fields[9]]
                strand = fields[11]
                if strand == '+':
                    genomicstart, genomicend = start + txstart - exonstart, start + txend - exonstart
                elif strand == '-':
                    genomicstart, genomicend = end - (txend - exonstart), end - (txstart - exonstart)
                else:
                    raise ValueError('Invalid exon strand: {!r}'.format(strand))
                if fields[5] == '-':
                    strand = '-' if strand == '+' else '+'
                fragments.add('\t'.join([chromosome, str(genomicstart), str(genomicend),
                                          description, '1', strand]) + '\n')
            if fragments:
                # Lexicographic tie-breaking matches the existing deterministic policy.
                output.write(min(fragments, key=lambda row: (-(int(row.split('\t')[2]) -
                                                               int(row.split('\t')[1])), row)))
    os.remove(overlapout)


def create_parser():

    parser = argparse.ArgumentParser(description='Chimeric Read Annotator: merge alignments and convert coordinates',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-bb', '--block_based', action='store_true', help='Use Blockbuster-based locus merging')

    parser.add_argument('-b', '--bed', required=True,
                        metavar='BED', help='Input BED file with alignments grouped by read ID')

    parser.add_argument('-o', '--outdir', required=True, metavar='DIR',
                        help='Output directory path for the whole analysis')

    parser.add_argument('-g', '--gtf', action='store', dest='gtf', required=False,
                        metavar='GTF', help='Annotation GTF/GFF file')

    parser.add_argument('-ao', '--alignment_overlap', action='store', type=chira_utilities.score_float, default=0.7,
                        dest='alignment_overlap_fraction',
                        help='Minimum overlap fraction of either BED interval to merge. [0-1.0]')

    parser.add_argument('-so', '--segment_overlap', action='store', type=chira_utilities.score_float, default=0.7,
                        dest='segment_overlap_fraction',
                        help='Matching read positions with greater than this %% overlap are merged into a segment')

    parser.add_argument('-lt', '--length_threshold', action='store', type=chira_utilities.score_float, default=0.9,
                        dest='length_threshold',
                        help='Minimum alignment length as a fraction of the longest alignment. [0-1.0]')

    parser.add_argument('-d', '--distance', action='store', type=int, default=30,
                        dest='distance',
                        help='Blockbuster parameter distance')

    parser.add_argument('-mc', '--min_cluster_height', action='store', type=int, default=10,
                        dest='min_cluster_height',
                        help='Blockbuster parameter minClusterHeight')

    parser.add_argument('-mb', '--min_block_height', action='store', type=int, default=10,
                        dest='min_block_height',
                        help='Blockbuster parameter minBlockHeight')

    parser.add_argument('-sc', '--scale', action='store', type=chira_utilities.score_float, default=0.1,
                        dest='scale',
                        help='Blockbuster parameter scale')

    parser.add_argument('-co', '--chimeric_overlap', action='store', type=int, default=2,
                        dest='chimeric_overlap',
                        help='Maximum number of overlapping bases allowed between chimeric segments')

    parser.add_argument('-f1', '--ref_fasta1', action='store', dest='ref_fasta1', required=False,
                        metavar='FASTA', help='First reference FASTA file')

    parser.add_argument('-f2', '--ref_fasta2', action='store', dest='ref_fasta2', required=False,
                        metavar='FASTA', help='Second reference FASTA file')

    parser.add_argument("-c", '--chimeric_only', action='store_true', dest='chimeric_only',
                        help="Consider chimeric reads only for merging")

    parser.add_argument('-ls', '--min_locus_size', action='store', type=int, default=1,
                        dest='min_locus_size',
                        help='Minimum number of alignments required per merged locus (overlap method)')

    parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.4.3')
    return parser


def validate_args(args):
    if args.ref_fasta2 and not args.ref_fasta1:
        raise ValueError('--ref_fasta2 requires --ref_fasta1')
    for name in ('min_locus_size', 'min_cluster_height', 'min_block_height'):
        if getattr(args, name) <= 0:
            raise ValueError('--{} must be positive'.format(name))
    for name in ('distance', 'chimeric_overlap'):
        if getattr(args, name) < 0:
            raise ValueError('--{} must be non-negative'.format(name))
    for path in (args.bed, args.gtf, args.ref_fasta1, args.ref_fasta2):
        if path and not os.path.isfile(path):
            raise ValueError('Input file does not exist or is not a regular file: {}'.format(path))


def print_args(args):
    print('Input BED file                       : ' + args.bed)
    print('Output directory                     : ' + args.outdir)
    print('Segment overlap fraction             : ' + str(args.segment_overlap_fraction))
    print('Alignment/Block overlap fraction     : ' + str(args.alignment_overlap_fraction))
    if args.gtf:
        print('Annotation file                      : ' + args.gtf)
    if args.block_based:
        print('Merge method                         : Blockbuster based')
        print('Blockbuster distance                 : ' + str(args.distance))
        print('Blockbuster minClusterHeight         : ' + str(args.min_cluster_height))
        print('Blockbuster minBlockHeight           : ' + str(args.min_block_height))
        print('Blockbuster scale                    : ' + str(args.scale))
    else:
        print('Merge method                         : overlap based')
        print('Minimum locus size                   : ' + str(args.min_locus_size))
    print('Minimum alignment length as % of longest : ' + str(args.length_threshold))
    print('Chimeric overlap                     : ' + str(args.chimeric_overlap))
    if args.ref_fasta1:
        print('1st priority reference fasta file    : ' + args.ref_fasta1)
    if args.ref_fasta2:
        print('2nd priority reference fasta file    : ' + args.ref_fasta2)
    print("===================================================================")


def run_merge(args):
    os.makedirs(args.outdir, exist_ok=True)
    d_reflen1 = {}
    d_reflen2 = {}
    if args.ref_fasta1 and args.ref_fasta2:
        chira_utilities.extract_reflengths(args.ref_fasta1, d_reflen1)
        chira_utilities.extract_reflengths(args.ref_fasta2, d_reflen2)

    chira_utilities.print_w_time("START: Reads to segments")
    reads_to_segments(args.bed, args.outdir, args.segment_overlap_fraction, args.chimeric_overlap,
                      d_reflen1.keys(), d_reflen2.keys(), args.chimeric_only, args.length_threshold)
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
        os.replace(os.path.join(args.outdir, 'segments.temp.bed'), os.path.join(args.outdir, 'segments.bed'))
    if args.block_based:
        chira_utilities.print_w_time("START: Blockbuster based merging")
        merge_loci_blockbuster(args.outdir, args.distance, args.min_cluster_height, args.min_block_height,
                               args.scale, args.alignment_overlap_fraction)
        chira_utilities.print_w_time("END: Blockbuster based merging")
    else:
        chira_utilities.print_w_time("START: overlap based merging")
        merge_loci_overlap(args.outdir, args.alignment_overlap_fraction, args.min_locus_size)
        chira_utilities.print_w_time("END: overlap based merging")


def main(argv=None):
    parser = create_parser()
    args = parser.parse_args(argv)
    try:
        validate_args(args)
        print_args(args)
        run_merge(args)
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        parser.exit(2, '{}: error: {}\n'.format(parser.prog, error))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
