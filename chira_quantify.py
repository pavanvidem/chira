#!/usr/bin/env python
import argparse
from collections import defaultdict
import math
import os
import subprocess

import chira_utilities


def validate_score(value, name):
    if not math.isfinite(value) or not 0 <= value <= 1:
        raise ValueError('{} must be finite and in [0, 1]'.format(name))


def bed_length(start, end):
    """Return the length of a non-empty, zero-based half-open interval."""
    start, end = int(start), int(end)
    if start < 0 or end <= start:
        raise ValueError('Invalid BED interval: {}-{}'.format(start, end))
    return end - start


def read_loci(merged_bed):
    """Read locus positions, segment sets and their transcript alignments."""
    positions, reads = [], []
    transcripts = defaultdict(set)
    with open(merged_bed) as handle:
        for locus_id, line in enumerate(handle):
            reference, start, end, strand, descriptions = line.rstrip('\n').split('\t')
            bed_length(start, end)
            positions.append(':'.join([reference, start, end, strand]))
            segments = set()
            for description in descriptions.split(';'):
                segment, transcript, tx_start, tx_end, tx_strand, cigar = description.split(',')
                segments.add(segment)
                transcripts[segment, locus_id].add('\t'.join([transcript, tx_start, tx_end, tx_strand, cigar]))
            reads.append(segments)
    return positions, reads, transcripts


def group_loci(locus_reads, build_crls_too, cutoff, min_locus_size):
    """Build CRLs by Jaccard similarity, then retain best locus shares.

    A locus may belong to multiple equally good CRLs. Below-threshold loci
    and loci without a qualifying membership get their own CRL.
    """
    validate_score(cutoff, 'CRL share')
    if min_locus_size <= 0:
        raise ValueError('Minimum locus size must be positive')
    if any(not reads for reads in locus_reads):
        raise ValueError('A locus must contain at least one read segment')
    ordered = sorted(enumerate(locus_reads), key=lambda item: len(item[1]), reverse=True)
    groups, group_reads = [], []
    remaining = []
    for locus_id, reads in ordered:
        if not build_crls_too or len(reads) < min_locus_size:
            remaining.append((locus_id, reads))
            continue
        matches = []
        for group_id in range(len(groups) - 1, -1, -1):
            members = group_reads[group_id]
            # Necessary size bounds for Jaccard >= cutoff, including 0 and 1.
            if cutoff and not len(reads) * cutoff <= len(members) <= len(reads) / cutoff:
                continue
            common = len(reads & members)
            if common and common / float(len(reads | members)) >= cutoff:
                matches.append(group_id)
        if not matches:
            matches = [len(groups)]
            groups.append({})
            group_reads.append(set())
        for group_id in matches:
            groups[group_id][locus_id] = reads
            group_reads[group_id].update(reads)
    groups.extend({locus_id: reads} for locus_id, reads in remaining)

    shares = defaultdict(dict)
    for group_id, group in enumerate(groups):
        union = set().union(*group.values())
        for locus_id, reads in group.items():
            shares[locus_id][group_id] = len(reads) / float(len(union))
    retained_loci = set()
    for group_id, group in enumerate(groups):
        for locus_id in list(group):
            share = shares[locus_id][group_id]
            if share < cutoff or share < max(shares[locus_id].values()):
                del group[locus_id]
            else:
                retained_loci.add(locus_id)
    groups = [group for group in groups if group]
    groups.extend({locus_id: locus_reads[locus_id]} for locus_id in sorted(shares)
                  if locus_id not in retained_loci)
    return groups


def read_segment_positions(bed):
    positions = {}
    with open(bed) as handle:
        for line in handle:
            reference, start, end, description, score, strand = line.rstrip('\n').split('\t')
            bed_length(start, end)
            segment, transcript, tx_start, tx_end, tx_strand, cigar = description.split(',')
            transcript_position = '\t'.join([transcript, tx_start, tx_end, tx_strand, cigar])
            positions[transcript_position, segment] = ':'.join([reference, start, end, strand])
    return positions


def build_crls(build_crls_too, bed, merged_bed, crl_file, crl_share_cutoff, min_locus_size):
    positions, locus_reads, transcripts = read_loci(merged_bed)
    groups = group_loci(locus_reads, build_crls_too, crl_share_cutoff, min_locus_size)
    segment_positions = read_segment_positions(bed)
    print('Number of loci: {}; number of CRLs: {}'.format(len(locus_reads), len(groups)))
    with open(crl_file, 'w') as output:
        for group_id, group in enumerate(groups):
            union = set().union(*group.values())
            for locus_id, reads in group.items():
                share = len(reads) / float(len(union))
                for segment in sorted(reads):
                    for transcript_position in sorted(transcripts[segment, locus_id]):
                        try:
                            genomic_position = segment_positions[transcript_position, segment]
                        except KeyError:
                            raise ValueError('No segment BED position for {!r} at locus {}'.format(segment, locus_id))
                        transcript, rest = transcript_position.split('\t', 1)
                        output.write('\t'.join([segment, transcript, str(locus_id), str(group_id), rest,
                                                genomic_position, positions[locus_id], '{:.2f}'.format(share)]) + '\n')


def em(d_alpha, d_rho, library_size, l_multimap_readids, em_threshold):
    """Update segment memberships and relative abundances until convergence.

    Each unique segment contributes one count, independently of its collapsed
    FASTA read count. Existing mapping arguments are updated in place.
    """
    validate_score(em_threshold, 'EM threshold')
    if not d_alpha:
        return d_alpha
    if not math.isfinite(library_size) or library_size <= 0:
        raise ValueError('Library size must be finite and positive')
    for iteration in range(1, 1001):
        old_rho = dict(d_rho)
        for read_id in l_multimap_readids:
            memberships = d_alpha[read_id]
            if len(memberships) <= 1:
                continue
            total = sum(old_rho[group_id] for group_id in memberships)
            if total <= 0 or not math.isfinite(total):
                raise ValueError('Invalid total abundance for segment {!r}'.format(read_id))
            for group_id in memberships:
                memberships[group_id] = old_rho[group_id] / total
        counts = defaultdict(float)
        for memberships in d_alpha.values():
            for group_id, fraction in memberships.items():
                counts[group_id] += fraction
        d_rho.clear()
        d_rho.update((group_id, count / library_size) for group_id, count in counts.items())
        difference = sum(abs(d_rho.get(group_id, 0) - old_rho.get(group_id, 0))
                         for group_id in old_rho.keys() | d_rho.keys())
        if difference <= em_threshold:
            return d_alpha
    raise ValueError('EM did not converge within 1000 iterations (difference: {})'.format(difference))


def tpm(d_crl_expression, d_crl_loci_len):
    """Normalize expression by median half-open locus length and total RPK."""
    rates = {}
    for group_id, expression in sorted(d_crl_expression.items()):
        lengths = list(d_crl_loci_len[group_id].values())
        if not lengths or any(length <= 0 or not math.isfinite(length) for length in lengths):
            raise ValueError('CRL {!r} must have positive, finite locus lengths'.format(group_id))
        if expression < 0 or not math.isfinite(expression):
            raise ValueError('CRL expression must be finite and non-negative')
        rates[group_id] = expression / (chira_utilities.median(sorted(lengths)) / 1000.0)
    total = sum(rates.values())
    if not total:
        return {group_id: 0.0 for group_id in rates}
    return {group_id: rate / total * 1000000.0 for group_id, rate in rates.items()}


def quantify_crls(crl_file, em_threshold):
    validate_score(em_threshold, 'EM threshold')
    lengths = defaultdict(dict)
    alpha = defaultdict(dict)
    rho = defaultdict(float)
    with open(crl_file) as handle:
        for line in handle:
            fields = line.rstrip('\n').split('\t')
            if len(fields) != 11:
                raise ValueError('Expected 11 columns in CRL file')
            segment, locus_id, group_id = fields[0], fields[2], fields[3]
            reference, start, end, strand = fields[9].rsplit(':', 3)
            alpha[segment][group_id] = 1.0
            lengths[group_id][locus_id] = bed_length(start, end)
    multimap = []
    for segment, memberships in alpha.items():
        for group_id in memberships:
            memberships[group_id] = 1.0 / len(memberships)
            rho[group_id] += memberships[group_id]
        if len(memberships) > 1:
            multimap.append(segment)
    library_size = len(alpha)
    for group_id in rho:
        rho[group_id] /= library_size
    em(alpha, rho, library_size, multimap, em_threshold)
    expression = defaultdict(float)
    for memberships in alpha.values():
        for group_id, fraction in memberships.items():
            expression[group_id] += fraction
    return alpha, tpm(expression, lengths)


def create_parser():
    parser = argparse.ArgumentParser(description='Chimeric Read Annotator: quantify mapped loci',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-b', '--bed', required=True, metavar='BED', help='Input segment BED file')
    parser.add_argument('-m', '--merged_bed', required=True, metavar='BED', help='Input merged locus file')
    parser.add_argument('-o', '--outdir', required=True, metavar='DIR', help='Output directory')
    parser.add_argument('-cs', '--crl_share', type=chira_utilities.score_float, default=0.7,
                        help='Minimum Jaccard similarity for grouping and minimum locus share within a CRL')
    parser.add_argument('-ls', '--min_locus_size', type=int, default=10,
                        help='Minimum number of read segments per locus to participate in CRL creation')
    parser.add_argument('-e', '--em_threshold', dest='em_thresh', type=chira_utilities.score_float, default=0.00001,
                        help='Maximum sum of absolute changes in relative CRL abundance for EM convergence')
    parser.add_argument('-crl', '--build_crls_too', action='store_true', help='Group loci into CRLs')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.4.3')
    return parser


def validate_args(args):
    if args.min_locus_size <= 0:
        raise ValueError('--min_locus_size must be positive')
    for path in (args.bed, args.merged_bed):
        if not os.path.isfile(path):
            raise ValueError('Input file does not exist or is not a regular file: {}'.format(path))


def write_counts(crl_file, counts_file, fractions, tpms):
    with open(crl_file) as source, open(counts_file, 'w') as output:
        for line in source:
            fields = line.rstrip('\n').split('\t')
            output.write('\t'.join([line.rstrip('\n'), '{:.2f}'.format(fractions[fields[0]][fields[3]]),
                                    '{:.4g}'.format(tpms[fields[3]])]) + '\n')


def run_quantification(args):
    os.makedirs(args.outdir, exist_ok=True)
    crl_file = os.path.join(args.outdir, 'loci.txt')
    temporary = os.path.join(args.outdir, 'loci.counts.temp')
    chira_utilities.print_w_time('START: Build CRLs')
    build_crls(args.build_crls_too, args.bed, args.merged_bed, crl_file, args.crl_share, args.min_locus_size)
    chira_utilities.print_w_time('END: Build CRLs')
    chira_utilities.print_w_time('START: Quantify CRLs')
    fractions, tpms = quantify_crls(crl_file, args.em_thresh)
    write_counts(crl_file, temporary, fractions, tpms)
    with open(os.path.join(args.outdir, 'loci.counts'), 'w') as output:
        subprocess.run(['sort', '-V', temporary], stdout=output, check=True)
    os.remove(crl_file)
    os.remove(temporary)
    chira_utilities.print_w_time('END: Quantify CRLs')


def main(argv=None):
    parser = create_parser()
    args = parser.parse_args(argv)
    try:
        validate_args(args)
        print('Input BED: {}; merged BED: {}; output: {}'.format(args.bed, args.merged_bed, args.outdir))
        print('Minimum locus size: {}; CRL share: {}; EM threshold: {}; build CRLs: {}'.format(
            args.min_locus_size, args.crl_share, args.em_thresh, args.build_crls_too))
        run_quantification(args)
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        parser.exit(2, '{}: error: {}\n'.format(parser.prog, error))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
