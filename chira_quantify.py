#!/usr/bin/env python
import os
import sys
from collections import defaultdict
import argparse
import chira_utilities
sys.setrecursionlimit(10000)


def build_crls(build_crls_too, bed, merged_bed, crl_file, crl_share_cutoff, min_locus_size):
    l_locilen = []
    l_locipos = []
    l_locireads = []
    d_readlocus_transcripts = defaultdict(list)
    with open(merged_bed) as fh_merged_bed:
        for n_locus, line in enumerate(fh_merged_bed):
            # chr14\t64814786\t64814804\t-\ttag_1308593|1|r,ENSMUST00000176386;tag_1308594|2|r,ENSMUST00000176386
            f = line.rstrip('\n').split('\t')
            locuslen = float(f[2])-float(f[1])+1
            l_locilen.append(locuslen)
            pos = ':'.join([f[0], f[1], f[2], f[3]])
            l_locipos.append(pos)
            alignments = f[4].split(';')  # in the description column of the bed file,alignments are seperated by ';'
            l_locusreads = set()
            for alignment in alignments:
                [segmentid, transcriptid, start, end, tx_strand, cigar] = alignment.split(',')
                transcriptid = alignment.split(',')[1]
                transcriptid_pos = '\t'.join([transcriptid, start, end, tx_strand, cigar])
                d_readlocus_transcripts[segmentid+str(n_locus)].append(transcriptid_pos)
                if segmentid not in l_locusreads:
                    l_locusreads.add(segmentid)
            l_locireads.append(set(l_locusreads))

    d_crl_reads = defaultdict(list)
    d_crl_locus_reads = defaultdict(lambda: defaultdict())
    l_remaining_locireads = defaultdict(list)
    chira_utilities.print_w_time("START: 1st iteration of CRLs")
    print("Number of loci: " + str(len(l_locireads)))
    n_crl = 0
    # create CRLs only if required
    if build_crls_too:
        l_qualified_locireads = sorted([(i, item) for i, item in enumerate(l_locireads) if len(item) >= min_locus_size],
                                       key=lambda x: len(x[1]), reverse=True)
        l_remaining_locireads = sorted([(i, item) for i, item in enumerate(l_locireads) if len(item) < min_locus_size],
                                       key=lambda x: len(x[1]), reverse=True)

        print("Number of qualified loci: ", len(l_qualified_locireads))

        for n_locus, l_locusreads in l_qualified_locireads:
            already_crl_member = False
            l_matched_crls = []
            # lower and uppder bounds for filtering crls based on their size
            lower_bound = len(l_locusreads) * (1 - crl_share_cutoff)
            upper_bound = len(l_locusreads) / (1 - crl_share_cutoff)
            # traverse in reverse order because the latest CRL is the last one
            for crlid in range(len(d_crl_reads) - 1, 0, -1):
                l_crlreads = set(d_crl_reads[crlid])
                # if the CRL has similar size
                if lower_bound <= len(l_crlreads) <= upper_bound:
                    # if there are significantly more reads in crl than in locus or otherway around
                    n_common_reads = len(l_locusreads.intersection(l_crlreads))
                    if n_common_reads == 0:
                        continue
                    n_union_reads = len(l_crlreads.union(l_locusreads))
                    # jaccard similarity score
                    if n_common_reads / float(n_union_reads) < crl_share_cutoff:
                        continue
                    # identical loci are multi-mapped loci with the same set of identical set of multi-mapped reads
                    l_matched_crls.append(crlid)
                    already_crl_member = True
                else:
                    break

            # n_locus is not a member of any crl, hence create a new crl
            if not already_crl_member:
                l_matched_crls.append(n_crl)
                n_crl += 1
            for matched_crl in l_matched_crls:
                l_reads_temp = d_crl_reads[matched_crl]
                d_crl_locus_reads[matched_crl][n_locus] = l_locusreads
                l_reads_temp.extend(l_locusreads)
                d_crl_reads[matched_crl] = l_reads_temp
        d_crl_reads.clear()
    else:
        l_remaining_locireads = sorted(enumerate(l_locireads), key=lambda x: len(x[1]), reverse=True)
    # every locus with less than min_locus_size reads gets its' own locus
    for n_locus, l_locusreads in l_remaining_locireads:
        d_crl_locus_reads[n_crl][n_locus] = l_locusreads
        n_crl += 1
    chira_utilities.print_w_time("END: 1st iteration of CRLs")

    d_locus_crl_share = defaultdict(lambda: defaultdict(float))
    d_highest_shares = defaultdict(float)
    for crlid, d_crlloci in d_crl_locus_reads.items():
        l_crlreads = []
        for l_locusreads in d_crlloci.values():
            l_crlreads.extend(l_locusreads)
        for locusid, l_locusreads in d_crlloci.items():
            locus_share = len(l_locusreads) / float(len(set(l_crlreads)))
            d_locus_crl_share[locusid][crlid] = locus_share
            if locus_share > d_highest_shares[locusid]:
                d_highest_shares[locusid] = locus_share

    d_isolated_loci = {}  # loci that are separated from crl because of not enough overall share
    chira_utilities.print_w_time("START: 2st iteration of CRLs")
    # in this iteration each locus is checked again for crl share and only one best crl for locus kept
    # if there are multiple eqaully good crls for a locus, all of the crls are considered for that locus
    for crlid in list(d_crl_locus_reads.keys()):
        d_crlloci = d_crl_locus_reads[crlid]
        for locusid in list(d_crlloci):
            if d_locus_crl_share[locusid][crlid] >= crl_share_cutoff:
                if d_locus_crl_share[locusid][crlid] < d_highest_shares[locusid]:
                    if locusid in d_crl_locus_reads[crlid]:
                        del d_crl_locus_reads[crlid][locusid]
                        # del d_locus_crl_share[locusid][crlid]
            else:
                # store isolated locus id and its reads
                d_isolated_loci[locusid] = d_crlloci[locusid]
                del d_crl_locus_reads[crlid][locusid]
                del d_locus_crl_share[locusid]

    chira_utilities.print_w_time("END: 2st iteration of CRLs")
    chira_utilities.print_w_time("START: creating CRLs with isolated loci")
    # every locus that is in isolated_loci makes its own crl
    crl_index = len(d_crl_locus_reads)
    for locusid, l_locusreads in d_isolated_loci.items():
        # already member of some crl, safely ignore it
        if locusid in d_locus_crl_share:
            continue
        d_crl_locus_reads[crl_index][locusid] = l_locusreads
        d_locus_crl_share[locusid][crl_index] = 1
        crl_index += 1
    d_locus_crl_share.clear()

    chira_utilities.print_w_time("END: creating CRLs with isolated loci")

    print("There are a total of " + str(len(d_crl_locus_reads)) + " uniq crls")

    # read the segments BED to get the genomic positions
    d_read_genomic_pos = defaultdict(str)
    with open(bed) as fh_bed:
        for line in fh_bed:
            b = line.rstrip('\n').split('\t')
            pos = ':'.join([b[0], b[1], b[2], b[5]])
            desc = b[3].split(',')  # in the description column of the bed file,alignments are seperated by ';'
            segmentid = desc[0]
            transcriptid_pos = '\t'.join(desc[1:])
            # at this level reads have unique ids preceeded by a serialnumber
            # each read can have multiple alignements on same transcript
            d_read_genomic_pos[transcriptid_pos+segmentid] = pos

    with open(crl_file, "w") as fh_crl_file:
        # enumerate again because of above processing some crlids might be missing
        for crlid, d_crlloci in enumerate(d_crl_locus_reads.values()):
            if len(d_crlloci) > 1:
                l_crlreads = set([val for sublist in d_crlloci.values() for val in sublist])
            for locusid, l_locusreads in d_crlloci.items():
                if len(d_crlloci) > 1:
                    locus_share = len(l_locusreads) / float(len(l_crlreads))
                else:
                    locus_share = 1.0
                for segmentid in sorted(l_locusreads):
                    for transcriptid_pos in sorted(d_readlocus_transcripts[segmentid+str(locusid)]):
                        entry = "\t".join([segmentid,
                                           transcriptid_pos.split('\t')[0],
                                           str(locusid),
                                           str(crlid),
                                           '\t'.join(transcriptid_pos.split('\t')[1:]),
                                           d_read_genomic_pos[transcriptid_pos+segmentid],
                                           l_locipos[locusid],
                                           "{:.2f}".format(locus_share)]
                                          )
                        fh_crl_file.write(entry + "\n")


def em(d_read_crl_shares, d_crl_counts, l_multimap_readids, em_threshold, i=1):
    print("iteration: " + str(i))
    d_crl_counts_new = {}

    for readid in l_multimap_readids:
        if len(d_read_crl_shares[readid]) == 1:
            continue
        total_crl_count = 0
        for crlid in d_read_crl_shares[readid]:
            total_crl_count += d_crl_counts[crlid]
        for crlid in d_read_crl_shares[readid]:
            d_read_crl_shares[readid][crlid] = d_crl_counts[crlid] / float(total_crl_count)

    for readid in d_read_crl_shares.keys():
        for crlid in d_read_crl_shares[readid]:
            if crlid not in d_crl_counts_new:
                d_crl_counts_new[crlid] = 0
            d_crl_counts_new[crlid] += d_read_crl_shares[readid][crlid]

    equal = True
    for crlid in d_crl_counts.keys():
        if abs(d_crl_counts[crlid] - d_crl_counts_new[crlid]) >= em_threshold:
            equal = False
            break
    if equal:
        return d_read_crl_shares
    else:
        d_crl_counts.clear()
        i += 1
        return em(d_read_crl_shares, d_crl_counts_new, l_multimap_readids, em_threshold, i)


def tpm(d_crl_expression, d_crl_loci_len):
    total_rpk = 0
    d_crl_tpm = defaultdict(float)
    for crlid in sorted(d_crl_expression.keys()):
        crl_expression = d_crl_expression[crlid]
        crl_len = chira_utilities.median(sorted(d_crl_loci_len[crlid].values())) / 1000.0  # length in kbs
        rpk = crl_expression / crl_len
        d_crl_tpm[crlid] = rpk
        total_rpk += rpk
    millions_of_rpk = total_rpk / 1000000.0
    for crlid in sorted(d_crl_expression.keys()):
        crl_tpm = d_crl_tpm[crlid] / millions_of_rpk
        d_crl_tpm[crlid] = crl_tpm
    return d_crl_tpm


def quantify_crls(crl_file, em_threshold):
    d_read_crl_shares = defaultdict(lambda: defaultdict(float))
    d_crl_loci_len = defaultdict(lambda: defaultdict(int))
    d_crl_expression = defaultdict(float)
    l_multimap_readids = []

    fh_crl_file = open(crl_file, "r")
    for line in fh_crl_file:
        f = line.rstrip("\n").split("\t")
        # consider the segment id and quantify individula segments than whole reads
        readid = f[0]
        locusid = f[2]
        crlid = f[3]
        pos = f[9].split(':')
        locuslength = int(pos[-2]) - int(pos[-3]) + 1
        # a single locus can belong to multiple crls
        # one read can be part of multiple crls
        d_read_crl_shares[readid][crlid] = 1
        d_crl_loci_len[crlid][locusid] = locuslength
    fh_crl_file.close()

    d_crl_counts = {}

    for readid in d_read_crl_shares.keys():
        for crlid in d_read_crl_shares[readid].keys():
            d_read_crl_shares[readid][crlid] = 1 / float(len(d_read_crl_shares[readid]))
            if crlid not in d_crl_counts:
                d_crl_counts[crlid] = 0
            d_crl_counts[crlid] += d_read_crl_shares[readid][crlid]

        if len(d_read_crl_shares[readid]) > 1:
            l_multimap_readids.append(readid)

    sys.setrecursionlimit(10000)
    d_res = em(d_read_crl_shares, d_crl_counts, l_multimap_readids, em_threshold)
    for readid in d_res.keys():
        for crlid in d_res[readid]:
            d_crl_expression[crlid] += d_res[readid][crlid]

    d_crl_tpm = tpm(d_crl_expression, d_crl_loci_len)

    return d_res, d_crl_tpm


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Chimeric Read Annotator: quantify mapped loci',
                                     usage='%(prog)s [-h] [-v,--version]',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-b', '--bed', action='store', dest='bed', required=True,
                        metavar='', help='Input BED file')

    parser.add_argument('-m', '--merged_bed', action='store', dest='merged_bed', required=True,
                        metavar='', help='Input merged BED file')

    parser.add_argument('-o', '--outdir', action='store', dest='outdir', required=True, metavar='',
                        help='Output file containing merged alignments')

    parser.add_argument('-cs', '--crl_share', action='store', type=chira_utilities.score_float, default=0.7, metavar='',
                        dest='crl_share',
                        help='Minimum fraction of reads of a locus that must overlap with all CRL loci '
                             'inorder to merge it into that CRL.')

    parser.add_argument('-ls', '--min_locus_size', action='store', type=int, default=10, metavar='',
                        dest='min_locus_size',
                        help='Minimum number of reads a locus should have in order to participate in CRL creation.'
                             'Always set this value relative to your sequencing depth. Setting this to lower leads'
                             'CRLs of random multimappings Also consider setting the --crl_share option '
                             'along with this')

    parser.add_argument('-e', '--em_threshold', action='store', type=chira_utilities.score_float, default=1, metavar='',
                        dest='em_thresh',
                        help='The maximum difference of transcripts expression between two consecutive iterations '
                             'of EM algorithm to converge.')

    parser.add_argument("-crl", '--build_crls_too', action='store_true', dest='build_crls_too',
                        help="Create CRLs too")

    parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.3.2')

    args = parser.parse_args()

    print('Input BED file                       : ' + args.bed)
    print('Input merged BED file                : ' + args.merged_bed)
    print('Output directory                     : ' + args.outdir)
    print('Minimum locus size                   : ' + str(args.min_locus_size))
    print('CRL share                            : ' + str(args.crl_share))
    print('EM threshold                         : ' + str(args.em_thresh))
    print('Create CRLs too                      : ' + str(args.build_crls_too))
    print("===================================================================")

    chira_utilities.print_w_time("START: Build CRLs")
    build_crls(args.build_crls_too, args.bed, args.merged_bed,
               os.path.join(args.outdir, 'loci.txt'), args.crl_share, args.min_locus_size)
    chira_utilities.print_w_time("END: Build CRLs")
    chira_utilities.print_w_time("START: Quantify CRLs")
    d_read_crl_fractions, d_crl_tpms = quantify_crls(os.path.join(args.outdir, 'loci.txt'), args.em_thresh)
    chira_utilities.print_w_time("END: Quantify CRLs")
    chira_utilities.print_w_time("START: Write CRLs")
    with open(os.path.join(args.outdir, 'loci.txt')) as fh_in:
        with open(os.path.join(args.outdir, 'loci.counts.temp'), "w") as fh_out:
            for l in fh_in:
                k = l.rstrip("\n").split("\t")
                read_id = k[0]
                crl_id = k[3]
                fh_out.write("\t".join([l.strip("\n"),
                                        "{:.2f}".format(d_read_crl_fractions[read_id][crl_id]),
                                        "{:.4g}".format(d_crl_tpms[crl_id])]) + "\n")
    chira_utilities.print_w_time("END: Write CRLs")
    os.remove(os.path.join(args.outdir, 'loci.txt'))
    chira_utilities.print_w_time("START: Sort CRLs file by read name")
    os.system("sort -V " + os.path.join(args.outdir, 'loci.counts.temp') + " > " + os.path.join(args.outdir, 'loci.counts'))
    os.remove(os.path.join(args.outdir, 'loci.counts.temp'))
    chira_utilities.print_w_time("END: Sort CRLs file by read name")
