#!/usr/bin/env python
import os
import sys
from collections import defaultdict
import datetime
import argparse
sys.setrecursionlimit(10000)


def create_crls(create_crls_too, bed, merged_bed, crl_file, crl_share, min_locus_size):
    l_locilen = []
    l_locipos = []
    l_locireads = []
    d_readlocus_transcripts = defaultdict(lambda: defaultdict(list))
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
                d_readlocus_transcripts[segmentid][n_locus].append(transcriptid_pos)
                if segmentid not in l_locusreads:
                    l_locusreads.add(segmentid)
            l_locireads.append(set(l_locusreads))

    d_crl_reads = defaultdict(list)
    d_crl_locus_reads = defaultdict(lambda: defaultdict())
    l_remaining_locireads = defaultdict(list)
    print(str(datetime.datetime.now()), "start 1st iteration of CRLs")
    print("Number of loci: " + str(len(l_locireads)))
    n_crl = 0
    # create CRLs only if required
    if create_crls_too:
        l_qualified_locireads = sorted([(i, item) for i, item in enumerate(l_locireads) if len(item) >= min_locus_size],
                                       key=lambda x: len(x[1]), reverse=True)
        l_remaining_locireads = sorted([(i, item) for i, item in enumerate(l_locireads) if len(item) < min_locus_size],
                                       key=lambda x: len(x[1]), reverse=True)

        print("Number of qualified loci: ", len(l_qualified_locireads))

        for n_locus, l_locusreads in l_qualified_locireads:
            already_crl_member = False
            l_matched_crls = []
            # substraction of 0.1 is a heuristic to loosly filter out the loci
            lower_bound = len(l_locusreads) * (crl_share - 0.1)
            upper_bound = len(l_locusreads) / (crl_share - 0.1)
            # traverse in reverse order because the latest CRL is the last one
            for crlid in range(len(d_crl_reads) - 1, 0, -1):
                l_crlreads = set(d_crl_reads[crlid])
                if lower_bound <= len(l_crlreads) <= upper_bound:
                    # if there are significantly more reads in crl than in locus or otherway around
                    n_common_reads = len(l_locusreads.intersection(l_crlreads))
                    if n_common_reads == 0:
                        continue
                    n_union_reads = len(l_crlreads.union(l_locusreads))
                    # jaccard similarity score, skip the whole crl even if one of the locus doesn't overlap enough
                    if n_common_reads / float(n_union_reads) < crl_share:
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
    print(str(datetime.datetime.now()), " done 1st iteration of CRLs")

    d_locus_crl_share = defaultdict(lambda: defaultdict(float))
    for crlid, d_crlloci in d_crl_locus_reads.items():
        l_crlreads = []
        for l_locusreads in d_crlloci.values():
            l_crlreads.extend(l_locusreads)
        for locusid, l_locusreads in d_crlloci.items():
            locus_share = len(l_locusreads) / float(len(set(l_crlreads)))
            d_locus_crl_share[locusid][crlid] = locus_share

    d_isolated_loci = {}  # loci that are separated from crl because of not enough overall share
    print(str(datetime.datetime.now()), "start 2st iteration of CRLs")
    # in this iteration each locus is checked again for crl share and only one best crl for locus kept
    # if there are multiple eqaully good crls for a locus, all of the crls are considered for that locus
    for crlid, d_crlloci in d_crl_locus_reads.items():
        for locusid in list(d_crlloci):
            if d_locus_crl_share[locusid][crlid] >= crl_share:
                for any_crlid in d_locus_crl_share[locusid]:
                    if d_locus_crl_share[locusid][crlid] < d_locus_crl_share[locusid][any_crlid]:
                        if locusid in d_crl_locus_reads[crlid]:
                            del d_crl_locus_reads[crlid][locusid]
            else:
                # store isolated locus id and its reads
                d_isolated_loci[locusid] = d_crlloci[locusid]
                del d_crl_locus_reads[crlid][locusid]
    print(str(datetime.datetime.now()), "End: finding isolated CRL loci")
    print(str(datetime.datetime.now()), "done 2st iteration of CRLs")
    print(str(datetime.datetime.now()), "Start: creating CRLs with isolated loci")
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

    print(str(datetime.datetime.now()), "End: creating CRLs with isolated loci")

    print(str(datetime.datetime.now()), "End: creating CRLs")
    print("There are a total of " + str(len(d_crl_locus_reads)) + " uniq crls")

    # read the segments BED once again to get the genomic positions
    # to reduce the memory usage
    d_read_genomic_pos = defaultdict(lambda: defaultdict(str))
    with open(bed) as fh_bed:
        for line in fh_bed:
            b = line.rstrip('\n').split('\t')
            pos = ':'.join([b[0], b[1], b[2], b[5]])
            desc = b[3].split(',')  # in the description column of the bed file,alignments are seperated by ';'
            readid = desc[0]
            transcriptid_pos = '\t'.join(desc[1:])
            # at this level reads have unique ids preceeded by a serialnumber
            # each read can have multiple alignements on same transcript
            d_read_genomic_pos[transcriptid_pos][readid] = pos

    print(str(datetime.datetime.now()), "Start: Writing CRLs")
    d_duplicate_entries = {}
    with open(crl_file, "w") as fh_groups_file:
        # enumerate again because of above processing some crlids might be missing
        for crlid, d_crlloci in enumerate(d_crl_locus_reads.values()):
            if len(d_crlloci) > 1:
                l_crlreads = set([val for sublist in d_crlloci.values() for val in sublist])
            for locusid, l_locusreads in d_crlloci.items():
                if len(d_crlloci) > 1:
                    locus_share = len(l_locusreads) / float(len(l_crlreads))
                else:
                    locus_share = 1.0
                for readid in sorted(l_locusreads):
                    for transcriptid_pos in sorted(d_readlocus_transcripts[readid][locusid]):
                        entry = "\t".join([readid,
                                           transcriptid_pos.split('\t')[0],
                                           str(locusid),
                                           str(crlid),
                                           '\t'.join(transcriptid_pos.split('\t')[1:]),
                                           d_read_genomic_pos[transcriptid_pos][readid],
                                           l_locipos[locusid],
                                           "{:.4g}".format(locus_share)]
                                          )
                        if entry in d_duplicate_entries:
                            continue
                        fh_groups_file.write(entry + "\n")
                        d_duplicate_entries[entry] = 1
    print(str(datetime.datetime.now()), "End: CRLs written")


def em(d_read_group_fractions, d_group_counts, em_threshold, i=1):
    print("iteration: " + str(i))
    d_read_group_fractions_new = {}
    d_group_counts_new = {}

    for readid in d_read_group_fractions.keys():
        total_group_count = 0
        d_read_group_fractions_new[readid] = {}
        for groupid in d_read_group_fractions[readid]:
            total_group_count += d_group_counts[groupid]
        for groupid in d_read_group_fractions[readid]:
            d_read_group_fractions_new[readid][groupid] = d_group_counts[groupid] / float(total_group_count)
    d_read_group_fractions.clear()

    for readid in d_read_group_fractions_new.keys():
        for groupid in d_read_group_fractions_new[readid]:
            if groupid not in d_group_counts_new:
                d_group_counts_new[groupid] = 0
            d_group_counts_new[groupid] += d_read_group_fractions_new[readid][groupid]

    equal = True
    for groupid in d_group_counts.keys():
        if abs(d_group_counts[groupid] - d_group_counts_new[groupid]) >= em_threshold:
            equal = False
            break
    if equal:
        return d_read_group_fractions_new
    else:
        d_group_counts.clear()
        i += 1
        return em(d_read_group_fractions_new, d_group_counts_new, em_threshold, i)


def median(x):
    n = len(x)
    mid = int(n/2)
    if not n % 2:
        return (x[mid-1] + x[mid]) / 2.0
    return x[mid]


def tpm(d_group_expression, d_group_locilen):
    total_rpk = 0
    d_group_tpm = defaultdict(float)
    for groupid in sorted(d_group_expression.keys()):
        crl_expression = d_group_expression[groupid]
        crl_len = median(sorted(d_group_locilen[groupid].values())) / 1000.0  # length in kbs
        rpk = crl_expression / crl_len
        d_group_tpm[groupid] = rpk
        total_rpk += rpk
    millions_of_rpk = total_rpk / 1000000.0
    for groupid in sorted(d_group_expression.keys()):
        group_tpm = d_group_tpm[groupid] / millions_of_rpk
        d_group_tpm[groupid] = group_tpm
    return d_group_tpm


def quantify_crls(crl_file, em_threshold):
    d_read_group_fractions = defaultdict(lambda: defaultdict(float))
    d_group_locilen = defaultdict(lambda: defaultdict(int))
    d_group_expression = defaultdict(float)

    print(str(datetime.datetime.now()), "Start: preparing EM")
    fh_crl_file = open(crl_file, "r")
    for line in fh_crl_file:
        f = line.rstrip("\n").split("\t")
        # consider the segment id and quantify individula segments than whole reads
        readid = f[0]
        locusid = f[2]
        groupid = f[3]
        pos = f[9].split(':')
        locuslength = int(pos[-2]) - int(pos[-3]) + 1
        # a single locus can belong to multiple crls
        # one read can be part of multiple crls
        d_read_group_fractions[readid][groupid] = 0
        d_group_locilen[groupid][locusid] = locuslength
    fh_crl_file.close()

    d_group_counts = {}
    for readid in d_read_group_fractions.keys():
        for groupid in d_read_group_fractions[readid].keys():
            d_read_group_fractions[readid][groupid] = 1 / float(len(d_read_group_fractions[readid]))
            if groupid not in d_group_counts:
                d_group_counts[groupid] = 0
            d_group_counts[groupid] += d_read_group_fractions[readid][groupid]

    print(str(datetime.datetime.now()), "End: preparing EM")
    print(str(datetime.datetime.now()), "Start: recursing EM")
    sys.setrecursionlimit(10000)
    d_res = em(d_read_group_fractions, d_group_counts, em_threshold)
    print(str(datetime.datetime.now()), "End: recursing EM")
    for readid in d_res.keys():
        # now d_read_groups should contain only multi mapped reads because because
        # uniquely mapped reads were already removed
        # set() because a read can occur 2 times at a locus in file with 2 tx ids
        # one read can be part of multiple crls
        for groupid in d_res[readid]:
            d_group_expression[groupid] += d_res[readid][groupid]

    d_group_tpm = tpm(d_group_expression, d_group_locilen)

    return d_res, d_group_tpm


def score_float(x):
    x = float(x)
    if x < 0.0 or x > 1.0:
        raise argparse.ArgumentTypeError("%r not in range [0.0, 1.0]" % (x,))
    return x


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

    parser.add_argument('-cs', '--crl_share', action='store', type=score_float, default=0.7, metavar='',
                        dest='crl_share',
                        help='Minimum fraction of reads of a locus that must overlap with all CRL loci '
                             'inorder to merge it into that CRL.')

    parser.add_argument('-ls', '--min_locus_size', action='store', type=int, default=5, metavar='',
                        dest='min_locus_size',
                        help='Minimum number of reads a locus should have in order to participate in CRL creation.'
                             'Always set this value relative to your sequencing depth. Setting this to lower leads'
                             'CRLs of random multimappings Also consider setting the --crl_share option '
                             'along with this')

    parser.add_argument('-e', '--em_threshold', action='store', type=score_float, default=1, metavar='',
                        dest='em_thresh',
                        help='The maximum difference of transcripts expression between two consecutive iterations '
                             'of EM algorithm to converge.')

    parser.add_argument("-crl", '--create_crls_too', action='store_true', dest='create_crls_too',
                        help="Create CRLs too")

    parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.0')

    args = parser.parse_args()

    print('Input BED file                       : ' + args.bed)
    print('Input merged BED file                : ' + args.merged_bed)
    print('Output directory                     : ' + args.outdir)
    print('Minimum locus size                   : ' + str(args.min_locus_size))
    print('CRL share                            : ' + str(args.crl_share))
    print('EM threshold                         : ' + str(args.em_thresh))
    print('Create CRLs todo                     : ' + str(args.create_crls_too))
    print("===================================================================")

    create_crls(args.create_crls_too, args.bed, args.merged_bed,
                os.path.join(args.outdir, 'loci.txt'), args.crl_share, args.min_locus_size)
    d_read_crl_fractions, d_crl_tpms = quantify_crls(os.path.join(args.outdir, 'loci.txt'), args.em_thresh)
    with open(os.path.join(args.outdir, 'loci.txt')) as fh_in:
        with open(os.path.join(args.outdir, 'loci.counts'), "w") as fh_out:
            for l in fh_in:
                k = l.rstrip("\n").split("\t")
                read_id = k[0]
                crl_id = k[3]
                fh_out.write("\t".join([l.strip("\n"),
                                        "{:.4g}".format(d_read_crl_fractions[read_id][crl_id]),
                                        "{:.4g}".format(d_crl_tpms[crl_id])]) + "\n")
