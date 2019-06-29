
import os
import sys
import re
import time
from collections import defaultdict
import pysam
import pkg_resources
import datetime
import functools
import operator
import multiprocessing
from multiprocessing import Process, Value, Lock, Pool, Manager, Array
import numpy as np
import itertools
from numpy import array
from functools import partial
import pprint
import utilities

sys.setrecursionlimit(10000)


def process_loci(i, d_sub_loci, d_uniq_crg_locus_reads, d_uniq_crg_locus_lengths, crg_share_threshold, lock):
    d_reads = {}
    n_crg = 0
    d_lengths = {}
    for n_locus, l_locusreads in d_sub_loci:
        already_crg_member = False
        l_matched_crgs = []
        for crgid, d_crgloci in d_reads.items():
            overlap_all_loci = True
            if d_lengths[crgid] * crg_share_threshold <= len(l_locusreads) and \
                    len(l_locusreads) * crg_share_threshold <= d_lengths[crgid]:
                for locusid, l_crg_locureads in d_crgloci.items():
                    # if there are significantly more reads in CRG than in locus or otherway around
                    if len(l_crg_locureads) * crg_share_threshold > len(l_locusreads) or \
                            len(l_locusreads) * crg_share_threshold > len(l_crg_locureads):
                        overlap_all_loci = False
                        break
                    n_common_reads = len(l_locusreads.intersection(l_crg_locureads))
                    if n_common_reads == 0:
                        overlap_all_loci = False
                        break
                    n_union_reads = len(l_crg_locureads.union(l_locusreads))
                    # jaccard similarity score, skip the whole crg even if one of the locus doesn't overlap enough
                    if n_common_reads / float(n_union_reads) < crg_share_threshold:
                        overlap_all_loci = False
                        break
                if overlap_all_loci:
                    # add to l_uniq_crgloci only if there was no 100% identical locus already present in list
                    # identical loci are multi-mapped loci with the same set of identical set of multi-mapped reads
                    # do not add them to l_uniq_crgloci so that next time skip scanning same read-set multiple times
                    # if not identical_locus_found:
                    l_matched_crgs.append(crgid)
                    already_crg_member = True
        # n_locus is not a member of any crg, hence create a new crg
        if not already_crg_member:
            l_matched_crgs.append((str(i) + '_' + str(n_crg)))
            n_crg += 1
        for matched_crg in l_matched_crgs:
            if matched_crg not in d_reads:
                d_reads[matched_crg] = {}
            d_reads[matched_crg][n_locus] = l_locusreads
            t = 0
            for l, r in d_reads[matched_crg].items():
                t += len(r)
            a = t/float(len(d_reads[matched_crg]))
            d_lengths[matched_crg] = a

    # lock.acquire()
    d_uniq_crg_locus_lengths.update(d_lengths)
    d_uniq_crg_locus_reads.update(d_reads)
    # lock.release()


def create_crgs(bed, merge_overlap, groups_file, crg_share_threshold, min_locus_size):
    d_read_genomic_pos = defaultdict(lambda: defaultdict(str))
    fh_bed = open(bed, "r")
    d_desc = defaultdict(lambda: defaultdict(list))
    for line in fh_bed:
        f = line.rstrip('\n').split('\t')
        d_desc[f[0] + "\t" + f[5]][(int(f[1]), int(f[2]))].append(f[3])
        pos = ':'.join([f[0], f[1], f[2], f[5]])
        desc = f[3].split(',')  # in the description column of the bed file,alignments are seperated by ';'
        readid = f[3].split(",")[0]
        transcriptid_pos = '\t'.join(desc[1:])
        # at this level reads have unique ids preceeded by a serialnumber
        # each read can have multiple alignements on same transcript
        d_read_genomic_pos[readid][transcriptid_pos] = pos
    fh_bed.close()
    print("size new:", sys.getsizeof(d_desc))
    print("Read bed file")

    d_readlocus_transcripts = defaultdict(lambda: defaultdict(list))
    l_locipos = []
    d_locireads = {}
    print(str(datetime.datetime.now()), "merging overlapping alignments")
    n_locus = 0
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
                    if overlap_len/float(currentpos[1]-currentpos[0]+1) >= merge_overlap \
                            or overlap_len/float(prevpos[1]-prevpos[0]+1) >= merge_overlap:
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

        for pos in t_merged:
            [chromid, strand] = chrom.split("\t")
            l_locipos.append(':'.join([chromid,
                                       str(pos[0]),
                                       str(pos[1]),
                                       strand]))
            d_longest_alignments = defaultdict(int)
            # choose per locus per transcript only one longest alignment
            for alignment in set(d_mergeddesc[pos[0]]):
                [readid, transcriptid, start, end, tx_strand, cigar] = alignment.split(',')
                if int(end) - int(start) >= d_longest_alignments[readid + "\t" + transcriptid]:
                    d_longest_alignments[readid+"\t"+transcriptid] = int(end) - int(start)
            l_alignments = []
            l_locusreads = set()
            for alignment in set(d_mergeddesc[pos[0]]):
                [readid, transcriptid, start, end, tx_strand, cigar] = alignment.split(',')
                if int(end) - int(start) >= d_longest_alignments[readid+"\t"+transcriptid]:
                    l_alignments.append(alignment)
                    transcriptid_pos = '\t'.join([transcriptid, start, end, tx_strand, cigar])
                    d_readlocus_transcripts[readid][n_locus].append(transcriptid_pos)
                    if readid not in l_locusreads:
                        l_locusreads.add(readid)
            d_locireads[n_locus] = l_locusreads
            # fh_mergedbed.write(chromid + "\t" +
            #                    str(pos[0]) + "\t" +
            #                    str(pos[1]) + "\t" +
            #                    strand + "\t" +
            #                    ";".join(l_alignments) + "\n")
            n_locus += 1


    print(str(datetime.datetime.now()), "start 1st iteration of CRGs")
    n_crg = 0
    d_uniq_crg_locus_reads = defaultdict(lambda: defaultdict())
    for n_locus, l_locusreads in sorted({k: v for k, v in d_locireads.items() if len(v) >= min_locus_size}.items(),
                                        key=lambda kv: (len(kv[1]), kv[0]), reverse=True):
        already_crg_member = False
        l_matched_crgs = []
        for crgid, d_crgloci in d_uniq_crg_locus_reads.items():
            overlap_all_loci = True
            for locusid, l_crg_locureads in d_crgloci.items():
                # if there are significantly more reads in CRG than in locus or otherway around
                if len(l_crg_locureads) * crg_share_threshold > len(l_locusreads) or \
                        len(l_locusreads) * crg_share_threshold > len(l_crg_locureads):
                    overlap_all_loci = False
                    break
                n_common_reads = len(l_locusreads.intersection(l_crg_locureads))
                if n_common_reads == 0:
                    overlap_all_loci = False
                    break
                n_union_reads = len(l_crg_locureads.union(l_locusreads))
                # jaccard similarity score, skip the whole crg even if one of the locus doesn't overlap enough
                if n_common_reads / float(n_union_reads) < crg_share_threshold:
                    overlap_all_loci = False
                    break
            if overlap_all_loci:
                # add to l_uniq_crgloci only if there was no 100% identical locus already present in list
                # identical loci are multi-mapped loci with the same set of identical set of multi-mapped reads
                # do not add them to l_uniq_crgloci so that next time skip scanning same read-set multiple times
                # d_uniq_crg_locus_reads[crgid][n_locus] = l_locusreads
                l_matched_crgs.append(crgid)
                already_crg_member = True
        # n_locus is not a member of any crg, hence create a new crg
        if not already_crg_member:
            # d_uniq_crg_locus_reads[n_crg][n_locus] = l_locusreads
            l_matched_crgs.append(n_crg)
            n_crg += 1
        for matched_crg in l_matched_crgs:
            d_uniq_crg_locus_reads[matched_crg][n_locus] = l_locusreads
    print(str(datetime.datetime.now()), "done 1st iteration of CRGs")

    d_isolated_loci = {}  # loci that are separated from crg because of not enough overall share
    integrated_loci = []  # loci already belong to any crg
    d_locus_crg_share = defaultdict(lambda: defaultdict(float))

    for n_locus, l_locusreads in sorted({k: v for k, v in d_locireads.items() if len(v) < min_locus_size}.items(),
                                        key=lambda kv: (len(kv[1]), kv[0]), reverse=True):
        d_uniq_crg_locus_reads[n_crg][n_locus] = l_locusreads
        n_crg += 1

    for crgid, d_crgloci in d_uniq_crg_locus_reads.items():
        l_crgreads = []
        for l_locusreads in d_crgloci.values():
            l_crgreads.extend(l_locusreads)
        for locusid, l_locusreads in d_crgloci.items():
            locus_share = len(l_locusreads) / float(len(set(l_crgreads)))
            d_locus_crg_share[locusid][crgid] = locus_share

    print(str(datetime.datetime.now()), "start 2st iteration of CRGs")
    # print(str(datetime.datetime.now()), "Start: finding isolated CRG loci")
    for crgid, d_crgloci in d_uniq_crg_locus_reads.items():
        for locusid in list(d_crgloci):
            if d_locus_crg_share[locusid][crgid] >= crg_share_threshold:
                for any_crgid in d_locus_crg_share[locusid]:
                    if d_locus_crg_share[locusid][crgid] < d_locus_crg_share[locusid][any_crgid]:
                        if locusid in d_uniq_crg_locus_reads[crgid]:
                            del d_uniq_crg_locus_reads[crgid][locusid]
                # there is at least a CRG associated with this locus
                if len(d_locus_crg_share[locusid]) > 0:
                    integrated_loci.append(locusid)
            else:
                # store isolated locus id and its reads
                d_isolated_loci[locusid] = d_crgloci[locusid]
                del d_uniq_crg_locus_reads[crgid][locusid]
    # print(str(datetime.datetime.now()), "End: finding isolated CRG loci")
    print(str(datetime.datetime.now()), "done 2st iteration of CRGs")
    crg_index = len(d_uniq_crg_locus_reads)
    # print(str(datetime.datetime.now()), "Start: creating CRGs with isolated loci")
    # every locus that is in isolated_loci makes its own crg
    for locusid, l_locusreads in d_isolated_loci.items():
        if locusid in integrated_loci:
            # already member of some crg, safely ignore it
            continue
        # d_uniq_crg_locus_reads[crg_index] = {locusid: l_locusreads}
        d_uniq_crg_locus_reads[crg_index][locusid] = l_locusreads
        d_locus_crg_share[locusid][crg_index] = 1
        crg_index += 1
        # print("Creating a new crg " + str(crg_index) + " with locus " + str(locusid))
    # print(str(datetime.datetime.now()), "End: creating CRGs with isolated loci")

    # print(str(datetime.datetime.now()), "End: creating CRGs")
    print("Number of loci: " + str(len(d_locireads.keys())))
    print("There are a total of " + str(len(d_uniq_crg_locus_reads)) + " uniq crgs")

    # print(str(datetime.datetime.now()), "Start: Writing CRGs")
    d_duplicate_entries = {}
    fh_groups_file = open(groups_file, "w")
    # enumerate again because of above processing some crgids may be missing
    for crgid, d_crgloci in enumerate(d_uniq_crg_locus_reads.values()):
        if len(d_crgloci) > 1:
            l_crgreads = set([val for sublist in d_crgloci.values() for val in sublist])
        for locusid, l_locusreads in d_crgloci.items():
            if len(d_crgloci) > 1:
                locus_share = len(l_locusreads) / float(len(l_crgreads))
            else:
                locus_share = 1.0
            for readid in l_locusreads:
                for transcriptid_pos in d_readlocus_transcripts[readid][locusid]:
                    entry = "\t".join([readid,
                                       transcriptid_pos.split('\t')[0],
                                       "locus_" + str(locusid),
                                       "group_" + str(crgid),
                                       '\t'.join(transcriptid_pos.split('\t')[1:]),
                                       d_read_genomic_pos[readid][transcriptid_pos],
                                       l_locipos[locusid],
                                       "{:.6g}".format(locus_share)]
                                      )
                    if entry in d_duplicate_entries:
                        # print("Duplicate entry found: " + entry)
                        continue
                    fh_groups_file.write(entry + "\n")
                    d_duplicate_entries[entry] = 1
    fh_groups_file.close()
    # print(str(datetime.datetime.now()), "End: CRGs written")


def em(d_read_group_fractions, em_threshold, i=1):
    d_read_group_fractions_new = defaultdict(lambda: defaultdict(float))
    print("iteration: " + str(i))
    d_group_counts = defaultdict(float)
    d_group_counts_new = defaultdict(float)
    for readid in d_read_group_fractions.keys():
        for groupid in d_read_group_fractions[readid]:
            d_group_counts[groupid] += d_read_group_fractions[readid][groupid]

    for readid in d_read_group_fractions.keys():
        total_group_count = 0
        for groupid in d_read_group_fractions[readid]:
            total_group_count += d_group_counts[groupid]
        for groupid in d_read_group_fractions[readid]:
            d_read_group_fractions_new[readid][groupid] = d_group_counts[groupid]/float(total_group_count)

    for readid in d_read_group_fractions_new.keys():
        for groupid in d_read_group_fractions_new[readid]:
            d_group_counts_new[groupid] += d_read_group_fractions_new[readid][groupid]

    equal = True
    for groupid in d_group_counts.keys():
        if abs(d_group_counts[groupid] - d_group_counts_new[groupid]) >= em_threshold:
            equal = False

    if equal:
        return d_read_group_fractions
    else:
        i += 1
        return em(d_read_group_fractions_new, em_threshold, i)


def tpm(d_group_expression, d_group_locilen):
    total_rpk = 0
    d_group_tpm = defaultdict(float)
    for groupid in sorted(d_group_expression.keys()):
        crg_expression = d_group_expression[groupid]
        crg_len = utilities.median(sorted(d_group_locilen[groupid].values()))/1000.0  # length in kbs
        # print(groupid, d_group_expression[groupid], sorted(d_group_locilen[groupid].values()), crg_len)
        rpk = crg_expression/crg_len
        d_group_tpm[groupid] = rpk
        total_rpk += rpk
    millions_of_rpk = total_rpk/1000000.0
    for groupid in sorted(d_group_expression.keys()):
        group_tpm = d_group_tpm[groupid]/millions_of_rpk
        d_group_tpm[groupid] = group_tpm
    return d_group_tpm


def quantify_crgs(loci_groups_file, em_threshold):
    d_read_group_fractions = defaultdict(lambda: defaultdict(float))
    d_group_locilen = defaultdict(lambda: defaultdict(int))
    d_group_expression = defaultdict(float)

    print(str(datetime.datetime.now()), "Start: preparing EM")
    fh_loci_groups_file = open(loci_groups_file, "r")
    for line in fh_loci_groups_file:
        f = line.rstrip("\n").split("\t")
        readid = f[0]
        locusid = f[2]
        groupid = f[3]
        pos = f[9].split(':')
        locuslength = int(pos[-2]) - int(pos[-3]) + 1
        # a single locus can belong to multiple crgs
        # one read can be part of multiple crgs
        d_read_group_fractions[readid][groupid] = 0
        d_group_locilen[groupid][locusid] = locuslength
    fh_loci_groups_file.close()

    for readid in d_read_group_fractions.keys():
        for groupid in set(d_read_group_fractions[readid].keys()):
            d_read_group_fractions[readid][groupid] = 1 / float(len(set(d_read_group_fractions[readid])))

    print(str(datetime.datetime.now()), "End: preparing EM")
    print(str(datetime.datetime.now()), "Start: recursing EM")
    sys.setrecursionlimit(10000)
    d_res = em(d_read_group_fractions, em_threshold)
    print(str(datetime.datetime.now()), "End: recursing EM")
    for readid in d_res.keys():
        # now d_read_groups should contain only multi mapped reads because because
        # uniquely mapped reads were already removed
        # set() because a read can occur 2 times at a locus in file with 2 tx ids
        # one read can be part of multiple crgs
        for groupid in d_res[readid]:
            d_group_expression[groupid] += d_res[readid][groupid]

    d_group_tpm = tpm(d_group_expression, d_group_locilen)

    return d_res, d_group_tpm
