import os
import sys
import re
from collections import defaultdict
import pkg_resources
import itertools
import pysam


def tagcounts(fasta, d_initial_tagcounts):
    fh_query = open(fasta, "r")
    for line in fh_query:
        if line.startswith('>'):
            tagid = line.rstrip('\n').replace('>', '')
            # here count whether a tag is seen each replicate
            d_initial_tagcounts[tagid] += 1
    fh_query.close()
    return d_initial_tagcounts


def readcounts(bam, d_readcounts):
    fh_bam = pysam.Samfile(bam, "rb")
    for alignment in fh_bam.fetch(until_eof=True):
        if alignment.is_unmapped:
            continue
        tagid = alignment.qname
        readcount = tagid.split('|')[1]
        d_readcounts[tagid] = int(readcount)
    fh_bam.close()
    return d_readcounts


def singletons(bam1, bam2):
    d_readcounts = readcounts(bam1)
    if bam2:
        d_readcounts2 = readcounts(bam2)
    d_readcounts.update(d_readcounts2)
    total_readcount = 0
    for c in d_readcounts.values():
        total_readcount += c
    return total_readcount


def hybrid_readcounts(pairs):
    fh_pairs = open(pairs, "r")
    d_readcounts = defaultdict()
    for line in fh_pairs:
        tagid = line.split('\t')[1].replace('|l','').replace('|r','')
        readcount = tagid.split('|')[1]
        d_readcounts[tagid] = readcount
    total_readcount = 0
    fh_pairs.close()
    for c in d_readcounts.values():
        total_readcount += c
    return total_readcount


def unmappedcounts(fasta, d_unmapped_tagcounts):
    fh_unmapped = open(fasta, "r")
    for line in fh_unmapped:
        if line.startswith('>'):
            tagid = line.rstrip('\n').replace('>','')
            d_unmapped_tagcounts[tagid] += 1
    fh_unmapped.close()
    return d_unmapped_tagcounts


def hybrid_abundances(pairs):
    fh_pairs = open(pairs, "r")
    d_tag_loci = defaultdict(str)
    for line in fh_pairs:
        f = line.rstrip('\n').split('\t')
        tagid = f[1].replace('|l','').replace('|r','')
        locus1 = f[20]
        locus2 = f[21]
        score = f[28]
        if score > d_tag_scores[tagid]:
            d_tag_scores[tagid] = score
            d_tag_loci[tagid] = locus1+locus2
    fh_pairs.close()

    d_hybrid_abundances = defaultdict()
    d_hybrid_read_abundances = defaultdict(float)
    for tagid in d_tag_loci.keys():
        loci = d_tag_loci[tagid]
        d_hybrid_abundances[loci] += 1
        readcount = tagid.split('|')[1]
        d_hybrid_read_abundances[loci] += readcount


def file_prefixes(outputdir, firstpart, secondpart):
    if firstpart == "miRNA":
        first_prefix = os.path.join(outputdir, 'vs_miRNA')
        pairs_prefix = os.path.join(outputdir, 'miRNA')
        if secondpart == "miRNA":
            pairs_prefix += '_vs_miRNA'
            # TODO: never used? check this
            second_prefix = os.path.join(outputdir, 'clipped_vs_mirna')
        elif secondpart == "transcriptome":
            second_prefix = os.path.join(outputdir, 'clipped_vs_transcriptome')
            pairs_prefix += '_vs_transcriptome'
        else:
            sys.exit("Unknown hybrid first part")
    elif firstpart == "transcriptome":
        first_prefix = os.path.join(outputdir, 'vs_transcriptome')
        pairs_prefix = os.path.join(outputdir, 'transcriptome')
        if secondpart == "transcriptome":
            pairs_prefix += '_vs_transcriptome'
            # TODO: never used? check this
            second_prefix = os.path.join(outputdir, 'clipped_vs_transcriptome')
        else:
            sys.exit("Unknown hybrid second part")
    else:
        sys.exit("Unknown hybrid part")

    return first_prefix, second_prefix, pairs_prefix
