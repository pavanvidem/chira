#!/usr/bin/env python
import re
import string


def overlap(f, s):
    return max(0, min(f[1], s[1]) - max(f[0], s[0]))


def median(x):
    n = len(x)
    mid = int(n/2)
    if not n % 2:
        return (x[mid-1] + x[mid]) / 2.0
    return x[mid]


def query_length(cigar, is_reverse):
    cigar_tup = re.findall(r'(\d+)([MISH=X])', cigar)
    read_length = 0
    if is_reverse:
        cigar_tup = reversed(cigar_tup)
    for c in cigar_tup:
        read_length += int(c[0])
    return read_length


def match_positions(cigar, is_reverse):
    cigar_tup = re.findall(r'(\d+)([MISH=X])', cigar)
    match_start = match_end = 0
    if is_reverse:
        cigar_tup = reversed(cigar_tup)
    for c in cigar_tup:
        if c[1] == "S" or c[1] == "H":
            if not match_start:
                match_start = match_end = int(c[0]) + 1
        elif c[1] == "M":
            if match_start:
                match_end = match_end + int(c[0]) - 1
            else:
                match_start = 1
                match_end = int(c[0])
    return match_start, match_end


def switch_alignments(cigar1, cigar2, is_reverse1, is_reverse2, max_allowed_overlap):
    match_start1, match_end1 = match_positions(cigar1, is_reverse1)
    match_start2, match_end2 = match_positions(cigar2, is_reverse2)
    if overlap([match_start1, match_end1], [match_start2, match_end2]) > max_allowed_overlap:
        return None
    else:
        if match_start1 < match_start2:
            return "n"
        else:
            return "y"


# TODO check if substracting the NM can give the no of matched bases
def alignment_score(cigar, nm):
    cigar_tup = re.findall(r'(\d+)(M)', cigar)
    n_matches = 0
    for c in cigar_tup:
        n_matches += int(c[0])
    # TODO be careful with the hardcoded mismatch score of 4
    return n_matches - int(nm) * 4


def alignment_end(start, cigar, is_reverse):
    cigar_tup = re.findall(r'(\d+)([MDN=X])', cigar)
    end = int(start)
    if is_reverse:
        cigar_tup = reversed(cigar_tup)
    for c in cigar_tup:
        end += int(c[0])
    return end - 1


def bedentry(referenceid, reference_start, reference_end, readid, strand, cigarstring):
    line = '\t'.join([referenceid,
                      reference_start,
                      reference_end,
                      ','.join([readid, referenceid, reference_start, reference_end, strand, cigarstring]),
                      "1",
                      strand])
    return line


def reverse_complement(seq):
    tab = str.maketrans("ACTGactg", "TGACtgac")
    return seq.translate(tab)[::-1]
