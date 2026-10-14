#!/usr/bin/env python
import argparse
import math
import re
from Bio import SeqIO
from datetime import datetime


_CIGAR_OPERATION = re.compile(r'([1-9][0-9]*)([MIDNSHP=X])')


def score_float(x):
    """Parse a finite command-line score in the closed interval [0, 1]."""
    x = float(x)
    if not math.isfinite(x) or not 0.0 <= x <= 1.0:
        raise argparse.ArgumentTypeError("%r not in range [0.0, 1.0]" % (x,))
    return x


def overlap(f, s):
    """Return overlap length for half-open intervals, such as BED coordinates."""
    return max(0, min(f[1], s[1]) - max(f[0], s[0]))


def inclusive_overlap(f, s):
    """Return overlap length for inclusive intervals, such as read positions."""
    return max(0, min(f[1], s[1]) - max(f[0], s[0]) + 1)


def median(x):
    """Return the median of a non-empty, already sorted sequence."""
    n = len(x)
    if not n:
        raise ValueError('Cannot calculate the median of an empty sequence')
    mid = n // 2
    if not n % 2:
        return (x[mid-1] + x[mid]) / 2.0
    return x[mid]


def parse_cigar(cigar):
    """Parse positive-length SAM CIGAR operations; reject gaps and bad clipping.

    Hard clips may appear only at the ends, and soft clips only at the ends
    inside any hard clips. The unmapped sentinel '*' is not a CIGAR here.
    """
    operations = []
    position = 0
    for match in _CIGAR_OPERATION.finditer(cigar):
        if match.start() != position:
            raise ValueError('Invalid CIGAR: {!r}'.format(cigar))
        operations.append((int(match.group(1)), match.group(2)))
        position = match.end()
    if not operations or position != len(cigar):
        raise ValueError('Invalid CIGAR: {!r}'.format(cigar))

    inner = operations[:]
    if inner and inner[0][1] == 'H':
        inner.pop(0)
    if inner and inner[-1][1] == 'H':
        inner.pop()
    if inner and inner[0][1] == 'S':
        inner.pop(0)
    if inner and inner[-1][1] == 'S':
        inner.pop()
    if any(operation in 'HS' for length, operation in inner):
        raise ValueError('Invalid clipping in CIGAR: {!r}'.format(cigar))
    return operations


def query_length(cigar, is_reverse):
    """Return original read span, including hard clips, on either strand.

    This intentionally includes bases absent from SAM SEQ after hard clipping.
    is_reverse is retained for compatibility; total length is strand-independent.
    """
    return sum(length for length, operation in parse_cigar(cigar) if operation in 'HIMS=X')


def aligned_query_length(cigar):
    """Return aligned query span (M/I/=/X), excluding both types of clipping."""
    return sum(length for length, operation in parse_cigar(cigar) if operation in 'MI=X')


def reference_length(cigar):
    """Return reference span (M/D/N/=/X), including deletions and skipped regions."""
    return sum(length for length, operation in parse_cigar(cigar) if operation in 'MDN=X')


def match_positions(cigar, is_reverse):
    """Return one-based, inclusive aligned positions on the original read.

    Reverse-strand CIGARs are traversed in original-read orientation. Insertions
    count towards the aligned span; deletions, skips and padding do not. At
    least one reference-aligned query operation (M/=/X) is required.
    """
    operations = parse_cigar(cigar)
    match_start = match_end = 0
    query_position = 0
    has_match = False
    if is_reverse:
        operations = reversed(operations)
    for length, operation in operations:
        if operation in 'MI=X':
            if not match_start:
                match_start = query_position + 1
            match_end = query_position + length
            has_match = has_match or operation in 'M=X'
        if operation in 'HIMS=X':
            query_position += length
    if not has_match:
        raise ValueError('CIGAR has no aligned query bases: {!r}'.format(cigar))
    return match_start, match_end


def is_chimeric(cigar1, cigar2, is_reverse1, is_reverse2, max_allowed_overlap):
    """Check whether inclusive read spans overlap by at most the given bases."""
    match_start1, match_end1 = match_positions(cigar1, is_reverse1)
    match_start2, match_end2 = match_positions(cigar2, is_reverse2)
    return inclusive_overlap((match_start1, match_end1), (match_start2, match_end2)) <= max_allowed_overlap


def alignment_length(cigar):
    """Return legacy alignment length (M/I/D/=/X), excluding skips and clips.

    This is neither query span nor reference span when indels are present.
    """
    return sum(length for length, operation in parse_cigar(cigar) if operation in 'DIM=X')


def alignment_end(start, cigar, is_reverse):
    """Return the inclusive reference end from a one-based SAM start.

    This numeric end is also the zero-based BED exclusive end. is_reverse is
    retained for compatibility; reference span is strand-independent.
    """
    return int(start) + reference_length(cigar) - 1


def bedentry(referenceid, reference_start, reference_end, readid, strand, cigarstring):
    """Format a BED row with zero-based, end-exclusive reference coordinates."""
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


def extract_reflengths(ref_fasta, d_reflen):
    """Populate the supplied mapping with FASTA record lengths."""
    with open(ref_fasta) as fasta_handle:
        for record in SeqIO.parse(fasta_handle, 'fasta'):
            d_reflen[record.id] = len(record)


def print_w_time(message):
    print("[" + datetime.now().strftime('%d-%m-%Y %H:%M:%S') + "] " + message)
