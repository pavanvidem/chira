#!/usr/bin/env python
import argparse
import os
import shlex
import subprocess
import pysam
import chira_utilities


def run_command(command, stdout=None):
    """Run an argument list without shell interpolation and check its status."""
    print(shlex.join(command))
    subprocess.run(command, stdout=stdout, check=True)


def run_pipeline(producer_command, consumer_command, output):
    """Stream between two commands, checking both and reaping the producer."""
    print('{} | {} > {}'.format(shlex.join(producer_command),
                               shlex.join(consumer_command), shlex.quote(output)))
    with open(output, 'wb') as output_handle:
        with subprocess.Popen(producer_command, stdout=subprocess.PIPE) as producer:
            try:
                subprocess.run(consumer_command, stdin=producer.stdout,
                               stdout=output_handle, check=True)
            except BaseException:
                if producer.poll() is None:
                    producer.terminate()
                raise
            finally:
                producer.stdout.close()
            returncode = producer.wait()
            if returncode:
                raise subprocess.CalledProcessError(returncode, producer_command)


def align_with_bwa(align_type, index_type, query_fasta, refindex, outdir, seed_length, align_score,
                   match_score, mismatch_score, gap_o, gap_e, n_aligns, processes):
    """Map the supplied reads using the requested iteration parameters."""
    bam = os.path.join(outdir, index_type + "." + align_type + ".bam")
    command = ['bwa', 'mem', '-r', '1', '-c', '1000',
               '-A', str(match_score), '-B', str(mismatch_score),
               '-O', str(gap_o), '-E', str(gap_e), '-L', '0', '-Y',
               '-h', str(n_aligns), '-k', str(seed_length),
               '-T', str(align_score), '-t', str(processes), refindex, query_fasta]
    run_pipeline(command, ['samtools', 'view', '-hb', '-'], bam)
    return bam


def strand_matches(strand, stranded):
    return stranded == 'both' or strand == ('-' if stranded == 'rc' else '+')


def parse_xa(tag):
    """Parse XA hits into BED coordinates and aligned query lengths."""
    hits = []
    for entry in tag.rstrip(';').split(';'):
        if not entry:
            continue
        reference, position, cigar, mismatches = entry.split(',')
        if not position or position[0] not in '+-' or int(position[1:]) < 1:
            raise ValueError('Invalid XA position: {!r}'.format(position))
        int(mismatches)
        start = int(position[1:]) - 1
        hits.append((reference, start, start + chira_utilities.reference_length(cigar),
                     position[0], cigar, chira_utilities.aligned_query_length(cigar)))
    return hits


def write_mapped_bed(bam, bed, fasta, stranded):
    """
        Extracts the mapped and unmapped reads from the BAM alignments and writes them to BED and fasta files
        Parameters:
            bam: BAM file containing all merged alignments
            bed: output BED file path
            fasta: output FASTA file path
            stranded: Strand specificity
    """
    prev_readid = None
    read_sequence = None
    prev_unmapped = True
    if stranded not in ('fw', 'rc', 'both'):
        raise ValueError('Unknown strand selection: {!r}'.format(stranded))
    with pysam.AlignmentFile(bam, "rb") as fh_bam, open(bed, "w") as fh_mapped_bed, open(fasta, "w") as fh_unmapped_fasta:
        def write_unmapped():
            if prev_readid is not None and prev_unmapped:
                if read_sequence is None:
                    raise ValueError('No sequence available for unmapped read {!r}'.format(prev_readid))
                fh_unmapped_fasta.write('>{}\n{}\n'.format(prev_readid, read_sequence))

        for alignment in fh_bam.fetch(until_eof=True):
            readid = alignment.query_name
            if readid != prev_readid:
                write_unmapped()
                prev_unmapped = True
                read_sequence = None
            prev_readid = readid
            # Secondary records can omit SEQ. Keep the longest available sequence.
            sequence = alignment.get_forward_sequence()
            if sequence is not None and (read_sequence is None or len(sequence) > len(read_sequence)):
                read_sequence = sequence

            if alignment.is_unmapped:
                continue
            optimal_alignment_len = 0
            # write the alignment only if it mapped on desired strand
            if strand_matches('-' if alignment.is_reverse else '+', stranded):
                fh_mapped_bed.write(chira_utilities.bedentry(alignment.reference_name,
                                                             str(alignment.reference_start),
                                                             str(alignment.reference_end),
                                                             readid,
                                                             "-" if alignment.is_reverse else "+",
                                                             alignment.cigarstring) + "\n")
                prev_unmapped = False
                optimal_alignment_len = chira_utilities.aligned_query_length(alignment.cigarstring)

            # XA tag present in bwa output only
            if not alignment.has_tag('XA'):
                continue

            hits = [hit for hit in parse_xa(alignment.get_tag('XA')) if strand_matches(hit[3], stranded)]
            optimal_alignment_len = max([optimal_alignment_len] + [hit[5] for hit in hits])
            for reference, start, end, strand, cigar, query_length in hits:
                if query_length < optimal_alignment_len:
                    continue
                fh_mapped_bed.write(chira_utilities.bedentry(reference, str(start), str(end),
                                                           readid, strand, cigar) + '\n')
                prev_unmapped = False

        write_unmapped()


def align_with_clan(query_fasta, outdir, ref_fasta1, ref_index1, ref_fasta2, ref_index2,
                    chimeric_overlap, align_score, stranded, n_aligns, processes):

    ref_index = ['-f', ref_fasta1, '-d', ref_index1]
    if ref_fasta2 and ref_index2:
        ref_index += ['-F', ref_fasta2, '-D', ref_index2]
    map_to_both_strands = "FALSE"
    if stranded == "both":
        map_to_both_strands = "TRUE"
    run_command(['clan_search', '-r', query_fasta, '-m', str(n_aligns),
                 '-l', str(align_score), '-s', map_to_both_strands,
                 '-t', str(processes), '-v', str(chimeric_overlap),
                 '-o', os.path.join(outdir, 'out.clan')] + ref_index)

    ref_fasta = ['-f', ref_fasta1]
    if ref_fasta2 and ref_index2:
        ref_fasta += ['-F', ref_fasta2]

    run_command(['clan_output', '-r', query_fasta, '-i', os.path.join(outdir, 'out.clan'),
                 '-o', os.path.join(outdir, 'out.map')] + ref_fasta)


def clan_to_bed(outdir):
    """Convert CLAN's one-based locations; its output lacks strand information."""
    with open(os.path.join(outdir, "out.map")) as fh_in, open(os.path.join(outdir, "mapped.bed"), "w") as fh_out:
        next(fh_in, None)
        for line in fh_in:
            [read_id,
             solution_id,
             read_mapped_begin,
             read_mapped_end,
             read_length,
             mapped_locations] = line.rstrip("\n").rstrip("\t").split("\t")
            lead_soft_clips = int(read_mapped_begin) - 1
            trail_soft_clips = int(read_length) - int(read_mapped_end)
            if not 1 <= int(read_mapped_begin) <= int(read_mapped_end) <= int(read_length):
                raise ValueError('Invalid CLAN read coordinates for {!r}'.format(read_id))
            cigar = ""
            if lead_soft_clips != 0:
                cigar += str(lead_soft_clips) + "S"
            cigar += str(int(read_mapped_end) - int(read_mapped_begin) + 1) + "M"
            if trail_soft_clips != 0:
                cigar += str(trail_soft_clips) + "S"

            for mapped_location in mapped_locations.rstrip(';').split(";"):
                d = mapped_location.split(":")
                # if header has spaces select the id only
                ref_id = ":".join(d[0:-1]).split(' ')[0]
                [ref_start, ref_end] = d[-1].split("-")
                if not 1 <= int(ref_start) <= int(ref_end):
                    raise ValueError('Invalid CLAN reference coordinates: {!r}'.format(mapped_location))
                # NOTE: at the moment there is no way to findout which strand it is mapping to
                fh_out.write("\t".join([ref_id, str(int(ref_start)-1), ref_end,
                                        ",".join([read_id, ref_id, str(int(ref_start)-1), ref_end, "+", cigar]),
                                        "1", "+"]) + "\n")


def create_parser():

    parser = argparse.ArgumentParser(description='Chimeric Read Annotator: map reads to the reference',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-a", '--aligner', type=str, choices=["bwa", "clan"], default='bwa', required=False,
                        dest='aligner', help='Alignment program to use, bwa or clan')

    parser.add_argument('-i', '--query_fasta', dest='fasta', required=True,
                        metavar='FASTA', help='Path to query FASTA file')

    parser.add_argument('-o', '--outdir', required=True, metavar='DIR',
                        help='Output directory path for the analysis')

    parser.add_argument('-x1', '--index1', action='store', dest='idx1', required=False,
                        metavar='INDEX', help='First reference index prefix')

    parser.add_argument('-x2', '--index2', action='store', dest='idx2', required=False,
                        metavar='INDEX', help='Second reference index prefix')

    parser.add_argument('-f1', '--ref_fasta1', action='store', dest='ref_fasta1', required=False,
                        metavar='FASTA', help='First reference FASTA file')

    parser.add_argument('-f2', '--ref_fasta2', action='store', dest='ref_fasta2', required=False,
                        metavar='FASTA', help='Second reference FASTA file')

    parser.add_argument("-b", '--build', action='store_true', dest='build_index',
                        help="Build indices from reference fasta files")

    parser.add_argument('-p', '--processes', action='store', type=int, default=1,
                        dest='processes',
                        help='Number of processes to use')

    parser.add_argument("-s", '--stranded', type=str, choices=["fw", "rc", "both"], default='fw',
                        dest='stranded',
                        help='''Strand-specificity of input samples.
                             fw = map to transcript strand;
                             rc = map to reverse complement of transcript strand;
                             both = try to map on both strands''')

    parser.add_argument("-l1", '--seed_length1', action='store', type=int, default=12,
                        dest='seed_length1',
                        help='''Seed length for 1st mapping iteration.
                                bwa-mem parameter "-k"''')

    parser.add_argument("-l2", '--seed_length2', action='store', type=int, default=16,
                        dest='seed_length2',
                        help='''Seed length for 2nd mapping iteration.
                                bwa-mem parameter "-k"''')

    parser.add_argument("-s1", '--align_score1', action='store', type=int, default=18,
                        dest='align_score1',
                        help='''Minimum alignment score in 1st mapping iteration.
                                bwa-mem parameter "-T"''')

    parser.add_argument("-s2", '--align_score2', action='store', type=int, default=16,
                        dest='align_score2',
                        help='''Minimum alignment score in 2nd mapping iteration.
                                Usually smaller than --align_score1.
                                bwa-mem parameter "-T" and clan_search parameter "-l"''')

    parser.add_argument("-ma1", '--match1', action='store', type=int, default=1,
                        dest='match1',
                        help='Matching score for 1st mapping iteration.')

    parser.add_argument("-mm1", '--mismatch1', action='store', type=int, default=4,
                        dest='mismatch1',
                        help='Mismatch penalty for 1st mapping iteration.')

    parser.add_argument("-ma2", '--match2', action='store', type=int, default=1,
                        dest='match2',
                        help='Matching score for 2nd mapping iteration.')

    parser.add_argument("-mm2", '--mismatch2', action='store', type=int, default=6,
                        dest='mismatch2',
                        help='Mismatch penalty for 2nd mapping iteration.')

    parser.add_argument("-go1", '--gapopen1', action='store', type=int, default=6,
                        dest='gapopen1',
                        help='Gap opening penalty for 1st mapping iteration.')

    parser.add_argument("-ge1", '--gapext1', action='store', type=int, default=1,
                        dest='gapext1',
                        help='Gap extension penalty for 1st mapping iteration.')

    parser.add_argument("-go2", '--gapopen2', action='store', type=int, default=100,
                        dest='gapopen2',
                        help='Gap opening penalty for 2nd mapping iteration.')

    parser.add_argument("-ge2", '--gapext2', action='store', type=int, default=100,
                        dest='gapext2',
                        help='Gap extension penalty for 2nd mapping iteration.')

    parser.add_argument("-h1", '--nhits1', action='store', type=int, default=50,
                        dest='nhits1',
                        help='Number of allowed multi hits per read')

    parser.add_argument("-h2", '--nhits2', action='store', type=int, default=100,
                        dest='nhits2',
                        help='Number of allowed multi hits per read in 2nd iteration')

    parser.add_argument('-co', '--chimeric_overlap', action='store', type=int, default=2,
                        dest='chimeric_overlap',
                        help='Maximum number of overlapping bases allowed between chimeric segments')

    parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.4.3')
    return parser


def validate_args(args):
    """Validate index/reference combinations before creating output files."""
    if args.build_index:
        if args.idx1 or args.idx2:
            raise ValueError('--build and existing index options are mutually exclusive')
        if not args.ref_fasta1:
            raise ValueError('--build requires --ref_fasta1')
    elif not args.idx1:
        raise ValueError('--index1 is required unless --build is used')
    if args.aligner == 'clan':
        if not args.ref_fasta1:
            raise ValueError('CLAN requires --ref_fasta1')
        if not args.build_index and bool(args.idx2) != bool(args.ref_fasta2):
            raise ValueError('CLAN requires --index2 and --ref_fasta2 together')
    elif not args.build_index and args.ref_fasta2 and not args.idx2:
        raise ValueError('--ref_fasta2 requires --index2 unless --build is used')

    for name in ('processes', 'seed_length1', 'seed_length2', 'match1', 'match2', 'nhits1', 'nhits2'):
        if getattr(args, name) <= 0:
            raise ValueError('--{} must be positive'.format(name))
    for name in ('align_score1', 'align_score2', 'mismatch1', 'mismatch2',
                 'gapopen1', 'gapopen2', 'gapext1', 'gapext2', 'chimeric_overlap'):
        if getattr(args, name) < 0:
            raise ValueError('--{} must be non-negative'.format(name))
    for path in (args.fasta, args.ref_fasta1, args.ref_fasta2):
        if path and not os.path.isfile(path):
            raise ValueError('Input file does not exist or is not a regular file: {}'.format(path))


def print_args(args):
    print('Query fasta                          : ' + args.fasta)
    print('Output directory                     : ' + args.outdir)
    print('Aligner                              : ' + args.aligner)
    print('Build index?                         : ' + str(args.build_index))
    if args.idx1:
        print('1st priority BWA/CLAN index          : ' + args.idx1)
    if args.idx2:
        print('2nd priority BWA/CLAN  index         : ' + args.idx2)
    if args.ref_fasta1:
        print('1st priority reference fasta file    : ' + args.ref_fasta1)
    if args.ref_fasta2:
        print('2nd priority reference fasta file    : ' + args.ref_fasta2)
    print('Number of processes                  : ' + str(args.processes))
    print('Stranded                             : ' + args.stranded)
    print('Seed length                          : ' + str(args.seed_length1))
    if args.seed_length2:
        print('Seed length for 2nd iteration        : ' + str(args.seed_length2))
    print('Alignment score                      : ' + str(args.align_score1))
    if args.align_score2:
        print('Alignment score for 2nd iteration    : ' + str(args.align_score2))
    print('Chimeric overlap                     : ' + str(args.chimeric_overlap))
    print("===================================================================")


def build_indices(args):
    if not args.build_index:
        return args.idx1, args.idx2
    indices = []
    for number, reference in enumerate((args.ref_fasta1, args.ref_fasta2), start=1):
        if not reference:
            indices.append(None)
            continue
        index = os.path.join(args.outdir, 'index{}'.format(number))
        if args.aligner == 'bwa':
            run_command(['bwa', 'index', '-p', index, reference])
        else:
            run_command(['clan_index', '-f', reference, '-d', index])
        indices.append(index)
    return tuple(indices)


def map_with_bwa(args, index1, index2):
    """Run both iterations against each reference, then merge and name-sort."""
    bams = []
    for index_type, index in (('index1', index1), ('index2', index2)):
        if not index:
            continue
        for iteration, align_type in ((1, 'long'), (2, 'short')):
            message = 'Map {} read segments to {} at {}'.format(align_type, index_type, index)
            chira_utilities.print_w_time('START: ' + message)
            bams.append(align_with_bwa(
                align_type, index_type, args.fasta, index, args.outdir,
                getattr(args, 'seed_length{}'.format(iteration)),
                getattr(args, 'align_score{}'.format(iteration)),
                getattr(args, 'match{}'.format(iteration)),
                getattr(args, 'mismatch{}'.format(iteration)),
                getattr(args, 'gapopen{}'.format(iteration)),
                getattr(args, 'gapext{}'.format(iteration)),
                getattr(args, 'nhits{}'.format(iteration)), args.processes))
            chira_utilities.print_w_time('END: ' + message)

    unsorted_bam = os.path.join(args.outdir, 'unsorted.bam')
    sorted_bam = os.path.join(args.outdir, 'sorted.bam')
    chira_utilities.print_w_time('START: Merge BAM files')
    pysam.merge('-f', '-@', str(args.processes), unsorted_bam, *bams)
    chira_utilities.print_w_time('END: Merge BAM files')
    chira_utilities.print_w_time('START: Sorting BAM file')
    pysam.sort('-m', '1G', '-@', str(args.processes), '-n', unsorted_bam,
               '-T', os.path.join(args.outdir, 'sorted'), '-o', sorted_bam)
    chira_utilities.print_w_time('END: Sorting BAM file')
    for bam in bams + [unsorted_bam]:
        os.remove(bam)
    chira_utilities.print_w_time('START: Write alignments to BED')
    write_mapped_bed(sorted_bam, os.path.join(args.outdir, 'mapped.bed'),
                     os.path.join(args.outdir, 'unmapped.fasta'), args.stranded)
    chira_utilities.print_w_time('END: Write alignments to BED')


def sort_bed(outdir):
    """Keep external sorting for large files and remove input only on success."""
    bed = os.path.join(outdir, 'mapped.bed')
    with open(os.path.join(outdir, 'sorted.bed'), 'w') as output:
        run_command(['sort', '-k', '4,4', '-u', bed], stdout=output)
    os.remove(bed)


def run_mapping(args):
    os.makedirs(args.outdir, exist_ok=True)
    index1, index2 = build_indices(args)
    if args.aligner == 'clan':
        chira_utilities.print_w_time('START: Map reads using CLAN')
        align_with_clan(args.fasta, args.outdir, args.ref_fasta1, index1, args.ref_fasta2, index2,
                        args.chimeric_overlap, args.align_score2, args.stranded, args.nhits2, args.processes)
        chira_utilities.print_w_time('END: Map reads using CLAN')
        clan_to_bed(args.outdir)
    else:
        map_with_bwa(args, index1, index2)
    chira_utilities.print_w_time('START: Sorting BED file')
    sort_bed(args.outdir)
    chira_utilities.print_w_time('END: Sorting BED file')


def main(argv=None):
    parser = create_parser()
    args = parser.parse_args(argv)
    try:
        validate_args(args)
        print_args(args)
        run_mapping(args)
    except (OSError, ValueError, subprocess.CalledProcessError, pysam.SamtoolsError) as error:
        parser.exit(2, '{}: error: {}\n'.format(parser.prog, error))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
