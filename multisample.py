from collections import defaultdict
import os

import utilities


def parse_pairs(locipairs, d_interaction_counts, d_interaction_scores):
    fh_locipairs = open(locipairs, "r")
    for line in fh_locipairs:
        f = line.rstrip('\n').split('\t')
        tagid = f[1]
        transcriptid1 = f[2]
        transcriptid2 = f[3]
        score = float(f[-1])
        if score <= 0:
            continue
        d_interaction_counts[transcriptid1 + "\t" + transcriptid2][self.sample].append(tagid)
        d_interaction_scores[transcriptid1 + "\t" + transcriptid2 + "\t" + tagid][self.sample] = score
    fh_locipairs.close()
    return


def create_loci_beds(locipairs, leftpart_bed, rightpart_bed):
    fh_locipairs = open(locipairs, "r")
    fh_leftpart_bed = open(leftpart_bed, "w")
    fh_rightpart_bed = open(rightpart_bed, "w")
    for line in fh_locipairs:
        f = line.rstrip('\n').split('\t')
        score = float(f[-1])
        if score <= 0:
            continue
        tagid = f[1]
        locusid1 = f[18]
        locusid2 = f[19]
        locus1 = f[20]
        locus2 = f[21]
        print(locus1)
        [chrom, start, end, strand] = locus1.split(':')
        fh_leftpart_bed.write("\t".join([chrom, start, end, locusid1, 0, strand]))
        [chrom, start, end, strand] = locus2.split(':')
        fh_rightpart_bed.write("\t".join([chrom, start, end, locusid2, 0, strand]))
    fh_readpairs.close()

    return


def reproducible_interactions(outputdir, samples, firstpart, secondpart):
    d_interaction_counts = defaultdict(lambda: defaultdict(list))
    d_interaction_scores = defaultdict(lambda: defaultdict(list))
    for sample in samples:
        sampledir = os.path.join(outputdir, sample)
        first_prefix, second_prefix, pairs_prefix = utilities_paired.file_prefixes(sampledir, firstpart, secondpart)
        locipairs = pairs_prefix + ".loci.pairs"
        parse_pairs(locipairs, sample, d_interaction_counts, d_interaction_scores)
        fh_sample_summary = open(pairs_prefix + ".summary", "w")
        for interaction in d_interaction_counts.keys():
            if sample not in d_interaction_counts[interaction]:
                continue
            abundance = len(set(d_interaction_counts[interaction][sample]))
            if abundance <= 0:
                continue
            [first_transcript, second_transcript] = interaction.split('\t')
            # TODO: averaging the probabilities?
            total_score = 0
            for tagid in set(d_interaction_counts[interaction][sample]):
                total_score += d_interaction_scores[interaction + "\t" + tagid][sample]
            # avg_score = total_score/float(len(set(d_interaction_counts[interaction][sample])))
            fh_sample_summary.write("\t".join([interaction,
                                               d_txid_geneid[first_transcript],
                                               d_geneid_genetype[d_txid_geneid[first_transcript]],
                                               d_txid_geneid[second_transcript],
                                               d_geneid_genetype[d_txid_geneid[second_transcript]],
                                               str(abundance),
                                               str(total_score) + "\n"]))
        fh_sample_summary.close()

    print(len(d_interaction_counts))
    fh_summary = open(os.path.join(outputdir, 'miRNA_vs_transcriptome.summary'), "w")
    for interaction in d_interaction_counts.keys():
        total_abundance = 0
        total_score = 0.0
        for sample in d_interaction_counts[interaction].keys():
            abundance = len(set(d_interaction_counts[interaction][sample]))
            for tagid in set(d_interaction_counts[interaction][sample]):
                total_score += d_interaction_scores[interaction + "\t" + tagid][sample]
            total_abundance += abundance
        n_samples = len(d_interaction_counts[interaction])
        avg_abundance = total_abundance / float(n_samples)
        avg_score = total_score / float(total_abundance)
        [first_transcript, second_transcript] = interaction.split('\t')
        fh_summary.write("\t".join([interaction,
                                    # d_txid_geneid[second_transcript],
                                    # d_geneid_genetype[d_txid_geneid[second_transcript]],
                                    str(n_samples),
                                    ','.join(d_interaction_counts[interaction].keys()),
                                    str(total_abundance),
                                    str(avg_abundance),
                                    str(total_score),
                                    str(avg_score) + "\n"]))
    fh_summary.close()

    return


def reproducible_interactions_at_loci_level(outputdir, samples, firstpart, secondpart):
    d_interaction_counts = defaultdict(lambda: defaultdict(list))
    d_interaction_samples = defaultdict(lambda: defaultdict(int))
    for sample in samples:
        sampledir = os.path.join(outputdir, sample)
        first_prefix, second_prefix, pairs_prefix = utilities_paired.file_prefixes(sampledir, firstpart, secondpart)
        locipairs = pairs_prefix + ".loci.pairs"
        leftpart_bed = pairs_prefix + ".leftpart.bed"
        rightpart_bed = pairs_prefix + ".rightpart.bed"
        create_loci_beds(locipairs, leftpart_bed, rightpart_bed)

    print(list(itertools.product(sample, sample)))
    exit()

    for pair in list(itertools.product(sample, sample)):
        sample1 = pair[0]
        sample2 = pair[1]
        sampledir = os.path.join(outputdir, sample)
        first_prefix1, second_prefix1, pairs_prefix1 = utilities_paired.file_prefixes(sampledir, firstpart, secondpart)
        first_prefix2, second_prefix2, pairs_prefix2 = utilities_paired.file_prefixes(sampledir, firstpart, secondpart)
        leftpart_bed1 = pairs_prefix1 + ".leftpart.bed"
        rightpart_bed1 = pairs_prefix1 + ".rightpart.bed"
        leftpart_bed2 = pairs_prefix2 + ".leftpart.bed"
        rightpart_bed2 = pairs_prefix2 + ".rightpart.bed"

        intersect_command = (INTERSECTBED + " -a " + leftpart_bed1 + " -b " + leftpart_bed2 + " -f 0.5 -r -s -wb > "
                             + overlapout)
        intersect_command = (INTERSECTBED + " -a " + rightpart_bed1 + " -b " + rightpart_bed2 + " -f 0.5 -r -s -wb > "
                             + overlapout)

        print(intersect_command)
        os.system(intersect_command)

        fh_overlapout = open(overlapout, "r")
        d_idsseen = {}
        d_idsdouble = {}
        d_withjunctions = defaultdict(list)

        # ENSMUST00000023614      5301    5319    2175|tag_1004650|1,ENSMUST00000023614,ENSMUSG00000022897,5301,5319 \
        #       0       +       ENSMUST00000023614      1817    5754    ENSMUST00000023614_e012 0       +
        for line in fh_overlapout:
            # txid, s, e, readid, exonstart, exonend, exonid, pol = line.rstrip('\n').split('\t')[0,1,2,3,7,8,9,11]
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
                n_plus += 1
            # Minus strand case.
            elif pol == "-":
                genomicstart = exid2end[exonid] - (txend - exonstart)
                genomicend = exid2end[exonid] - (txstart - exonstart)
                n_minus += 1
            else:
                sys.exit("ERROR: Check polarity " + pol + "in line: " + line)

            # Store ID.
            if readid in d_idsseen:
                d_idsdouble[readid] = d_idsdouble.get(readid, 0) + 1
            else:
                d_idsseen[readid] = d_idsseen.get(readid, 0) + 1

            n_out += 1
            d_withjunctions[readid].append(id2chr[transcriptid] + "\t" +
                                           str(genomicstart) + "\t" +
                                           str(genomicend) + "\t" +
                                           readid +
                                           "\t0\t" +
                                           pol + '\n')
            if str(id2pol[transcriptid]) != pol:
                sys.exit("Different strands!!")
        fh_overlapout.close()

        parse_pairs(locipairs, sample, d_interaction_counts, d_interaction_samples)
        fh_sample_summary = open(pairs_prefix + ".summary", "w")
        for interaction in d_interaction_counts.keys():
            if sample not in d_interaction_counts[interaction]:
                continue
            abundance = d_interaction_counts[interaction][sample]
            if abundance > 0:
                [first_transcript, second_transcript] = interaction.split('\t')
                score = d_interaction_samples[interaction][sample]
                fh_sample_summary.write("\t".join([interaction,
                                                   d_txid_geneid[second_transcript],
                                                   d_geneid_genetype[d_txid_geneid[second_transcript]],
                                                   str(abundance),
                                                   str(score) + "\n"]))
        fh_sample_summary.close()

    print(len(d_interaction_counts))
    fh_summary = open(os.path.join(outputdir, 'miRNA_vs_transcriptome.summary'), "w")
    for interaction in d_interaction_counts.keys():
        total_abundance = 0
        weighted_score = 0.0
        for sample in d_interaction_counts[interaction].keys():
            abundance = int(d_interaction_counts[interaction][sample])
            score = float(d_interaction_samples[interaction][sample])
            total_abundance += abundance
            weighted_score += abundance * score
        n_samples = len(d_interaction_samples[interaction])
        avg_abundance = total_abundance / float(n_samples)
        weighted_avg_score = weighted_score / float(total_abundance)
        [first_transcript, second_transcript] = interaction.split('\t')
        fh_summary.write("\t".join([interaction,
                                    d_txid_geneid[second_transcript],
                                    d_geneid_genetype[d_txid_geneid[second_transcript]],
                                    str(n_samples),
                                    ','.join(d_interaction_samples[interaction].keys()),
                                    str(total_abundance),
                                    str(avg_abundance),
                                    str(weighted_avg_score) + "\n"]))
    fh_summary.close()

    return
