import os
import sys
import re
from collections import defaultdict
import pysam
import pkg_resources
from BCBio import GFF
from BCBio.GFF import GFFExaminer


def transcript_to_genomic_pos(transcripts_bed, f_geneexonbed, f_txexonbed):
        id2pol = {}
        id2chr = {}
        exid2start = {}
        exid2end = {}
        print("Read in genomic exon features .. ")
        fh_genomic_exonsne_exon_bed = open(f_geneexonbed, "r")
        for line in fh_genomic_exonsne_exon_bed:
            # chrName, s, e, exid, pol = line.rstrip('\n').split('\t')[0,1,2,3,5]
            f = line.rstrip('\n').split('\t')
            chrname = f[0]
            s = int(f[1])
            e = int(f[2])
            exid = f[3]
            pol = f[5]
            match = re.match("(.+?)_e", exid)
            transcriptid = match.group(1)
            id2pol[transcriptid] = pol
            id2chr[transcriptid] = chrname
            exid2start[exid] = s
            exid2end[exid] = e
        fh_genomic_exonsne_exon_bed.close()
        print("done")
        print("Transcripts with exons:   " + str(len(id2pol)))

        overlapout = transcripts_bed.replace(".bed", ".overlap.txt")
        wojunctions = transcripts_bed.replace(".bed", ".genomic.bed")

        print("Calculating bed overlap .. ")
        intersect_command = ("intersectBed -a " + transcripts_bed + " -b " + f_txexonbed + " -wb > " + overlapout)
        print(intersect_command)
        os.system(intersect_command)

        fh_overlapout = open(overlapout, "r")

        d_idsseen = {}
        d_idsdouble = {}
        d_ids_plus = {}
        d_ids_minus = {}
        d_withjunctions = defaultdict(list)

        # ENSMUST00000023614      5301    5319    2175|tag_1004650|1,ENSMUST00000023614,ENSMUSG00000022897,5301,5319 \
        #       0       +       ENSMUST00000023614      1817    5754    ENSMUST00000023614_e012 0       +
        for line in fh_overlapout:
            # transcriptid, s, e, readid, exonstart, exonend, exonid, pol
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
                d_ids_plus[readid] = d_ids_plus.get(readid, 0) + 1
            # Minus strand case.
            elif pol == "-":
                genomicstart = exid2end[exonid] - (txend - exonstart)
                genomicend = exid2end[exonid] - (txstart - exonstart)
                d_ids_minus[readid] = d_ids_minus.get(readid, 0) + 1
            else:
                sys.exit("ERROR: Check polarity " + pol + "in line: " + line)

            # Store ID.
            if readid in d_idsseen:
                d_idsdouble[readid] = d_idsdouble.get(readid, 0) + 1
            else:
                d_idsseen[readid] = d_idsseen.get(readid, 0) + 1

            d_withjunctions[readid].append(id2chr[transcriptid] + "\t" +
                                           str(genomicstart) + "\t" +
                                           str(genomicend) + "\t" +
                                           readid +
                                           "\t0\t" +
                                           pol + '\n')
            if str(id2pol[transcriptid]) != pol:
                sys.exit("Different strands!!")
        fh_overlapout.close()

        print("Total alignments:        " + str(len(d_ids_plus)+len(d_ids_minus)))
        print("Hits on plus strand:     " + str(len(d_ids_plus)))
        print("Hits on minus strand:    " + str(len(d_ids_minus)))
        print("Hits over exon borders:  " + str(sum(1 for i in {**d_ids_plus, **d_ids_minus}.values() if i >= 2)))

        # Write alignments spanning exon borders to a file
        fh_wojunctions = open(wojunctions, "w")
        n_final = 0
        print("Filtering out exon border reads .. ")
        for readid, lines in d_withjunctions.items():
            # Alignments over exon borders have more than 2 lines in hash
            # This way we can identify them and filter them out.
            # TODO: Considering them as 2 seperate alignments at genomic level doesn't harm. Just need to make sure that
            # they should not be treated as multi-mapped read. Give them a special ID while quantification
            if len(lines) > 1:
                # print readid
                continue
            n_final += 1
            for line in lines:
                fh_wojunctions.write(line)
        fh_wojunctions.close()

        print("done")
        print("Alignments after filtering junction reads: " + str(n_final))
        os.system("rm " + overlapout)

        return


def overlap(f, s):
    return max(0, min(f[1], s[1]) - max(f[0], s[0]))


def strandardize(strand):
    if strand == '-1':
        strand = '-'
    elif strand == '1':
        strand = '+'
    return strand


def median(x):
    n = len(x)
    mid = int(n/2)
    if not n%2:
        return (x[mid-1] + x[mid]) / 2.0
    return x[mid]


def chunks(ids, n):
    return [ids[i:i+n] for i in range(0, len(ids), n)]


def fasta_to_dict(fasta):
    seqid = None
    seq = ""
    fasta_dict = {}
    with open(fasta) as fh_fasta:
        for line in fh_fasta:
            if line.startswith(">"):
                fasta_dict[seqid] = seq
                seqid = line
                seq = ""
            else:
                seq += trim(line)
    return fasta_dict


def parse_counts_file(loci_groups_counts_file, tpm_cutoff, score_cutoff):
    d_group_tpm = defaultdict(float)
    d_read_lines = defaultdict(list)

    fh_loci_groups_counts_file = open(loci_groups_counts_file, "r")
    for line in fh_loci_groups_counts_file:
        f = line.rstrip('\n').split('\t')
        groupid = f[3]
        group_tpm = f[12]
        d_group_tpm[groupid] = float(group_tpm)
    fh_loci_groups_counts_file.close()

    uniq_tpms = sorted(list(set(d_group_tpm.values())))
    tpm_threshold = uniq_tpms[int(tpm_cutoff * len(uniq_tpms))]

    fh_loci_groups_counts_file = open(loci_groups_counts_file, "r")
    for line in fh_loci_groups_counts_file:
        f = line.rstrip('\n').split('\t')
        readid = '|'.join(f[0].split("|")[:-1])
        group_tpm = f[12]
        prob = f[10]
        locusshare = f[11]
        if float(group_tpm) < tpm_threshold:
            continue
        locus_score = float(prob) * float(locusshare)
        if locus_score < score_cutoff:
            continue
        d_read_lines[readid].append(line)
    fh_loci_groups_counts_file.close()

    return d_read_lines


def pairspos_to_bed(pairs_prefix, bed1, bed2):
    fh_readpairs = open(pairs_prefix + ".loci.pairs", "r")
    fh_bed1 = open(bed1, "w")
    fh_bed2 = open(bed2, "w")
    for line in fh_readpairs:
        f = line.rstrip('\n').split('\t')
        chimeraid = f[0]
        pos1 = f[16]
        f1 = pos1.split(':')
        # TODO : add 20nt flanking
        f1[1] = str(int(f1[1]) - 0)
        f1[2] = str(int(f1[2]) + 0)
        f1.insert(3, chimeraid)
        f1.insert(4, "0")
        fh_bed1.write('\t'.join(f1)+"\n")

        pos2 = f[17]
        f2 = pos2.split(':')
        # TODO : add 20nt flanking
        f2[1] = str(int(f2[1]) - 0)
        f2[2] = str(int(f2[2]) + 0)
        f2.insert(3, chimeraid)
        f2.insert(4, "0")
        fh_bed2.write('\t'.join(f2)+"\n")

    fh_readpairs.close()
    fh_bed1.close()
    fh_bed2.close()
    return


def run_intarna(d_seq_records1, d_seq_records2, ids, matrix_file):
    os.system("[ -e " + matrix_file + "] && rm " + matrix_file + ";  touch " + matrix_file)
    for seqId in ids:
        record1 = d_seq_records1[seqId]
        record2 = d_seq_records2[seqId]
        os.system("python " +
                  os.path.dirname(os.path.realpath(__file__)) + "/intarna_predict_hybrid_return_base_pairs.py " +
                  " -m " + str(record1.seq) +
                  " -t " + str(record2.seq) +
                  " >> " + matrix_file)
    return


def run_hybridmin(d_seq_records1, d_seq_records2, ids, hybrid_file, pairs_prefix, p_hybmin):
    fh_hybrids = open(hybrid_file, "w")
    for seqId in ids:
        record1 = d_seq_records1[seqId]
        record2 = d_seq_records2[seqId]
        hybrid_prefix = pairs_prefix + "." + seqId
        with open(hybrid_prefix + ".s1", "w") as f1:
            f1.write(str(record1.seq))
        with open(hybrid_prefix + ".s2", "w") as f2:
            f2.write(str(record2.seq))
        os.system(p_hybmin + " -o " + hybrid_prefix + " " + hybrid_prefix + ".s1 " + hybrid_prefix + ".s2 > /dev/null") # + hybrid_prefix + ".log"

        mfe = None
        with open(hybrid_prefix + ".dG") as fh_dg:
            next(fh_dg)
            for line in fh_dg:
                mfe = line.split("\t")[1]
        i = []
        j = []
        with open(hybrid_prefix + ".37.plot") as fh_plot:
            next(fh_plot)
            for line in fh_plot:
                f = line.split("\t")
                i.append(int(f[0]))
                j.append(int(f[1]))
        dotbracket = ""
        for index in range(1, len(str(record1.seq))+1):
            if index in i:
                dotbracket += "("
            else:
                dotbracket += "."
        dotbracket += "&"
        for index in range(1, len(str(record2.seq))+1):
            if index in j:
                dotbracket += ")"
            else:
                dotbracket += "."
        if not i:
            i.append(-1)
            j.append(-1)
        fh_hybrids.write("\t".join([seqId,
                                   str(record1.seq) + "&" + str(record2.seq),
                                   dotbracket,
                                   str(i[0]) + "&" + str(j[-1]),
                                   mfe,
                                   "\n"]))
        os.system("rm " + hybrid_prefix + ".*")
    fh_hybrids.close()
    return


def parse_annotations(d_gene_annotations,
                      d_transcript_annotations,
                      f_goterms,
                      f_mirfammap,
                      f_gff,
                      f_bed,
                      outputdir):
    n_exon = 1
    exon_rel_start = 0
    exon_rel_end = 0
    exon_len = 0
    prev_transcript_id = None

    limit_info = dict(gff_type=["exon", "UTR", "CDS", "miRNA_primary_transcript", "miRNA", "tRNA"])

    d_attributes = defaultdict(list)
    d_attributes['tid'] = ['transcript_id', 'Name']
    d_attributes['gid'] = ['gene_id', 'ID']
    d_attributes['name'] = ['gene_name', 'Name']
    d_attributes['type'] = ['gene_biotype', 'Type']

    fh_genomic_exons = open(os.path.join(outputdir, "genomic_exons.bed"), "w")
    fh_transcriptomic_exons = open(os.path.join(outputdir, "transcriptomic_exons.bed"), "w")

    with open(f_gff) as gff_handle:
        # Chromosome seq level
        for rec in GFF.parse(gff_handle, limit_info=limit_info, target_lines=1):
            # for each selected sub_feature
            for sub_feature in rec.features:
                # for some reason each qualifier is a list! take the first element
                for i in d_attributes['gid']:
                    if i in sub_feature.qualifiers:
                        gene_id = sub_feature.qualifiers[i][0]
                        break
                for i in d_attributes['tid']:
                    if i in sub_feature.qualifiers:
                        transcript_id = sub_feature.qualifiers[i][0]
                        break
                if sub_feature.type == 'UTR':
                    if transcript_id not in d_transcript_annotations['UTR']:
                        d_transcript_annotations['UTR'][transcript_id] = []
                    d_transcript_annotations['UTR'][transcript_id].append([rec.id,
                                                                           str(sub_feature.location.start),
                                                                           str(sub_feature.location.end),
                                                                           str(sub_feature.location.strand)])
                elif sub_feature.type == 'CDS':
                    if transcript_id not in d_transcript_annotations['CDS']:
                        d_transcript_annotations['CDS'][transcript_id] = []
                    d_transcript_annotations['CDS'][transcript_id].append([rec.id,
                                                                           str(sub_feature.location.start),
                                                                           str(sub_feature.location.end),
                                                                           str(sub_feature.location.strand)])
                # mature miRNA annotaion
                elif sub_feature.type == 'miRNA':
                    if 'Derives_from' in sub_feature.qualifiers:
                        gene_id = sub_feature.qualifiers['Derives_from'][0]
                    for i in d_attributes['name']:
                        if i in sub_feature.qualifiers:
                            mature_mirna_name = sub_feature.qualifiers[i][0]
                            if gene_id not in d_transcript_annotations['mature']:
                                # TODO: a quick fix used gene_id in transcript annotations dict
                                d_transcript_annotations['mature'][gene_id] = []
                            d_transcript_annotations['mature'][gene_id].append([rec.id,
                                                                                str(sub_feature.location.start),
                                                                                str(sub_feature.location.end),
                                                                                str(sub_feature.location.strand),
                                                                                mature_mirna_name])
                # remaining are exon and miRNA_primary_transcript lines
                else:
                    # reset the variables if it is a new transcript
                    if prev_transcript_id != transcript_id:
                        d_transcript_annotations['len'][transcript_id] = 0
                        n_exon = 1
                        exon_rel_start = 0
                        exon_rel_end = 0
                        exon_len = 0
                    d_transcript_annotations['gid'][transcript_id] = gene_id
                    # biotype
                    bio_type = 'NA'
                    for i in d_attributes['type']:
                        if i in sub_feature.qualifiers:
                            bio_type = sub_feature.qualifiers[i][0]
                    if sub_feature.type == "miRNA_primary_transcript":
                        bio_type = "miRNA"
                    d_gene_annotations['type'][gene_id] = bio_type

                    # name
                    try:
                        for i in d_attributes['name']:
                            if i in sub_feature.qualifiers:
                                d_gene_annotations['name'][gene_id] = sub_feature.qualifiers[i][0]
                    except KeyError:
                        d_gene_annotations['name'][gene_id] = "NA"

                    gene_exon_bed_entry = '\t'.join([rec.id,
                                                     str(sub_feature.location.start),
                                                     str(sub_feature.location.end),
                                                     transcript_id + "_e" + str(n_exon).zfill(3),
                                                     '0',
                                                     '-' if str(sub_feature.location.strand) == '-1' else '+'
                                                     ])
                    exon_len = sub_feature.location.end - sub_feature.location.start
                    exon_rel_end = exon_rel_start + exon_len
                    tx_exon_bed_entry = '\t'.join([transcript_id,
                                                   str(exon_rel_start),
                                                   str(exon_rel_end),
                                                   transcript_id + "_e" + str(n_exon).zfill(3),
                                                   '0',
                                                   '-' if str(sub_feature.location.strand) == '-1' else '+'
                                                   ])
                    fh_genomic_exons.write(gene_exon_bed_entry + "\n")
                    fh_transcriptomic_exons.write(tx_exon_bed_entry + "\n")

                    # TODO: check +1 to length ?
                    d_transcript_annotations['len'][transcript_id] += exon_len + 1
                    exon_rel_start = exon_rel_end
                    prev_transcript_id = transcript_id
                    n_exon += 1
                # print(gene_id, d_gene_annotations['type'][gene_id])
    gff_handle.close()
    fh_genomic_exons.close()
    fh_transcriptomic_exons.close()

    if f_bed:
        # TODO: this is only an example for parsing gtRNAdb bed file
        fh_bed = open(f_bed, "r")
        fh_genomic_exons = open(os.path.join(outputdir, "genomic_exons.bed"), "a")
        fh_transcriptomic_exons = open(os.path.join(outputdir, "transcriptomic_exons.bed"), "a")
        for line in fh_bed:
            f = line.split('\t')
            family = re.sub(r'-\d.+', '', f[3])
            d_transcript_annotations['gid'][f[3]] = f[3]
            # TODO: check no +1 to length?
            d_transcript_annotations['len'][f[3]] = int(f[2])-int(f[1])
            d_gene_annotations['type'][f[3]] = "tRNA"
            d_gene_annotations['name'][f[3]] = f[3]
            d_gene_annotations['family'][f[3]] = family
            d_transcript_annotations['family'][f[3]] = family
            fh_genomic_exons.write('\t'.join([f[0], f[1], f[2], f[3]+'_e_001', '0', f[5], '\n']))
            fh_transcriptomic_exons.write(
                '\t'.join([f[3], '0', str(int(f[2]) - int(f[1])), f[3] + '_e_001', '0', f[5], '\n']))
        fh_genomic_exons.close()
        fh_transcriptomic_exons.close()

    if f_goterms:
        fh_go_terms = open(f_goterms, "r")
        for line in fh_go_terms:
            f = line.rstrip('\n').split('\t')
            d_transcript_annotations['go'][f[1]].append(f[2])
        fh_go_terms.close()

    if f_mirfammap:
        fh_mirna_family = open(f_mirfammap, "r")
        for l in fh_mirna_family:
            f= l.rstrip('\n').split('\t')
            d_gene_annotations['family'][f[0]] = f[3]
            d_transcript_annotations['family'][f[1]] = f[3]
        fh_mirna_family.close()


def guess_region(transcriptid, read_pos, d_transcript_annotations):
    [read_chr, read_start, read_end, read_strand] = read_pos.split(':')[-4:]
    region = 'NA'
    overlap_length = 0
    utr_start = utr_end = first_cds_start = last_cds_end = None
    read_strand = strandardize(read_strand)
    if transcriptid in d_transcript_annotations['UTR']:
        for pos in d_transcript_annotations['UTR'][transcriptid]:
            chrom = pos[0]
            start = int(pos[1])
            end = int(pos[2])
            strand = strandardize(pos[3])
            if read_chr != chrom or read_strand != strand:
                continue

            if overlap([start, end], [int(read_start), int(read_end)]) > overlap_length:
                overlap_length = overlap([start, end], [int(read_start), int(read_end)])
                region = 'UTR'
                utr_start = start
                utr_end = end

    if transcriptid in d_transcript_annotations['CDS']:
        for pos in d_transcript_annotations['CDS'][transcriptid]:
            chrom = pos[0]
            start = int(pos[1])
            end = int(pos[2])
            strand = strandardize(pos[3])
            if read_chr != chrom or read_strand != strand:
                continue
            if overlap([start, end], [int(read_start), int(read_end)]) > overlap_length:
                overlap_length = overlap([start, end], [int(read_start), int(read_end)])
                region = 'CDS'
            if not first_cds_start or start < first_cds_start:
                first_cds_start = start
            if not last_cds_end or end > last_cds_end:
                last_cds_end = end

    # if region is still a UTR
    if region == 'UTR' and utr_end <= first_cds_start:
        region = '5_prime_UTR'
        if strand == '-':
            region = '3_prime_UTR'
    elif region == 'UTR' and utr_start >= last_cds_end:
        region = '3_prime_UTR'
        if strand == '-':
            region = '5_prime_UTR'

    # used geneid in transcript annotation dict
    if transcriptid in d_transcript_annotations['gid']:
        geneid = d_transcript_annotations['gid'][transcriptid]
        if geneid in d_transcript_annotations['mature']:
            for pos in d_transcript_annotations['mature'][geneid]:
                chrom = pos[0]
                start = int(pos[1])
                end = int(pos[2])
                strand = strandardize(pos[3])
                name = pos[4]
                # print('mature', ':'.join(pos))
                if read_chr != chrom and read_strand != strand:
                    continue
                if overlap([start, end], [int(read_start), int(read_end)]) > overlap_length:
                    overlap_length = overlap([start, end], [int(read_start), int(read_end)])
                    region = name
    return region


def match_positions(cigar):
    matches = re.findall(r'(\d+)([MINSH=X]{1})', cigar)
    match_start = match_end = 0
    for m in matches:
        if not match_start:
            if m[1] != "M" or m[1] != "=":
                match_start = match_end
        elif m[1] == "S" and match_end:
            break
        match_end += int(m[0])
    return match_start, match_end


def cigars_overlapping(cigar1, cigar2):
    match_start1, match_end1 = match_positions(cigar1)
    match_start2, match_end2 = match_positions(cigar2)
    if max(0, min(match_end1, match_end2) - max(match_start1, match_start2)) > 0:
        return True
    else:
        return False