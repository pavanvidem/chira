# ChiRA - Chimeric Read Analyzer

ChiRA is a set of tools to analyze RNA-RNA interactome experimental data such as CLASH, CLEAR-CLIP, PARIS, SPLASH etc. Following are the descriptions of the each tool. Here we provide descriptions about the input and ouptput files. For the detailed description of the other parameters please look at the help texts of tools.

## Development environment

Create the environment containing the Python dependencies and external pipeline tools:

```sh
conda env create -f environment.yml
conda activate chira-refactor
python -B -m unittest discover -s tests -v
```

The test suite includes real BEDTools and Blockbuster integration checks. These checks are skipped if their executables are unavailable.

## chira_collapse.py
* Description: This tool deduplicates reads from a FASTQ file and writes each unique sequence/UMI pair once to FASTA with its read count. Output is sorted by trimmed sequence and UMI for reproducible sequence IDs.
* Inputs
  * Quality and adapter trimmed FASTQ file
  * Optional `-u/--umi_len`: number of bases to trim from the 5' end as a UMI (default: 0). The length must be non-negative and strictly less than the length of every read. Reads that would have an empty sequence after trimming cause an error identifying the read, before the output file is opened.
* Outputs
  * FASTA file with unique sequences. Headers use `>sequence_id|read_count` without UMIs, or `>sequence_id|UMI|read_count` with UMIs. Identical trimmed sequences with different UMIs are kept as separate records.

## chira_map.py
* Description: This tool handles the mapping of the reads to reference transcriptome. User can choose between the bwa-mem and CLAN alignment tools.
* Inputs
  * A fasta file containing reads
  * A FASTA file containing reference transcript sequences
  * An optional second reference FASTA file, if split-reference is used
* Outputs
  * BED file containing the alignments
  * Optional unmapped FASTA file
* Mapping behavior
  * BWA runs both long- and short-segment iterations on the input reads for each supplied reference index. `--nhits1` and `--nhits2` are passed to BWA's `-h` option (defaults: 50 and 100).
  * Primary alignments on the requested strand are retained. XA alternatives on that strand are compared by aligned query span (`M/I/=/X`), not by reference span; only alternatives attaining the maximum span for that BAM record are retained.
  * External commands use literal argument lists and failures stop processing, including failures on either side of the BWA/SAMtools pipeline. Paths containing spaces are supported.
  * Existing indices require `--index1`. Building indices requires `--build` and `--ref_fasta1`. CLAN also needs reference FASTA files when using existing indices, with the second reference and index provided together.
  * CLAN output lacks strand information, so converted CLAN BED records retain the existing `+` strand representation.

## chira_merge.py
* Description: This tool merges the overlapping aligned positions to define the read concentrated loci. If an annotation GTF file produced, the transcriptomic alignment positions are first converted to their corresponding genomic positions.
* Inputs
  * Alignments in BED format
  * An annotation GFF file contaning reference genomic positions.
  * A FASTA file containing reference transcript sequences
  * An optional second reference FASTA file, if split-reference is used
* Outputs
  * BED file containing the alignments with reads categorized into segments depending on which part of the read is aligned.
  * Five-column merged locus file: reference, start, end, strand, and semicolon-separated alignment descriptions.
* Merging behavior
  * Input alignments must be grouped by read ID, as in the mapping stage's `sorted.bed`. Segment positions use inclusive read coordinates; locus merging uses half-open BED coordinates. Adjacent BED intervals do not overlap.
  * The overlap method greedily merges intervals when their shared bases meet the overlap fraction of either interval. It retains all ties for the longest reference span per read/transcript pair within a locus.
  * Annotation exons are deduplicated by position and sorted in transcript orientation, independently of input feature order. Transcript identifiers come from `transcript_id`, falling back to `Name`.
  * Transcript-to-genome conversion clips alignments to each intersected exon and keeps the longest genomic fragment for spliced alignments, including the final alignment. Ties are resolved lexicographically. Genomic strand combines transcript strand and alignment orientation.
  * External sorting, Blockbuster and BEDTools commands use checked argument lists; failures stop processing. The output directory is created after argument validation.

## chira_quantify.py
* Description: This tool first creates CRLS from merged BED file and quantifies them based on the mapped reads.
* Inputs
  * A BED file containing alignment segments
  * A tabular  file containing merged alignments
* Outputs
  * Tabular file containing the reads and their CRLs with TPM values.
* Quantification behavior
  * With `--build_crls_too`, eligible loci are grouped using Jaccard similarity of their read-segment sets. All existing CRLs are considered, including CRL 0. `--crl_share 1` requires identical sets; loci with no shared segments remain separate even at cutoff 0.
  * A second pass retains each locus's highest qualifying CRL shares, including ties. Loci without a qualifying membership become standalone CRLs; empty CRLs are removed before numbering.
  * Each unique read segment contributes one count to EM, independently of the collapsed FASTA read count. Duplicate alignment rows do not add counts. EM stops when the sum of absolute abundance changes meets `--em_threshold`, or reports an error after 1000 unconverged iterations.
  * TPM uses the median half-open BED locus length (`end - start`) for each CRL. Nonzero expression is normalized to a total of one million TPM; empty input produces an empty counts file.
  * Output retains the existing 13-column format and GNU version-sort order. Sorting failures stop processing and preserve intermediate files.

## chira_extract.py
* Description: This tool extracts the best chimeric alignments for each read. User can optionally hybridize the loci where the chimeric arms are mapping to. 
* Inputs
  * Tabular file containg CRLs information
  * Annotation GFF file
  * A FASTA file containing reference transcript sequences
  * An optional second reference FASTA file, if split-reference is used  
  * If your alignments are merged at genomic level in previous step (chira_merge.py), then provide a reference genomic FASTA fille.
* Outputs
  * Tabular file of chimeric reads information
  * Tabular file of detected interactions
* Extraction behavior
  * Counts must be grouped by read ID. Worker ranges partition whole read groups without gaps or duplication; worker errors stop processing before results are merged. Annotation data is explicitly passed to workers.
  * Empty counts produce header-only chimera/singleton files. The TPM cutoff uses the existing percentile of distinct CRL TPM values; cutoff 1 selects the maximum value.
  * Annotation exon lengths use deduplicated half-open intervals. Missing annotation lengths fall back independently to reference FASTA lengths for both chimeras and singletons. UTRs without CDS information remain labeled `UTR`.
  * BEDTools, IntaRNA and sorting commands are checked. Hybrid coordinates are projected back from oriented sequences, including reverse-strand loci. Sequences are retained when IntaRNA predicts no interaction.
  * Existing output columns and scoring are preserved: filtering and singleton scores use locus share multiplied by segment probability; chimera arm scores use segment probability and the combined score is their product. Best hits are selected by TPM, retaining ties.
  * Both `--summarize` and the original `--summerize` spelling are accepted.

## Shared coordinate conventions

* BED and parsed annotation intervals are zero-based and end-exclusive. Raw GTF/GFF coordinates are one-based and inclusive before parsing. `chira_utilities.overlap()` uses the half-open convention.
* CIGAR-derived positions on the original read are one-based and inclusive, including offsets from hard-clipped bases. `match_positions()` returns these positions in original-read orientation on either strand; `inclusive_overlap()` counts their shared bases. Chimera overlap limits and segment overlap fractions use this inclusive convention.
* `query_length()` includes hard and soft clips to represent the original read span. `aligned_query_length()` counts `M/I/=/X`; `reference_length()` counts `M/D/N/=/X`. The legacy `alignment_length()` counts `M/I/D/=/X` and is retained for existing callers.
* `alignment_end()` accepts a one-based SAM reference start and returns its inclusive end, which is numerically equal to the zero-based BED exclusive end.
* CIGAR helpers reject malformed operations, zero-length operations and misplaced clipping. Read match positions require at least one `M`, `=` or `X` operation.
