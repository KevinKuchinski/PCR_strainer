# PCR_strainer
PCR_strainer is a tool for assessing the inclusivity of primer and probe oligonucleotides from diagnostic qPCR assays and amplicon sequencing schemes. It depends on thermonucleotideBLAST (TNTBLAST), which conducts local alignments between query oligonucleotides and subject sequences that include a thermodynamic assessment of the alignment. PCR_strainer parses and tabulates the TNTBLAST output to generate reports on assay performance and sequence variants in oligo sites. TNTBLAST will need to be installed separately from: https://github.com/jgans/thermonucleotideBLAST

# PCR_strainer Usage
<b>Usage example:</b>
```
$ pcr_strainer -a <assay CSV file> -g <genomes FASTA file> -o <output dir>/<output name> [<optional args>]
```
<b>Required arguments:</b>

	-a : path to assay details in CSV file
	-g : path to target genomes in FASTA file
	-o : path to output directory and name to append to output files

<b>Optional arguments:</b>

	-t : minimum prevalence (%) of primer site variants reported in reports (default = 0 (reports everything), min > 0, max < 100)
	-m : minimum Tm (degrees C) for primers and probes (default = 45)
	-p : molar concentration of primer oligos (uM) (default = 1, min > 0)
	-P : molar concentration of probe oligos (uM) (default = 1, min > 0)

<b>The assay CSV file:</b>

PCR_strainer expects a csv file where each line describes a PCR assay using the following format: 
assay_name, forward_primer_name, forward_primer_seq, reverse_primer_name, reverse_primer_seq, probe_name, probe_seq
Example assay file entry:
```
BCCDC_SARS2_RdRP,BCCDC_RdRP_Fwd,TGCCGATAAGTATGTCCGCA,BCCDC_RdRP_Rev,CAGCATCGTCAGAGAGTATCATCATT,BCCDC_RdRP_Probe,TTGACACAGACTTTGTGAATG
```
  * all oligo sequences should be writen in the 5' to 3' orientation
  * degenerate nucleotides are permitted in the assay oligo sequences
  * the probe name and probe sequence can be omitted for a conventional PCR
  * for amplicon sequencing schemes, enter primer pairs as lines in the same file

<b>The reference genomes</b>: 

PCR_strainer expects DNA sequences in FASTA format without spaces in the header. For single-stranded genomes, ensure all sequences represent the same sense (e.g. all coding strand). We recommend you filter your reference genomes to remove sequences containing degenerate nucleotides in target locations to limit false negatives; thermonucleotideBLAST does not expand degenerate nucleotide possibilities for the subject sequences.

<b>The name of the output</b>: 

PCR_strainer generates four TSV files. The output name will be appended to these file names (no spaces). Including a file path before the output name will write output files to that directory.

# PCR_Strainer Reports
PCR_strainer generates three report files and a table of raw results from TNTBLAST.

## assay_report
The assay_report indicates how many reference sequences are impacted by nucleotide mismatches and gaps acrosss all oligos for each assay. Filter this table for rows with 0 in the <b>errors</b> columns for quick overview of assay inclusivity; this will quickly show what percentage are the provided reference sequences had no gaps or mismatches against the provided assays.

<b>COLUMN : DESCRIPTION

  assay_name</b> : The name of the assay from the assay CSV file

  <b>total_targets</b> : The total number of reference sequences in the genomes file

  <b>detected_targets</b> : The number of reference sequences in the genomes file in which thermonucleotideBLAST was able to identify all oligo sites and generate an amplicon

  <b>perc_detected</b> : detected_targets as a percentage of total_targets

  <b>total_errors</b> : The number of nucleotide errors across all of the assay's oligonucleotides; this includes gaps and the total number of unannealed nucleotides (including those impacted by nearby mismatches despite having complementary base pairing)

  <b>target_count</b>: The number of reference sequences with the indicated number of errors for this assay

  <b>perc_of_detected</b> : target_count as a percentage of detected_targets
  
  <b>perc_of_total</b> : target_count as a percentage of total_targets

## variant_report
The variant_report provides information about locations in the provided reference sequences that are targeted by assay oligos, but contain gaps and mismatches. This report identifies common sequence variants in oligo sites, facilitating oligo re-design. In oligo site variant sequences, mismatched bases are written in lower case, deletions are indicated with dashes, and insertions are surrounded by parentheses. 

<b>COLUMN : DESCRIPTION

  assay_name</b> : The name of the assay from the assay file
  
  <b>oligo</b> : Forward primer, reverse primer, or probe
  
  <b>oligo_name</b> : Name of the oligo from the assay file
  
  <b>oligo_seq</b> : The sequence of the oligo provided in the assay file

  <b>total_targets</b> : The total number of reference sequences in the genomes file

  <b>detected_targets</b> : The number of reference sequences in the genomes file in which thermonucleotideBLAST was able to identify all oligo sites and generate an amplicon

  <b>perc_detected</b> : detected_targets as a percentage of total_targets

  <b>oligo_site_variant</b> : The variant sequence at the oligo site, written in 'oligo sense', ie the same sense as the PCR oligo
  
  <b>oligo_errors</b> : The number of nucleotide errors present at this variant site; this includes gaps and the total number of unannealed nucleotides (including those impacted by nearby mismatches despite having complementary base pairing)

  <b>target_count</b>: The number of reference sequences with the indicated oligo site variant

  <b>perc_of_detected</b> : target_count as a percentage of detected_targets
  
  <b>perc_of_total</b> : target_count as a percentage of total_targets

## missed_seqs_report
The missed_report provides the name of target reference sequences in the genomes files that were not aligned by thermonucleotideBLAST. PCR_strainer provides the headers of these missed targets for trouble-shooting assays with high percentages of missed targets (i.e. low perc_detected values). These targets are generally either a) poor quality and contain too many Ns in/around the oligo target sites, or b) too divergent from the oligos.

<b>COLUMN : DESCRIPTION

  assay_name</b> : The name of the assay from the assay file
  
  <b>target</b> : The FASTA header of the missed reference sequence in the genomes file
  
  <b>target_length</b> : The length of the target sequence in nucleotides

  <b>total_Ns</b> : The number of nucleotide positions in the target sequence represented by ambiguous N bases
  
  <b>perc_Ns</b> : The percentage of the target sequence represented by ambiguous N bases

Questions, feedback, and bug reports are welcome! kevin.kuchinski@bccdc.ca
