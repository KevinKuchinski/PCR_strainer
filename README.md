# PCR_strainer
A tool for assessing the inclusivity of primer and probe oligonucleotides from diagnostic qPCR assays and amplicon sequencing schemes.

# PCR_strainer Setup
Before you use PCR_strainer, you will have to install thermonucleotideBLAST. This extremely helpful tool conducts local alignments between query oligonucleotides and subject sequences that include a thermodynamic assessment of the alignment. Installation instructions are available from the author's Github:
https://github.com/jgans/thermonucleotideBLAST

Extensive documentation can be found at:
https://public.lanl.gov/jgans/tntblast/tntblast_doc.html

After installing thermonucleotideBLAST, you will have to create an environment for Python 3 with the modules Numpy and Pandas:
```
conda create -n PCR_strainer python=3.7.1 pandas=0.20.3
```
After creating your environment for PCR_strainer, you can download the PCR_strainer script:
```
git clone https://github.com/KevinKuchinski/PCR_strainer.git
```

# PCR_strainer Usage
Running PCR_strainer requires four arguments:

  -a : the assay file
  
  -g : the reference genomes
  
  -t : the variant prevlance threshold for reporting (0-100)
  
  -o : the name for the PCR_strainer output (and a path to the output directory)
  
  eg. python PCR_strainer.py -a BCCDC_SARS-CoV-2_RdRP.csv -g SARS-CoV-2_genomes.fasta -t 1 -o BCCDC_SARS-CoV-2_RdRP_results

<b>The assay file</b>: PCR_strainer expects a csv file where each line describes a PCR assay using the following format:
  assay_name,forward_primer_name,forward_primer_seq,reverse_primer_name,reverse_primer_seq,probe_name,probe_seq
  
  Example assay file entry:
  ```
  BCCDC_SARS2_RdRP,BCCDC_RdRP_Fwd,TGCCGATAAGTATGTCCGCA,BCCDC_RdRP_Rev,CAGCATCGTCAGAGAGTATCATCATT,BCCDC_RdRP_Probe,TTGACACAGACTTTGTGAATG
  ```
  * all oligo sequences should be writen in the 5' to 3' orientation
  * degenerate nucleotides are permitted in the assay oligo sequences
  * the probe name and probe sequence can be omitted for a conventional PCR
  * for amplicon sequencing schemes, enter primer pairs as lines in the same file
 
<b>The reference genomes</b>: PCR_strainer expects DNA sequences in FASTA format without spaces in the header. For single-stranded genomes, ensure all sequences represent the same sense (e.g. all coding strand). We recommend you filter your reference genomes to remove sequences containing degenerate nucleotides in target locations to limit false negatives; thermonucleotideBLAST does not expand degenerate nucleotide possibilities for the subject sequences.

<b>The name of the output</b>: PCR_strainer generates four TSV files. The output name will be appended to these file names (no spaces). Including a file path before the output name will write output files to that directory.

# PCR_Strainer Reports
PCR_strainer generates three report files and a table of raw results from TNTBLAST.

## assay_report
The assay_report indicates how many reference sequences are impacted by nucleotide mismatches and gaps acrosss all oligos for each assay. Filter this table for rows with 0 in the <b>errors</b> columns for quick overview of assay inclusivity; this will quickly show what percentage are the provided reference sequences had no gaps or mismatches against the provided assays.

<b>COLUMN : DESCRIPTION

  assay</b> : The name of the assay from the -assay file
  
  <b>total_errors</b> : The number of nucleotide errors among all of the assay's oligonucleotides; this includes gaps and the total number of unannealed nucleotides (including those with complementary base pairing impacted by nearby mismatches)
  
  <b>target_count</b>: The number of reference sequences with the indicated number of errors for the indicated assay
  
  <b>detected_targets</b> : The number of reference sequences in the -genomes file that were analyzed by thermonucleotideBLAST
  
  <b>total_targets</b> : The total number of reference sequences in the -genomes file
  
  <b>perc_detected</b> : detected_targets as a percentage of total_targets
  
  <b>perc_of_detected</b> : target_count as a percentage of detected_targets
  
  <b>perc_of_total</b> : target_count as a percentage of total_targets

## variant_report
The variant_report provides information about locations in the provided reference sequences that are targeted by assay oligos, but contain gaps and mismatches. This report identifies common oligo site variants, facilitating oligo re-design.

<b>COLUMN : DESCRIPTION

  assay</b> : The name of the assay from the -assay file
  
  <b>oligo</b> : Forward primer, reverse primer, or probe
  
  <b>oligo_seq</b> : The sequence of the oligo provided in the -assay file
  
  <b>oligo_site_variant</b> : The sequence of the hypothetical oligo that would be a perfect match for the target location in the reference sequence
  
  <b>oligo_errors</b> : The number of nucleotide errors between this row's oligo and this row's target site variant ; this includes gaps and the total number of unannealed nucleotides (including those with complementary base pairing impacted by nearby mismatches)
  
  <b>target_count</b>: The number of reference sequences with the indicated oligo site variant
  
  <b>detected_targets</b> : The number of reference sequences in the -genomes file that were detected and analyzed by thermonucleotideBLAST
  
  <b>total_targets</b> : The total number of reference sequences in the -genomes file
  
  <b>perc_detected</b> : detected_targets as a percentage of total_targets
  
  <b>perc_of_detected</b> : target_count as a percentage of detected_targets
  
  <b>perc_of_total</b> : target_count as a percentage of total_targets

## missed_report
The missed_report provides the name of target reference sequences in the -genomes files that were not aligned by thermonucleotideBLAST. PCR_strainer provides the headers of these missed targets for trouble-shooting assays with high percentages of missed targets. These targets are generally either a) poor quality and contain too many Ns in/around the oligo target sites, or b) too divergent from the oligos.

<b>COLUMN : DESCRIPTION

  assay</b> : The name of the assay from the -assay file
  
  <b>target</b> : The FASTA header of the missed reference sequence in the -genomes file
  
  <b>target_length</b> : The length of the target sequence
  
  <b>perc_Ns</b> : The percentage of the target sequence represented by ambiguous N bases

Questions, feedback, and bug reports are welcome! kevin.kuchinski@bccdc.ca
