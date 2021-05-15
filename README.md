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
  -o : the name for the PCR_strainer output
  
  eg. python PCR_strainer.py -a BCCDC_SARS-CoV-2_RdRP.csv -g SARS-CoV-2_genomes.fasta -t 1 -o BCCDC_SARS-CoV-2_RdRP_results

The assay file: PCR_strainer expects a csv file where each line describes a PCR assay using the following format:
  assay_name,forward_primer_name,reverse_primer_name,probe_name,forward_primer_seq,reverse_primer_seq,probe_seq
  
  eg. BCCDC_SARS-CoV-2_RdRP,BCCDC_SARS2_RdRP_Fwd,BCCDC_SARS2_RdRP_Rev,BCCDC_SARS2_RdRP_Probe,TGCCGATAAGTATGTCCGCA,CAGCATCGTCAGAGAGTATCATCATT,TTGACACAGACTTTGTGAATG
  * all oligo sequences should be writen in the 5' to 3' orientation
  ** degenerate nucleotides are permitted in the assay oligo sequences
  *** the probe name and probe sequence can be omitted for a conventional PCR
  **** for amplicon sequencing schemes, enter primer pairs as lines in the same file
 
The reference genomes: PCR_strainer expects DNA sequences in FASTA format without spaces in the header. For single-stranded genomes, ensure all sequences represent the same sense (ie. all coding strand). We recommend you filter your reference genomes to remove sequences containing degenerate nucleotides in target locations to limit false negatives; thermonucleotideBLAST does not expand degenerate nucleotide possibilities for the subject sequences.

The name of the output: PCR_strainer generates two TSV files. The output name will be appended to these file names (no spaces).

# PCR_Strainer Reports
PCR_strainer generates two report files. One provides metrics for the overall performance of each assay. The other provides metrics for the performance of individual oligonucleotides.

## assay_report
The assay_report indicates how many reference sequences are impacted by nucleotide mismatches and gaps acrosss all oligos for each assay. Filter this table for rows with 0 in the <b>errors</b> columns for quick overview of assay inclusivity; this will quickly show what percentage are the provided reference sequences had no gaps or mismatches against the provided assays.

<b>COLUMN : DESCRIPTION<b/>

  <b>assay</b> : The name of the assay from the -assay file
  
  <b>targets_detected</b> : The number of reference sequences in the -genomes file that were analyzed by tntblast
  
  <b>targets_total</b> : The total number of reference sequences in the -genomes file
  
  <b>perc_detected</b> : Percentage of reference sequences analyzed by tntblast; missed reference sequences either indicate poor-quality sequences in the -genomes file or sequences with extreme divergence from the assay oligonucleotides
  
  <b>errors</b> : The number of nucleotide errors between all of the assay's oligonucleotides and oligo target sites; this includes gaps and the total number of unannealed nucleotides (including those with complementary base pairing impacted by nearby mismatches)
  
  <b>count</b>: The number of reference sequences with that row's number of errors against that row's assay
  
  <b>perc_of_detected</b> : Percentage of detected reference sequences with that row's number of errors against that row's assay
  
  <b>perc_of_total</b> : Percentage of total reference sequences with that row's number of errors against that row's assay

## variant_report
The variant_report provides information about locations in the provided reference sequences that are targeted by assay oligos, but contain gaps and mismatches. This report identifies common oligo site variants, facilitating oligo re-design.

<b>COLUMN : DESCRIPTION</b>

  <b>assay</b> : The name of the assay from the -assay file
  
  <b>targets_detected</b> : The number of reference sequences in the -genomes file that were analyzed by tntblast
  
  <b>targets_total</b> : The total number of reference sequences in the -genomes file
  
  <b>perc_detected</b> : Percentage of reference sequences analyzed by tntblast; missed reference sequences either indicate poor-quality sequences in the -genomes file or sequences with extreme divergence from the assay oligonucleotides
  
  <b>oligo</b> : Forward primer, reverse primer, or probe
  
  <b>oligo_name</b> : The name of the oligo provided in the -assay file
  
  <b>oligo_seq</b> : The sequence of the oligo provided in the -assay file
  
  <b>oligo_site_variant</b> : The sequence of the hypothetical oligo that would be a perfect match for this oligo target site variant
  
  <b>errors</b> : The number of nucleotide errors between this row's oligo and this row's target site variant ; this includes gaps and the total number of unannealed nucleotides (including those with complementary base pairing impacted by nearby mismatches)
  
  <b>count</b>: The number of reference sequences with this row's oligo site variant
  
  <b>perc_of_detected</b> : Percentage of detected reference sequences with this row's oligo site variant; this row's value is used for by -threshold option
  
  <b>perc_of_total</b> : Percentage of total reference sequences with this row's oligo site variant

Questions, feedback, and bug reports are welcome! kevin.kuchinski@bccdc.ca
