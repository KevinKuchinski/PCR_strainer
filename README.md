# PCR_strainer
A tool for assessing the inclusivity of primer and probe oligonucleotides from Taqman qPCR assays against numerous reference genome sequences

Before you use PCR_strainer, you will have to install Thermonucleotide BLAST: https://public.lanl.gov/jgans/tntblast/tntblast_doc.html

Running PCR_strainer requires three arguments:
  -a : the assay file
  -g : the reference genomes
  -o : the name for the PCR_strainer output
  
  eg. python PCR_strainer.py -a BCCDC_SARS-CoV-2_RdRP.csv -g SARS-CoV-2_genomes.fasta -o BCCDC_SARS-CoV-2_RdRP_results

The assay file: PCR_strainer expects a csv file where each line describes a Taqman qPCR assay using the following format:
  assay_name,forward_primer_seq,reverse_primer_seq,probe_seq
  
  eg. BCCDC_SARS-CoV-2_RdRP,TGCCGATAAGTATGTCCGCA,CAGCATCGTCAGAGAGTATCATCATT,TTGACACAGACTTTGTGAATG
  * all oligo sequences should be writen in the 5' to 3' orientation
  ** degenerate nucleotides are permitted in the assay oligo sequences

The reference genomes: PCR_strainer expects DNA sequences in FASTA format without spaces in the header. For single-stranded genomes, ensure all sequences represent the same sense (ie. all coding strand). We recommend you filter your reference genomes to remove sequences containing degenerate nucleotides in target locations to limit false negatives; Thermonucleotide BLAST does not expand degenerate nucleotide possibilities for the alignment subject sequences.

The name of the output: PCR_strainer generates several TSV and FASTA files. The output name will be appended to these file names (no spaces).

Questions, feedback, and bug reports are welcome! kevin.kuchinski@bccdc.ca
