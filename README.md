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
  
  eg. BCCDC_SARS-CoV-2_RdRP,TGCCGATAAGTATGTCCGCA,CAGCATCGTCAGAGAGTATCATCATT,TTGACACAGACTTTGTGAATG
  * all oligo sequences should be writen in the 5' to 3' orientation
  ** degenerate nucleotides are permitted in the assay oligo sequences
  *** the probe name and probe sequence can be omitted for a conventional PCR
  **** for amplicon sequencing schemes, enter primer pairs as lines in the same file
 
The reference genomes: PCR_strainer expects DNA sequences in FASTA format without spaces in the header. For single-stranded genomes, ensure all sequences represent the same sense (ie. all coding strand). We recommend you filter your reference genomes to remove sequences containing degenerate nucleotides in target locations to limit false negatives; thermonucleotideBLAST does not expand degenerate nucleotide possibilities for the subject sequences.

The name of the output: PCR_strainer generates several TSV files. The output name will be appended to these file names (no spaces).

Questions, feedback, and bug reports are welcome! kevin.kuchinski@bccdc.ca
