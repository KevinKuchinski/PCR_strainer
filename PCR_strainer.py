import argparse
import subprocess
import os
import numpy as np
import pandas as pd

parser=argparse.ArgumentParser()
parser.add_argument("-a", "--assay", type=str, required=True, help="input TSV file containing assay name, primer and probe oligo sequences, and (optionally) reaction molarities")
parser.add_argument("-g", "--genomes", type=str, required=True, help="input FASTA file containing sequences assay is meant to detect")
parser.add_argument("-t", "--threshold", type=float, required=True, help="reporting threshold, report oligo if more than this percentage of target sequences are not a perfect match (0-100)")
parser.add_argument("-o", "--output", type=str, required=True, help="brief descriptor to append to output file names")
args=parser.parse_args()

print()
print('PCR_strainer - A tool for assessing inclusivity of PCR primers')
print()

## Prepare/overwrite existing output files for tabulated PCR results
outputFile=open(args.output+'_PCR_results.tsv','w')
oligos=['fwd_primer','rev_primer','probe']
outputFields=['assay']+['target']+[oligo+field for oligo in oligos for field in ['_name','_seq', '_loc_seq','_tm','_mismatches','_gaps','_errors']]
outputFields+=['amplicon'+field for field in ['_coords','_length','_loc_seq']]
outputFileHeader='\t'.join(outputFields)
outputFile.write(outputFileHeader+'\n')
outputFile.close()

## Parse comma-separated line from input file containing assay details and create separate file for each assay formatted in TSV for TNTBLAST
def createAssayFile(assayDetails,args):
    if len(assayDetails)==7: ## Taqman assay
        print('Creating Taqman-based assay file for '+assayDetails[0]+'...')
        assayFile=open(args.output+'_'+assayDetails[0]+'_TNTBLAST_assay.tsv','w') ## Create TSV file containing each separate assay's oligos for use by TNTBLAST                                           
        assayFile.write(assayDetails[0]+'\t'+assayDetails[4]+'\t'+assayDetails[5]+'\t'+assayDetails[6]+'\n') ## Write assay name, fwd primer seq, rev primer seq, probe seq
        assayFile.close()
    elif len(assayDetails)==5: ## Standard assay
        print('Creating standard assay file for '+assayDetails[0]+'...')
        assayFile=open(args.output+'_'+assayDetails[0]+'_TNTBLAST_assay.tsv','w') ## Create TSV file containing each separate assay's oligos for use by TNTBLAST                                           
        assayFile.write(assayDetails[0]+'\t'+assayDetails[3]+'\t'+assayDetails[4]+'\n') ## Write assay name, fwd primer seq, rev primer seq
        assayFile.close()

## Run TNTBLAST from user-provided assay details
def runTNTBLAST(assayDetails,args):
    print('Running TNTBLAST for',assayDetails[0],'against',args.genomes,'...')
    targets=args.genomes
    tm=50
    command='tntblast -i '+ args.output+'_'+assayDetails[0] +'_TNTBLAST_assay.tsv' +' -d '+ targets ## Create command for running TNTBLAST
    command+=' -o '+ args.output + '_' + assayDetails[0] +'_TNTBLAST_output.txt ' ## Create command for running TNTBLAST (cont)   
    command+='-e 30 -E 30 --best-match -m 0 -v F' ## Create command for running TNTBLAST (cont)    
    subprocess.call(command,shell=True) ## Run TNTBLAST command
    command='rm '+ args.output+'_'+assayDetails[0] +'_TNTBLAST_assay.tsv'
    subprocess.call(command,shell=True)

## Parse TNTBLAST output into TSV file
def parseOutput(assayDetails,args):
    complementaryBases={'A':'T','T':'A','G':'C','C':'G','W':'W','S':'S','M':'K','K':'M','R':'Y','Y':'R','B':'V','V':'B','D':'H','H':'D','N':'N','-':'-'}
    reverseComplement=lambda seq: ''.join([complementaryBases[base] for base in seq[::-1]])
    print('Parsing TNTBLAST output from',assayDetails[0],'...')
    ## List of relevant fields in TNTBLAST output
    oligos=['forward primer','reverse primer','probe']
    tntFields=[oligo+field for oligo in oligos for field in ['',' tm',' mismatches',' gaps']]
    tntFields+=['amplicon'+field for field in [' range',' length']]+['probe range']
    outputFields=['assay']+['target']+[oligo+field for oligo in oligos for field in [' name','', ' loc seq',' tm',' mismatches',' gaps',' errors']]
    outputFields+=['amplicon'+field for field in [' range',' length',' seq']]
    ## Open PCR_strainer output file
    outputFile=open(args.output+'_PCR_results.tsv','a')
    ## Open TNTBLAST output file and read into memory
    with open(args.output + '_' + assayDetails[0] +'_TNTBLAST_output.txt','r') as inputFile:
        inputLines=inputFile.readlines()
    for i in range(len(inputLines)):
        if ' = ' in inputLines[i]:
            key=inputLines[i].split(' = ')[0]
            value=inputLines[i].split(' = ')[1].rstrip()
            if key=='name':
                entry={}
                for key in outputFields:
                    entry[key]='---'
                entry['assay']=value
            elif key in tntFields:
                entry[key]=value
        elif inputLines[i][0]=='>':
            entry['target']=inputLines[i].rstrip()
            entry['amplicon seq']=inputLines[i+1].rstrip()
            for oligo in oligos:
                if entry[oligo]!='---':
                    entry[oligo+' errors']=0
                    if entry[oligo+' mismatches']=='---':
                        entry[oligo+' mismatches']=0
                    entry[oligo+' errors']+=int(entry[oligo+' mismatches'])
                    if entry[oligo+' gaps']=='---':
                        entry[oligo+' gaps']=0
                    entry[oligo+' errors']+=int(entry[oligo+' gaps'])
            entry['forward primer']=entry['forward primer'].lstrip(" -'53").rstrip(" -'53")
            entry['reverse primer']=entry['reverse primer'].lstrip(" -'53").rstrip(" -'53")
            entry['probe']=entry['probe'].lstrip(" -'53").rstrip(" -'53") if len(assayDetails)==7 else '---'
            entry['forward primer name']=assayDetails[1]
            entry['reverse primer name']=assayDetails[2]
            entry['probe name']=assayDetails[3] if len(assayDetails)==7 else '---'
            entry['forward primer loc seq']=entry['amplicon seq'][:len(entry['forward primer'])]
            entry['reverse primer loc seq']=reverseComplement(entry['amplicon seq'][-len(entry['reverse primer']):])
            entry['forward primer loc seq']=''.join([locBase if locBase==oligoBase else locBase.lower() for locBase,oligoBase in zip(entry['forward primer loc seq'],entry['forward primer'])])
            entry['reverse primer loc seq']=''.join([locBase if locBase==oligoBase else locBase.lower() for locBase,oligoBase in zip(entry['reverse primer loc seq'],entry['reverse primer'])])
            if len(assayDetails)==7:
                start=int(entry['probe range'].split(' .. ')[0])-int(entry['amplicon range'].split(' .. ')[0])
                end=int(entry['probe range'].split(' .. ')[1])-int(entry['amplicon range'].split(' .. ')[1])
                entry['probe loc seq']=entry['amplicon seq'][start:end]
                entry['probe loc seq']=''.join([locBase if locBase==oligoBase else locBase.lower() for locBase,oligoBase in zip(entry['probe loc seq'],entry['probe'])])
            outputFile.write('\t'.join([str(entry[field]) for field in outputFields])+'\n')
    outputFile.close()
    command='rm '+args.output + '_' + assayDetails[0] +'_TNTBLAST_output.txt'
    subprocess.call(command,shell=True)

