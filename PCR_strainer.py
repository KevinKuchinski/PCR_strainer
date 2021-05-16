import argparse
import subprocess
import os
import numpy as np
import pandas as pd

parser=argparse.ArgumentParser()
parser.add_argument("-a", "--assay", type=str, required=True, help="input TSV file containing assay name, primer and probe oligo sequences, and (optionally) reaction molarities")
parser.add_argument("-g", "--genomes", type=str, required=True, help="input FASTA file containing sequences assay is meant to detect")
parser.add_argument("-t", "--threshold", type=float, required=True, help="reporting threshold, report oligo if more than this percentage of target sequences are not a perfect match (0-100)")
parser.add_argument("-m", "--melting", type=float, required=False, default=45, help="melting temperature cutoff in degrees Celcius (default=45)")
parser.add_argument("-o", "--output", type=str, required=True, help="brief descriptor to append to output file names")
args=parser.parse_args()

## Parse comma-separated line from input file containing assay details and create separate file for each assay formatted in TSV for TNTBLAST
def createAssayFile(assayDetails,args):
    if len(assayDetails)==7: ## Taqman assay
        print('Creating Taqman PCR assay file for '+assayDetails[0]+'...')
        assayFile=open(args.output+'_'+assayDetails[0]+'_TNTBLAST_assay.tsv','w') ## Create TSV file containing each separate assay's oligos for use by TNTBLAST                                           
        assayFile.write(assayDetails[0]+'\t'+assayDetails[4]+'\t'+assayDetails[5]+'\t'+assayDetails[6]+'\n') ## Write assay name, fwd primer seq, rev primer seq, probe seq
        assayFile.close()
    elif len(assayDetails)==5: ## Standard assay
        print('Creating standard PCR assay file for '+assayDetails[0]+'...')
        assayFile=open(args.output+'_'+assayDetails[0]+'_TNTBLAST_assay.tsv','w') ## Create TSV file containing each separate assay's oligos for use by TNTBLAST                                           
        assayFile.write(assayDetails[0]+'\t'+assayDetails[3]+'\t'+assayDetails[4]+'\n') ## Write assay name, fwd primer seq, rev primer seq
        assayFile.close()

## Run TNTBLAST from user-provided assay details
def runTNTBLAST(assayDetails,args):
    print('Running TNTBLAST for',assayDetails[0],'against',args.genomes,'...')
    targets=args.genomes
    tm=args.melting
    command='tntblast -i '+ args.output+'_'+assayDetails[0] +'_TNTBLAST_assay.tsv' +' -d '+ targets ## Create command for running TNTBLAST
    command+=' -o '+ args.output + '_' + assayDetails[0] +'_TNTBLAST_output.txt ' ## Create command for running TNTBLAST (cont)   
    command+='-e '+str(tm)+' -E '+str(tm)+' --best-match -m 0 -v F' ## Create command for running TNTBLAST (cont)    
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
    outputFields+=['total errors']
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
                entry['total errors']=sum([entry[oligo+' errors'] for oligo in ['forward primer','reverse primer','probe']])
            else:
                entry['total errors']=sum([entry[oligo+' errors'] for oligo in ['forward primer','reverse primer']])
            outputFile.write('\t'.join([str(entry[field]) for field in outputFields])+'\n')
    outputFile.close()
    command='rm '+args.output + '_' + assayDetails[0] +'_TNTBLAST_output.txt'
    subprocess.call(command,shell=True)

## MAIN
print()
print('PCR_strainer - A tool for assessing inclusivity of PCR primers and probes against reference nucleotide sequences')

## Prepare/overwrite existing output files for tabulated PCR results                                                                                                                                      
outputFile=open(args.output+'_PCR_results.tsv','w')
oligos=['fwd_primer','rev_primer','probe']
outputFields=['assay']+['target']+[oligo+field for oligo in oligos for field in ['_name','_seq', '_loc_seq','_tm','_mismatches','_gaps','_errors']]
outputFields+=['total_errors']
outputFields+=['amplicon'+field for field in ['_coords','_length','_loc_seq']]
outputFileHeader='\t'.join(outputFields)
outputFile.write(outputFileHeader+'\n')
outputFile.close()

## Read assay file
with open(args.assay) as assayFile:
    assays=[assay.rstrip().split(',') for assay in assayFile.readlines()]

assayNames=[assay[0] for assay in assays]
if len(assayNames)==len(set(assayNames)) and len(assayNames)>0:
    for assay in assays:
        print() 
        if len(assay) in [5,7]:
            oligoNames=[assay[1],assay[2]] if len(assay)==5 else [assay[1],assay[2],assay[3]]
            if len(oligoNames)==len(set(oligoNames)):
                createAssayFile(assay,args)
                runTNTBLAST(assay,args)
                parseOutput(assay,args)
            else:
                print('ERROR: Assay file entry '+assay[0]+' does not have unique oligo names.')
        else:
            print('ERROR: Assay file entry '+assay[0]+' is not in a recognizable format.')
else:
    print('ERROR: Each assay file entry  must have a unique name.')

## WRITE REPORTS
print()
print('Writing reports...')
with open(args.genomes,'r') as genomeFile:
    genomes=genomeFile.readlines()
totalGenomes=[line[0]=='>' for line in genomes].count(True)

dataFile=pd.read_csv(args.output+'_PCR_results.tsv',sep='\t').replace('---',np.nan)

### Write assay report
assayReport=pd.DataFrame(index=assayNames,columns=range(dataFile.total_errors.max()+1)).reset_index()
assayReport.columns=['assay']+[i for i in range(dataFile.total_errors.max()+1)]
assayReport=pd.melt(assayReport,'assay',[i for i in range(dataFile.total_errors.max()+1)],'errors','count')[['assay','errors']]
errorCounts=pd.DataFrame()
for assay in assayNames:
    assayData=dataFile[dataFile.assay==assay]
    assayData=assayData.total_errors.value_counts().reset_index()
    assayData.columns=['errors','count']
    assayData['assay']=[assay]*len(assayData)
    errorCounts=pd.concat([errorCounts,assayData],sort=True)
assayReport=pd.merge(assayReport,errorCounts,on=['assay','errors'],how='outer')
assayReport=assayReport[(assayReport['errors']==0)|(assayReport['count']>0)]
assayReport['count']=assayReport['count'].replace(np.nan,0)
assayReport['targets_detected']=assayReport.apply(lambda row: len(dataFile[dataFile.assay==row.assay]['target'].unique()),axis=1)
assayReport['targets_total']=totalGenomes
assayReport['perc_detected']=assayReport.apply(lambda row: round(row['targets_detected']*100/row['targets_total'],2), axis=1)
assayReport['perc_of_detected']=assayReport.apply(lambda row: round(row['count']*100/row['targets_detected'],2) if row.targets_detected>0 else np.nan, axis=1)
assayReport['perc_of_total']=assayReport.apply(lambda row: round(row['count']*100/row['targets_total'],2), axis=1)
cols=['assay','targets_detected','targets_total','perc_detected','errors','count','perc_of_detected','perc_of_total']
assayReport=assayReport[cols].sort_values(by=['assay','errors'])
assayReport.to_csv(args.output+'_assay_report.tsv',sep='\t',index=False)

## Write oligo site variant report
variantReport=pd.DataFrame()
for assay in assayNames:
    assayData=dataFile[dataFile.assay==assay]
    for oligo in ['fwd_primer','rev_primer','probe']:
        variantData=assayData[oligo+'_loc_seq'].value_counts().reset_index()
        variantData.columns=['oligo_site_variant','count']
        variantData['assay']=assay
        variantData['oligo']=oligo
        variantData['oligo_name']=assayData[oligo+'_name'].value_counts().idxmax() if len(assayData[oligo+'_name'].value_counts())>0 else np.nan
        variantData['oligo_seq']=assayData[oligo+'_seq'].value_counts().idxmax() if len(assayData[oligo+'_seq'].value_counts())>0 else np.nan
        errors=assayData[[oligo+'_loc_seq',oligo+'_errors']].drop_duplicates()
        errors.columns=['oligo_site_variant','errors']
        variantData=pd.merge(variantData,errors,on='oligo_site_variant')
        variantReport=pd.concat([variantReport,variantData],sort=True)
variantReport=pd.merge(variantReport,assayReport[['assay','targets_total','targets_detected']].drop_duplicates(),on='assay')
variantReport['perc_detected']=variantReport.apply(lambda row: round(row['targets_detected']*100/row['targets_total'],2), axis=1)
variantReport['perc_of_detected']=variantReport.apply(lambda row: round(row['count']*100/row['targets_detected'],2), axis=1)
variantReport['perc_of_total']=variantReport.apply(lambda row: round(row['count']*100/row['targets_total'],2), axis=1)
cols=['assay','targets_detected','targets_total','perc_detected','oligo','oligo_name','oligo_seq','oligo_site_variant','errors','count','perc_of_detected','perc_of_total']
variantReport=variantReport[cols].sort_values(by=['assay','oligo','errors'])
variantReport=variantReport[(variantReport['perc_of_detected']>args.threshold)&(variantReport['errors']>0)]
variantReport.to_csv(args.output+'_variant_report.tsv',sep='\t',index=False)

## Write missed targets report
missedReport=pd.DataFrame()
for assay in assayNames:
    assayData=dataFile[dataFile.assay==assay]
    foundTargets=set(assayData.target.unique())
    missedTargets=allTargets-foundTargets
    assayData=pd.DataFrame()
    assayData['target']=list(missedTargets)
    assayData['assay']=assay
    missedReport=pd.concat([missedReport,assayData],sort=True)
missedReport[['assay','target']].to_csv(args.output+'_missed_report.tsv',sep='\t',index=False)

command='rm '+args.output + '_PCR_results.tsv'
subprocess.call(command,shell=True)

print()
print('Done.')
print()
