import argparse
import subprocess
import os
import pandas as pd

parser=argparse.ArgumentParser()
parser.add_argument("-a", "--assay", type=str, required=True, help="input TSV file containing assay name, primer and probe oligo sequences, and (optionally) reaction molarities")
parser.add_argument("-g", "--genomes", type=str, required=True, help="input FASTA file containing sequences assay is meant to detect")
parser.add_argument("-o", "--output", type=str, required=True, help="brief descriptor to append to output file names")
args=parser.parse_args()

print()
print('PCR_STRAINER - A tool for assessing inclusivity of oligonucleotides from Taqman PCR assays')
print()

## Prepare/overwrite existing output files for tabulated PCR results
outputFile=open(args.output+'_PCR_results.tsv','w')
recordedValues=['name','sequence header',
            'forward primer tm','reverse primer tm','probe tm',
            'forward primer mismatches','reverse primer mismatches','probe mismatches','total mismatches',
            'forward primer gaps','reverse primer gaps','probe gaps','total gaps',
            'forward primer errors','reverse primer errors','probe errors',
            'total errors','amplicon location','amplicon range',
            'forward primer location','reverse primer location','probe location',
            'concatenated oligo locations','degenerate bases']
outputFileHeader='\t'.join(recordedValues)
outputFile.write(outputFileHeader+'\n')
outputFile.close()

## Open user-provided file containing assay details and reformat as input TSV for TNTBLAST
def createAssayFile(args,assayDetailsList):
    print('Creating',assayDetailsList[0],'assay file...')
    assayFile=open(args.output+'_'+assayDetailsList[0]+'_TNTBLAST_assay.tsv','w') ## Create TSV file containing each separate assay's oligos for use by TNTBLAST
    assayFile.write(assayDetailsList[0]+'\t'+assayDetailsList[1]+'\t'+assayDetailsList[2]+'\t'+assayDetailsList[3]+'\n') ## Write assay name, for primer seq, rev primer seq, probe seq
    assayFile.close()

## Run TNTBLAST from user-provided assay details
def runTNTBLAST(args,assayDetailsList):
    print('Running TNTBLAST with',assayDetailsList[0],'assay against',args.genomes,'...')
    targets=args.genomes ## direct command at genomes
    if len(assayDetailsList)>4:
        primerMol=assayDetailsList[4] ## Set primer molarity
        probeMol=assayDetailsList[5] ## Set probe molarity
    else:
        primerMol=9.0*(10**-7) ## Use TNTBLAST default if primer molarity not provided
        probeMol=primerMol ## Use TNTBLAST default if probe molarity not provided
    command='tntblast -i '+ args.output+'_'+assayDetailsList[0] +'_TNTBLAST_assay.tsv' +' -d '+ targets ## Create command for running TNTBLAST
    command+=' -o '+ args.output + '_' + assayDetailsList[0] +'_TNTBLAST_output.txt ' ## Create command for running TNTBLAST (cont)
    command+='-e 30 -E 30 --best-match -m 0 -v F' ## Create command for running TNTBLAST (cont)
    #print('TNTBLAST command:',command)
    subprocess.call(command,shell=True) ## Run TNTBLAST command

## Parse TNTBLAST output into TSV file
def parseTNTBLASTOutput(args,assayDetailsList,recordedValues): ## Function to parse output file from TNTBLAST then tabulated and write relevant values into TSV
    print('Parsing TNTBLAST output from',assayDetailsList[0],'assay...')
    def updateAmplicon(amplicon,currentLine):
        if ' = ' in currentLine:
            key=currentLine.split(' = ')[0]
            value=currentLine.split(' = ')[1].rstrip()
            amplicon[key]=value
        elif currentLine[0]=='>':
            amplicon['sequence header']=currentLine.rstrip()
        elif 'sequence header' in amplicon.keys() and 'amplicon location' not in amplicon.keys():
            amplicon['total mismatches']=sum([int(amplicon[oligo+' mismatches']) for oligo in ['forward primer','reverse primer','probe']])
            amplicon['total gaps']=sum([int(amplicon[oligo+' gaps']) for oligo in ['forward primer','reverse primer','probe']])
            for oligo in ['forward primer','reverse primer','probe']:
                amplicon[oligo+' errors']=int(amplicon[oligo+' gaps'])+int(amplicon[oligo+' mismatches'])
            amplicon['total errors']=amplicon['total mismatches']+amplicon['total gaps']
            amplicon['amplicon location']=currentLine.rstrip()
            amplicon['forward primer location']=amplicon['amplicon location'][0:len(amplicon['forward primer'].lstrip("5'-").rstrip("-3'"))]
            amplicon['reverse primer location']=amplicon['amplicon location'][-len(amplicon['reverse primer'].lstrip("5'-").rstrip("-3'")):]
            probeStart=int(amplicon['probe range'].split(' .. ')[0])-int(amplicon['amplicon range'].split(' .. ')[0])
            probeEnd=int(amplicon['probe range'].split(' .. ')[1])-int(amplicon['amplicon range'].split(' .. ')[0])
            amplicon['probe location']=amplicon['amplicon location'][probeStart:probeEnd+1]
            amplicon['concatenated oligo locations']=amplicon['forward primer location']+'-'+amplicon['probe location']+'-'+amplicon['reverse primer location']
            totalCanonicalBases=sum([amplicon['concatenated oligo locations'].count(symbol) for symbol in ['A','T','G','C','-']])
            totalBases=len(amplicon['concatenated oligo locations'])
            amplicon['degenerate bases']=str(totalCanonicalBases!=totalBases)
            complementaryBases={'A':'T','T':'A','G':'C','C':'G','W':'W','S':'S','M':'K','K':'M','R':'Y','Y':'R','B':'V','V':'B','D':'H','H':'D','N':'N'}
            amplicon['reverse primer location']=''.join([complementaryBases[base] for base in amplicon['reverse primer location'][::-1]])
        return amplicon
    inputFile=open(args.output + '_' + assayDetailsList[0] +'_TNTBLAST_output.txt','r')
    outputFile=open(args.output+'_PCR_results.tsv','a')
    currentAmplicon={}
    currentLine=inputFile.readline()
    while currentLine!='':
        if 'name = ' in currentLine:
            if currentAmplicon!={}:
                outputFile.write('\t'.join([str(currentAmplicon[value]) for value in recordedValues])+'\n')
            currentAmplicon={}
            currentAmplicon=updateAmplicon(currentAmplicon,currentLine)
        else:
            currentAmplicon=updateAmplicon(currentAmplicon,currentLine)
        currentLine=inputFile.readline()
    if currentAmplicon!={}:
        outputFile.write('\t'.join([str(currentAmplicon[value]) for value in recordedValues])+'\n')
    inputFile.close()
    outputFile.close()

inputFile=open(args.assay,'r') ## Open user-provided file containing assay details
currentLine=inputFile.readline() ## Get first assay from user-provided assay file
while currentLine!='':
    assayDetailsList=currentLine.split(',')
    createAssayFile(args,assayDetailsList)
    runTNTBLAST(args,assayDetailsList)
    parseTNTBLASTOutput(args,assayDetailsList,recordedValues) ## Recorded values defined above; the data fields to tabluate into the TSV file
    #os.remove(args.output + '_' + assayDetailsList[0] +'_TNTBLAST_output.txt')
    os.remove(args.output+'_'+assayDetailsList[0]+'_TNTBLAST_assay.tsv')
    print()
    currentLine=inputFile.readline() ## Get next assay from user-provided assay file

### Get names of all assays
inputFile=open(args.assay,'r')
assayNames=set()
currentLine=inputFile.readline()
while currentLine!='':
    assayNames.add(currentLine.split(',')[0])
    currentLine=inputFile.readline()
inputFile.close()
assayNames=sorted(assayNames)

### Get headers of all inclusion target genome sequences and store a separate copy for each assay
inputFile=open(args.genomes,'r')
genomes=set()
currentLine=inputFile.readline()
while currentLine!='':
    if currentLine[0]=='>':
        genomes.add(currentLine.rstrip())
    currentLine=inputFile.readline()
inputFile.close()
genomes={assayName:genomes for assayName in assayNames}

### Load TSV data for PCR_strainer
PCRData=pd.read_csv(args.output+'_PCR_results.tsv', sep='\t')

### Find target genomes that failed to generate amplicons and write headers to TXT file
### Store missed targets in dictionary where keys are assay names and values are the targets for which that assay failed to generate amplicons
missedGenomes={}
for assayName in assayNames:
    missedGenomes[assayName]=genomes[assayName]-set(PCRData[PCRData['name']==assayName]['sequence header'].unique())
    if len(missedGenomes[assayName])>0:
        outputFile=open(args.output+'_'+assayName+'_missed_genomes.txt','w')
        for genome in missedGenomes[assayName]:
            outputFile.write(genome+'\n')
        outputFile.close()

### Find target genomes that with degenerate bases in locations targeted by assays and write headers to TXT file
### Store degenerate targets in dictionary where keys are assay names and values are targets for which there are degenerate bases in that assay's target locations
degenerateGenomes={}
for assayName in assayNames:
    degenerateGenomes[assayName]=set(PCRData[(PCRData['name']==assayName)&(PCRData['degenerate bases']==True)]['sequence header'].unique())
    if len(degenerateGenomes[assayName])>0:
        outputFile=open(args.output+'_'+assayName+'_degenerate_genomes.txt','w')
        for genome in degenerateGenomes[assayName]:
            outputFile.write(genome+'\n')
        outputFile.close()

### Write out TSV containing frequency of gaps/mismatches in combined asasy oligo target locations
combinedOligoErrors=pd.DataFrame()
for assay in assayNames:
    totalGenomes=len(genomes[assay])-len(degenerateGenomes[assayName])
    assayData=PCRData[(PCRData['name']==assay)&(PCRData['degenerate bases']==False)]
    errorFrequencies=assayData['total errors'].value_counts()
    errorFrequencies=errorFrequencies.reset_index()
    errorFrequencies.columns=['gaps/mismatches','count (#)']
    errorFrequencies['assay']=[assay]*len(errorFrequencies)
    errorFrequencies['frequency (%)']=errorFrequencies.apply(lambda row:round(row['count (#)']/totalGenomes*100,1),axis=1)
    combinedOligoErrors=pd.concat([combinedOligoErrors,errorFrequencies])
combinedOligoErrors=combinedOligoErrors[['assay','gaps/mismatches','count (#)','frequency (%)']].sort_values(by=['assay','gaps/mismatches'])
combinedOligoErrors=combinedOligoErrors.reset_index(drop=True)
combinedOligoErrors.to_csv(args.output+'_combined_oligo_errors.tsv',sep='\t')

### Write out TSV containing frequency of gaps/mismatches in separate asasy oligo target locations
separatedOligoErrors=pd.DataFrame()
for assay in assayNames:
    totalGenomes=len(genomes[assay])-len(degenerateGenomes[assayName])
    assayData=PCRData[(PCRData['name']==assay)&(PCRData['degenerate bases']==False)]
    for oligo in ['forward primer','reverse primer','probe']:
        errorFrequencies=assayData[oligo+' errors'].value_counts()
        errorFrequencies=errorFrequencies.reset_index()
        errorFrequencies.columns=['gaps/mismatches','count (#)']
        errorFrequencies['assay']=[assay]*len(errorFrequencies)
        errorFrequencies['oligo']=[oligo]*len(errorFrequencies)
        errorFrequencies['frequency (%)']=errorFrequencies.apply(lambda row:round(row['count (#)']/totalGenomes*100,1),axis=1)
        separatedOligoErrors=pd.concat([separatedOligoErrors,errorFrequencies])
separatedOligoErrors=separatedOligoErrors[['assay','oligo','gaps/mismatches','count (#)','frequency (%)']].sort_values(by=['assay','oligo','gaps/mismatches'])
separatedOligoErrors=separatedOligoErrors.reset_index(drop=True)
separatedOligoErrors.to_csv(args.output+'_separated_oligo_errors.tsv',sep='\t')

### Write out TSV containing frequency of combined assay oligo location variants
combinedOligoVariants=pd.DataFrame()
for assay in assayNames:
    totalGenomes=len(genomes[assay])-len(degenerateGenomes[assayName])
    assayData=PCRData[(PCRData['name']==assay)&(PCRData['degenerate bases']==False)]
    variantFrequencies=assayData['concatenated oligo locations'].value_counts()
    variantFrequencies=variantFrequencies.reset_index()
    variantFrequencies.columns=['variant','count (#)']
    variantFrequencies['assay']=[assay]*len(variantFrequencies)
    variantFrequencies['frequency (%)']=variantFrequencies.apply(lambda row:round(row['count (#)']/totalGenomes*100,1),axis=1)
    variantFrequencies['gaps/mismatches']=variantFrequencies.apply(lambda row: assayData[assayData['concatenated oligo locations']==row['variant']]['total errors'].mode(),axis=1)
    combinedOligoVariants=pd.concat([combinedOligoVariants,variantFrequencies])
combinedOligoVariants=combinedOligoVariants[['assay','variant','count (#)','frequency (%)','gaps/mismatches']].sort_values(by=['assay','count (#)'],ascending=[True,False])
combinedOligoVariants=combinedOligoVariants.reset_index(drop=True)
combinedOligoVariants.to_csv(args.output+'_combined_oligo_variants.tsv',sep='\t')

### Write out TSV containing frequency of separated assay oligo location variants
separatedOligoVariants=pd.DataFrame()
for assay in assayNames:
    totalGenomes=len(genomes[assay])-len(degenerateGenomes[assayName])
    assayData=PCRData[(PCRData['name']==assay)&(PCRData['degenerate bases']==False)]
    for oligo in ['forward primer','reverse primer','probe']:
        variantFrequencies=assayData[oligo+' location'].value_counts()
        variantFrequencies=variantFrequencies.reset_index()
        variantFrequencies.columns=['variant','count (#)']
        variantFrequencies['assay']=[assay]*len(variantFrequencies)
        variantFrequencies['oligo']=[oligo]*len(variantFrequencies)
        variantFrequencies['frequency (%)']=variantFrequencies.apply(lambda row:round(row['count (#)']/totalGenomes*100,1),axis=1)
        variantFrequencies['gaps/mismatches']=variantFrequencies.apply(lambda row: assayData[assayData[oligo+' location']==row['variant']][oligo+' errors'].mode(),axis=1)
        variantFrequencies['predicted Tm (C)']=variantFrequencies.apply(lambda row: round(assayData[assayData[oligo+' location']==row['variant']][oligo+' tm'].mode(),1),axis=1)
        separatedOligoVariants=pd.concat([separatedOligoVariants,variantFrequencies])
separatedOligoVariants=separatedOligoVariants[['assay','oligo','variant','count (#)','frequency (%)','gaps/mismatches','predicted Tm (C)']].sort_values(by=['assay','oligo','count (#)'],ascending=[True,True,False])
separatedOligoVariants=separatedOligoVariants.reset_index(drop=True)
separatedOligoVariants.to_csv(args.output+'_separated_oligo_variants.tsv',sep='\t')

### Write out FASTA files containing genome sequence variants for each assay's target locations
for assay in assayNames:
    assayData=PCRData[PCRData['name']==assay]
    for oligo in ['forward primer','reverse primer','probe']:
        outputFile=open(args.output+'_'+assay+'_'+oligo.replace(' ','_')+'_location_seqs.fa','w')
        for header,variant in zip(assayData['sequence header'],assayData[oligo+' location']):
            outputFile.write(header+'\n')
            outputFile.write(variant+'\n')
        outputFile.close()
        
