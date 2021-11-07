# PCR_strainer - A tool for assessing inclusivity of PCR primers and probes against reference nucleotide sequences
#
# Kevin Kuchinski
# British Columbia Centre for Disease Control, Public Health Laboratory
# University of British Columbia, Department of Pathology and Laboratory Medicine
# kevin.kuchinski@bccdc.ca


import argparse
import os
import subprocess
import numpy as np
import pandas as pd


def main():
    """Main function called when script is run."""
    version = '0.1.0'
    print()
    print('PCR_strainer v' + version)
    print('https://github.com/KevinKuchinski/PCR_strainer')
    print()
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--assay', type = str, required = True, help = '')
    parser.add_argument('-g', '--genomes', type = str, required = True, help = '')
    parser.add_argument('-o', '--output', type = str, required = True, help = '')
    parser.add_argument('-t', '--threshold', type = float, required = True, help = '')
    parser.add_argument('-m', '--melting', type = float, required = False, default=45, help = '')
    args=parser.parse_args()
    assays = read_assay_file(args.assay)
    job_name, path_to_output = get_job_name_and_output_path(args.output)
    final_tntblast_results = pd.DataFrame()
    for assay_details in assays:
        run_TNTBLAST(assay_details, args.genomes, path_to_output, args.melting)
        tntblast_results = parse_tntblast_output(assay_details, job_name, path_to_output)
        final_tntblast_results = pd.concat([final_tntblast_results, tntblast_results], sort=True)
    write_assay_report(job_name, path_to_output, final_tntblast_results, args.genomes, args.threshold)
    write_variant_report(job_name, path_to_output, final_tntblast_results, args.genomes, args.threshold)
    write_missed_seqs_report(job_name, path_to_output, tntblast_results, args.genomes)
    write_tntblast_results(job_name, path_to_output, tntblast_results)


def read_assay_file(path_to_file):
    """Read input assay file and parse each list into a tuple:
    (assay_name,fwd_oligo_name,fwd_oligo_seq,rev_oligo_name,rev_oligo_seq,probe_oligo_name,probe_oligo_seq)
    Return a tuple of these tuples."""
    assays = []
    with open(path_to_file, 'r') as input_file:
        for line in input_file:
            line = line.rstrip().split(',')
            # Check if line has the accepted number of comma-separated values
            if len(line) not in [5, 7]:
                print('ERROR: Line in assay file not properly formatted:')
                print(','.join(line))
                print()
                print('Assay file lines must be a comma-separated list:')
                print('assay_name,fwd_primer_name,fwd_primer_seq,rev_primer_name,rev_primer_seq,probe_name,probe_seq') 
                exit()
            else:
                assay_name = line[0]
                fwd_primer_name = line[1]
                fwd_primer_seq = line[2]
                rev_primer_name = line[3]
                rev_primer_seq = line[4]
                if len(line) == 7:
                    probe_name = line[5]
                    probe_seq = line[6]
                else:
                    probe_name = ''
                    probe_seq = ''
                # Check if any assay or primer names are empty
                names = [assay_name, fwd_primer_name, rev_primer_name]
                names += [probe_name] if probe_seq != '' else []
                if any([name == '' for name in names]): 
                    print('ERROR: Assay and oligo names cannot be empty:')
                    print(','.join(line))
                    print()
                    exit()
                # Create assay info tuple and append to assays list
                assay = (assay_name, fwd_primer_name, fwd_primer_seq, rev_primer_name, rev_primer_seq, probe_name, probe_seq)
                assays.append(assay)
    # Check that names used for assays and oligos are unique
    for name_index in [0, 1, 3, 5]:
        names = [assay[name_index] for assay in assays if assay[name_index] != '']
        if len(names) != len(set(names)):
            print('ERROR: Assay and oligo names must be unique!')
            print()
            exit()
    return tuple(assays)


def get_job_name_and_output_path(output_arg):
    """Get job name and path for output files."""
    output_arg = output_arg.lstrip().rstrip().rstrip('/')
    if output_arg == '':
        print('ERROR: Must define path to output and job name!')
        exit()
    else:
        job_name = output_arg.split('/')[-1]
        path_to_output = '/'.join(output_arg.split('/')[:-1])
        path_to_output = '.' if path_to_output == '' else path_to_output
        path_to_output += '/'
        return job_name, path_to_output


def run_TNTBLAST(assay_details, path_to_genomes, path_to_output, melting_temp):
    """Run TNTBLAST with provided assay details."""
    # Create TNTBLAST assay file from assay details (TSV file of assay name and oligo seqs)
    assay_name = assay_details[0]
    fwd_primer_seq = assay_details[2]
    rev_primer_seq = assay_details[4]
    probe_seq = assay_details[6]
    print('Running TNTBLAST for', assay_name, 'against', path_to_genomes, '...')
    assay = [assay_name, fwd_primer_seq, rev_primer_seq]
    assay += [probe_seq] if probe_seq != '' else []
    assay = '\t'.join(assay)
    path_to_tntblast_assay = path_to_output + assay_name + '_tntblast_assay.tsv'
    with open(path_to_tntblast_assay, 'w') as output_file:
        output_file.write(assay + '\n')    
    # Create terminal command for TNTBLAST and run
    path_to_tntblast_txt_output = path_to_output + assay_name + '_tntblast_output.txt'
    terminal_command = 'tntblast -i' + path_to_tntblast_assay + ' -d ' + path_to_genomes
    terminal_command += ' -o ' + path_to_tntblast_txt_output + ' -e ' + str(melting_temp) + ' -E ' + str(melting_temp)
    terminal_command += ' --best-match -m 0 -v F'
    subprocess.run(terminal_command, shell=True)
    # Garbage collection (remove TNTBLAST assay TSV file)
    os.remove(path_to_tntblast_assay)


def parse_tntblast_output(assay_details, job_name, path_to_output):
    """Parses txt format output from TNTBLAST and tabulates relevant data into Pandas dataframes.
    Return the dataframe."""
    # Values for getting complement of base
    comp_bases = {'A': 'T', 'T': 'A',
                  'G': 'C', 'C': 'G',
                  'W': 'W', 'S': 'S',
                  'M': 'K', 'K': 'M',
                  'R': 'Y', 'Y': 'R',
                  'B': 'V', 'V': 'B',
                  'D': 'H', 'H': 'D',
                  'N': 'N', '-': '-'}
    # Function for getting reverse complement of sequence
    rev_comp = lambda seq: ''.join([comp_bases[base] for base in seq[::-1]])
    # Create list of relevant fields in TNTBLAST output
    oligo_names = ['forward primer', 'reverse primer', 'probe']
    fields = ['', ' tm', ' mismatches', ' gaps']
    input_fields = [oligo + field for oligo in oligo_names for field in fields]
    input_fields += ['amplicon' + field for field in [' range', ' length']] + ['probe range']
    # Create list of fields to store in tntblast results
    oligo_names = ['forward primer', 'reverse primer', 'probe']
    fields = [' name', ' seq', ' mismatches', ' gaps', ' errors', ' site seq']
    output_fields = ['assay name', 'target', 'total errors']
    output_fields += [oligo + field for oligo in oligo_names for field in fields]
    output_fields += ['amplicon seq']
    # Create dict of lists for storing results from parsing
    tntblast_results = {field: [] for field in output_fields}
    # Open TNTBLAST output and parse into dataframe
    assay_name = assay_details[0]
    print('Parsing TNTBLAST output from', assay_name, '...')
    path_to_tntblast_txt_input = path_to_output + assay_name + '_tntblast_output.txt'
    with open(path_to_tntblast_txt_input, 'r') as input_file:
        # Write header to output file
        input_lines = input_file.readlines()
        for line_index, line in enumerate(input_lines):
            if ' = ' in line: 
                key, value = line.rstrip().split(' = ')
                if key == 'name': # Initialize new entry if key is 'name'
                    entry = {}
                    entry['assay name'] = assay_name
                    for key in input_fields:
                        entry[key] = 'none'
                elif key in input_fields:
                    entry[key] = value
            elif line[0] == '>':
                # If line starts with '>', it signals 3 things:
                # 1. the name of the target is on this line
                # 2. the seq of the amplicon is on the following line
                # 3. the end of the current entry, so it can be added to the overall results
                entry['target'] = line.lstrip('>').rstrip()
                entry['amplicon seq'] = input_lines[line_index+1].rstrip()
                # Count number of mismatches and gaps for each oligo and add together for errors
                oligos = ['forward primer', 'reverse primer', 'probe']
                oligo_names = [assay_details[i] for i in [1, 3, 5]]
                oligo_seqs = [assay_details[i] for i in [2, 4, 6]]
                for oligo, oligo_name, oligo_seq in zip(oligos, oligo_names, oligo_seqs):
                    if entry[oligo] != 'none':
                        if entry[oligo + ' mismatches'] == 'none':
                            entry[oligo + ' mismatches'] = 0
                        if entry[oligo + ' gaps'] == 'none':
                            entry[oligo+' gaps'] = 0
                        mismatches = int(entry[oligo + ' mismatches'])
                        gaps = int(entry[oligo + ' gaps'])
                        entry[oligo + ' errors'] =  mismatches + gaps
                    # Add name and sequence of oligo to entry
                    entry[oligo + ' name'] = oligo_name
                    entry[oligo + ' seq'] = oligo_seq  
                # Extract oligo site sequences from amplicon seq
                entry['forward primer site seq'] = entry['amplicon seq'][:len(entry['forward primer seq'])]
                entry['reverse primer site seq'] = rev_comp(entry['amplicon seq'][-len(entry['reverse primer seq']):])
                if entry['probe range'] != 'none':
                    start = int(entry['probe range'].split(' .. ')[0]) - int(entry['amplicon range'].split(' .. ')[0])
                    end = int(entry['probe range'].split(' .. ')[1]) - int(entry['amplicon range'].split(' .. ')[1])
                    entry['probe site seq'] = entry['amplicon seq'][start:end]
                entry['total errors'] = sum([entry[oligo + ' errors'] for oligo in ['forward primer', 'reverse primer', 'probe']])
                # Add values from this entry to overall results
                for field in output_fields:
                    tntblast_results[field].append(entry[field])
    # Replace spaces in column names with underscores
    tntblast_results = {field.replace(' ', '_'): tntblast_results[field] for field in output_fields}
    # Replace full spelling of primers with shortened versions (forward to fwd, and reverse to rev)
    for oligo, short_oligo in zip(['forward_primer', 'reverse_primer'], ['fwd_primer', 'rev_primer']):
        tntblast_results = {field.replace(oligo, short_oligo): tntblast_results[field] for field in tntblast_results.keys()}
    # Create dataframe from dictionary of parsed results
    tntblast_results = pd.DataFrame.from_dict(tntblast_results)
    # Garbage collection
    os.remove(path_to_tntblast_txt_input)
    return tntblast_results


def write_tntblast_results(job_name, path_to_output, tntblast_results):
    """Writes tabulated results in Pandas dataframe to a TSV file."""
    # Re-order columns and save TNTBLAST results to TSV
    oligo_names = ['fwd_primer', 'rev_primer', 'probe']
    fields = ['_name', '_seq', '_mismatches', '_gaps', '_errors', '_site_seq']
    cols = ['assay_name', 'target', 'total_errors']
    cols += [oligo + field for oligo in oligo_names for field in fields]
    cols += ['amplicon_seq']
    path_to_tntblast_tsv = path_to_output + job_name + '_PCR_results.tsv'
    tntblast_results[cols].to_csv(path_to_tntblast_tsv, sep='\t', index=False)


def count_targets(path_to_genomes):
    """Counts the number of target sequences in the provided genomes FASTA file.
    Returns an int of the count."""
    with open(path_to_genomes, 'r') as input_file:
        seq_counter = 0
        for line in input_file:
            if line[0] == '>':
                seq_counter += 1
    return int(seq_counter)


def write_assay_report(job_name, path_to_output, tntblast_results, path_to_genomes, threshold):
    """Takes TNTBLAST results in Pandas dataframe and writes the assay report, which describes the number of
    target sequences with 0, 1, 2, 3... total errors for each assay. Writes report to TSV file."""
    print('Write assay report...')
    # Count number of targets with each level of total errors for all assays
    cols = ['assay_name', 'total_errors']
    assay_report = tntblast_results[cols].groupby(cols).size().reset_index()
    assay_report.columns = ['assay_name', 'total_errors', 'target_count']
    # Count targets detected by each assay
    targets_detected = tntblast_results['assay_name'].value_counts().reset_index()
    targets_detected.columns = ['assay_name', 'detected_targets']
    assay_report = pd.merge(assay_report, targets_detected, on='assay_name')
    # Count total targets
    assay_report['total_targets'] = count_targets(path_to_genomes)
    # Calculate percentage of detected targets representing each level of total errors
    assay_report['perc_of_detected'] = round(assay_report['target_count'] * 100 / assay_report['detected_targets'],2)
    # Calculate percentage of total targets representing each level of total errors
    assay_report['perc_of_total'] = round(assay_report['target_count'] * 100 / assay_report['total_targets'],2)
    # Drop report rows with percent of detected below threshold
    assay_report = assay_report[assay_report['perc_of_detected']>=float(threshold)]
    # Write assay report to TSV file
    path_to_report_file = path_to_output + job_name + '_assay_report.tsv'
    assay_report = assay_report.sort_values(['assay_name', 'target_count'], ascending=[True, False])
    assay_report.to_csv(path_to_report_file, sep='\t', index=False)


def write_variant_report(job_name, path_to_output, tntblast_results, path_to_genomes, threshold):
    """Takes TNTBLAST results in Pandas dataframe and writes the variant report, which describes
    oligo site variants whose prevalence exceeded the provided threshold. Writes report to TSV file."""
    print('Write variant report...')
    # Create variant report
    variant_report = pd.DataFrame()
    for oligo in ['fwd_primer', 'rev_primer', 'probe']:
        cols = ['assay_name', oligo + '_name', oligo + '_seq', oligo + '_site_seq', oligo + '_errors']
        oligo_data = tntblast_results[cols]
        cols = ['assay_name', 'oligo_name', 'oligo_seq', 'oligo_site_seq', 'oligo_errors']
        oligo_data.columns = cols
        oligo_data = oligo_data[cols].groupby(cols).size().reset_index()
        oligo_data.columns = cols + ['target_count']
        oligo_data['oligo'] = oligo
        variant_report = pd.concat([variant_report, oligo_data], sort=True)
    # Add lower case bases to oligo site seqs
    def compare_seqs(row):
        a, b = row['oligo_site_seq'], row['oligo_seq']
        seq = ''.join([site_base if site_base == oligo_base else site_base.lower() for site_base, oligo_base in zip(a, b)])
        return seq
    variant_report['oligo_site_seq'] = variant_report.apply(compare_seqs, axis=1)
    # Count targets detected by each assay
    targets_detected = tntblast_results['assay_name'].value_counts().reset_index()
    targets_detected.columns = ['assay_name', 'detected_targets']
    variant_report = pd.merge(variant_report, targets_detected, on='assay_name')
    # Count total targets
    variant_report['total_targets'] = count_targets(path_to_genomes)
    # Calculate percentage of detected targets representing each level of total errors
    variant_report['perc_of_detected'] = round(variant_report['target_count'] * 100 / variant_report['detected_targets'],2)
    # Calculate percentage of total targets representing each level of total errors
    variant_report['perc_of_total'] = round(variant_report['target_count'] * 100 / variant_report['total_targets'],2)
    # Re-order columns and sort rows
    cols = ['assay_name', 'oligo', 'oligo_seq', 'oligo_site_seq', 'oligo_errors', 'target_count', 'detected_targets', 'total_targets']
    cols += ['perc_of_detected', 'perc_of_total']
    variant_report = variant_report[cols]
    # Drop rows for variants with no errors
    variant_report = variant_report[variant_report['oligo_errors']!=0]
    # Drop report rows with percent of detected below threshold
    variant_report = variant_report[variant_report['perc_of_detected']>=float(threshold)]
    # Write assay report to TSV file
    path_to_report_file = path_to_output + job_name + '_variant_report.tsv'
    variant_report = variant_report.sort_values(['assay_name', 'oligo', 'target_count'], ascending=[True, True, False])
    variant_report.to_csv(path_to_report_file, sep='\t', index=False)


def get_targets(path_to_genomes):
    """Reads all sequence in the provided genomes FASTA file and returns them in
    a dict where keys are FASTA headers and values are the nucleotide sequences."""
    with open(path_to_genomes, 'r') as input_file:
        seqs = {}
        for line in input_file:
            if line[0] == '>':
                header = line.lstrip('>').rstrip()
                seqs[header] = ''
            else:
                seqs[header] += line.rstrip()
    return seqs


def write_missed_seqs_report(job_name, path_to_output, tntblast_results, path_to_genomes):
    """Takes TNTBLAST results in Pandas dataframe and the provided genomes FASTA file and writes
    a report describing sequences in the provided genome targets that were not detected by TNTBLAST.
    Writes report to TSV file."""
    missed_seqs_report = pd.DataFrame()
    seqs = get_targets(path_to_genomes)
    all_seqs = set(seqs.keys())
    for assay in tntblast_results['assay_name'].unique():
        found_seqs = set(tntblast_results[tntblast_results['assay_name']==assay]['target'].unique())
        missed_seqs = tuple(all_seqs - found_seqs)
        assay_data = pd.DataFrame()
        assay_data['target'] = missed_seqs
        assay_data['assay_name'] = assay
        missed_seqs_report = pd.concat([missed_seqs_report, assay_data], sort=True)
    missed_seqs_report['target_length'] = missed_seqs_report.apply(lambda row: len(seqs[row['target']]), axis=1)
    get_perc_Ns = lambda row: round(seqs[row['target']].count('N') * 100 / row['target_length'], 2)
    missed_seqs_report['perc_Ns'] = missed_seqs_report.apply(get_perc_Ns, axis=1)
    path_to_report_file = path_to_output + job_name + '_missed_seqs_report.tsv'
    missed_seqs_report.to_csv(path_to_report_file, sep='\t', index=False)


if __name__ == '__main__':
    report = main()
