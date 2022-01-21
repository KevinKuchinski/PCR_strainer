#!/usr/bin/env python3

# PCR_strainer - A tool for assessing inclusivity of PCR primers and probes against reference nucleotide sequences
#
# Kevin Kuchinski
# British Columbia Centre for Disease Control, Public Health Laboratory
# University of British Columbia, Department of Pathology and Laboratory Medicine
# kevin.kuchinski@bccdc.ca

import sys
import subprocess
import os
import numpy as np
import pandas as pd


def main():
    version = '0.2.0'
    # Parse command line arguments
    args = parse_args(sys.argv, version)
    print(f'\nPCR_strainer v{version}')
    print('https://github.com/KevinKuchinski/PCR_strainer\n')
    out_path, name = os.path.split(args['-o'])
    out_path = '.' if out_path == '' else out_path
    # Make sure output directory exists
    if os.path.exists(out_path) == False or os.path.isdir(out_path) == False:
        print(f'\nERROR: Output path {out_path} does not exist.\n')
        exit(1)
    # Check input genomes file
    check_genomes_file(args['-g'])
    # Parse assay file to get oligo names and seqs
    assays = read_assay_file(args['-a'])
    # Create empty dataframe for tabulated tntblast results
    final_tntblast_results = pd.DataFrame()
    # Get tntblast results for each assay
    for assay_details in assays:
        run_TNTBLAST(assay_details, args['-g'], out_path, args['-m'], args['-p'], args['-P'])
        tntblast_results = parse_tntblast_output(assay_details, name, out_path)
        final_tntblast_results = pd.concat([final_tntblast_results, tntblast_results], sort=True)
    write_assay_report(name, out_path, final_tntblast_results, args['-g'], args['-t'])
    write_variant_report(name, out_path, final_tntblast_results, args['-g'], args['-t'])
    write_missed_seqs_report(name, out_path, final_tntblast_results, args['-g'])
    write_tntblast_results(name, out_path, final_tntblast_results)
    print('\nDone.\n')


def parse_args(args, version):
    arg_values = {}
    for arg_1, arg_2 in zip(args[1:-1], args[2:]):
        if arg_1[0] == '-':
            if arg_2[0] != '-':
                arg_values[arg_1] = arg_2
            else:
                arg_values[arg_1] = ''
    if args[-1][0] == '-':
        arg_values[args[-1]] = ''
    # Set defaults, mins, and maxs
    required_args = {'-a', '-g', '-o'}
    arg_value_types = {'-a': str, '-g':str, '-o': str, '-t': float, '-m': float, '-p': float, '-P': float}
    min_arg_values = {'-t': 0, '-p': 0, '-P': 0}
    max_arg_values = {'-t': 100}
    default_arg_values = {'-t': 0, '-m': 45, '-p': 1, '-P': 1}
    # Check if all required arguments were provided
    missing_args = set()
    for required_arg in required_args:
        if required_arg not in arg_values.keys() or arg_values[required_arg] == '':
            missing_args = missing_args | {required_arg}
    if missing_args != set():
        print(f'\nERROR: Values must be provided for the argument following arguments: {", ".join(sorted(missing_args))}')
        print_usage(version)
        exit(1)
    # Check if unrecognized arguments were provided
    recognized_args = required_args | set(arg_value_types.keys())
    unrecognized_args = set()
    for provided_arg in arg_values.keys():
        if provided_arg not in recognized_args:
            unrecognized_args = unrecognized_args | {provided_arg}
    if unrecognized_args != set():
        print(f'\nERROR: The following arguments are not recognized: {", ".join(sorted(unrecognized_args))}')
        print_usage(version)
        exit(1)
    # Check if arguments were provided without values
    empty_args = set()
    for arg, value in arg_values.items():
        if value == '':
            empty_args = empty_args | {arg}
    if empty_args != set():
        print(f'\nERROR: The following arguments were provided without values: {", ".join(sorted(empty_args))}')
        print_usage(version)
        exit(1)
    # Check if provided values are of the correct type
    for arg, value in arg_values.items():
        try:
            arg_values[arg] = arg_value_types[arg](value)
        except ValueError:
            print(f'\nERROR: Value for argument {arg} must be of type {str(arg_value_types[arg].__name__)}')
            print_usage(version)
            exit(1)
    # Check if provided values are within the correct range
    for arg, value in arg_values.items():
        if arg in min_arg_values.keys() and value <= min_arg_values[arg]:
            print(f'\nERROR: Value for argument {arg} must be greater than {min_arg_values[arg]}')
            print_usage(version)
            exit(1)
        if arg in max_arg_values.keys() and value >= max_arg_values[arg]:
            print(f'\nERROR: Value for argument {arg} must be less than {max_arg_values[arg]}')
            print_usage(version)
            exit(1)        
    # Assign default values to unspecified arguments
    for arg, value in default_arg_values.items():
        if arg not in arg_values.keys():
            arg_values[arg] = value
    # Return keyword args and their values
    return arg_values


def print_usage(version):
    print(f'\nPCR_strainer v{version}')
    print('https://github.com/KevinKuchinski/PCR_strainer\n')
    print('Usage: pcr_strainer -a <assay CSV> -g <genomes FASTA> -o <path to output directory>/<output name> [<optional args>]\n')
    print('Required arguments:')
    print(' -a : path to assay details in CSV file')
    print(' -g : path to target genomes in FASTA file')
    print(' -o : path to output directory and name to append to output files')
    print('Optional arguments:')
    print(' -t : minimum prevalence (%) of primer site variants reported in reports (default=5, min=0, max=100)')
    print(' -m : minimum Tm (degrees C) for primers and probes (default=45, min=0, max=100)')
    print(' -p : molar concentration of primer oligos (uM) (default=1, min=0)')
    print(' -P : molar concentration of probe oligos (uM) (default=1, min=0)\n')


def check_genomes_file(path_to_file):
    '''Checks that file exists, is not empty, and headers are unique and hashable.'''
    if os.path.exists(path_to_file) == False or os.path.isfile(path_to_file) == False:
        print(f'\nERROR: Genomes FASTA file {path_to_file} does not exist.\n')
        exit(1)
    if os.path.getsize(path_to_file) == 0:
        print(f'\nERROR: Genomes FASTA file {path_to_file} is empty.\n')
        exit(1)
    with open(path_to_file, 'r') as input_file:
        headers = []
        for line in input_file:
            if line[0] == '>':
                headers.append(line.strip().lstrip('>'))
        if len(headers) != len(set(headers)):
            print('\nERROR: FASTA headers in input genomes file must be unique.\n')
            exit(1)


def read_assay_file(path_to_file):
    """Read input assay file and parse each list into a tuple:
    (assay_name,fwd_oligo_name,fwd_oligo_seq,rev_oligo_name,rev_oligo_seq,probe_oligo_name,probe_oligo_seq)
    Return a tuple of these tuples."""
    if os.path.exists(path_to_file) == False or os.path.isfile(path_to_file) == False:
        print(f'\nERROR: Assay CSV file {path_to_file} does not exist.\n')
        exit(1)
    if os.path.getsize(path_to_file) == 0:
        print(f'\nERROR: Assay CSV file {path_to_file} is empty.\n')
        exit(1)
    assays = []
    with open(path_to_file, 'r') as input_file:
        for line in input_file:
            line = line.rstrip().split(',')
            # Check if line has the accepted number of comma-separated values
            if len(line) not in [5, 7]:
                print('\nERROR: Line in assay file not properly formatted:')
                print(','.join(line))
                print('\nAssay file lines must be a comma-separated list:')
                print('assay_name,fwd_primer_name,fwd_primer_seq,rev_primer_name,rev_primer_seq,probe_name,probe_seq') 
                exit(1)
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
                    print('\nERROR: Assay and oligo names cannot be empty:')
                    print(','.join(line))
                    print()
                    exit(1)
                # Check if oligo seqs contain invalid bases
                bases = 'ATGCWSMKRYBVDHN'
                for oligo, seq in zip(['fwd primer', 'rev primer', 'probe'], [fwd_primer_seq, rev_primer_seq, probe_seq]):
                    if len(set(seq.upper()) - set(bases)) > 0:
                        print(f'\nERROR: Oligo seq for {oligo} contains invalid bases:')
                        print(seq)
                        print()
                        exit(1)
                # Create assay info tuple and append to assays list
                assay = (assay_name, fwd_primer_name, fwd_primer_seq, rev_primer_name, rev_primer_seq, probe_name, probe_seq)
                assays.append(assay)
    # Check that names used for assays and oligos are unique
    for name_index in [0, 1, 3, 5]:
        names = [assay[name_index] for assay in assays if assay[name_index] != '']
        if len(names) != len(set(names)):
            print('\nERROR: Assay and oligo names must be unique!\n')
            exit()
    return tuple(assays)


def run_TNTBLAST(assay_details, path_to_genomes, path_to_output, melting_temp, primer_molarity, probe_molarity):
    """Run TNTBLAST with provided assay details."""
    # Create TNTBLAST assay file from assay details (TSV file of assay name and oligo seqs)
    assay_name = assay_details[0]
    fwd_primer_seq = assay_details[2]
    rev_primer_seq = assay_details[4]
    probe_seq = assay_details[6]
    print(f'Running TNTBLAST for {assay_name} against {path_to_genomes}...')
    print(f'Fwd primer seq: {fwd_primer_seq}\nRev primer seq: {rev_primer_seq}')
    if probe_seq != '':
        print(f'Probe seq: {probe_seq}')
    assay = [assay_name, fwd_primer_seq, rev_primer_seq]
    assay += [probe_seq] if probe_seq != '' else []
    assay = '\t'.join(assay)
    path_to_tntblast_assay = os.path.join(path_to_output, assay_name + '_tntblast_assay.tsv')
    with open(path_to_tntblast_assay, 'w') as output_file:
        output_file.write(assay + '\n')    
    # Create terminal command for TNTBLAST and run
    path_to_tntblast_txt_output = os.path.join(path_to_output, assay_name + '_tntblast_output.txt')
    terminal_command = (f'tntblast -i {path_to_tntblast_assay} -d {path_to_genomes} -o {path_to_tntblast_txt_output}'
                        f' -e {melting_temp} -E {melting_temp} -t {primer_molarity / 1000000} -T {probe_molarity / 1000000}'
                        f' --best-match -m 0 -v F')
    subprocess.run(terminal_command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, shell=True)
    # GARBAGE COLLECTION (remove TNTBLAST assay TSV file)
    os.remove(path_to_tntblast_assay)
    print()


def parse_tntblast_output(assay_details, job_name, path_to_output):
    """Parses txt format output from TNTBLAST and tabulates relevant data into Pandas dataframes.
    Returns the dataframe."""
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
    path_to_tntblast_txt_input = os.path.join(path_to_output, assay_name + '_tntblast_output.txt')
    if os.path.exists(path_to_tntblast_txt_input) == False:
        print('\nERROR: Expected TNTBLAST output does not exist.\n')
        exit(1)
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
                        else:
                            entry[oligo + ' mismatches'] = int(entry[oligo + ' mismatches'])
                        if entry[oligo + ' gaps'] == 'none':
                            entry[oligo + ' gaps'] = 0
                        else:
                            entry[oligo + ' gaps'] = int(entry[oligo + ' gaps'])
                        mismatches = int(entry[oligo + ' mismatches'])
                        gaps = int(entry[oligo + ' gaps'])
                        entry[oligo + ' errors'] =  mismatches + gaps
                    # Add name and sequence of oligo to entry
                    entry[oligo + ' name'] = oligo_name
                    entry[oligo + ' seq'] = oligo_seq  
                # Extract oligo site sequences from amplicon seq
                primer_length = len(entry['forward primer seq']) + entry['forward primer gaps']
                entry['forward primer site seq'] = entry['amplicon seq'][:primer_length]
                primer_length = len(entry['reverse primer seq']) + entry['reverse primer gaps']
                entry['reverse primer site seq'] = rev_comp(entry['amplicon seq'][-primer_length:])
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
    print()
    return tntblast_results


def write_tntblast_results(job_name, path_to_output, tntblast_results):
    """Writes tabulated results in Pandas dataframe to a TSV file."""
    print('Writing PCR results...')
    # Create formatted oligo site seqs
    oligos = ['fwd_primer', 'rev_primer', 'probe']
    for oligo in oligos:
        cols = [oligo + '_name', oligo + '_seq', oligo + '_site_seq']
        oligo_site_seqs = tntblast_results[cols].drop_duplicates()
        write_seq = lambda row: write_oligo_site_variant(row[oligo + '_seq'], row[oligo + '_site_seq'])
        oligo_site_seqs[oligo + '_site_seq_new'] = oligo_site_seqs.apply(write_seq, axis=1)
        tntblast_results = pd.merge(tntblast_results, oligo_site_seqs, on=cols)
        tntblast_results[oligo + '_site_seq'] = tntblast_results[oligo + '_site_seq_new']
    # Re-order columns and save TNTBLAST results to TSV
    oligos = ['fwd_primer', 'rev_primer', 'probe']
    fields = ['_name', '_seq', '_mismatches', '_gaps', '_errors', '_site_seq']
    cols = ['assay_name', 'target', 'total_errors']
    cols += [oligo + field for oligo in oligos for field in fields]
    cols += ['amplicon_seq']
    path_to_tntblast_tsv = os.path.join(path_to_output, job_name + '_PCR_results.tsv')
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
    print('Writing assay report...')
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
    # Calculate percentage of total targets detected by TNTBLAST
    assay_report['perc_detected'] = round(assay_report['detected_targets'] * 100 / assay_report['total_targets'],2)
    # Calculate percentage of detected targets representing each level of total errors
    assay_report['perc_of_detected'] = round(assay_report['target_count'] * 100 / assay_report['detected_targets'],2)
    # Calculate percentage of total targets representing each level of total errors
    assay_report['perc_of_total'] = round(assay_report['target_count'] * 100 / assay_report['total_targets'],2)
    # Drop report rows with percent of detected below threshold
    assay_report = assay_report[assay_report['perc_of_detected']>=float(threshold)]
    # Reorder columns and write assay report to TSV file
    path_to_report_file = os.path.join(path_to_output, job_name + '_assay_report.tsv')
    assay_report = assay_report.sort_values(['assay_name', 'target_count'], ascending=[True, False])
    cols = ['assay_name', 'total_targets', 'detected_targets', 'perc_detected', 'total_errors',
            'target_count', 'perc_of_detected', 'perc_of_total']
    assay_report[cols].to_csv(path_to_report_file, sep='\t', index=False)


def write_variant_report(job_name, path_to_output, tntblast_results, path_to_genomes, threshold):
    """Takes TNTBLAST results in Pandas dataframe and writes the variant report, which describes
    oligo site variants whose prevalence exceeded the provided threshold. Writes report to TSV file."""
    print('Writing variant report...')
    # Create variant report
    variant_report = pd.DataFrame()
    for oligo in ['fwd_primer', 'rev_primer', 'probe']:
        cols = ['assay_name', oligo + '_name', oligo + '_seq', oligo + '_site_seq', oligo + '_errors']
        oligo_data = tntblast_results[cols]
        cols = ['assay_name', 'oligo_name', 'oligo_seq', 'oligo_site_variant', 'oligo_errors']
        oligo_data.columns = cols
        oligo_data = oligo_data[cols].groupby(cols).size().reset_index()
        oligo_data.columns = cols + ['target_count']
        oligo_data['oligo'] = oligo
        variant_report = pd.concat([variant_report, oligo_data], sort=True)
    # Add lower case bases to oligo site seqs
    write_seq = lambda row: write_oligo_site_variant(row['oligo_seq'], row['oligo_site_variant'])
    variant_report['oligo_site_variant'] = variant_report.apply(write_seq, axis=1)
    # Count targets detected by each assay
    targets_detected = tntblast_results['assay_name'].value_counts().reset_index()
    targets_detected.columns = ['assay_name', 'detected_targets']
    variant_report = pd.merge(variant_report, targets_detected, on='assay_name')
    # Count total targets
    variant_report['total_targets'] = count_targets(path_to_genomes)
    # Calculate percentage of total targets detected by TNTBLAST
    variant_report['perc_detected'] = round(variant_report['detected_targets'] * 100 / variant_report['total_targets'],2)
    # Calculate percentage of detected targets representing each level of total errors
    variant_report['perc_of_detected'] = round(variant_report['target_count'] * 100 / variant_report['detected_targets'],2)
    # Calculate percentage of total targets representing each level of total errors
    variant_report['perc_of_total'] = round(variant_report['target_count'] * 100 / variant_report['total_targets'],2)
    # Drop rows for variants with no errors
    variant_report = variant_report[variant_report['oligo_errors']!=0]
    # Drop report rows with percent of detected below threshold
    variant_report = variant_report[variant_report['perc_of_detected']>=float(threshold)]
    # Reorder columns an write assay report to TSV file
    path_to_report_file = os.path.join(path_to_output, job_name + '_variant_report.tsv')
    variant_report = variant_report.sort_values(['assay_name', 'oligo', 'target_count'], ascending=[True, True, False])
    cols = ['assay_name', 'oligo', 'oligo_name', 'oligo_seq', 'total_targets', 'detected_targets', 'perc_detected',
            'oligo_site_variant', 'oligo_errors', 'target_count', 'perc_of_detected', 'perc_of_total']
    variant_report[cols].to_csv(path_to_report_file, sep='\t', index=False)


def longest_common_kmer(seq_1, seq_2):
    enum_kmers = lambda seq, k: set(seq[i:i+k] for i in range(len(seq) - k + 1))
    common_kmers = lambda seq_1, seq_2, k: enum_kmers(seq_1, k) & enum_kmers(seq_2, k)
    uniq_kmer = lambda seq_1, seq_2, kmer: True if seq_1.count(kmer) == 1 and seq_2.count(kmer) == 1 else False
    l = min(len(seq_1), len(seq_2))
    k = [kmer for kmer in common_kmers(seq_1, seq_2, l) if uniq_kmer(seq_1, seq_2, kmer) == True]
    while l > 0 and len(k) == 0:
        l -= 1
        k = [kmer for kmer in common_kmers(seq_1, seq_2, l) if uniq_kmer(seq_1, seq_2, kmer) == True]
    if len(k) > 0:
        return k[0]
    else:
        return None
    

def align_seqs(seq_1, seq_2):
    if type(seq_1) != list:
        seq_1 = [seq_1]
    if type(seq_2) != list:
        seq_2 = [seq_2]
    new_seq_1 = []
    new_seq_2 = []
    for sub_seq_1, sub_seq_2 in zip(seq_1, seq_2):
        if sub_seq_1 == sub_seq_2:
            k = None
        else:
            k = longest_common_kmer(sub_seq_1, sub_seq_2)
        if k != None:
            sub_seq_1 = sub_seq_1.split(k)
            sub_seq_2 = sub_seq_2.split(k)
            sub_seq_1 = [sub_seq_1[0], k, sub_seq_1[1]]
            sub_seq_2 = [sub_seq_2[0], k, sub_seq_2[1]]
            sub_seq_1, sub_seq_2 = align_seqs(sub_seq_1, sub_seq_2)
        else:
            sub_seq_1 = [sub_seq_1]
            sub_seq_2 = [sub_seq_2]
        new_seq_1 += sub_seq_1
        new_seq_2 += sub_seq_2
    return new_seq_1, new_seq_2


def write_oligo_site_variant(oligo_seq, oligo_site_variant):
    seq_1_alignment, seq_2_alignment = align_seqs(oligo_seq, oligo_site_variant)
    formatted_oligo_site_variant = []
    for sub_seq_1, sub_seq_2 in zip(seq_1_alignment, seq_2_alignment):
        if sub_seq_1 == sub_seq_2:
            formatted_oligo_site_variant += [sub_seq_2.upper()]
        elif len(sub_seq_1) == len(sub_seq_2):
            formatted_oligo_site_variant += [sub_seq_2.lower()]
        elif len(sub_seq_1) > len(sub_seq_2):
            formatted_oligo_site_variant += ['-' * (len(sub_seq_1) - len(sub_seq_2)) + sub_seq_2.lower()]
        elif len(sub_seq_2) > len(sub_seq_1):
            formatted_oligo_site_variant += ['(' + sub_seq_2.upper() + ')']
    if seq_1_alignment[-1] == '':
        formatted_oligo_site_variant = formatted_oligo_site_variant[:-1]
    formatted_oligo_site_variant = ''.join(formatted_oligo_site_variant)
    return formatted_oligo_site_variant


def get_targets(path_to_genomes):
    """Reads all sequence in the provided genomes FASTA file and returns them in
    a dict where keys are FASTA headers and values are the nucleotide sequences."""
    with open(path_to_genomes, 'r') as input_file:
        seqs = {}
        for line in input_file:
            if line[0] == '>':
                header = line.strip().lstrip('>')
                seqs[header] = ''
            else:
                seqs[header] += line.strip()
    return seqs


def write_missed_seqs_report(job_name, path_to_output, tntblast_results, path_to_genomes):
    """Takes TNTBLAST results in Pandas dataframe and the provided genomes FASTA file and writes
    a report describing sequences in the provided genome targets that were not detected by TNTBLAST.
    Writes report to TSV file."""
    print('Writing missed sequences report...')
    missed_seqs_report = pd.DataFrame()
    seqs = get_targets(path_to_genomes)
    all_seqs = set(seqs.keys())
    for assay in tntblast_results['assay_name'].unique():
        found_seqs = set(tntblast_results[tntblast_results['assay_name']==assay]['target'].unique())
        missed_seqs = sorted(all_seqs - found_seqs)
        assay_data = pd.DataFrame()
        assay_data['target'] = missed_seqs
        assay_data['assay_name'] = assay
        if len(missed_seqs) > 0:
            assay_data['target_length'] = assay_data.apply(lambda row: len(seqs[row['target']]), axis=1)
            assay_data['total_Ns'] = assay_data.apply(lambda row: seqs[row['target']].count('N'), axis=1)
            assay_data['perc_Ns'] = round(assay_data['total_Ns'] * 100 / assay_data['target_length'], 2)
        else:
            assay_data['target_length'] = []
            assay_data['total_Ns'] = []
            assay_data['perc_Ns'] = []
        missed_seqs_report = pd.concat([missed_seqs_report, assay_data], sort=True)
    # Reorder columns and write report
    path_to_report_file = os.path.join(path_to_output, job_name + '_missed_seqs_report.tsv')
    cols = ['assay_name', 'target', 'target_length', 'total_Ns', 'perc_Ns']
    missed_seqs_report[cols].to_csv(path_to_report_file, sep='\t', index=False)


if __name__ == '__main__':
    main()
