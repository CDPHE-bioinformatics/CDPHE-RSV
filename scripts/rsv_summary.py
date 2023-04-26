#! /usr/bin/env python

import argparse
import sys
import pandas as pd
from datetime import date
import re

#### FUNCTIONS #####
def getOptions(args=sys.argv[1:]):
    parser = argparse.ArgumentParser(description="Parses command.")
    parser.add_argument('--sample_name_array')
    parser.add_argument('--workbook_path')
    parser.add_argument("--cov_out_files",  help= "txt file with list of bam file paths")
    parser.add_argument('--percent_cvg_files', help = 'txt file with list of percent cvg file paths')
    parser.add_argument('--assembler_version')
    parser.add_argument('--project_name')

    options = parser.parse_args(args)
    return options

def create_list_from_write_lines_input(write_lines_input):
    list = []
    with open(write_lines_input, 'r') as f:
        for line in f:
            list.append(line.strip())
    return list

def concat_cov_out(cov_out_file_list):
     # initiate dataframe for concatenation
    df = pd.DataFrame()
    sample_name_list = []
    samtools_mapped_reads_list = []
    samtools_depth_list = []
    samtools_baseq_list = []
    samtools_mapq_list = []

    #loop through bam file stats files and pull data
    for file in cov_out_file_list:
        d = pd.read_csv(file, sep = '\t')
        if re.search('barcode', file):
            # for nanopore runs
            sample_name = re.findall('/([0-9a-zA-Z_\-\.]+)_barcode', file)[0]
        else:
            # for illumina runs
            sample_name = re.findall('/([0-9a-zA-Z_\-\.]+)_coverage.txt', file)[0]

        # pull data from samtools output
        num_reads = d.numreads[0]
        depth = d.meandepth[0]
        baseq = d.meanbaseq[0]
        mapq = d.meanmapq[0]

        sample_name_list.append(sample_name)
        samtools_mapped_reads_list.append(num_reads)
        samtools_depth_list.append(depth)
        samtools_baseq_list.append(baseq)
        samtools_mapq_list.append(mapq)

    df['sample_name'] = sample_name_list
    df['mapped_reads'] = samtools_mapped_reads_list
    df['mean_depth'] = samtools_depth_list
    df['mean_base_quality'] = samtools_baseq_list
    df['mean_map_quality'] = samtools_mapq_list

    return df

def concat_percent_cvg(percent_cvg_file_list):
    df_list = []
    for file in percent_cvg_file_list:
        d = pd.read_csv(file, dtype = {'sample_name' : object})
        df_list.append(d)

    df = pd.concat(df_list)

    return df

def get_df_spike_mutations(variants_csv):
    def get_sample_name(fasta_header):
        sample_name = str(re.findall('CO-CDPHE-([0-9a-zA-Z_\-\.]+)', fasta_header)[0])
        return sample_name 
    
    # read in variants file
    variants = pd.read_csv(variants_csv, dtype = {'sample_name' : object})
    variants = variants.drop(columns = 'fasta_header')

    #### filter variants for spike protein varaints in rbd and pbcs #####
    crit = variants.gene == 'S'
    critRBD = (variants.codon_position >= 461) & (variants.codon_position <= 509)
    critPBCS = (variants.codon_position >= 677) & (variants.codon_position <= 694)
    crit732 = variants.codon_position == 732
    crit452 = variants.codon_position == 452
    crit253 = variants.codon_position == 253
    crit13 = variants.codon_position == 13
    crit145 = variants.codon_position == 145 # delta plus AY.4.2
    crit222 = variants.codon_position == 222 # delta plus AY.4.2
    critdel = variants.variant_name.str.contains('del') # mainly for 69/70 del

    spike_variants_df = variants[crit & (critRBD | critPBCS | crit732 | critdel | crit452 | crit253 | crit13 | crit145 | crit222)]
    spike_variants_df = spike_variants_df.reset_index(drop = True)

    # generate a df of the sample and their spike variants
    sample_name_list_variants = spike_variants_df.sample_name.unique().tolist()

    df = pd.DataFrame()
    sample_name_list = []
    variant_name_list = []

    seperator = '; '

    for sample_name in sample_name_list_variants:
        sample_name_list.append(sample_name)

        crit = spike_variants_df.sample_name == sample_name
        f = spike_variants_df[crit]
        f = f.reset_index()

        mutations = []
        for row in range(f.shape[0]):
            mutations.append(f.variant_name[row])
        mutations_string = seperator.join(mutations)
        variant_name_list.append(mutations_string)

    df['sample_name'] = sample_name_list
    df['spike_mutations'] = variant_name_list

    return df

def concat_results(sample_name_list, workbook_path, project_name, 
                   assembler_version,
                   cov_out_df, percent_cvg_df, spike_variants_df):

    # set some functions for getting data formatted
    def get_sample_name_from_fasta_header(fasta_header):
        sample_name = str(re.findall('CO-CDPHE-([0-9a-zA-Z_\-\.]+)', fasta_header)[0])
        return sample_name

    def create_fasta_header(sample_name):
        return 'CO-CDPHE-%s' % sample_name
    
    # create dataframe and fill with constant strings
    df = pd.DataFrame()
    df['sample_name'] = sample_name_list
    df = df.set_index('sample_name')
    df['analysis_date'] = str(date.today())
    df['assembler_version'] = assembler_version
    print(df)

    # read in workbook
    workbook = pd.read_csv(workbook_path, sep = '\t', dtype = {'sample_name' : object, 'hsn' : object})
    workbook = workbook.set_index('sample_name')

    # read in nextclade csv
    nextclade = nextclade.drop(columns = 'fasta_header')
    nextclade = nextclade.set_index('sample_name')

    # set index on the samtools_df and percent_cvg_df and variants_df to prepare for joining
    cov_out_df = cov_out_df.set_index('sample_name')
    percent_cvg_df = percent_cvg_df.set_index('sample_name')
    spike_variants_df = spike_variants_df.set_index('sample_name')

    # join
    j = df.join(workbook, how = 'left')
    j = j.join(percent_cvg_df, how = 'left')
    j = j.join(cov_out_df, how = 'left')
    j = j.join(nextclade, how = 'left')
    j = j.join(pangolin, how = 'left')
    j = j.join(spike_variants_df, how = 'left')
    j = j.reset_index()

    # add fasta header
    j['fasta_header'] = j.apply(lambda x:create_fasta_header(x.sample_name), axis=1)

    # add assembled column and fill in failed assembles with 0% coveage
    j.percent_coverage = j.percent_coverage.fillna(value = 0)

    def get_assembly_pass(percent_coverage):
        if percent_coverage == 0:
            return False
        if percent_coverage > 0:
            return True      
    j['assembly_pass'] = j.apply(lambda x:get_assembly_pass(x.percent_coverage), axis = 1)

    # order columns
    columns = j.columns.tolist()
    columns.sort()
    primary_columns = ['hsn', 'sample_name', 'project_name', 'plate_name', 
                       'run_name', 'analysis_date', 'run_date', 'assembly_pass', 
                 'percent_coverage', 'nextclade', 'pangolin_lineage', 
                 'expanded_lineage', 'spike_mutations']
    for column in columns:
         if column not in primary_columns:
              primary_columns.append(column)
    
    j = j[primary_columns]

    outfile = '%s_sequencing_results.csv' % project_name
    j.to_csv(outfile, index = False)

    return j

def make_wgs_horizon_output (results_df, project_name):
    results_df['report_to_epi'] = ''
    results_df['Run_Date'] = str(date.today())

    # rename columns 
    results_df = results_df.rename(columns = {'hsn' : 'accession_id', 'pango_designation_version' : 'pangoLEARN_version'})

    col_order = ['accession_id', 'percent_coverage', 'pangolin_lineage', 'pangolin_version',
                 'report_to_epi', 'Run_Date', 'pangoLEARN_version']
    
    results_df = results_df[col_order]

    outfile = "%s_wgs_horizon_report.csv" % project_name
    results_df.to_csv(outfile, index = False)

if __name__ == '__main__':
    options = getOptions()

    sample_name_array = options.sample_name_array
    workbook_path = options.workbook_path
    cov_out_files = options.cov_out_files
    percent_cvg_files = options.percent_cvg_files
    assembler_version = options.assembler_version
    project_name = options.project_name

    # create lists from the column table txt file input
    sample_name_list = create_list_from_write_lines_input(write_lines_input=sample_name_array)
    cov_out_file_list = create_list_from_write_lines_input(write_lines_input = cov_out_files)
    percent_cvg_file_list = create_list_from_write_lines_input(write_lines_input=percent_cvg_files)
    
    # concat cov_out files and percent_cvg files
    cov_out_df = concat_cov_out(cov_out_file_list=cov_out_file_list)
    percent_cvg_df = concat_percent_cvg(percent_cvg_file_list=percent_cvg_file_list)

    # create results file
    results_df = concat_results(sample_name_list = sample_name_list,
                                workbook_path = workbook_path,
                                project_name = project_name,
                                assembler_version = assembler_version,
                                cov_out_df=cov_out_df,
                                percent_cvg_df=percent_cvg_df, 
                                spike_variants_df = spike_variants_df)
    
    # create wgs horizon output
    make_wgs_horizon_output(project_name = project_name,
                            results_df=results_df)

    print('DONE!')