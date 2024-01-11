'''Module containing general utility functions'''

import os
import shlex
import subprocess

import numpy as np
import pandas as pd


def remove_folder(folder, is_on_s3=False):
    if is_on_s3:
        cmd = f'aws s3 rm --recursive {folder}'
    else:
        cmd = f'rm -rf {folder}'
    return cmd.split()


def download_gencode_and_chess(data_folder):
    """Deprecated."""
    des_folder = os.path.join(data_folder, 'hg38_CHESS_supplemented_GENCODE')
    cmd = (f'aws s3 cp --recursive '
           f's3://velia-annotation-dev/hg38_CHESS_supplemented_GENCODE/ '
           f'{des_folder} ')
    return cmd.split()

    
def download_contaminant_fasta(data_folder):
    """Deprecated."""
    des_folder = os.path.join(data_folder, 'contaminant')
    cmd = (f'aws s3 cp --recursive '
           f's3://velia-annotation-dev/contaminant_fasta/ '
           f'{des_folder} ')
    return cmd.split()
    

def download_gene_names(data_folder):
    """Deprecated."""
    des_folder = os.path.join(data_folder, 'gene_names')
    cmd = (f'aws s3 cp --recursive '
           f's3://velia-annotation-dev/gene_names/ '
           f'{des_folder} ')
    return cmd.split()
    

def download_genomes(data_folder, genome_name):
    """Deprecated."""
    des_file = os.path.join(data_folder, 'genomes', genome_name)
    cmd = (f'aws s3 cp '
           f's3://velia-annotation-dev/genomes/{genome_name} '
           f'{des_file} ')
    return cmd.split()
    
    
def download_riboseq_file(riboseq_filepaths_s3, data_folder):
    command_list = []
    for riboseq in riboseq_filepaths_s3:
        cmd = (f'aws s3 cp {riboseq} {os.path.join(data_folder, os.path.basename(riboseq))} ')
        command_list.append(cmd.split())
    return command_list


# def collect_input_data_for_orfrater(riboseq_filepaths_s3, reference_genome, data_folder):
def collect_input_data_for_orfrater(riboseq_filepaths_s3, data_folder):
    """Returns a list of commands for data preparation.

    The commands include (in order):
        (No longer) Download genome annotation file.
        (No longer) Download reference genome sequences for contaminants.
        (No longer) Download genome names.
        (No longer) Download the hg38 reference genome.
        Download riboseq reads.
    """
    
    command_list = []
    #command_list.append(download_gencode_and_chess(data_folder))
    #command_list.append(download_contaminant_fasta(data_folder))
    #command_list.append(download_gene_names(data_folder))
    #command_list.append(download_genomes(data_folder, reference_genome))
    command_list += download_riboseq_file(riboseq_filepaths_s3, data_folder)
    
    return command_list
    

def copy_folder_to_or_from_s3(src_folder, des_folder):
    """Returns a command that downloads or uploads data folder from/to AWS S3."""

    cmd = (f'aws s3 cp --recursive '
           f'{src_folder} '
           f'{des_folder} ')
    return cmd.split()


def upload_parameter_file_to_s3(experiment_name, parameter_filepath):
    """Uploads the parameter file from local machine to AWS S3.

    Note: this function executes immediately and does not submit any AWS jobs.
    """

    upload_parameters_to_db = (
            f'aws s3 cp '
            f'{parameter_filepath} '
            f's3://velia-piperuns-dev/{experiment_name}/')
    subprocess.run(shlex.split(upload_parameters_to_db))


def process_data_frame(df):
    df = df.dropna(axis=0, how="all", ignore_index=True)
    df = df.T
    header = df.iloc[0]
    df = df[1:]
    df.columns = header
    df.fillna('', inplace=True)
    return df


def parse_samplesheet(sample_csv_path):
    """
    """
    sample_df = pd.DataFrame()
    job_df = pd.DataFrame()
    reading_status = {"sample": False, "job": False}
    for index, row in pd.read_csv(sample_csv_path, header=None).iterrows():
        if np.all([str(val).strip() in ["", "nan"] for _, val in row.items()]):
            continue

        if str(row.iloc[0]).strip().startswith("["):
            for content_type in reading_status:
                reading_status[content_type] = False
            if str(row.iloc[0]).strip() == "[Samples]":
                reading_status["sample"] = True
            if str(row.iloc[0]).strip() == "[Riboseq_jobs]":
                reading_status["job"] = True
            continue

        if reading_status["sample"]:
            sample_df = pd.concat([sample_df, row], ignore_index = True, axis=1)
        if reading_status["job"]:
            job_df = pd.concat([job_df, row], ignore_index = True, axis=1)

    sample_df = process_data_frame(sample_df)
    job_df = process_data_frame(job_df)
    
    sample_df.set_index('sample_id', inplace=True)
    
    return sample_df, job_df


def build_param_dicts(sample_df, jobs_df, default_params):
    """Constructs full pipeline run parameter dictionaries"""

    piperun_dicts = {}

    for index, row in jobs_df.iterrows():

        piperun_dict = default_params.copy()
        
        job_name = f"VPR_orfcalling_{default_params['timestamp']}"
        if default_params["skip_orfcalling"]:
            job_name += "_mapping_only"
        
        chx_samples = row["CHX"].split(';')
        chx_sample_df = sample_df.loc[chx_samples]
        job_name = f"{job_name}_{'_'.join(chx_sample_df.index)}"

        sample_paths_s3 = chx_sample_df['containing_folder'] + chx_sample_df['R1_fastq_file'].to_list()
        
        piperun_dict['experiment_name'] = job_name
        piperun_dict['sample_paths_s3'] = ",".join(sample_paths_s3)
        piperun_dict['umi'] = ','.join(chx_sample_df['umi'])
        piperun_dict['adapter_sequence'] = ','.join(chx_sample_df['adaptor_sequence'])
        piperun_dict['input_dir'] = default_params['input_dir'] + f'_{index}'
        piperun_dict['output_dir'] = default_params['output_dir'] + f'_{index}'
        
        if row["TIS"] != '':
            tis_samples = row["TIS"].split(';')
            tis_sample_df = sample_df.loc[tis_samples]
            tis_paths = tis_sample_df['containing_folder'] + tis_sample_df['R1_fastq_file'].to_list()
            piperun_dict['umi_tis'] = ','.join(tis_sample_df['umi'])
            piperun_dict['adapter_sequence_tis']: ','.join(tis_sample_df['adapter_sequence'])

        piperun_dicts[job_name] = piperun_dict

    return piperun_dicts