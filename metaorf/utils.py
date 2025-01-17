'''Module containing general utility functions'''
import boto3
import os
import shlex
import subprocess

import numpy as np
import pandas as pd

 
def create_piperun_folders(experiment_name, params):
    """
    Create S3 directory for piperun output

    Parameters:
    -----------
    experiment_name: str
        Name of job being run
    params: dict
        All parameters for pipeline run

    """
    s3 = boto3.client('s3')

    input_path = f'{experiment_name}/input/'
    output_path = f'{experiment_name}/output/'

    s3.put_object(Bucket=params['bucket_name'], Key=input_path)
    s3.put_object(Bucket=params['bucket_name'], Key=output_path)


def remove_folder(folder, is_on_s3=False):
    if is_on_s3:
        cmd = f'aws s3 rm --recursive {folder}'
    else:
        cmd = f'rm -rf {folder}'
    return cmd.split()
    
    
def download_riboseq_file(riboseq_filepaths_s3, data_folder):
    """
    """
    command_list = []
    for riboseq in riboseq_filepaths_s3:
        cmd = (f'aws s3 cp {riboseq} {os.path.join(data_folder, os.path.basename(riboseq))} ')
        command_list.append(cmd.split())
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
    Parse CSV samplesheet to extract samples and jobs dataframes
    
    Parameters:
    -----------
    sample_csv_path : pathlib.Path
        Absolute path to a CSV sample sheet file

    Returns:
    --------
    sample_df: pandas.DataFrame
        Table representing sample information
    job_df:
        Table representing job information

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


def get_s3_paths(sample_df):
    sample_path_list = []
    for idx, row in sample_df.iterrows():
        containing_folder = row["containing_folder"]
        R1_fastq_files = row["R1_fastq_file"].split(";")
        if containing_folder.startswith("s3"):
            for fastq_file in R1_fastq_files:
                sample_path_list.append(os.path.join(
                    containing_folder, fastq_file))
        else:
            for fastq_file in R1_fastq_files:
                sample_path_list.append(os.path.join(
                    "s3://velia-data-dev/VDC_003_ngs/primary/raw_data/",
                    containing_folder,
                    "Analysis/1/Data/fastq",
                    fastq_file))
    return sample_path_list


def get_uniform_param(params):
    current_param = None
    for _, param in params.items():
        if current_param is None:
            current_param = param
        elif current_param != param:
            error_msg = f"Each experiment can only have one type of UMI and adapter. But we got: {param_df}"
            exit(error_msg)
    return current_param


def build_param_dicts(sample_df, jobs_df, params):
    """
    Build a dictionary of parameter dictionaries, one for each individual job
    
    Parameters:
    -----------
    sample_df: pandas.DataFrame
        Table representing sample information
    job_df: pandas.DataFrame
        Table representing job information
    params: dict
        All parameters for pipeline run

    Returns:
    --------
    piperun_dicts: dict
        Map of experiment_name -> parameter dictionary

    """
    piperun_dicts = {}

    for index, row in jobs_df.iterrows():

        piperun_dict = params.copy()
        
        job_name = f"VPR_orfcalling_{params['timestamp']}"
        if params["skip_orfcalling"]:
            job_name += "_mapping_only"
        
        chx_samples = row["CHX"].split(';')
        chx_sample_df = sample_df.loc[chx_samples]
        job_name = f"{job_name}_{'_'.join(chx_sample_df.index)}"


        sample_paths_s3 = get_s3_paths(chx_sample_df)
        piperun_dict['experiment_name'] = job_name
        piperun_dict['sample_paths_s3'] = ",".join(sample_paths_s3)
        piperun_dict['umi'] = get_uniform_param(chx_sample_df['umi'])
        piperun_dict['adapter_sequence'] = get_uniform_param(chx_sample_df['adaptor_sequence'])
        piperun_dict['organism'] = get_uniform_param(chx_sample_df['organism'])
        piperun_dict['rna_seq_name'] = get_uniform_param(chx_sample_df['rna_seq_name'])
        piperun_dict['input_dir'] = params['input_dir'] + f'_{index}'
        piperun_dict['output_dir'] = params['output_dir'] + f'_{index}'
        
        if row["TIS"] != '':
            tis_samples = row["TIS"].split(';')
            tis_sample_df = sample_df.loc[tis_samples]
            tis_paths = get_s3_paths(tis_sample_df)
            piperun_dict['tis_paths'] = ",".join(tis_paths)
            piperun_dict['umi_tis'] = get_uniform_param(tis_sample_df['umi'])
            piperun_dict['adapter_sequence_tis'] = get_uniform_param(tis_sample_df['adaptor_sequence'])

        piperun_dicts[job_name] = piperun_dict

    return piperun_dicts
