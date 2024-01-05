'''Module containing job definitions'''
import boto3
import click
import json
import os
import subprocess

from metaorf import utils


def submit_list_job(experiment_name, parameter_filepath, dependencies):
    """Submits jobs to clean up the temporary folders on the AWS EFS drive."""

    response = boto3.client('batch', 'us-west-2').submit_job(
        jobName=f'{experiment_name}_list',
        jobQueue='bfx-jq-general',
        jobDefinition='arn:aws:batch:us-west-2:328315166908:job-definition/ribo_mapping:3',
        containerOverrides = {
            'command': ['ls', '-R', '/'],
        },
        dependsOn=[{'jobId': job_id, 'type': 'N_TO_N'} for job_id in dependencies],
    )
    return [response['jobId']]


def submit_cleaning_job(experiment_name, parameter_filepath, dependencies):
    """Submits jobs to clean up the temporary folders on the AWS EFS drive."""

    with open(parameter_filepath) as parameter_file:
        parameters = json.load(parameter_file)

    command_list = []
    command_list.append(utils.remove_folder(parameters['input_dir']))
    command_list.append(utils.remove_folder(parameters['output_dir']))

    job_id_list = []
    for command in command_list:
        response = boto3.client('batch', 'us-west-2').submit_job(
            jobName=f'{experiment_name}_clearning',
            jobQueue='bfx-jq-general',
            jobDefinition='arn:aws:batch:us-west-2:328315166908:job-definition/ribo_mapping:3',
            containerOverrides = {
                'command': command,
            },
            dependsOn=[{'jobId': job_id, 'type': 'N_TO_N'} for job_id in dependencies],
        )
        job_id_list.append(response['jobId'])
    return job_id_list


def submit_s3_cleaning_job(experiment_name, parameter_filepath, dependencies):
    """Submits jobs to clean up the folders on the AWS S3 drive."""

    command_list = []
    command_list.append(utils.remove_folder(
        os.path.join('s3://velia-piperuns-dev', f'{experiment_name}'), is_on_s3=True))

    job_id_list = []
    for command in command_list:
        response = boto3.client('batch', 'us-west-2').submit_job(
            jobName=f'{experiment_name}_clearning',
            jobQueue='bfx-jq-general',
            jobDefinition='arn:aws:batch:us-west-2:328315166908:job-definition/ribo_mapping:3',
            containerOverrides = {
                'command': command,
            },
            dependsOn=[{'jobId': job_id, 'type': 'N_TO_N'} for job_id in dependencies],
        )
        job_id_list.append(response['jobId'])
    return job_id_list

        
def submit_data_preparation_job(experiment_name, parameter_filepath, dependencies):
    """Submits jobs to download data from AWS S3 to EFS."""

    with open(parameter_filepath) as parameter_file:
        parameters = json.load(parameter_file)
        
    riboseq_filepaths_s3 = parameters['sample_paths_s3'].split(",") 
    if 'tis_paths' in parameters:
        riboseq_filepaths_s3 += parameters['tis_paths'].split(",")
    command_list = utils.collect_input_data_for_orfrater(
        riboseq_filepaths_s3, parameters['input_dir'])
    #        riboseq_filepaths_s3, parameters['reference_genomes'], parameters['input_dir'])

    job_id_list = []
    for command in command_list:
        response = boto3.client('batch', 'us-west-2').submit_job(
            jobName=f'{experiment_name}_data_prep',
            jobQueue='bfx-jq-general',
            jobDefinition='arn:aws:batch:us-west-2:328315166908:job-definition/ribo_mapping:3',
            containerOverrides = {
                'command': command,
            },
            dependsOn=[{'jobId': job_id, 'type': 'N_TO_N'} for job_id in dependencies],
        )
        job_id_list.append(response['jobId'])
    return job_id_list
    

def submit_ribo_code_job(experiment_name, parameter_filepath, dependencies):
    """Submits a job to execute orfrater for ribo-seq calling."""

    with open(parameter_filepath) as parameter_file:
        parameters = json.load(parameter_file)

    riboseq_filenames = [os.path.basename(riboseq)
                         for riboseq in parameters['sample_paths_s3'].split(",")]
    orf_calling_cmd = [
        'python3', 'src/main_run_ribocode.py',
        '--experiment_name', parameters['experiment_name'],
        '--annotation_dir', parameters['annotation_dir'],
        '--output_dir', parameters['output_dir'],
        '--reference_genomes', parameters['reference_genomes'],
        '--genome_annotation_prefix', parameters['genome_annotation_prefix'],]

    response = boto3.client('batch', 'us-west-2').submit_job(
        jobName=f'{experiment_name}_ribocode',
        jobQueue='bfx-jq-general',
        jobDefinition='arn:aws:batch:us-west-2:328315166908:job-definition/ribocode:2',
        containerOverrides = {
            'command': orf_calling_cmd,
        },
        dependsOn=[{'jobId': job_id, 'type': 'N_TO_N'} for job_id in dependencies],
    )
    
    return [response['jobId']]


def submit_data_preprocessing_job(experiment_name, parameter_filepath, dependencies):
    """Submits a job to execute orfrater for ribo-seq calling."""

    with open(parameter_filepath) as parameter_file:
        parameters = json.load(parameter_file)

    riboseq_filenames = [os.path.basename(riboseq)
                         for riboseq in parameters['sample_paths_s3'].split(",")]
    data_prep_cmd = [
        'python3', 'src/main_run_data_prep.py',
        '--experiment_name', parameters['experiment_name'],
        '--annotation_dir', parameters['annotation_dir'],
        '--input_dir', parameters['input_dir'],
        '--output_dir', parameters['output_dir'],
        '--contaminant_genomes', parameters['contaminant_genomes'],
        '--reference_genomes', parameters['reference_genomes'],
        '--genome_annotation_prefix', parameters['genome_annotation_prefix'],
        '--contaminant_genome_index', parameters['contaminant_genome_index'],
        '--reference_genome_index', parameters['reference_genome_index']]
    
    # CHX
    if riboseq_filenames[0].endswith("fastq.gz") or riboseq_filenames[0].endswith("fastq"):
        data_prep_cmd += ['--samples', ','.join(riboseq_filenames)]
        if 'umi' in parameters and parameters['umi']:
            data_prep_cmd += ['--umi_pattern', parameters['umi']]
        if 'adapter_sequence' in parameters and parameters['adapter_sequence']:
            data_prep_cmd += ['--adapter_sequence', parameters['adapter_sequence']]
    elif riboseq_filenames[0].endswith("bam"):
        data_prep_cmd += ['--sample_bam_files', ','.join(riboseq_filenames)]
    else:
        exit(f'Unexpected type of input {riboseq_filenames[0]}')
        
    # TIS
    if 'tis_paths' in parameters:
        tis_filenames = [os.path.basename(riboseq)
                         for riboseq in parameters['tis_paths'].split(",")]
        if tis_filenames[0].endswith("fastq.gz") or tis_filenames[0].endswith("fastq"):
            data_prep_cmd += ['--fastq_tis_files', ','.join(tis_filenames)]
            if 'umi_pattern_tis' in parameters and parameters['umi_pattern_tis']:
                data_prep_cmd += ['--umi_pattern_tis', parameters['umi_pattern_tis']]
            if 'adapter_sequence_tis' in parameters and parameters['adapter_sequence_tis']:
                data_prep_cmd += ['--adapter_sequence_tis', parameters['adapter_sequence_tis']]
        elif tis_filenames[0].endswith("bam"):
            data_prep_cmd += ['--bam_tis_files', ','.join(tis_filenames)]
        else:
            exit(f'Unexpected type of input {tis_filenames[0]}')

    if 'multimap' in parameters:
        data_prep_cmd += ['--multimap', parameters['multimap']]

    response = boto3.client('batch', 'us-west-2').submit_job(
        jobName=f'{experiment_name}_data_preprocessing',
        jobQueue='bfx-jq-general',
        jobDefinition='arn:aws:batch:us-west-2:328315166908:job-definition/ribo_mapping:3',
        containerOverrides = {
            'command': data_prep_cmd,
        },
        dependsOn=[{'jobId': job_id, 'type': 'N_TO_N'} for job_id in dependencies],
    )
    
    return [response['jobId']]


def submit_data_upload_job(experiment_name, parameter_filepath, dependencies):
    """Submits jobs to upload data from AWS EFS to S3."""

    with open(parameter_filepath) as parameter_file:
        parameters = json.load(parameter_file)

    command_list = []
    command_list.append(utils.copy_folder_to_or_from_s3(
        src_folder=parameters['input_dir'],
        des_folder=os.path.join('s3://velia-piperuns-dev', f'{experiment_name}', 'input')))
    command_list.append(utils.copy_folder_to_or_from_s3(
        src_folder=parameters['output_dir'],
        des_folder=os.path.join('s3://velia-piperuns-dev', f'{experiment_name}', 'output')))
        
    job_id_list = []
    for command in command_list:
        response = boto3.client('batch', 'us-west-2').submit_job(
            jobName=f'{experiment_name}_data_upload',
            jobQueue='bfx-jq-general',
            jobDefinition='arn:aws:batch:us-west-2:328315166908:job-definition/ribo_mapping:3',
            containerOverrides = {
                'command': command,
            },
            dependsOn=[{'jobId': job_id, 'type': 'N_TO_N'} for job_id in dependencies],
        )
        job_id_list.append(response['jobId'])
    return job_id_list
