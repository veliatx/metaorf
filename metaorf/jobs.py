'''Module containing job definitions'''
import boto3
import json
import os

from metaorf import utils
from pathlib import Path
import subprocess
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


JOB_DEFINITION_TO_CONTAINER_NAME = {
    'arn:aws:batch:us-west-2:328315166908:job-definition/ribo_mapping:4': 'ribo_mapping',
    'arn:aws:batch:us-west-2:328315166908:job-definition/ribocode:2': 'ribocode',
    'arn:aws:batch:us-west-2:328315166908:job-definition/ribotish:1': 'ribotish',
    'arn:aws:batch:us-west-2:328315166908:job-definition/price_docker_test:7': 'price',
    'arn:aws:batch:us-west-2:328315166908:job-definition/orfquant:3': 'orfquant',
    'arn:aws:batch:us-west-2:328315166908:job-definition/post_processing:2': 'post-processing',
}


class Job:
    """
    A class to represent an AWS Batch job.

    Attributes:
    -----------
    params : dict
            All parameters for the full pipeline run containing this job
    command_list : list(str)
            Shell compatible list of commands to be run


    Methods:
    --------
    submit():
        Submits this job
    stop():
        Cancel or terminate job dependent on queue status

    """

    def __init__(self, params, command_list):
        """

        Parameters:
        -----------
        params : dict
            All parameters for the full pipeline run containing this job
        command_list : list(str)
            Shell compatible list of commands to be run

        """
        self.params = params
        self.command_list = command_list


    def submit(self, dependencies, run_locally=False):
        """"""
        job_id_list = []
        for command in self.command_list:
            command = [str(c) for c in command]
            if run_locally:
                docker_command = ' '.join(command)
                container_name = f'328315166908.dkr.ecr.us-west-2.amazonaws.com/{JOB_DEFINITION_TO_CONTAINER_NAME[self.params["jobDefinition"]]}:latest'
                docker_run_command = f"docker run {self.params['docker_mount_flag']} {container_name} {docker_command}"
                logging.info("Docker command to run inside container: %s", docker_run_command)
                try:
                    result = subprocess.run(docker_run_command, shell=True, capture_output=True, text=True, check=True)
                    logging.info("Docker run output:\n%s", result.stdout)
                except subprocess.CalledProcessError as error:
                    logging.error("Docker run failed: \n %s", error.stderr)
                    raise subprocess.CalledProcessError
            else:
                response = boto3.client('batch', 'us-west-2').submit_job(
                    jobName = self.params['jobName'],
                    jobQueue = self.params['jobQueue'],
                    jobDefinition = self.params['jobDefinition'],
                    containerOverrides = {
                        'command': command,
                    },
                    dependsOn=[{'jobId': job_id, 'type': 'N_TO_N'} for job_id in dependencies],
                )
                job_id_list.append(response['jobId'])

        return job_id_list


    def stop(self):
        """Cancel or terminate job dependent on queue status """
        # TODO
        # boto3.client('batch', 'us-west-2').describe_job()
        # if status == 'queue':
        #    boto3.client('batch', 'us-west-2').cancel_job()
        # else:
        #    boto3.client('batch', 'us-west-2').terminate_job()


class PrepareData(Job):
    """
    A class to facilitate data preparation.
    """

    def __init__(self, experiment_name, params):

        job_type = 'data_prep'
        params['jobName'] = f'{experiment_name}_{job_type}'
        params['jobDefinition'] = 'arn:aws:batch:us-west-2:328315166908:job-definition/ribo_mapping:4'

        riboseq_filepaths_s3 = params['sample_paths_s3'].split(",") 
        if 'tis_paths' in params:
            riboseq_filepaths_s3 += params['tis_paths'].split(",")

        command_list = utils.download_riboseq_file(riboseq_filepaths_s3, params['input_dir'])

        super().__init__(params, command_list)


class PreprocessData(Job): 
    """
    A class to facilitate preprocessing job.
    """

    def __init__(self, experiment_name, params):

        job_type = 'data_preprocessing'
        params['jobName'] = f'{experiment_name}_{job_type}'
        params['jobDefinition'] = 'arn:aws:batch:us-west-2:328315166908:job-definition/ribo_mapping:4'

        riboseq_filenames = [os.path.basename(riboseq)
                            for riboseq in params['sample_paths_s3'].split(",")]
        data_prep_cmd = [
            'python3', 'src/main_run_data_prep.py',
            '--experiment_name', experiment_name,
            '--annotation_dir', params['annotation_dir'],
            '--input_dir', params['input_dir'],
            '--output_dir', params['output_dir'],
            '--contaminant_genomes', params['contaminant_genomes'],
            '--reference_genomes', params['reference_genomes'],
            '--genome_annotation_prefix', params['genome_annotation_prefix'],
            '--contaminant_genome_index', params['contaminant_genome_index'],
            '--reference_genome_index', params['reference_genome_index']]
        
        # CHX
        if riboseq_filenames[0].endswith("fastq.gz") or riboseq_filenames[0].endswith("fastq"):
            data_prep_cmd += ['--samples', ','.join(riboseq_filenames)]
            if 'umi' in params and params['umi']:
                data_prep_cmd += ['--umi_pattern', params['umi']]
            if 'adapter_sequence' in params and params['adapter_sequence']:
                data_prep_cmd += ['--adapter_sequence', params['adapter_sequence']]
        elif riboseq_filenames[0].endswith("bam"):
            data_prep_cmd += ['--sample_bam_files', ','.join(riboseq_filenames)]
        else:
            exit(f'Unexpected type of input {riboseq_filenames[0]}')
            
        # TIS
        if 'tis_paths' in params:
            tis_filenames = [os.path.basename(riboseq)
                            for riboseq in params['tis_paths'].split(",")]
            if tis_filenames[0].endswith("fastq.gz") or tis_filenames[0].endswith("fastq"):
                data_prep_cmd += ['--fastq_tis_files', ','.join(tis_filenames)]
                if 'umi_pattern_tis' in params and params['umi_pattern_tis']:
                    data_prep_cmd += ['--umi_pattern_tis', params['umi_pattern_tis']]
                if 'adapter_sequence_tis' in params and params['adapter_sequence_tis']:
                    data_prep_cmd += ['--adapter_sequence_tis', params['adapter_sequence_tis']]
            elif tis_filenames[0].endswith("bam"):
                data_prep_cmd += ['--bam_tis_files', ','.join(tis_filenames)]
            else:
                exit(f'Unexpected type of input {tis_filenames[0]}')

        if 'multimap' in params:
            data_prep_cmd += ['--multimap', params['multimap']]

        super().__init__(params, [data_prep_cmd])


class Ribocode(Job):
    """
    A class to represent a RiboCode job.
    """

    def __init__(self, experiment_name, params):

        job_type = 'ribocode'
        params['jobName'] = f'{experiment_name}_{job_type}'
        params['jobDefinition'] = 'arn:aws:batch:us-west-2:328315166908:job-definition/ribocode:2'

        orf_calling_cmd = [
            'python3', 'src/main_run_ribocode.py',
            '--experiment_name', experiment_name,
            '--annotation_dir', params['annotation_dir'],
            '--output_dir', params['output_dir'],
            '--reference_genomes', params['reference_genomes'],
            '--genome_annotation_prefix', params['genome_annotation_prefix'],]
    
        super().__init__(params, [orf_calling_cmd])


class RiboTish(Job):
    """
    A class to represent a Ribo-TISH job.
    """

    def __init__(self, experiment_name, params):

        job_type = 'ribotish'
        params['jobName'] = f'{experiment_name}_{job_type}'
        params['jobDefinition'] = 'arn:aws:batch:us-west-2:328315166908:job-definition/ribotish:1'

        orf_calling_cmd = [
            'python3', 'src/main_run_ribotish.py',
            '--experiment_name', experiment_name,
            '--annotation_dir', params['annotation_dir'],
            '--output_dir', params['output_dir'],
            '--reference_genomes', params['reference_genomes'],
            '--genome_annotation_prefix', params['genome_annotation_prefix'],]
    
        super().__init__(params, [orf_calling_cmd])


class Price(Job):
    """
    A class to represent a Ribo-TISH job.
    """

    def __init__(self, experiment_name, params):

        job_type = 'price'
        params['jobName'] = f'{experiment_name}_{job_type}'
        params['jobDefinition'] = 'arn:aws:batch:us-west-2:328315166908:job-definition/price_docker_test:7'
        
        orf_calling_cmd = [
            'python3', 'src/main_run_price.py',
            '--experiment_name', experiment_name,
            '--annotation_dir', params['annotation_dir'],
            '--output_dir', params['output_dir'],
            '--reference_genomes', params['reference_genomes'],
            '--genome_annotation_prefix', params['genome_annotation_prefix'],]
    
        super().__init__(params, [orf_calling_cmd])


class Orfquant(Job):
    """
    A class to represent a Orfquant job.
    """

    def __init__(self, experiment_name, params):

        job_type = 'orfquant'
        params['jobName'] = f'{experiment_name}_{job_type}'
        params['jobDefinition'] = 'arn:aws:batch:us-west-2:328315166908:job-definition/orfquant:3'
        
        orf_calling_cmd = [
            'python3', 'src/main_run_orfquant.py',
            '--experiment_name', experiment_name,
            '--annotation_dir', params['annotation_dir'],
            '--output_dir', params['output_dir'],
            '--reference_genomes', params['reference_genomes'],
            '--genome_annotation_prefix', params['genome_annotation_prefix'],]
    
        params['jobQueue'] = params['jobQueue_large_instance']
        super().__init__(params, [orf_calling_cmd])


class PostprocessData(Job):
    """
    A class to process ORF-calls into formatted datasets.
    """

    def __init__(self, experiment_name, params):

        job_type = 'data_post_process'
        params['jobName'] = f'{experiment_name}_{job_type}'
        params['jobDefinition'] = 'arn:aws:batch:us-west-2:328315166908:job-definition/post_processing:2'
        
        if 'transcript_list_file' in params:
            data_process_cmd = [
                'python3', 'src/main_run_post_processing.py',
                '--experiment_name', experiment_name,
                '--annotation_dir', params['annotation_dir'],
                '--output_dir', params['output_dir'],
                '--reference_genomes', params['reference_genomes'],
                '--genome_annotation_prefix', params['genome_annotation_prefix'],
                '--callers', params['callers'],
                '--transcript_list_file', params['transcript_list_file'],
                '--rna_seq_name', params['rna_seq_name'],
                '--organism', params['organism'],
                '--transcript_tpm_threshold', 0.5]
        else:
            data_process_cmd = [
                'python3', 'src/main_run_post_processing.py',
                '--experiment_name', experiment_name,
                '--annotation_dir', params['annotation_dir'],
                '--output_dir', params['output_dir'],
                '--reference_genomes', params['reference_genomes'],
                '--genome_annotation_prefix', params['genome_annotation_prefix'],
                '--organism', params['organism'],
                '--callers', params['callers']]
    
        super().__init__(params, [data_process_cmd])

 
class UploadData(Job):
    """
    A class to facilitate data upload to S3 upon job completion.
    """

    def __init__(self, experiment_name, params):

        job_type = 'data_upload'
        params['jobName'] = f'{experiment_name}_{job_type}'
        params['jobDefinition'] = 'arn:aws:batch:us-west-2:328315166908:job-definition/ribo_mapping:4'

        command_list = []
        
        command_list.append(
            utils.copy_folder_to_or_from_s3(
                src_folder = params['input_dir'],
                des_folder = '/'.join(['s3://velia-piperuns-dev', f'{experiment_name}', 'input'])))
        
        command_list.append(
            utils.copy_folder_to_or_from_s3(
                src_folder = params['output_dir'],
                des_folder = '/'.join(['s3://velia-piperuns-dev', f'{experiment_name}', 'output'])))
            
        super().__init__(params, command_list)


class CleanDirectories(Job):
    """
    A class to remove the local copy of the data upon job completion.
    """

    def __init__(self, experiment_name, params):

        job_type = 'clean_directories'
        params['jobName'] = f'{experiment_name}_{job_type}'
        params['jobDefinition'] = 'arn:aws:batch:us-west-2:328315166908:job-definition/ribo_mapping:4'

        command_list = []
        command_list.append(utils.remove_folder(params['input_dir']))
        command_list.append(utils.remove_folder(params['output_dir']))
            
        super().__init__(params, command_list)

