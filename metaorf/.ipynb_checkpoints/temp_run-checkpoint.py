import boto3
import json
import os

from pathlib import Path

import click
import pathlib

from datetime import datetime
from metaorf import utils, jobs

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


    def submit(self, dependencies):
        """"""
        job_id_list = []
        for command in self.command_list:
            command = [str(c) for c in command]
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

class Orfquant(Job):
    """
    A class to represent a Orfquant job.
    """

    def __init__(self, experiment_name, params):

        job_type = 'orfquant'
        params['jobName'] = f'{experiment_name}_{job_type}'
        params['jobDefinition'] = 'arn:aws:batch:us-west-2:328315166908:job-definition/orfquant:2'
        
        orf_calling_cmd = [
            'python3', 'src/main_run_orfquant.py',
            '--experiment_name', experiment_name,
            '--annotation_dir', params['annotation_dir'],
            '--output_dir', params['output_dir'],
            '--reference_genomes', params['reference_genomes'],
            '--genome_annotation_prefix', params['genome_annotation_prefix'],]
    
        super().__init__(params, [orf_calling_cmd])
        
def define_jobs(sample_sheet, params):
    """
    Setup job definition params 

    Parameters:
    -----------
    sample_sheet : pathlib.Path
        Absolute path to a CSV sample sheet file
    params : dict
        Any non-default parameters for pipeline run

    Returns:
    --------
    sample_df: pandas.DataFrame
        Table representing sample information
    job_df: pandas.DataFrame
        Table representing job information
    params: dict
        All parameters for pipeline run

    """
    sample_df, jobs_df = utils.parse_samplesheet(sample_sheet)

    timestamp = datetime.today().strftime('%Y%m%d%H%M%S')

    default_params = {
        'orfcallers': 'ribocode',
        'multimap': 10,
        'skip_orfcalling': False,
        'timestamp': 20240228215916,
        'input_dir': f'/mount/efs/riboseq_callers/data/ORFrater/input_batch_tmp_{timestamp}',
        'output_dir': f'/mount/efs/riboseq_callers/data/ORFrater/output_batch_tmp_{timestamp}',
        'annotation_dir': '/mount/efs/riboseq_callers/data/ORFrater/human_GRCh38_p14',
        'contaminant_genomes': 'hg38_rRNA.fa,hg38_tRNA.fa',
        'reference_genomes': 'GRCh38.p14.genome.fa',
        'genome_annotation_prefix': 'veliadb_v1_1.fixed',
        'contaminant_genome_index': 'star',
        'reference_genome_index': 'star',
        'callers': 'price,ribotish,ribocode,orfquant',
        'jobQueue': 'bfx-jq-metaorf',
        'bucket_name': 'velia-piperuns-dev'
    }

    default_params.update(params)

    return sample_df, jobs_df, default_params


def submit_jobs(experiment_name, params, job_list):
    """
    Submit batch jobs

    Parameters:
    -----------
    experiment_name : str
        Absolute path to a CSV sample sheet file
    params : dict
        All parameters for pipeline run
    job_list: list
        List of jobs to run in DAG format

    """
    prev_job_ids = []
    curr_job_ids = []

    for job in job_list:
        if type(job) == tuple:
            for subjob in job:
                subjob = subjob(experiment_name, params)
                curr_job_ids.append(subjob.submit(dependencies=prev_job_ids))
            curr_job_ids = [job_id for sublist in curr_job_ids for job_id in sublist]
        else:
            job = job(experiment_name, params)
            curr_job_ids = job.submit(dependencies=prev_job_ids)

        prev_job_ids = curr_job_ids
        curr_job_ids = []


def main():
    """
    SAMPLE_SHEET is a conforming CSV file
    """
    options = {"skip_orfcalling": False}
    sample_df, jobs_df, params = define_jobs(
        sample_sheet="/home/ec2-user/bfx-containers/metaorf/SampleSheet_founders_data_orf_calls_test.csv",
        params=options)
    piperun_dicts = utils.build_param_dicts(sample_df, jobs_df, params)

    for experiment_name, params in piperun_dicts.items():
        submit_jobs(experiment_name,
                    params=params,
                    job_list=[jobs.Orfquant])
        
main()