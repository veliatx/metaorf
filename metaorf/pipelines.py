'''Module containing pipeline entrypoints and definitions'''
import boto3
import click
import json
import os
import subprocess

from metaorf import utils, jobs


@click.command()
@click.argument('experiment_parameter_folder', type=str,
              help='Path to the folder containing all parameter files')
def main(experiment_parameter_folder):
    """
    """
    for parameter_filename in experiment_parameter_folder.glob('*.json'):
        experiment_name = parameter_filename.stem
        parameter_filepath = experiment_parameter_folder.joinpath(parameter_filename)

        subprocess.run(utils.remove_folder(f's3://velia-piperuns-dev/{experiment_name}', is_on_s3=True))

        job_ids = []
        #job_list = [submit_ribo_code_job]
        job_list = [jobs.submit_cleaning_job,
                    jobs.submit_data_preparation_job, jobs.submit_data_preprocessing_job,
                    jobs.submit_ribo_code_job,
                    jobs.submit_data_upload_job, jobs.submit_cleaning_job]
        for job in job_list:
            job_ids = job(experiment_name, parameter_filepath, dependencies=job_ids)

        utils.upload_parameter_file_to_s3(experiment_name, parameter_filepath)

