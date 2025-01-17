import os
import json
import argparse
import boto3
import shutil


def remove_folder(folder, is_on_s3=False):
    if is_on_s3:
        cmd = f'aws s3 rm --recursive {folder}'
    else:
        cmd = f'rm -rf {folder}'
    return cmd.split()


def download_gencode_and_chess(data_folder):
    des_folder = os.path.join(data_folder, 'hg38_CHESS_supplemented_GENCODE')
    cmd = (f'aws s3 cp --recursive '
           f's3://velia-annotation-dev/hg38_CHESS_supplemented_GENCODE/ '
           f'{des_folder} ')
    return cmd.split()

    
def download_contaminant_fasta(data_folder):
    des_folder = os.path.join(data_folder, 'contaminant')
    cmd = (f'aws s3 cp --recursive '
           f's3://velia-annotation-dev/contaminant_fasta/ '
           f'{des_folder} ')
    return cmd.split()
    

def download_genomes(data_folder, genome_name):
    des_file = os.path.join(data_folder, 'genomes', genome_name)
    cmd = (f'aws s3 cp '
           f's3://velia-annotation-dev/genomes/{genome_name} '
           f'{des_file} ')
    return cmd.split()
    
    
def download_riboseq_file(riboseq_filepaths_s3, data_folder):
    command_list = []
    for riboseq in riboseq_filepaths_s3:
        cmd = (f'aws s3 cp '
               f'{os.path.dirname(riboseq)} '
               f'{data_folder} '
               f'--recursive '
               f'--exclude "*" '
               f'--include "{os.path.basename(riboseq)}*" ')
        command_list.append(cmd.split())
    return command_list


def collect_input_data_for_price(riboseq_filepaths_s3, reference_genome, data_folder):
    """Returns a list of commands for data preparation.

    The commands include (in order):
        Download genome annotation file.
        Download reference genome sequences for contaminants.
        Download the hg38 reference genome.
        Download riboseq reads.
    """
    
    command_list = []
    command_list.append(download_gencode_and_chess(data_folder))
    command_list.append(download_contaminant_fasta(data_folder))
    command_list.append(download_genomes(data_folder, reference_genome))
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
    os.system(upload_parameters_to_db)


def submit_cleaning_job(experiment_name, parameter_filepath, dependencies):
    """Submits jobs to clean up the temporary folders on the AWS EFS drive."""

    with open(parameter_filepath) as parameter_file:
        parameters = json.load(parameter_file)

    command_list = []
    command_list.append(remove_folder(parameters['input_dir']))
    command_list.append(remove_folder(parameters['output_dir']))

    job_id_list = []
    for command in command_list:
        response = boto3.client('batch', 'us-west-2').submit_job(
            jobName=f'{experiment_name}_clearning',
            jobQueue='bfx-jq-general',
            jobDefinition='arn:aws:batch:us-west-2:328315166908:job-definition/price_docker_test:6',
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
    command_list.append(remove_folder(
        os.path.join('s3://velia-piperuns-dev', f'{experiment_name}'), is_on_s3=True))

    job_id_list = []
    for command in command_list:
        response = boto3.client('batch', 'us-west-2').submit_job(
            jobName=f'{experiment_name}_clearning',
            jobQueue='bfx-jq-general',
            jobDefinition='arn:aws:batch:us-west-2:328315166908:job-definition/price_docker_test:6',
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
    command_list = collect_input_data_for_price(
            riboseq_filepaths_s3, parameters['reference_genomes'], parameters['input_dir'])

    job_id_list = []
    for command in command_list:
        response = boto3.client('batch', 'us-west-2').submit_job(
            jobName=f'{experiment_name}_data_prep',
            jobQueue='bfx-jq-general',
            jobDefinition='arn:aws:batch:us-west-2:328315166908:job-definition/price_docker_test:6',
            containerOverrides = {
                'command': command,
            },
            dependsOn=[{'jobId': job_id, 'type': 'N_TO_N'} for job_id in dependencies],
        )
        job_id_list.append(response['jobId'])
    return job_id_list
    

def submit_orf_caller_job(experiment_name, parameter_filepath, dependencies):
    """Submits a job to execute price for ribo-seq calling."""

    with open(parameter_filepath) as parameter_file:
        parameters = json.load(parameter_file)

    riboseq_filenames = [os.path.basename(riboseq)
                         for riboseq in parameters['sample_paths_s3'].split(",")]
    price_cmd = [
        'python3', 'src/main_run_price.py',
        '--experiment_name', parameters['experiment_name'],
        '--input_dir', parameters['input_dir'],
        '--output_dir', parameters['output_dir'],
        '--contaminant_genomes', parameters['contaminant_genomes'],
        '--reference_genomes', parameters['reference_genomes'],
        '--genome_annotation_prefix', parameters['genome_annotation_prefix'],
        '--contaminant_genome_index', parameters['contaminant_genome_index'],
        '--reference_genome_index', parameters['reference_genome_index'],
        '--multimap', str(parameters['multimap'])]

    # '--samples', ','.join(riboseq_filenames),
    if riboseq_filenames[0].endswith("fastq.gz") or riboseq_filenames[0].endswith("fastq"):
        price_cmd += ['--samples', ','.join(riboseq_filenames)]
    elif riboseq_filenames[0].endswith("bam"):
        price_cmd += ['--sample_bam_files', ','.join(riboseq_filenames)]
    else:
        exit(f'Unexpected type of input {riboseq_filenames[0]}')

    response = boto3.client('batch', 'us-west-2').submit_job(
        jobName=f'{experiment_name}_price',
        jobQueue='bfx-jq-general',
        jobDefinition='arn:aws:batch:us-west-2:328315166908:job-definition/price_docker_test:6',
        containerOverrides = {
            'command': price_cmd,
        },
        dependsOn=[{'jobId': job_id, 'type': 'N_TO_N'} for job_id in dependencies],
    )
    
    return [response['jobId']]


def submit_data_upload_job(experiment_name, parameter_filepath, dependencies):
    """Submits jobs to upload data from AWS EFS to S3."""

    with open(parameter_filepath) as parameter_file:
        parameters = json.load(parameter_file)

    command_list = []
    command_list.append(copy_folder_to_or_from_s3(
        src_folder=parameters['input_dir'],
        des_folder=os.path.join('s3://velia-piperuns-dev', f'{experiment_name}', 'input')))
    command_list.append(copy_folder_to_or_from_s3(
        src_folder=parameters['output_dir'],
        des_folder=os.path.join('s3://velia-piperuns-dev', f'{experiment_name}', 'output')))
        
    job_id_list = []
    for command in command_list:
        response = boto3.client('batch', 'us-west-2').submit_job(
            jobName=f'{experiment_name}_data_upload',
            jobQueue='bfx-jq-general',
            jobDefinition='arn:aws:batch:us-west-2:328315166908:job-definition/price_docker_test:6',
            containerOverrides = {
                'command': command,
            },
            dependsOn=[{'jobId': job_id, 'type': 'N_TO_N'} for job_id in dependencies],
        )
        job_id_list.append(response['jobId'])
    return job_id_list


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Submit AWS batch jobs.')
    parser.add_argument('--experiment_parameter_folder', type=str,
                        help=('Path to the folder containing all the parameter files.'
                             'This program will execute exper'))
    args = parser.parse_args()
    
    for parameter_filename in os.listdir(args.experiment_parameter_folder):
        experiment_name = parameter_filename[:-5] # to exlude .json suffix
        parameter_filepath = os.path.join(args.experiment_parameter_folder, parameter_filename)

        os.system(" ".join(remove_folder(
            os.path.join('s3://velia-piperuns-dev', f'{experiment_name}'), is_on_s3=True)))

        job_ids = []
        job_list = [submit_cleaning_job,
                    submit_data_preparation_job, submit_orf_caller_job,
                    submit_data_upload_job, submit_cleaning_job]
        for job in job_list:
            job_ids = job(experiment_name, parameter_filepath, dependencies=job_ids)

        upload_parameter_file_to_s3(experiment_name, parameter_filepath)

