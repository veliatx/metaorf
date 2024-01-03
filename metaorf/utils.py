'''Module containing general utility functions'''

import boto3
import os
import shlex
import subprocess


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
