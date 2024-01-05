'''Module containing pipeline entrypoints and definitions'''
import click
import subprocess

from datetime import datetime
from metaorf import utils, jobs


@click.command()
@click.argument('param_folder', type=str)
def main(param_folder):
    """
    PARAM_FOLDER is a path to a folder containing all json parameter files
    """
    timestamp = datetime.today().strftime('%Y%m%d%H%M%S')

    default_params = {
        'orfcallers': 'ribocode',
        'multimap': 3,
        'skip_orfcalling': False,
        'timestamp': timestamp,
        'input_dir': f'/mount/efs/riboseq_callers/data/ORFrater/input_batch_tmp_{timestamp}',
        'output_dir': f'/mount/efs/riboseq_callers/data/ORFrater/output_batch_tmp_{timestamp}',
        'annotation_dir': '/mount/efs/riboseq_callers/data/ORFrater/human_GRCh38_p14',
        'contaminant_genomes': 'hg38_rRNA.fa,hg38_tRNA.fa',
        'reference_genomes': 'GRCh38.p14.genome.fa',
        'genome_annotation_prefix': 'veliadb_v1',
        'contaminant_genome_index': 'star',
        'reference_genome_index': 'star',
        'transcriptome_bed_file': 'veliadb_v1.bed',
        'pseudogenes': 'transcript_ids_on_pseudogenes.txt',
        'gene_names': 'genename_mapping.txt'
    }


    for parameter_filename in param_folder.glob('*.json'):
        experiment_name = parameter_filename.stem
        parameter_filepath = param_folder.joinpath(parameter_filename)

        subprocess.run(utils.remove_folder(f's3://velia-piperuns-dev/{experiment_name}', is_on_s3=True))

        job_ids = []
        job_list = [jobs.submit_cleaning_job,
                    jobs.submit_data_preparation_job, 
                    jobs.submit_data_preprocessing_job,
                    jobs.submit_ribo_code_job,
                    jobs.submit_data_upload_job, 
                    jobs.submit_cleaning_job]

        for job in job_list:
            job_ids = job(experiment_name, parameter_filepath, dependencies=job_ids)

        utils.upload_parameter_file_to_s3(experiment_name, parameter_filepath)

