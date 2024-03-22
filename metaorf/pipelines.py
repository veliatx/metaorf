'''Module containing pipeline entrypoints and definitions'''
import click
import pathlib
import copy

from datetime import datetime
from metaorf import utils, jobs


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

    # default_params = {
    #     'multimap': 10,
    #     'skip_orfcalling': False,
    #     'timestamp': timestamp,
    #     'input_dir': f'/mount/efs/riboseq_callers/data/ORFrater/input_batch_tmp_{timestamp}',
    #     'output_dir': f'/mount/efs/riboseq_callers/data/ORFrater/output_batch_tmp_{timestamp}',
    #     'annotation_dir': '/mount/efs/riboseq_callers/data/ORFrater/human_GRCh38_p14',
    #     'contaminant_genomes': 'hg38_rRNA.fa,hg38_tRNA.fa',
    #     'reference_genomes': 'GRCh38.p14.genome.fa',
    #     'genome_annotation_prefix': 'veliadb_v1_1.fixed',
    #     'contaminant_genome_index': 'star',
    #     'reference_genome_index': 'star',
    #     'callers': 'price,ribotish,ribocode,orfquant',
    #     'jobQueue': 'bfx-jq-metaorf',
    #     'jobQueue_large_instance': 'bfx-jq-general',
    #     'bucket_name': 'velia-piperuns-dev'
    # }
    default_params = {
        'multimap': 10,
        'skip_orfcalling': False,
        'timestamp': timestamp,
        'input_dir': f'/mount/efs/riboseq_callers/data/ORFrater/input_batch_tmp_{timestamp}',
        'output_dir': f'/mount/efs/riboseq_callers/data/ORFrater/output_batch_tmp_{timestamp}',
        'annotation_dir': '/mount/efs/riboseq_callers/data/ORFrater/human_GRCh38_p14',
        'contaminant_genomes': 'hg38_rRNA.fa,hg38_tRNA.fa',
        'reference_genomes': 'GRCh38.p14.genome.fa', # 'GRCh38.p13.genome.fa',
        'genome_annotation_prefix': 'veliadb_v1_1.fixed', # 'veliadb_v2.fixed',
        'contaminant_genome_index': 'star_veliadb_v1_1',
        'reference_genome_index': 'star_veliadb_v1_1',
        'callers': 'price,ribotish,ribocode,orfquant',
        'jobQueue': 'bfx-jq-metaorf',
        'jobQueue_large_instance': 'bfx-jq-general',
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
    utils.create_piperun_folders(experiment_name, params)

    prev_job_ids = []
    curr_job_ids = []

    for job in job_list:
        if type(job) == tuple:
            for subjob in job:
                subjob = subjob(experiment_name, copy.deepcopy(params))
                curr_job_ids.append(subjob.submit(dependencies=prev_job_ids))
            curr_job_ids = [job_id for sublist in curr_job_ids for job_id in sublist]
        else:
            job = job(experiment_name, copy.deepcopy(params))
            curr_job_ids = job.submit(dependencies=prev_job_ids)

        prev_job_ids = curr_job_ids
        curr_job_ids = []


@click.command()
@click.argument('sample_sheet', type=pathlib.Path)
@click.option('--skip_orfcalling', default=False, help='Only run alignment tasks')
def main(sample_sheet, skip_orfcalling):
    """
    SAMPLE_SHEET is a conforming CSV file
    """
    options = {"skip_orfcalling": skip_orfcalling}
    sample_df, jobs_df, params = define_jobs(sample_sheet, options)
    piperun_dicts = utils.build_param_dicts(sample_df, jobs_df, params)

    orf_callers = {
        "price": jobs.Price,
        "ribotish": jobs.RiboTish,
        "ribocode": jobs.Ribocode,
        "orfquant": jobs.Orfquant}
    orf_call_jobs = tuple([orf_callers[caller] for caller in params["callers"].split(",")])
    for experiment_name, params in piperun_dicts.items():
        job_list = [jobs.PrepareData, 
                    jobs.PreprocessData,
                    orf_call_jobs,
                    jobs.PostprocessData,
                    jobs.UploadData,
                    jobs.CleanDirectories]
        submit_jobs(experiment_name, params, job_list) 


