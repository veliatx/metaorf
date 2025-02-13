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
    params['timestamp'] = timestamp
    return sample_df, jobs_df, params


def submit_jobs(experiment_name, params, job_list, run_locally=False):
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
            curr_job_ids = job.submit(dependencies=prev_job_ids, run_locally=run_locally)

        prev_job_ids = curr_job_ids
        curr_job_ids = []


@click.command()
@click.argument('sample_sheet', type=pathlib.Path)
@click.option('--multimap', default=10, help='Maximum number of multimapping reads')
@click.option('--skip_orfcalling', is_flag=True, help='Only run alignment tasks')
@click.option('--run_locally', is_flag=True, help='Run jobs locally')
@click.option('--docker_mount_flag', default='-v /efs:/mount/efs', help='Docker flag to mount additional drives within containers. E.g. -v /src:/container/dest')
@click.option('--input_dir', default=None, help='Input directory path')
@click.option('--output_dir', default=None, help='Output directory path')
@click.option('--annotation_dir', default='/mount/efs/riboseq_callers/data/ORFrater/human_GRCh38_p14', help='Annotation directory path (within containers)')
@click.option('--contaminant_genomes', default='hg38_rRNA.fa,hg38_tRNA.fa', help='Comma-separated list of contaminant genome files')
@click.option('--reference_genomes', default='GRCh38.p14.genome.fa', help='Reference genome file')
@click.option('--genome_annotation_prefix', default='veliadb_v1_1.fixed', help='Genome annotation prefix')
@click.option('--contaminant_genome_index', default='star_veliadb_v1_1', help='Contaminant genome index')
@click.option('--reference_genome_index', default='star_veliadb_v1_1', help='Reference genome index')
@click.option('--callers', default='price,ribotish,ribocode', help='Comma-separated list of ORF callers')
@click.option('--jobQueue', default='bfx-jq-metaorf', help='AWS Batch job queue')
@click.option('--jobQueue_large_instance', default='bfx-jq-general', help='AWS Batch job queue for large instances')
@click.option('--bucket_name', default='velia-piperuns-dev', help='S3 bucket name')
@click.option('--transcript_list_file', default='transcript_TPM_240419.tsv', help='Transcript list file')
def main(sample_sheet, multimap, skip_orfcalling, run_locally, docker_mount_flag, input_dir, output_dir, annotation_dir, 
         contaminant_genomes, reference_genomes, genome_annotation_prefix, contaminant_genome_index,
         reference_genome_index, callers, jobqueue, jobqueue_large_instance, bucket_name, transcript_list_file):
    """
    SAMPLE_SHEET is a conforming CSV file
    """
    timestamp = datetime.today().strftime('%Y%m%d%H%M%S')
    if input_dir is None:
        input_dir = f'/mount/efs/riboseq_callers/data/ORFrater/input_batch_tmp_{timestamp}'
    if output_dir is None:
        output_dir = f'/mount/efs/riboseq_callers/data/ORFrater/output_batch_tmp_{timestamp}'

    params = {
        'multimap': multimap,
        'skip_orfcalling': skip_orfcalling,
        'docker_mount_flag': docker_mount_flag,
        'input_dir': input_dir,
        'output_dir': output_dir,
        'annotation_dir': annotation_dir,
        'contaminant_genomes': contaminant_genomes,
        'reference_genomes': reference_genomes,
        'genome_annotation_prefix': genome_annotation_prefix,
        'contaminant_genome_index': contaminant_genome_index,
        'reference_genome_index': reference_genome_index,
        'callers': callers,
        'jobQueue': jobqueue,
        'jobQueue_large_instance': jobqueue_large_instance,
        'bucket_name': bucket_name,
        'transcript_list_file': transcript_list_file
    }

    sample_df, jobs_df, params = define_jobs(sample_sheet, params)
    piperun_dicts = utils.build_param_dicts(sample_df, jobs_df, params)

    orf_callers = {
        "price": jobs.Price,
        "ribotish": jobs.RiboTish,
        "ribocode": jobs.Ribocode,
        "orfquant": jobs.Orfquant}
    orf_call_jobs = tuple([orf_callers[caller] for caller in params["callers"].split(",")])
    for experiment_name, params in piperun_dicts.items():
        if run_locally:
            job_list = [jobs.PrepareData, 
                        jobs.PreprocessData,
                        *orf_call_jobs,
                        jobs.PostprocessData,
                        jobs.UploadData,
                        jobs.CleanDirectories]
        else:
            job_list = [jobs.PrepareData, 
                        jobs.PreprocessData,
                        orf_call_jobs,
                        jobs.PostprocessData,
                        jobs.UploadData,
                        jobs.CleanDirectories]
        submit_jobs(experiment_name, params, job_list, run_locally=run_locally)