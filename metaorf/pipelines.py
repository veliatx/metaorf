'''Module containing pipeline entrypoints and definitions'''
import click

from datetime import datetime
from metaorf import utils, jobs

@click.command()
@click.argument('sample_sheet', type=str)
@click.option('--skip_orfcalling', default=False, help='Only run alignment tasks')
def main(sample_sheet, skip_orfcalling):
    """
    SAMPLE_SHEET is a conforming CSV file
    """

    sample_df, jobs_df = utils.parse_samplesheet(sample_sheet)

    timestamp = datetime.today().strftime('%Y%m%d%H%M%S')

    default_params = {
        'orfcallers': 'ribocode',
        'multimap': 3,
        'skip_orfcalling': skip_orfcalling,
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

    piperun_dicts = utils.build_param_dicts(sample_df, jobs_df, default_params)

    for experiment_name, param_dict in piperun_dicts.items():

        utils.create_piperun_folder(experiment_name, param_dict)

        prev_job_ids = []
        curr_job_ids = []
        job_list = [jobs.PrepareData, 
                    jobs.PreprocessData,
                    (jobs.Ribocode, jobs.Price, jobs.RiboTish),
                    jobs.UploadData]

        for job in job_list:
            if len(job) > 1:
                for subjob in job:
                    curr_job_ids.append(subjob.submit(experiment_name, param_dict, dependencies=prev_job_ids))
            else:
                curr_job_ids = job.submit(experiment_name, param_dict, dependencies=prev_job_ids)
                
            prev_job_ids = curr_job_ids

        #utils.upload_parameter_file_to_s3(experiment_name, parameter_filepath)

