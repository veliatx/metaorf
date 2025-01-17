"""
This python program processes Ribo-seq data and uses orfquant to identify ORFs.
"""

import os
from pathlib import Path
import argparse
import pandas as pd

import orfquant_functions


def parse_arg_array(array_in_string_form):
    """Parses arguments and restores their array form."""

    if array_in_string_form is None:
        return []
    else:
        return [elem.strip() for elem in array_in_string_form.split(',')]


def init_folder(folder_path, remove_exisiting_content=True):
    """Initializes a folder and sets its I/O access to public."""

    if os.path.exists(folder_path):
        if remove_exisiting_content:
            os.system(f'rm -rf {folder_path}/*') 
    else:
        os.makedirs(folder_path)

    os.system(f'chmod -R 777 {folder_path}')
    return folder_path


def get_data_from_log(logfile, keywords):
    """Retrieves data from a log file."""

    if not os.path.exists(logfile):
        return 'Not available'

    if not keywords:
        exit('Received no keyword when getting data from a log file.')

    count = 0
    log_data = None
    for line in open(logfile):
        elements = line.split("|")
        if elements and elements[0].strip() == keywords:
            count += 1
            log_data = elements[1].strip()

    if count == 0:
        exit("Received incorrect keywords when getting data from a log file.")
    elif count > 1:
        exit("Received keywords matching with multiple lines.")
    else:
        return log_data
            

def is_validated(qc_file, log_file):
    """Filter requirement:
        a. Footprint size within 24 to 34
	    b. ≥ 5 million usable reads (multimappers + uniquely mapped reads)
	    c. Orfik periodicity score ≥ 0.5
    """
    # read Orfik QC tsv
    try:
        with open(qc_file) as tsv_file:
            tsv = pd.read_table(tsv_file)
            if tsv["size_peak"][0] < 24 or tsv["size_peak"][0] > 34:
                print(f"Peak size {tsv['size_peak'][0]}")
                return False
            if tsv["periodicity_score"][0] < 0.5:
                print(f"Peak size {tsv['periodicity_score'][0]}")
                return False
    except:
        print(f"QC file not found: {qc_file}")
        return False

    # read log file
    reads_aligned_unique_locus = get_data_from_log(log_file, 'Uniquely mapped reads number')
    reads_aligned_multiple_loci = get_data_from_log(log_file, 'Number of reads mapped to multiple loci')
    reads_useful = int(reads_aligned_unique_locus) + int(reads_aligned_multiple_loci)
    if reads_useful < 2*(10**6):
        print(f"Too few reads: {reads_useful}")
        return False
    return True


def main(annotation_dir, output_dir, experiment_name, reference_genomes, genome_annotation_prefix):
    reference_fasta = [Path(f'{annotation_dir}/genomes/{genome_file}').as_posix()
                       for genome_file in reference_genomes]
    gtf_file = Path(f'{annotation_dir}/annotations/{genome_annotation_prefix}.gtf')
    aligned_reads_prefix = Path(f'{output_dir}/aligned/{experiment_name}')

    init_folder(Path(f'{output_dir}/orfquant_results/'))
    if is_validated(qc_file=f'{aligned_reads_prefix}_qc.tsv', log_file=f'{aligned_reads_prefix}_Log.final.out'):
        orfquant_functions.run_orfquant_pipeline(
                experiment_name=experiment_name,
                bam_file= Path(f'{aligned_reads_prefix}_Aligned.filtered.sorted.bam'),
                annotation_gtf_file=gtf_file,
                genome_annotation_prefix=genome_annotation_prefix,
                orfquant_data_folder=Path(f'{output_dir}/orfquant_results'),
                genome_fasta_file=" ".join(reference_fasta))
    else:
        print("Orfquant did not run.")
    print("Orfquant has completed the task.") 


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='My program')
    parser.add_argument('--experiment_name', type=str, 
                        default='VPR_orfcalling_20240109220623_iPSC-rep1_SRR9113064',
                        help=('Name of this experiment. This name will be used as '
                              'the prefix of result files.'))
    parser.add_argument('--annotation_dir', type=str,
                        default='/usr/src/data/human_GRCh38_p14/',
                        help='Path to the annotation directory.')
    parser.add_argument('--output_dir', type=str,
                        default='/usr/src/data/output_test/',
                        help='Path to the output directory.')
    parser.add_argument('--reference_genomes', type=str,
                        default='GRCh38.p14.genome.fa',
                        help=('Nafme of the tRNA genome file in the input_dir/annotation folder, '
                              'seperated by comma.'))
    parser.add_argument('--genome_annotation_prefix', type=str,
                        default ='veliadb_v1_1',
                        help=('Prefix of the genome annotation files in the input_dir/annotation folder.'))
    args = parser.parse_args()

    if args.experiment_name is None:
        exit('Please specify --experiment_name')
    if args.output_dir is None:
        exit('Please specify --output_dir')
    if args.annotation_dir is None:
        exit('Please specify --annotation_dir')
    if args.reference_genomes is None:
        exit('Please specify --reference_genomes')
    if args.genome_annotation_prefix is None:
        exit('Please specify --genome_annotation_prefix')

    main(annotation_dir=Path(args.annotation_dir),
         output_dir=Path(args.output_dir),
         experiment_name=args.experiment_name,
         reference_genomes=parse_arg_array(args.reference_genomes),
         genome_annotation_prefix=args.genome_annotation_prefix)

