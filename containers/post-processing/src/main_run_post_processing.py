"""
This python program processes Ribo-seq data and maps them to genomes and transcriptomes.
"""

from pathlib import Path
import argparse

import trim_reads_to_psites
import find_common_calls
import generate_features


def parse_arg_array(array_in_string_form):
    """Parses arguments and restores their array form."""

    if array_in_string_form is None:
        return []
    else:
        return [elem.strip() for elem in array_in_string_form.split(',')]


def main(annotation_dir, output_dir, experiment_name, genome_annotation_prefix,
         reference_genomes, fastq_tis_files, bam_tis_files, orf_callers, organism,
         transcript_list_file, rna_seq_name, transcript_tpm_threshold):
    reference_fasta = [Path(f'{annotation_dir}/genomes/{genome_file}').as_posix()
                       for genome_file in reference_genomes]
    trim_reads_to_psites.trim_reads_to_psites_main(
        data_path=Path(f'{output_dir}/aligned/'),
        experiment_name=experiment_name,
        reference_fasta=reference_fasta)
    
    if fastq_tis_files or bam_tis_files:
        trim_reads_to_psites.trim_reads_to_psites_main(
            data_path=Path(f'{output_dir}/aligned_tis/'),
            experiment_name=experiment_name,
            reference_fasta=reference_fasta)
        
    find_common_calls.find_common_calls_main(
        data_path=output_dir,
        experiment_name=experiment_name,
        transcript_bed_file=Path(f'{annotation_dir}/annotations/{genome_annotation_prefix}.bed'), 
        reference_fasta=reference_fasta,
        orf_callers=orf_callers,
        transcript_list_file=Path(f'{annotation_dir}/annotations/{transcript_list_file}'),
        rna_seq_name=rna_seq_name,
        transcript_tpm_threshold=transcript_tpm_threshold)

    generate_features.generate_features_main(
        experiment_name=experiment_name,
        data_dir=output_dir,
        transcript_bed_file=Path(f'{annotation_dir}/annotations/{genome_annotation_prefix}.bed'), 
        tis_transformer_file=Path(f'{annotation_dir}/tis_transformer/tis_transformer_predictions.{genome_annotation_prefix}.txt'),
        callers=orf_callers,
        organism=organism,
        qc_file=Path(f'{output_dir}/aligned/{experiment_name}_qc.tsv'))
    print("Data post-processing has completed.")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='My program')
    parser.add_argument('--experiment_name', type=str,
                        help=('Name of this experiment. This name will be used as '
                              'the prefix of result files.'))
    parser.add_argument('--annotation_dir', type=str, help='Path to the annotation directory.')
    parser.add_argument('--output_dir', type=str, help='Path to the input directory.')
    parser.add_argument('--genome_annotation_prefix', type=str,
                        help=('Prefix of the genome annotation files in the input_dir/annotation folder.'))
    parser.add_argument('--reference_genomes', type=str,
                        help=('Nafme of the tRNA genome file in the input_dir/annotation folder, '
                              'seperated by comma.'))
    parser.add_argument('--bam_tis_files', type=str, default=None,
                        help=('Path to the transcriptome-aligned TIS BAM file. Set this as None '
                              'if no TIS file is used. Default: %(default)s'))
    parser.add_argument('--fastq_tis_files', type=str, default=None,
                        help=('Path to the transcriptome-aligned TIS fastq file. Set this as None '
                              'if no TIS file is used. Default: %(default)s'))
    parser.add_argument('--callers', type=str, default=None,
                        help=('Names of the callers being used. Default: %(default)s'))
    parser.add_argument('--organism', type=str, default="human",
                        help=('Organism name of the samples being used. Default: %(default)s'))
    parser.add_argument('--transcript_list_file', type=str, default=None,
                        help=('List of transcript activities in diffrent RNA-seq experiments. Default: %(default)s'))
    parser.add_argument('--rna_seq_name', type=str, default=None,
                        help=('Name of the RNA-seq experiment correspounding to the ribo-seq experiment. '
                              'Default: %(default)s'))
    parser.add_argument('--transcript_tpm_threshold', type=float, default=0.1,
                        help=('The TPM threshold for qualified transcripts. '
                              'Default: %(default)s'))
    args = parser.parse_args()

    main(annotation_dir=Path(args.annotation_dir),
         output_dir=Path(args.output_dir),
         experiment_name=args.experiment_name,
         genome_annotation_prefix=args.genome_annotation_prefix,
         reference_genomes=parse_arg_array(args.reference_genomes),
         fastq_tis_files=parse_arg_array(args.fastq_tis_files),
         bam_tis_files=parse_arg_array(args.bam_tis_files),
         orf_callers=parse_arg_array(args.callers),
         organism=args.organism,
         transcript_list_file=args.transcript_list_file,
         rna_seq_name=args.rna_seq_name,
         transcript_tpm_threshold=args.transcript_tpm_threshold)
