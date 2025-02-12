"""
This python program processes Ribo-seq data and maps them to genomes and transcriptomes.
"""

import logging
import os
from pathlib import Path
import argparse
import json
import shutil
import subprocess
import pandas as pd

import star_functions
import fastp_functions
import samtools_functions
import umitools_functions


def index_genomes(reference_fasta, contaminant_fasta,
                  reference_genome_index, contaminant_genome_index,
                  gtf_file, n_threads):
    """Generates indices for hg38 and contaminants."""
    
    if not os.path.exists(contaminant_genome_index):
        os.mkdir(contaminant_genome_index)
        star_functions.generate_genome_indices(
            genome_filepaths=" ".join(contaminant_fasta),
            dir_genome_index=contaminant_genome_index,
            n_threads=n_threads,
            genome_SA_index_Nbases=5,
            genome_chr_bin_Nbits=11,
            sjdb_overhang=100,
            limit_genome_generate_RAM=31000000000)

    if not os.path.exists(reference_genome_index):
        os.mkdir(reference_genome_index)
        star_functions.generate_genome_indices(
            genome_filepaths=" ".join(reference_fasta),
            dir_genome_index=reference_genome_index,
            sjdb_GTFfile=gtf_file,
            n_threads=n_threads,
            genome_SA_index_Nbases=14,
            genome_chr_bin_Nbits=18,
            sjdb_overhang=30,
            limit_genome_generate_RAM=31000000000)


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
        os.mkdir(folder_path)

    os.system(f'chmod -R 777 {folder_path}')
    return folder_path

    
def compress_fastq(sample_list, input_dir):
    """Compresses fastq files into fastq.gz."""
    compressed_samples = []

    for sample in sample_list:
        if sample.endswith('.gz'):
            compressed_samples.append(sample)
        else:
            fastq_path = Path(f'{input_dir}/{sample}')
            cmd_compress = f'gzip -k {fastq_path}'
            try:
                os.system(cmd_compress)
            except:
                raise RuntimeError(
                        f'A error occurs when compressing fastq files.')
            compressed_samples.append(sample+'.gz')
    return compressed_samples


def merge_fastq_files(input_dir, samples, output_file):
    """Merges multiple fastq files into one."""

    compressed_samples = compress_fastq(samples, input_dir)
    input_files = ' '.join([Path(f'{input_dir}/{sample}').as_posix()
                             for sample in compressed_samples])
    cmd_combine = f'cat {input_files} > {output_file}'
    try:
        os.system(cmd_combine)
    except:
        raise RuntimeError(f'A error occurs when combining fastq files.')

    
def map_reads_to_reference(reads_file, contaminant_index, contaminants_prefix,
                           reference_index, aligned_reads_prefix,
                           reference_out_filter_multimap_Nmax, n_threads):
    """Maps reads to the contaminant and reference genomes."""
    
    # map contaminant reads
    star_functions.map_reads_to_reference(
        reads_file=reads_file, 
        index_dir=contaminant_index, 
        output_prefix=Path(f'{contaminants_prefix}_'), 
        n_threads=n_threads,
        align_ends_type=f'Local', 
        out_SAM_strand_field=f'intronMotif',
        out_SAM_attributes=f'Standard',
        out_SAM_primary_flag=f'AllBestScore',    # f'OneBestScore', 
        chim_score_separation=10,
        chim_score_min=0, 
        chim_segment_min=0, 
        out_filter_multimap_Nmax=50000,
        out_filter_mismatch_Nmax=2,
        extra_filters=False)
    
    # map the remaining reads to the reference genome
    star_functions.map_reads_to_reference(
        reads_file=Path(f'{contaminants_prefix}_Unmapped.out.mate1'), 
        index_dir=reference_index,
        output_prefix=Path(f'{aligned_reads_prefix}_'), 
        n_threads=n_threads,
        align_ends_type=f'EndToEnd',
        out_SAM_strand_field=f'None',
        out_SAM_attributes=f'All',
        out_SAM_primary_flag=f'AllBestScore', # f'OneBestScore', 
        chim_score_separation=10,
        chim_score_min=20, 
        chim_segment_min=15,   ### remove all?
        out_filter_multimap_Nmax=reference_out_filter_multimap_Nmax,
        out_filter_mismatch_Nmax=2,
        extra_filters=True)

    #samtools_functions.convert_sam_to_bam(
    #    sam_file=Path(f'{aligned_reads_prefix}_Aligned.out.sam'),
    #    bam_file=Path(f'{aligned_reads_prefix}_Aligned.out.bam'))
    
    
def prepare_bam_files(aligned_reads_prefix, umi_pattern):
    """Sorts, filters, and indexes a reads file in the bam format."""

    samtools_functions.index_reads(
        bam_file=Path(f'{aligned_reads_prefix}_Aligned.sortedByCoord.out.bam'))

    if umi_pattern:
        umitools_functions.dedup(
                input_bam=Path(f'{aligned_reads_prefix}_Aligned.sortedByCoord.out.bam'),
                output_bam=Path(f'{aligned_reads_prefix}_Aligned.filtered.sorted.bam'),
                stats_prefix=Path(f'{aligned_reads_prefix}_umi_stats'))
    else:
        os.rename(Path(f'{aligned_reads_prefix}_Aligned.sortedByCoord.out.bam'),
                    Path(f'{aligned_reads_prefix}_Aligned.filtered.sorted.bam'))
        
    samtools_functions.index_reads(
        bam_file=Path(f'{aligned_reads_prefix}_Aligned.filtered.sorted.bam'))
    samtools_functions.filter_out_multimapper(
        raw_bam_file=Path(f'{aligned_reads_prefix}_Aligned.filtered.sorted.bam'),
        filtered_bam_file=Path(f'{aligned_reads_prefix}_Aligned.unique.sorted.bam'))
    samtools_functions.index_reads(
        bam_file=Path(f'{aligned_reads_prefix}_Aligned.unique.sorted.bam'))


def conduct_qc(qc_source_folder, gtf_file, bam_file, fasta_file, experiment_name, output_file_prefix):
    """Executes QC pipeline for processed reads."""

    cmd = (f"Rscript {qc_source_folder}/RiboSeqQC.r "
           f"-g {gtf_file} "
           f"-b {bam_file} "
           f"-f {fasta_file} "
           f"-n {experiment_name} "
           f"-o {output_file_prefix}")
    os.system(cmd)


def get_data_from_json(jsonfile, keywords):
    """Retreves data from a json file."""

    if not os.path.exists(jsonfile):
        return 'Not available'

    if not keywords:
        exit('Received no keyword when getting data from a json file.')

    with open(jsonfile) as f_json:
        log_data = json.load(f_json)
        for kw in keywords:
            if kw in log_data:
                log_data = log_data[kw]
            else:
                exit("Received incorrect keywords when getting data from a json file.")
    return str(log_data)


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


def create_bam_file(samples, sample_bam_files, annotation_dir, input_dir, output_dir, suffix, experiment_name,
                    umi_pattern, adapter_sequence, contaminant_genome_index, reference_genome_index,
                    multimap_Nmax):
    """Converts fastq read files to a merged BAM file."""

    processed_reads_folder = init_folder(Path(f'{output_dir}/processed_read_files{suffix}'))
    processed_reads_prefix = Path(f'{processed_reads_folder}/{experiment_name}')
    contaminant_reads_folder = init_folder(Path(f'{output_dir}/contaminants_depletion{suffix}'))
    contaminant_reads_prefix = Path(f'{contaminant_reads_folder}/contaminants_{experiment_name}')
    aligned_reads_folder = init_folder(Path(f'{output_dir}/aligned{suffix}'))
    aligned_reads_prefix = Path(f'{aligned_reads_folder}/{experiment_name}')

    if not sample_bam_files:
        merged_fastq_file = Path(f'{processed_reads_prefix}.fastq.gz')
        merge_fastq_files(input_dir, samples, merged_fastq_file) 
        if umi_pattern:
            umi_extracted = Path(f'{processed_reads_prefix}_umi_extracted.fastq.gz')
            umi_log = Path(f'{processed_reads_prefix}_extract_umi.log')
            umitools_functions.extract_umi(
                    input_fastq=merged_fastq_file,
                    output_fastq=umi_extracted,
                    pattern=umi_pattern,
                    log_file=umi_log)
            merged_fastq_file = umi_extracted
            
        trimmed_fastq = Path(f'{processed_reads_prefix}.trimmed.fastq')
        if not adapter_sequence:
            os.system(f'gunzip {merged_fastq_file} -c > {trimmed_fastq}')
        else:
            fastp_functions.trim_adapter(
                fastq_file=merged_fastq_file,
                trimming_output_prefix=processed_reads_prefix,
                adapter_sequence=adapter_sequence,
                min_length=20,
                n_threads=16)
        
        map_reads_to_reference(
            reads_file=trimmed_fastq,
            contaminant_index=Path(f'{annotation_dir}/contaminants/{contaminant_genome_index}/'),
            contaminants_prefix=contaminant_reads_prefix,
            reference_index=Path(f'{annotation_dir}/genomes/{reference_genome_index}/'),
            aligned_reads_prefix = aligned_reads_prefix,
            reference_out_filter_multimap_Nmax=multimap_Nmax,
            n_threads=16)
    else:
        input_bam_files = [Path(f'{input_dir}/{filename}').as_posix()
                           for filename in sample_bam_files]
        output_bam_file = Path(f'{aligned_reads_prefix}_Aligned.out.bam')
        if not os.path.exists(os.path.dirname(output_bam_file)):
            os.mkdir(os.path.dirname(output_bam_file))
        samtools_functions.merge(input_bam_files, output_bam_file) 
    prepare_bam_files(aligned_reads_prefix=aligned_reads_prefix, umi_pattern=umi_pattern)
    return contaminant_reads_prefix, aligned_reads_prefix


def compute_read_number_in_bam(bam_file):
    return int(subprocess.check_output(f"samtools view -c {bam_file}", shell=True))


def create_log_file(output_log, read_trimming_log, contaminant_log, aligned_log):
    """Logs essential statistics for data produced by the pipeline."""

    with open(output_log, 'w') as output_file:
        # total input reads
        total_reads = get_data_from_json(
                read_trimming_log, ['summary', 'before_filtering', 'total_reads'])
        output_file.write(f'Number of input reads\t|\t{total_reads}\n')
        
        # reads after trimming and filtering
        reads_after_qc = get_data_from_json(
                read_trimming_log, ['summary', 'after_filtering', 'total_reads'])
        output_file.write(f'Number of reads after adapter trimming and filtering\t|\t{reads_after_qc}\n')

        # reads mapped to contaminant genomes
        reads_contaminant_unique_locus = get_data_from_log(
                contaminant_log, 'Uniquely mapped reads number')
        output_file.write(
                f'Number of uniquely mapped reads from the contaminant log\t'
                f'|\t{reads_contaminant_unique_locus}\n')

        reads_contaminant_multiple_loci = get_data_from_log(
                contaminant_log, 'Number of reads mapped to multiple loci')
        output_file.write(
                f'Number of reads from the contaminant log with multiple mapped loci\t'
                f'|\t{reads_contaminant_multiple_loci}\n')

        reads_contaminant_too_many_mapped_loci = get_data_from_log(
                contaminant_log, 'Number of reads mapped to too many loci')
        output_file.write(
                f'Number of reads from the contaminant log with too many mapped loci\t'
                f'|\t{reads_contaminant_too_many_mapped_loci}\n')

        # reads mapped to the transcriptome genomes
        reads_input_transcriptome_alignment = get_data_from_log(
                aligned_log, 'Number of input reads')
        output_file.write(
                f'Total number of reads for transcriptome alignment\t'
                f'|\t{reads_input_transcriptome_alignment}\n')

        reads_aligned_unique_locus = get_data_from_log(
                aligned_log, 'Uniquely mapped reads number')
        output_file.write(
                f'Number of uniquely mapped reads from the aligned-read log\t'
                f'|\t{reads_aligned_unique_locus}\n')

        reads_aligned_multiple_loci = get_data_from_log(
                aligned_log, 'Number of reads mapped to multiple loci')
        output_file.write(
                f'Number of reads from the aligned-read log with multiple mapped loci\t'
                f'|\t{reads_aligned_multiple_loci}\n')

        reads_aligned_too_many_mapped_loci = get_data_from_log(
                aligned_log, 'Number of reads mapped to too many loci')
        output_file.write(
                f'Number of reads from the aligned-read log with too many mapped loci\t'
                f'|\t{reads_aligned_too_many_mapped_loci}\n')

        reads_unmapped_too_many_mismatches = get_data_from_log(
                aligned_log, 'Number of reads unmapped: too many mismatches')
        reads_unmapped_too_short = get_data_from_log(
                aligned_log, 'Number of reads unmapped: too short')
        reads_unmapped_other_reasons = get_data_from_log(
                aligned_log, 'Number of reads unmapped: other')
        if 'Not available' in [reads_unmapped_too_many_mismatches,
                               reads_unmapped_too_short, reads_unmapped_other_reasons]:
            reads_unmapped = 'Not available\n'
        else:
            reads_unmapped = str(int(reads_unmapped_too_many_mismatches)
                                 + int(reads_unmapped_too_short)
                                 + int(reads_unmapped_other_reasons))
        output_file.write(
                f'Number of unmapped reads from the aligned-read log\t'
                f'|\t%s\n' %(reads_unmapped))


def create_log_for_dup(output_log, reads_prefix, umi_pattern):
    """Logs essential statistics for data produced by the pipeline."""

    if umi_pattern:
        with open(output_log, 'w') as output_file:
            n_reads_no_filtering = compute_read_number_in_bam(Path(f'{reads_prefix}_Aligned.sortedByCoord.out.bam'))
            output_file.write(
                    f'Number of reads before filtering\t'
                    f'|\t{n_reads_no_filtering}\n')
        
            n_reads_after_umi = compute_read_number_in_bam(Path(f'{reads_prefix}_Aligned.filtered.sorted.bam'))
            output_file.write(
                    f'Number of reads after deduplicated using UMIs\t'
                    f'|\t{n_reads_after_umi}\n')
            output_file.write(
                    f'% reads are removed due to duplicating UMIs\t'
                    f'|\t{1 - n_reads_after_umi/float(n_reads_no_filtering)}\n')


def main(annotation_dir, input_dir, output_dir, samples, sample_bam_files,
         experiment_name, reference_genomes, contaminant_genomes,
         genome_annotation_prefix, reference_genome_index, contaminant_genome_index,
         umi_pattern, adapter_sequence, umi_pattern_tis, adapter_sequence_tis,
         qc_source_folder, fastq_tis_files, bam_tis_files, multimap_Nmax):

    # prepare indices for tRNA and rRNA reference genomes.
    reference_fasta = [Path(f'{annotation_dir}/genomes/{genome_file}').as_posix()
                       for genome_file in reference_genomes]
    contaminant_fasta=[Path(f'{annotation_dir}/contaminants/{genome_file}').as_posix()
                       for genome_file in contaminant_genomes]
    gtf_file = Path(f'{annotation_dir}/annotations/{genome_annotation_prefix}.gtf')
    index_genomes(
        reference_fasta=reference_fasta,
        contaminant_fasta=contaminant_fasta,
        reference_genome_index=Path(f'{annotation_dir}/genomes/{reference_genome_index}/'),
        contaminant_genome_index=Path(f'{annotation_dir}/contaminants/{contaminant_genome_index}/'),
        gtf_file=gtf_file,
        n_threads=16)

    # map reads to tRNA and rRNA reference genomes. 
    contaminant_reads_prefix, aligned_reads_prefix = create_bam_file(
            samples, sample_bam_files, annotation_dir, input_dir, output_dir, '', experiment_name,
            umi_pattern, adapter_sequence, contaminant_genome_index, reference_genome_index,
            multimap_Nmax)
    conduct_qc(qc_source_folder=qc_source_folder,
               gtf_file=gtf_file,
               bam_file=f'{aligned_reads_prefix}_Aligned.filtered.sorted.bam',
               fasta_file=" ".join(reference_fasta),
               experiment_name=experiment_name,
               output_file_prefix=f'{aligned_reads_prefix}_qc')

    if fastq_tis_files or bam_tis_files:
        contaminant_tis_reads_prefix, aligned_tis_reads_prefix = create_bam_file(
                fastq_tis_files, bam_tis_files, annotation_dir, input_dir, output_dir, "_tis", experiment_name,
                umi_pattern_tis, adapter_sequence_tis, contaminant_genome_index, reference_genome_index,
                multimap_Nmax)
        tis_bam_merged = Path(f'{aligned_tis_reads_prefix}_Aligned.filtered.sorted.bam')
        conduct_qc(qc_source_folder=qc_source_folder,
                   gtf_file=gtf_file,
                   bam_file=f'{aligned_tis_reads_prefix}_Aligned.filtered.sorted.bam',
                   fasta_file=" ".join(reference_fasta),
                   experiment_name=experiment_name+'_tis',
                   output_file_prefix=f'{aligned_tis_reads_prefix}_tis_qc')
    else:
        tis_bam_merged = None

    # Logging stats
    if samples and (not sample_bam_files):
        create_log_file(output_log=Path(f'{output_dir}/pipeline_stats.log'),
                        read_trimming_log=Path(f'{output_dir}/processed_read_files/{experiment_name}_report.json'),
                        contaminant_log=Path(f'{contaminant_reads_prefix}_Log.final.out'),
                        aligned_log=Path(f'{aligned_reads_prefix}_Log.final.out'))
        create_log_for_dup(output_log=Path(f'{output_dir}/dup_stats.log'),
                           reads_prefix=aligned_reads_prefix,
                           umi_pattern=umi_pattern)
    if fastq_tis_files and (not bam_tis_files):
        create_log_file(output_log=Path(f'{output_dir}/pipeline_stats_tis.log'),
                        read_trimming_log=Path(f'{output_dir}/processed_read_files_tis/{experiment_name}_report.json'),
                        contaminant_log=Path(f'{contaminant_tis_reads_prefix}_Log.final.out'),
                        aligned_log=Path(f'{aligned_tis_reads_prefix}_Log.final.out'))
        create_log_for_dup(output_log=Path(f'{output_dir}/dup_stats_tis.log'),
                           reads_prefix=aligned_tis_reads_prefix,
                           umi_pattern=umi_pattern_tis)
    print("Data-process has completed.")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='My program')
    parser.add_argument('--experiment_name', type=str,
                        help=('Name of this experiment. This name will be used as '
                              'the prefix of result files.'))
    parser.add_argument('--samples', type=str,
                        help=('Names of the samples to process, seperated by comma. '
                              'These files need to be in the `input_dir` and '
                              'in the format of `.fastq` or `fastq.gz`.'))
    parser.add_argument('--sample_bam_files', type=str,
                        help=('Names of the bam files for the samples to process, '
                              'seperated by comma. If specified, these files will directly be '
                              'used as the input for Price. Any read preparation steps will '
                              'be skipped.'))
    parser.add_argument('--annotation_dir', type=str, help='Path to the annotation directory.')
    parser.add_argument('--input_dir', type=str, help='Path to the input directory.')
    parser.add_argument('--output_dir', type=str, help='Path to the input directory.')
    parser.add_argument('--contaminant_genomes', type=str,
                        help=('Names of the contaminant genome files in the input_dir/annotation folder, '
                              'seperated by comma.'))
    parser.add_argument('--reference_genomes', type=str,
                        help=('Nafme of the tRNA genome file in the input_dir/annotation folder, '
                              'seperated by comma.'))
    parser.add_argument('--genome_annotation_prefix', type=str,
                        help=('Prefix of the genome annotation files in the input_dir/annotation folder.'))
    parser.add_argument('--contaminant_genome_index', type=str,
                        help=('Name of the folder for the contaminant genome indices. '
                              'If this folder does not exist under input_dir/annotation, '
                              'this pipeline will call STAR and generate riboseq indices.'))
    parser.add_argument('--reference_genome_index', type=str,
                        help=('Name of the folder of the reference genome indices. '
                              'If this folder does not exist under input_dir/annotation, '
                              'this pipeline will call STAR and generate riboseq indices.'))
    parser.add_argument('--umi_pattern', type=str, default='',
                        help=('UMI pattern. Leave this empty if no UMI is used.'))
    parser.add_argument('--adapter_sequence', type=str, default='auto',
                        help=('Adapter sequence for reads. If not specified, '
                              'the adapter will be auto-detected. Default: %(default)s'))
    parser.add_argument('--umi_pattern_tis', type=str, default='',
                        help=('UMI pattern for the TIS reads. Leave this empty if no UMI is used.'))
    parser.add_argument('--adapter_sequence_tis', type=str, default='auto',
                        help=('Adapter sequence for TIS reads. If not specified, '
                              'the adapter will be auto-detected. Default: %(default)s'))
    parser.add_argument('--qc_source_folder', type=str, default='src',
                        help=('Path to the folder where the QC code is. Default: %(default)s'))
    parser.add_argument('--bam_tis_files', type=str, default=None,
                        help=('Path to the transcriptome-aligned TIS BAM file. Set this as None '
                              'if no TIS file is used. Default: %(default)s'))
    parser.add_argument('--fastq_tis_files', type=str, default=None,
                        help=('Path to the transcriptome-aligned TIS fastq file. Set this as None '
                              'if no TIS file is used. Default: %(default)s'))
    parser.add_argument('--multimap', type=int, default=10,
                        help=('The multimap value for the reads mapping to the reference genomes. '
                              'Default: %(default)s'))
    args = parser.parse_args()

    if args.samples is None and args.sample_bam_files is None:
        exit('Please specify --samples or --sample_bam_files')
    if args.experiment_name is None:
        exit('Please specify --experiment_name')
    if args.input_dir is None:
        exit('Please specify --input_dir')
    if args.output_dir is None:
        exit('Please specify --output_dir')
    if args.annotation_dir is None:
        exit('Please specify --annotation_dir')
    if args.contaminant_genome_index is None:
        exit('Please specify --contaminant_genome_index')
    if args.reference_genome_index is None:
        exit('Please specify --reference_genome_index')
    if args.contaminant_genomes is None:
        exit('Please specify --contaminant_genomes')
    if args.reference_genomes is None:
        exit('Please specify --reference_genomes')

    main(annotation_dir=Path(args.annotation_dir),
         input_dir=Path(args.input_dir),
         output_dir=init_folder(Path(args.output_dir)),
         samples=parse_arg_array(args.samples),
         sample_bam_files=parse_arg_array(args.sample_bam_files),
         experiment_name=args.experiment_name,
         reference_genomes=parse_arg_array(args.reference_genomes),
         contaminant_genomes=parse_arg_array(args.contaminant_genomes),
         genome_annotation_prefix=args.genome_annotation_prefix,
         reference_genome_index=args.reference_genome_index,
         contaminant_genome_index=args.contaminant_genome_index,
         umi_pattern=args.umi_pattern,
         adapter_sequence=args.adapter_sequence,
         umi_pattern_tis=args.umi_pattern_tis,
         adapter_sequence_tis=args.adapter_sequence_tis,
         qc_source_folder=Path(args.qc_source_folder),
         fastq_tis_files=parse_arg_array(args.fastq_tis_files),
         bam_tis_files=parse_arg_array(args.bam_tis_files),
         multimap_Nmax=args.multimap)

