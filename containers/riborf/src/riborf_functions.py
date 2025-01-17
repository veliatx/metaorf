"""This package includes functions to execute components in RibORF.2.0.

Use the run_riborf_pipeline function if you need to execute the whole
RibORF pipeline.
"""


import os
from pathlib import Path
import pandas as pd


def run_riborf_pipeline(transcriptome_bam_file, genome_annotation_prefix, offset_file,
                        riborf_data_folder, genome_fasta_file):
    """Executes all the components in riborf for Ribo-seq calling.""" 

    # run_ORFannotate(genome_fasta_file, annotation_genePred_file, riborf_data_folder)
    annotation_genePred_file = Path(f'{genome_annotation_prefix}.genePred.txt')
    #run_readDist(transcriptome_bam_file, annotation_genePred_file, riborf_data_folder)
    create_offset_file(offset_file, riborf_data_folder)
    run_offsetCorrect(transcriptome_bam_file, riborf_data_folder, offset_file)

    candidate_orfs_file = Path(f'{genome_annotation_prefix}.candidateORF.genepred.txt')
    run_ribORF(transcriptome_bam_file, candidate_orfs_file, riborf_data_folder)


def run_ORFannotate(genome_fasta_file, annotation_genePred_file, riborf_data_folder):
    cmd = (
        f"perl RibORF/RibORF.2.0/ORFannotate.pl "
        f"-g {genome_fasta_file} "
        f"-t {annotation_genePred_file} "
        f"-o {riborf_data_folder}")
    os.system(cmd)


def run_readDist(transcriptome_bam_file, annotation_genePred_file, riborf_data_folder):
    cmd = (
        f"perl RibORF/RibORF.2.0/readDist.pl "
        f"-f {transcriptome_bam_file} "
        f"-g {annotation_genePred_file} "
        f"-o {riborf_data_folder} "
        f"-d 24,25,26,27,28,29,30,31,32,33,34, ")
    os.system(cmd)


def create_offset_file(offset_file, riborf_data_folder):
    df = pd.read_csv(offset_file, sep='\t')
    with open(f"{riborf_data_folder}/offset.corretion.parameters.txt", "w") as ofile:
        for col in df.columns:
            if col.startswith("Offset"):
                fragment_size = int(col.split(".")[1])
                offset = abs(int(df[col][0]))
                ofile.write(f"{fragment_size}\t{offset}\n")


def run_offsetCorrect(transcriptome_bam_file, riborf_data_folder, offset_file):
    transcriptome_bam_file = str(transcriptome_bam_file)
    cmd = (
        f"perl RibORF/RibORF.2.0/offsetCorrect.pl "
        f"-r {transcriptome_bam_file} "
        f"-p {offset_file} "
        f"-o {transcriptome_bam_file[:-3]}corrected.{transcriptome_bam_file[-3:]} "
    )
    os.system(cmd)


def run_ribORF(transcriptome_bam_file, candidate_orfs_file, riborf_data_folder):
    cmd = (
        f"perl RibORF/RibORF.2.0/ribORF.pl "
        f"-f {transcriptome_bam_file} "
        f"-c {candidate_orfs_file} "
        f"-o {riborf_data_folder}")
    os.system(cmd)