"""This package includes functions to execute components in RiboTISH.

Use the run_ribocode_pipeline function if you need to execute the whole
RiboTISH pipeline.
"""


import os


def run_ribotish_pipeline(bam_file, bam_tis_file, annotation_gtf_file,
                          ribotish_data_folder, genome_fasta_file):
    """Executes all the components in ribocode for Ribo-seq calling.""" 

    check_riboseq_quality(bam_file, annotation_gtf_file, is_tis_enriched=False)
    if bam_tis_file is not None:
        check_riboseq_quality(bam_tis_file, annotation_gtf_file, is_tis_enriched=True)
    run_orf_calling(bam_file, bam_tis_file, annotation_gtf_file, genome_fasta_file, ribotish_data_folder)


def check_riboseq_quality(bam_file, gtf_file, is_tis_enriched):
    """Run the quality control component of RiboTISH for riboseq bam data."""

    cmd = (f"ribotish quality "
           f"-b {bam_file} "
           f"-g {gtf_file} "
           f"-p 16 ")
    if is_tis_enriched:
        cmd += "-t "
    os.system(cmd)


def run_orf_calling(bam_file, bam_tis_file, gtf_file, genome_fasta_file, ribotish_data_folder):
    """Run RiboTISH to predict ORF/TIS with riboseq bam files."""

    cmd = (f"ribotish predict "
           f"-b {bam_file} "
           f"-g {gtf_file} "
           f"-f {genome_fasta_file} "
           f"-o {ribotish_data_folder}/pred.txt "
           f"--alt --altcodons CTG,GTG,TTG "
           f"--longest "
           f"--minaalen 5 "
           f"--seq --aaseq --blocks ")
    if bam_tis_file is not None:
        cmd += f"-t {bam_tis_file} --harr "
    os.system(cmd)

