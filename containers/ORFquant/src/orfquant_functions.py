"""This package includes functions to execute components in ORFquant.

Use the run_orfquant_pipeline function if you need to execute the whole
ORFquant pipeline.
"""


import os

def run_orfquant_pipeline(experiment_name, bam_file, annotation_gtf_file,
                          genome_annotation_prefix, orfquant_data_folder, genome_fasta_file):
    """Executes all the components in ORFquant for Ribo-seq calling.""" 

    fasta_two_bit = genome_fasta_file[:-2] + "2bit"
    cmd = (f"Rscript {os.getcwd()}/src/orfquant_functions.r "
           f"-g {annotation_gtf_file} "
           f"-b {bam_file} "
           f"-f {fasta_two_bit} "
           f"-n {genome_annotation_prefix} "
           f"-o {orfquant_data_folder} "
           f"-x {experiment_name} ")
    os.system(cmd)
