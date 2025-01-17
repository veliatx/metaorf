"""This package includes functions to execute components in RiboCode.

Use the run_ribocode_pipeline function if you need to execute the whole
RiboCode pipeline.
"""


import os


def run_ribocode_pipeline(transcriptome_bam_file, annotation_gtf_file,
                          ribocode_data_folder, genome_fasta_file):
    """Executes all the components in ribocode for Ribo-seq calling.""" 

    cmd = (f"RiboCode_onestep "
        f"-g {annotation_gtf_file} "
        f"-f {genome_fasta_file} "
        f"-r {transcriptome_bam_file} "
        f"-l no "
        f"-mA 5 "
        f"-A CTG,GTG,TTG "
        f"-outbed "
        f"-o {ribocode_data_folder}")
    print(cmd)
    os.system(cmd)
