import os
from pathlib import Path

 
def generate_index_for_genome(fasta_file, gtf_file, output_prefix):
    """Generates indices for genomes using gedi."""
    
    gedi_cmd = (
        f'gedi -e IndexGenome '
        f'-s {fasta_file} '
        f'-a {gtf_file} '
        f'-o {output_prefix} '
        f'-nobowtie -nostar -nokallisto '
    )
    try:
        os.system(gedi_cmd)
    except:
        raise RuntimeError(f'An error occurs when Gedi is indexing reference genomes.')

    
def call_ribo_seq_peaks(reads_file, genomic_meta_file, price_result):
    """Identifies ORFs using Price with Ribo-seq data as input."""
    
    gedi_cmd = (
        f'gedi -e Price '
        f'-reads {reads_file} '
        f'-genomic {genomic_meta_file} '
        f'-prefix {price_result} '
        f'-progress '
        f'-plot '
        f'-fdr 0.20 ')
    try:
        os.system(gedi_cmd)
    except:
        raise RuntimeError(f'An error occurs when Price is identifying peaks.')

def view_cit(prefix):
    """Views .cit file in the output of Price."""

    cit_file = Path(f'{prefix}.orfs.cit')
    cmd_view_in_bed_file = (
        f'gedi -e ViewCIT '
        f'-m bed '
        f'{cit_file} > '
        f'{prefix}.orfs.FDR20pc.bed')

    try:
        os.system(cmd_view_in_bed_file)
    except:
        raise RuntimeError(f'An error occurs when viewing .cit file in the BED format.')

    cmd_view_in_location_file = (
        f'gedi -e ViewCIT '
        f'-m Location '
        f'{cit_file} > '
        f'{prefix}.orfs.location.txt')

    try:
        os.system(cmd_view_in_location_file)
    except:
        raise RuntimeError(f'An error occurs when viewing .cit locations.')

