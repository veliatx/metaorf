"""This package includes functions to execute components in ORF-RATER.

Use the run_orfrater_pipeline function if you need to execute the whole
ORF-RATER pipeline.
"""


import os
import shutil
from pathlib import Path


def run_orfrater_pipeline(orfrater_source_folder, bam_files, bam_tis_file,
                          orfrater_data_folder, min_read_length, max_read_length,
                          genome_fasta_file, pseudogenes, transcriptome_bed_file,
                          gene_names, exclude_genes, init_codon):
    """Executes all the components in ORF-Rater for Ribo-seq calling.""" 
    
    if bam_tis_file is None:
        all_bam_files = bam_files
    else:
        all_bam_files = bam_files + [bam_tis_file]

    for bam in all_bam_files:
        bam_basename = bam.as_posix().split("/")[-1].split(".")[0]
        tally_filename = f'{bam_basename}.tallies.txt'
        offset_filename = f'{bam_basename}.offsets_auto.txt'
        find_psite(orfrater_source_folder=orfrater_source_folder,
                   subdir=orfrater_data_folder,
                   offset_file=offset_filename,
                   tally_file=tally_filename,
                   cds_bed=transcriptome_bed_file,
                   bam_file=bam)
        offset_filtered = f'{bam_basename}.offsets.txt'
        filter_offset_files(orfrater_data_folder, offset_filename, offset_filtered,
                            min_read_length, max_read_length)

    prune_transcripts(orfrater_source_folder, genome_fasta_file, bam_files,
                      transcriptome_bed_file, orfrater_data_folder, pseudogenes)
    create_transcript_families(orfrater_source_folder, gene_names, orfrater_data_folder)
    find_orfs_and_types(orfrater_source_folder, orfrater_data_folder, genome_fasta_file,
                        " ".join(init_codon.split(",")))

    if bam_tis_file is None:
        print("No TIS provided. Running regress_orfs.py on regular dataset")
        offset_filename = move_offset_file(bam_files[0], orfrater_data_folder, "ND")
        regress_orfs(orfrater_source_folder=orfrater_source_folder,
                     bamfiles=bam_files,
                     orfrater_data_folder=orfrater_data_folder,
                     subdir_name="ND",
                     offset_file=offset_filename,
                     exclude_genes=exclude_genes)
        rate_regression_output(orfrater_source_folder, orfrater_data_folder,
                               folder_names=["ND"])
    else:
        print("Running regress_orfs.py with initiation site enriched dataset")
        offset_filename_tis = move_offset_file(bam_tis_file, orfrater_data_folder, "TIS")
        regress_orfs(orfrater_source_folder=orfrater_source_folder,
                     bamfiles=[bam_tis_file],
                     orfrater_data_folder=orfrater_data_folder,
                     subdir_name="TIS",
                     offset_file=offset_filename_tis,
                     exclude_genes=exclude_genes,
                     startonly=True,
                     startcount=1)

        offset_filename = move_offset_file(bam_files[0], orfrater_data_folder, "ND")
        regress_orfs(orfrater_source_folder=orfrater_source_folder,
                     bamfiles=bam_files,
                     orfrater_data_folder=orfrater_data_folder,
                     subdir_name="ND",
                     offset_file=offset_filename,
                     exclude_genes=exclude_genes,
                     restrictbystarts=True)
        rate_regression_output(orfrater_source_folder, orfrater_data_folder,
                               folder_names=["ND", "TIS"])

    make_orf_bed(orfrater_source_folder, orfrater_data_folder)
    quantify_orfs(orfrater_source_folder, orfrater_data_folder, bam_files, offset_filename)
    return


def convert_filename_list_to_string(filenames):
    return " ".join([name.as_posix() for name in filenames])

    
def find_psite(orfrater_source_folder, subdir, offset_file, tally_file,
               cds_bed, bam_file):
    """Finds most common P-site offset for each read length in a ribosome profiling experiment."""

    cmd = (
        f'python {orfrater_source_folder}/psite_trimmed.py '
        f'--force '
        f'--subdir {subdir} '
        f'--offsetfile {offset_file} '
        f'--tallyfile {tally_file} '
        f'-p 3 '
        f'--cdsbed {cds_bed} '
        f'--minrdlen 18 '
        f'--maxrdlen 40 '
        f'{bam_file} '
    )
    try:
        os.system(cmd)
    except:
        raise RuntimeError(f'An error occurs when using ORF-RATER to find p-sites.')


def filter_offset_files(data_folder, offset_filename, filtered_offset_filename,
                        min_length, max_length):
    """Filters out reads with length below and above thresholds."""

    with open(Path(data_folder, filtered_offset_filename), "w") as filtered_offset_file:
        for line in open(Path(data_folder, offset_filename)):
            elements = line.split()
            footprint_length = int(elements[0])
            if footprint_length >= min_length and footprint_length <= max_length:
                filtered_offset_file.write(line)


def prune_transcripts(orfrater_source_folder, genome_fasta_file, bam_files,
                      transcriptome_bed_file, orfrater_data_folder, pseudogenes):
    """Uses ribosome profiling data to remove unwanted transcripts from a transcriptome."""

    cmd = (
        f'python {orfrater_source_folder}/prune_transcripts.py '
        f'{genome_fasta_file} '
        f'{convert_filename_list_to_string(bam_files)} '
        f'--force '
        f'--peakfrac 1 '
        f'--minreads 20 '
        f'--minlen 28 '
        f'--maxlen 31 '
        f'--inbed {transcriptome_bed_file} '
        f'--outbed {orfrater_data_folder}/transcripts.bed '
        f'--summarytable {orfrater_data_folder}/tidsummary.txt '
        f'-v -p 8 '
    )
    
    if pseudogenes is not None:
        cmd += f'--pseudogenes {pseudogenes} '
        
    try:
        print(cmd)
        os.system(cmd)
    except:
        raise RuntimeError(f'An error occurs when using ORF-RATER to prune transcripts.')


def create_transcript_families(orfrater_source_folder, gene_names, orfrater_data_folder):
    """Identifies overlapping transcripts (i.e. transcript families) from a bed file."""

    cmd = (
        f'python {orfrater_source_folder}/make_tfams.py '
        f'-v '
        f'--force '
        f'--tfamstem {orfrater_data_folder}/tfams '
        f'--inbed {orfrater_data_folder}/transcripts.bed '
    )
    if gene_names is not None:
        cmd += f'-g {gene_names} '

    try:
        print(cmd)
        os.system(cmd)
    except:
        raise RuntimeError(f'An error occurs when using ORF-RATER to create '
                           f'transcript families.')


def find_orfs_and_types(orfrater_source_folder, orfrater_data_folder, genome_fasta_file,
                        init_codon):
    """Identifies all possible ORFs in a transcriptome."""

    cmd = (
        f'python {orfrater_source_folder}/find_orfs_and_types.py '
        f'{genome_fasta_file} '
        f'--tfamstem {orfrater_data_folder}/tfams '
        f'--inbed {orfrater_data_folder}/transcripts.bed '
        f'--orfstore {orfrater_data_folder}/orf.h5 '
        f'--codons {init_codon} '
        f'-v '
        f'-p 24 '
        f'2> {orfrater_data_folder}/find_orfs_and_types.errors.txt'
    )
    try:
        print(cmd)
        os.system(cmd)
    except:
        raise RuntimeError(f'An error occurs when using ORF-RATER to '
                           f'identify all possible ORFs in a transcriptome.')


def move_offset_file(bam_file, orfrater_data_folder, subdir_name):
    """Moves the offset file into a new folder."""

    subdir = os.path.join(orfrater_data_folder, subdir_name)
    if os.path.exists(subdir):
        os.system(f'rm -rf {subdir}/*')
    else:
        os.mkdir(subdir)

    bam_basename = bam_file.as_posix().split("/")[-1].split(".")[0]
    offset_file = f'{bam_basename}.offsets.txt'
    shutil.copyfile(os.path.join(orfrater_data_folder, offset_file),
                    os.path.join(subdir, offset_file))

    return offset_file


def regress_orfs(orfrater_source_folder, bamfiles,
                 orfrater_data_folder, subdir_name, offset_file,
                 exclude_genes, startonly=False,
                 startcount=None, restrictbystarts=False):
    """Uses linear regression to identify likely sites of translation."""

    subdir = os.path.join(orfrater_data_folder, subdir_name)
    if not os.path.exists(subdir):
        os.mkdir(subdir)

    cmd = (
        f'{orfrater_source_folder}/regress_orfs.py '
        f'{convert_filename_list_to_string(bamfiles)} '
        f'--subdir {subdir} '
        f'-p 8 -v '
        f'--orfstore {orfrater_data_folder}/orf.h5 '
        f'--inbed {orfrater_data_folder}/transcripts.bed '
        f'--offsetfile {offset_file} '
        f'--exclude {" ".join(exclude_genes)} '
    )
    if startonly:
        cmd += f'--startonly '

    if startcount is not None:
        cmd += f'--startcount {startcount} '

    if restrictbystarts:
        cmd += f'--restrictbystarts {orfrater_data_folder}/TIS '

    cmd += f'2> {subdir}/regress_orfs.errors.txt '

    try:
        print(cmd)
        os.system(cmd)
    except:
        raise RuntimeError(f'An error occurs when using ORF-RATER to '
                           f'identify likely sites of translation.')

    if not os.path.exists(os.path.join(subdir, 'regression.h5')):
        raise FileNotFoundError(f'The regression.h5 was not generated in '
                                f'the regress_orfs function of ORF-rater.')

def rate_regression_output(orfrater_source_folder, orfrater_data_folder,
                           folder_names):
    """Creates a translation rating for each ORF."""

    regressfiles = " ".join([os.path.join(orfrater_data_folder, name)
                             for name in folder_names])

    cmd = (
        f'{orfrater_source_folder}/rate_regression_output.py '
        f'{regressfiles} '
        f'--minperleaf 1 '
        f'--minforestscore 0.5 '
        f'-v -p 16 '
        f'--ratingsfile {orfrater_data_folder}/orfratings.h5 '
        f'2> {orfrater_data_folder}/rate_regression_output.errors.txt'
    )
    try:
        print(cmd)
        os.system(cmd)
    except:
        raise RuntimeError(f'An error occurs when using ORF-RATER to '
                           f'create a translation rating for each ORF.')

def make_orf_bed(orfrater_source_folder, orfrater_data_folder):
    """Generates a BED file of the final ORF ratings."""

    cmd = (
        f'python {orfrater_source_folder}/make_orf_bed.py '
        f'--ratingsfile {orfrater_data_folder}/orfratings.h5 '
        f'--minrating 0.1 '
        f'-c '
        f'--inbed {orfrater_data_folder}/transcripts.bed '
        f'--outbed {orfrater_data_folder}/orfratings.bed '
    )
    try:
        print(cmd)
        os.system(cmd)
    except:
        raise RuntimeError(f'An error occurs when using ORF-RATER to '
                           f'generate a BED file of the final ORF ratings.')

def quantify_orfs(orfrater_source_folder, orfrater_data_folder, bam_files, offset_file):
    """Uses linear regression to quantify expression of the ORFs identified by ORF-RATER."""

    base_names = ""
    for bam in bam_files:
        base_names += (" " + os.path.basename(bam.as_posix()).split(".")[0])
    base_names = base_names.strip()

    cmd = (
        f'python {orfrater_source_folder}/quantify_orfs.py '
        f'{convert_filename_list_to_string(bam_files)} '
        f'--names {base_names} '
        f'--force '
        f'-v -p 16 '
        f'--subdir {orfrater_data_folder} '
        f'--offsetfile ND/{offset_file} '
        f'--metagenefile ND/metagene.txt '
        f'--quantfile quant.h5 '
        f'--ratingsfile {orfrater_data_folder}/orfratings.h5 '
        f'--inbed  {orfrater_data_folder}/transcripts.bed ' 
        f'--minrating 0.4 '
    )
    try:
        print(cmd)
        os.system(cmd)
    except:
        raise RuntimeError(f'An error occurs when using ORF-RATER to '
                           f'quantify expression of the ORFs identified by ORF-RATER.')
