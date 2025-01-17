import os


def generate_genome_indices(genome_filepaths, dir_genome_index, n_threads,
                            genome_SA_index_Nbases=14, genome_chr_bin_Nbits=18,
                            sjdb_GTFfile=None, sjdb_overhang=100,
                            limit_genome_generate_RAM=31000000000):
    """Configures and runs STAR for genome indexing."""
    
    star_cmd = (
        f'STAR '
        f'--runMode genomeGenerate '
        f'--genomeFastaFiles {genome_filepaths} '
        f'--genomeDir {dir_genome_index} '
        f'--runThreadN {n_threads} '
        f'--genomeSAindexNbases {genome_SA_index_Nbases} '
        f'--genomeChrBinNbits {genome_chr_bin_Nbits} '
        f'--outFileNamePrefix {dir_genome_index} '
        f'--limitGenomeGenerateRAM {limit_genome_generate_RAM} '
    )
    if sjdb_GTFfile is not None:
        star_cmd += (
                f'--sjdbGTFfile {sjdb_GTFfile} '
                f'--sjdbOverhang {sjdb_overhang} ')

    print(star_cmd) 
    try:
        os.system(star_cmd)
    except:
        raise RuntimeError(f'A error occurs when indexing reference genomes.')
        

def map_reads_to_reference(reads_file, index_dir, output_prefix, n_threads,
                           align_ends_type='Local', out_SAM_attributes='Standard',
                           out_SAM_strand_field='None', out_SAM_primary_flag='OneBestScore', 
                           chim_score_separation=10, chim_score_min=0, chim_segment_min=0,
                           out_filter_multimap_Nmax=10, out_filter_mismatch_Nmax=999,
                           transcriptome_SAM=False, random_multimapper_order=False):
    """Maps reads to the reference genomes using STAR."""
    
    star_cmd = ( 
        f'STAR '
        f'--seedSearchStartLmaxOverLread 0.5 '
        f'--readFilesIn {reads_file} '
        f'--genomeDir {index_dir} '
        f'--outFileNamePrefix {output_prefix} '
        f'--outReadsUnmapped Fastx '
        f'--runThreadN {n_threads} '
        f'--alignEndsType {align_ends_type} '
        f'--outSAMattributes {out_SAM_attributes} '
        f'--outSAMstrandField {out_SAM_strand_field} '
        f'--outSAMprimaryFlag {out_SAM_primary_flag} '
        f'--chimScoreSeparation {chim_score_separation} '
        f'--chimScoreMin {chim_score_min} '
        f'--chimSegmentMin {chim_segment_min} '
        f'--outFilterMultimapNmax {out_filter_multimap_Nmax} '
        f'--outFilterMismatchNmax {out_filter_mismatch_Nmax} '
    )

    if transcriptome_SAM:
        star_cmd += '--quantMode TranscriptomeSAM '

    if random_multimapper_order:
        star_cmd += '--outMultimapperOrder Random '
    print(star_cmd)
    try:
        os.system(star_cmd)
    except:
        raise RuntimeError(f'A error occurs when mapping reads to reference genomes.')
