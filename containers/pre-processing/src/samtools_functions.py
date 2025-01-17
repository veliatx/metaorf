import os


def convert_sam_to_bam(sam_file, bam_file):
    samtools_s2b_cmd = f'samtools view -@ 16 -bS {sam_file} -o {bam_file}'
    os.system(samtools_s2b_cmd)

    
def sort_bam_file(unsorted_bam_file, sorted_bam_file):
    samtools_sort_cmd = f'samtools sort -@ 16 -l 0 -o {sorted_bam_file} {unsorted_bam_file} '
    os.system(samtools_sort_cmd)
    os.remove(unsorted_bam_file)

    
def remove_duplicates(raw_bam_file, filtered_bam_file):
    samtools_remove_dup_cmd = f'samtools view -@ 16 -bS -F 0X100 {raw_bam_file} -o {filtered_bam_file} '
    os.system(samtools_remove_dup_cmd)

    
def index_reads(bam_file):
    samtool_indexing_cmd = f'samtools index -@ 16 {bam_file} '
    os.system(samtool_indexing_cmd)


def merge(input_bam_files, output_bam_file):
    samtool_merge_cmd = f'samtools merge -@ 16 {output_bam_file} %s' %(' '.join(input_bam_files))
    os.system(samtool_merge_cmd)


def filter_out_multimapper(raw_bam_file, filtered_bam_file):
    samtool_remove_multimapper = f'samtools view -d NH:1 -@ 16 -b {raw_bam_file} -o {filtered_bam_file}'
    print(samtool_remove_multimapper)
    os.system(samtool_remove_multimapper)
