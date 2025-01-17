import os


def extract_umi(input_fastq, output_fastq, pattern, log_file):
    cmd = (f'umi_tools extract --stdin={input_fastq} --bc-pattern={pattern} '
           f'--log={log_file} --stdout {output_fastq}')
    os.system(cmd)

    
def dedup(input_bam, output_bam, stats_prefix):
    cmd = f'umi_tools dedup -I {input_bam} --output-stats={stats_prefix} -S {output_bam}'
    os.system(cmd)

