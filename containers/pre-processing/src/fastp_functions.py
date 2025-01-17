import os


def trim_adapter(fastq_file, trimming_output_prefix, adapter_sequence,
                 min_length, n_threads):
    """Configures and runs fastp for adapter trimming."""
    
    fastp_cmd = (
        f'fastp '
        f'--in1={fastq_file} '
        f'--out1={trimming_output_prefix}.trimmed.fastq '
        f'--json={trimming_output_prefix}_report.json '
        f'--html={trimming_output_prefix}_report.html '
        f'--adapter_sequence {adapter_sequence} '
        f'--length_required={min_length} '
        f'--thread {n_threads}'
    )
    
    try:
        os.system(fastp_cmd)
    except:
        raise RuntimeError(f'A error occurs when trimming adapters.')
