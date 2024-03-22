"""This script trims ribo-seq reads to their p-sites and computes the psite coverage.

The p-sites offsets are computed by ORFik. The p-site coverage values are computed
per base pair and stored in the bigwig format.
"""

import os
import pandas as pd


def remove_reads_with_mismatches(data_path, experiment_name):
    cigar_limit = '(cigar =~ "^[0-9]*M$" || cigar =~ "^[0-9]*M[0-9]*N[0-9]*M$")'
    cmd_samtools = (f"samtools view -e '{cigar_limit}' "
                    f"{data_path}/{experiment_name}_Aligned.filtered.sorted.bam "
                    f"-o {data_path}/{experiment_name}_all_match.bam")
    os.system(cmd_samtools)

    
def convert_bam_to_bed(data_path, experiment_name):    
    cmd_bedtools = (f"bedtools bamtobed -bed12 "
                    f"-i {data_path}/{experiment_name}_all_match.bam "
                    f"> {data_path}/{experiment_name}.bed")
    os.system(cmd_bedtools)
    

def get_psite_offset(data_path, experiment_name):
    df = pd.read_csv(f"{data_path}/{experiment_name}_qc.tsv", sep='\t')
    offset = {}
    for col in df.columns:
        if col.startswith("Offset"):
            fragment_size = int(col.split(".")[1])
            offset[fragment_size] = abs(int(df[col][0]))
    return offset


def write_count_val(data_path, experiment_name, count):
    with open(f"{data_path}/{experiment_name}_count.txt", "w") as ofile:
        ofile.write(f"{count}\n")


def trim_reads(data_path, experiment_name, offsets, group):
    """Trims reads to their psites.
    
    This function trims each read based on their length. This function can properly
    process reads mapped to mulitple exon blocks. For each psite, only the first base
    is reported and stored in the BED format.
    """

    count_valid_reads = 0
    with open(f"{data_path}/{experiment_name}_{group}.psite.bed", "w") as ofile:
        for line in open(f"{data_path}/{experiment_name}.bed"):
            elements = line.strip().split()
            start = int(elements[1])
            block_sizes = [int(elem) for elem in elements[10].split(",")]
            block_starts = [int(elem) for elem in elements[11].split(",")]

            fragment_size = sum(block_sizes)
            if fragment_size in offsets:
                offset = offsets[fragment_size]
                count_valid_reads += 1
            else:
                continue

            strand = elements[5]
            if strand == "+":
                for block_idx, block_size in enumerate(block_sizes):
                    if offset >= block_size:
                        offset -= block_size
                    else:
                        new_start = start + block_starts[block_idx] + offset
                        break
            elif strand == "-":
                for block_idx in range(len(block_sizes)-1, -1, -1):
                    if offset >= block_sizes[block_idx]:
                        offset -= block_sizes[block_idx]
                    else:
                        new_start = start + block_starts[block_idx] + block_sizes[block_idx] - offset - 1
                        break
            else:
                print(line)

            elements[1] = str(new_start)
            elements[2] = str(new_start+1)
            new_line = "%s\n" %("\t".join(elements[:6]))
            ofile.write(new_line)


    write_count_val(data_path, experiment_name, count_valid_reads)


def normalized(ifilename, ofilename, total_reads):
    """Normalizes coverage values of each base pair.
    
    normalized_value = raw_coverage * 10 ** 6 / total_reads
    """

    with open(ofilename, "w") as ofile:
        for line in open(ifilename):
            elements = line.strip().split()
            elements[3] = str(float(elements[3])*(10**6)/(float(total_reads)))
            ofile.write("\t".join(elements)+"\n")

            
def convert_bg_to_bw(bg_file, bw_file, fai_file):
    cmd = (f"/usr/src/tools/bedGraphToBigWig "
           f"{bg_file} "
           f"{fai_file} "
           f"{bw_file} ")
    os.system(cmd)


def compute_coverage(data_path, experiment, fai_file, bam_suffix="_all.psite.bed"):
    """Computes psite converage.
    
    The input psites are stored in the BED format. We computes their convergae on each strands
    and report them seperately. The psite coverage is reported in the bigwig format.
    """


    bed_to_bam_cmd = (f"bedtools bedtobam -i {data_path}/{experiment}{bam_suffix} -g {fai_file} " 
                      f"> {data_path}/{experiment}{bam_suffix}.bam")
    os.system(bed_to_bam_cmd)
    sort_cmd = f"samtools sort {data_path}/{experiment}{bam_suffix}.bam > {data_path}/{experiment}{bam_suffix}.sorted.bam"
    os.system(sort_cmd)
    os.remove(f"{data_path}/{experiment}{bam_suffix}")
    os.remove(f"{data_path}/{experiment}{bam_suffix}.bam")
    bam_suffix = bam_suffix + ".sorted.bam"

    for symbol, suffix in [("+", "pos"), ("-", "neg")]:
        convert_to_bg_cmd = (f"bedtools genomecov -ibam {data_path}/{experiment}{bam_suffix} "
                             f"-bg -split -strand {symbol} > {data_path}/{experiment}{bam_suffix}.{suffix}.raw.bg")
        os.system(convert_to_bg_cmd)
            
        sort_cmd = (f"sort -k 1,1 -k2,2n {data_path}/{experiment}{bam_suffix}.{suffix}.raw.bg"
                    f" > {data_path}/{experiment}{bam_suffix}.{suffix}.bg")
        os.system(sort_cmd)
            
        total_reads = int(open(f"{data_path}/{experiment}_count.txt").readline())
        normalized(f"{data_path}/{experiment}{bam_suffix}.{suffix}.bg",
                   f"{data_path}/{experiment}{bam_suffix}.{suffix}.norm.bg",
                   total_reads)
        convert_bg_to_bw(bg_file=f"{data_path}/{experiment}{bam_suffix}.{suffix}.norm.bg",
                         bw_file=f"{data_path}/{experiment}{bam_suffix}.{suffix}.norm.bw",
                         fai_file=fai_file)
            
        os.remove(f"{data_path}/{experiment}{bam_suffix}.{suffix}.raw.bg")
        os.remove(f"{data_path}/{experiment}{bam_suffix}.{suffix}.bg")
        os.remove(f"{data_path}/{experiment}{bam_suffix}.{suffix}.norm.bg")
    os.remove(f"{data_path}/{experiment}{bam_suffix}")
        

def get_fai(reference_fasta):
    if len(reference_fasta) == 1:
        if os.path.exists(f"{reference_fasta[0]}.fai"):
            faidx_cmd = f"samtools faidx {reference_fasta[0]}"
            os.system(faidx_cmd)
        return f"{reference_fasta[0]}.fai"
    else:
        merge_cmd = "cat"
        for fasta_file in reference_fasta:
            merge_cmd += f" {fasta_file}"
        merge_cmd += f" > {reference_fasta[0]}.merged.fa"
        os.system(merge_cmd)

        faidx_cmd = f"samtools faidx {reference_fasta[0]}.merged.fa"
        os.system(faidx_cmd)
        return f"{reference_fasta[0]}.merged.fa.fai"


def trim_reads_to_psites_main(data_path, experiment_name, reference_fasta):
    """This is the main function to run functions trimming reads to psites and computing the psite coverage."""

    remove_reads_with_mismatches(data_path, experiment_name)
    convert_bam_to_bed(data_path, experiment_name)
    trim_reads(
        data_path, experiment_name,
        offsets=get_psite_offset(data_path, experiment_name),
        group="all")

    compute_coverage(
        data_path, experiment_name,
        fai_file=get_fai(reference_fasta),
        bam_suffix="_all.psite.bed")

