"""
This script converts coordinates to sequences in the fasta format.
"""

import os
from Bio import SeqIO

chrom_groups = [
        [1, 7, 13, 19],
        [2, 8, 14, 20],
        [3, 9, 15, 21],
        [4, 10, 16, 22],
        [5, 11, 17, "X"],
        [6, 12, 18, "Y"],
        ]

num_groups = len(chrom_groups)
chrom_to_group = {}
for group_idx in range(num_groups):
    for chrom_id in chrom_groups[group_idx]:
        chrom_to_group["chr"+str(chrom_id)] = group_idx
print(chrom_to_group)


def divide_bed_file_by_chrom(bed_file, output_bed, chrom_to_group, num_groups):
    output_bedfiles = [open(output_bed %(group_idx+1), "w") for group_idx in range(num_groups)]

    for line in open(bed_file):
       chrom_id = line.strip().split()[0]
       if chrom_id not in chrom_to_group:
           print(chrom_id)
           continue
       output_bedfiles[chrom_to_group[chrom_id]].write(line)

    for output_bed in output_bedfiles:
        output_bed.close()

def remove_long_sequences(input_file, output_file, limit):
    with open(output_file, "w") as output_handle:
        for record in SeqIO.parse(input_file, "fasta"):
            if len(record.seq) <= limit:
                SeqIO.write(record, output_handle, "fasta")
            else:
                print(len(record.seq))

        os.remove(input_file)


if __name__ == "__main__":
    data_folder = "/home/ec2-user/bfx-containers/riboseq_callers/TIS_transformer/data/"
    ref_genome = os.path.join(data_folder, "GRCh38.p14.genome.fa")
    bed_file = os.path.join(data_folder, "veliadb_v1_1.bed")
    output_bed_file_format = os.path.join(data_folder, "veliadb_v1_1.%d.bed")
    divide_bed_file_by_chrom(bed_file, output_bed_file_format, chrom_to_group, num_groups)

    for group_idx in range(num_groups):
        print(f"group {group_idx}")
        current_bed_file = output_bed_file_format %(group_idx+1)
        fasta_output = os.path.join(data_folder, "veliadb_v1_1.transcript.%d.fasta" %(group_idx+1))
        cmd = (f"bedtools getfasta "
               f"-fi {ref_genome} -bed {current_bed_file} -split -s -name > {fasta_output}.tmp ")
        os.system(cmd)
        remove_long_sequences(input_file=fasta_output+".tmp", output_file=fasta_output, limit=49990)
