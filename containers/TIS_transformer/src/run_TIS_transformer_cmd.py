"""This program executes TIS transformer through its command line function.
"""

import os

if __name__ == "__main__":
    fasta_input = "/home/ec2-user/bfx-containers/riboseq_callers/TIS_transformer/data/veliadb_v1_1.transcript.%d.fasta"
    model_path = "/home/ec2-user/bfx-containers/riboseq_callers/TIS_transformer/TIS_transformer/models/proteome/TIS_transformer_L_%d.ckpt"
    output_prediction = "/home/ec2-user/bfx-containers/riboseq_callers/TIS_transformer/data/TIS_predictions_%d"

    for group in range(1, 7):
        cmd = (f"transcript_transformer predict "
               f"{fasta_input %group} fa {model_path %group} "
               f"--out_prefix {output_prediction %group} "
               f"--num_workers 32 "
               f"--prob_th 0.001 "
               f"--accelerator cpu ")
        print(cmd)
        os.system(cmd)

    print("Done!")
