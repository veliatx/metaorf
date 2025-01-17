import os
import shutil


cmd = (f"docker run "
           f"-v /home/ec2-user/efs-mount-point/mnt/efs0/riboseq_callers/data/ORFrater/human_GRCh38_p14:/usr/src/human_GRCh38_p14 "
           f"-v /home/ec2-user/efs-mount-point/mnt/efs0/riboseq_callers/data/ORFrater/output_test/:/usr/src/output riborf "
           f"python3 src/main_run_riborf.py "
           f"--experiment_name VPR_orfcalling_20240109220623_iPSC-rep1_SRR9113064 "
           f"--annotation_dir /usr/src/human_GRCh38_p14 "
           f"--output_dir /usr/src/output "
           f"--reference_genomes GRCh38.p13.genome.fa "
           f"--genome_annotation_prefix veliadb_v2.fixed ")
os.system(cmd)