import os
import shutil

cmd = (f"docker run "
       f"-v /home/ec2-user/efs-mount-point/mnt/efs0/riboseq_callers/data/ORFrater/mouse_GRCm39:/usr/src/mouse_GRCm39 "
       f"-v /home/ec2-user/efs-mount-point/mnt/efs0/riboseq_callers/data/ORFrater/output_batch_tmp_20240610181528_1:/usr/src/output post-processing "
       f"python3 src/main_run_post_processing.py "
       f"--experiment_name VPR_orfcalling_20240610181528_SRR5262890 "
       f"--annotation_dir /usr/src/mouse_GRCm39 "
       f"--output_dir /usr/src/output "
       f"--reference_genomes GRCm39.primary_assembly.genome.fa "
       f"--genome_annotation_prefix gencode.vM35.primary_assembly.annotation "
       f"--organism mouse "
       f"--callers price,ribotish,ribocode ")
os.system(cmd)    
