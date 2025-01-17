import os
import shutil

for aws_project in open("experiments_to_rerun_test.txt"):
    aws_project_name, rna_seq_name = aws_project.strip().split("\t")
    rna_seq_name = rna_seq_name
    aws_project_name = aws_project_name.strip("/")
    
    cmd = f"aws s3 cp --recursive --quiet s3://velia-piperuns-dev/{aws_project_name} ./data/{aws_project_name}"
    os.system(cmd)
    
    if not os.path.exists(f"./data/{aws_project_name}/output/pipeline_stats.log"):
        continue

    cmd = (f"docker run "
           f"-v /home/ec2-user/efs-mount-point/mnt/efs0/riboseq_callers/data/ORFrater/human_GRCh38_p14:/usr/src/human_GRCh38_p14 "
           f"-v /home/ec2-user/bfx-containers/riboseq_callers/post-processing/data/{aws_project_name}/output:/usr/src/output post-processing "
           f"python3 src/main_run_post_processing.py "
           f"--experiment_name {aws_project_name} "
           f"--annotation_dir /usr/src/human_GRCh38_p14 "
           f"--output_dir /usr/src/output "
           f"--reference_genomes GRCh38.p13.genome.fa "
           f"--genome_annotation_prefix veliadb_v2.fixed "
           f"--callers price,ribotish,ribocode "
           f"--transcript_list_file transcript_TPM_240419.tsv "
           f"--rna_seq_name {rna_seq_name} "
           f"--transcript_tpm_threshold 0.5 ")
    os.system(cmd)

    for filename in [f"{aws_project_name}_found_by_all_caller.bed",
                     f"{aws_project_name}_found_by_any_caller.csv",
                     f"{aws_project_name}_found_by_any_caller_all_transcripts.csv",
                     f"{aws_project_name}_orf_features.csv"]:
        cmd = f"aws s3 cp ./data/{aws_project_name}/output/{filename} s3://velia-piperuns-dev/{aws_project_name}/output/higher_thres_0_5_{filename}"
        os.system(cmd)
    shutil.rmtree(f"./data/{aws_project_name}")
    
