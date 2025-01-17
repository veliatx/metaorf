import pandas as pd
import os

def download_feature_set(experiment_name, local_path):
    orf_path = f"s3://velia-piperuns-dev/{experiment_name}/output/{experiment_name}_found_by_any_caller.csv"
    local_path = f"{local_folder}/{experiment_name}_found_by_any_caller.csv"
    cmd = f"aws s3 cp {orf_path} {local_path}"
    print(cmd)
    os.system(cmd)
    return local_path

def load_orfs(experiment_name, local_folder, merged_orfs_info, orf_callers):
    data_path = download_feature_set(experiment_name, local_folder)
    if not os.path.exists(data_path):
        return 0
    
    orfs_info = {}
    for index, row in pd.read_csv(data_path, sep='\t').iterrows():
        orf_key = (row["chrom_id"], row["orf_start"], row["orf_end"], row["strand"], row["exon_blocks"])
        if orf_key not in merged_orfs_info:
            merged_orfs_info[orf_key] = row
            merged_orfs_info[orf_key] = pd.concat([merged_orfs_info[orf_key],
                                                   pd.Series({"#called": 1,
                                                              "cell_types": experiment_name})])
        else:
            merged_orfs_info[orf_key]["#called"] += 1
            merged_orfs_info[orf_key]["cell_types"] += ("," + experiment_name)
            for orf_caller in orf_callers:
                if (type(row[f"orf_score_{orf_caller}"]) == type("")
                    and (not row[f"orf_score_{orf_caller}"].strip())):
                    continue
                else:
                    if ((type(merged_orfs_info[orf_key][f"orf_score_{orf_caller}"]) == type("") and
                         (not merged_orfs_info[orf_key][f"orf_score_{orf_caller}"].strip()))
                        or (merged_orfs_info[orf_key][f"orf_score_{orf_caller}"] is None)
                        or (float(merged_orfs_info[orf_key][f"orf_score_{orf_caller}"]) <= float(row[f"orf_score_{orf_caller}"]))):
                        merged_orfs_info[orf_key][f"orf_score_{orf_caller}"] = row[f"orf_score_{orf_caller}"]
                        merged_orfs_info[orf_key][f"ORF_type_{orf_caller}"] = row[f"ORF_type_{orf_caller}"]
                        merged_orfs_info[orf_key][f"gene_id_{orf_caller}"] = row[f"gene_id_{orf_caller}"]
                        merged_orfs_info[orf_key][f"transcript_id_{orf_caller}"] = row[f"transcript_id_{orf_caller}"]    
    os.remove(data_path)
    return 1

if __name__ == "__main__":
    local_folder= "./data"
    orf_callers = ["price", "ribotish", "ribocode"]

    merged_orfs_info = {}
    count = 0
    for experiment_name in open(f"{local_folder}/experiments.txt"):
        experiment_name = experiment_name.split("PRE")[1].strip()[:-1]
        count += load_orfs(experiment_name, local_folder, merged_orfs_info, orf_callers)
    print(count) 
    pd.DataFrame(merged_orfs_info.values()).reset_index().to_csv(
        f"{local_folder}/merged_orfs_found_by_any_caller.csv",
        sep="\t",
        index=False) 