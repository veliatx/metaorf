import numpy as np
import pandas as pd
import math
import os
import pickle


from metaorf.modeling.ensemble import Dataset


def download_feature_set(experiment_name, local_folder):
    feature_path = f"s3://velia-piperuns-dev/{experiment_name}/output/{experiment_name}_orf_features.csv"
    cmd = f"aws s3 cp {feature_path} {local_folder}"
    os.system(cmd)


def remove_local_copy(experiment_name, local_folder):    
    os.remove(f"{local_folder}/{experiment_name}_orf_features.csv")


def read_labels(data_path, labels, label_columns):
    for index, row in pd.read_csv(data_path, sep='\t').iterrows():
        block_sizes = [int(val) for val in row["blockSizes"].split(",")[:-1]]
        block_starts = [int(val) for val in row["blockStarts"].split(",")[:-1]]
        blocks = []
        for block_idx in range(row["blockCount"]):
            block_start = row["chromStart"] + block_starts[block_idx]
            block_end = row["chromStart"] + block_starts[block_idx] + block_sizes[block_idx]
            blocks.append(f"{block_start}-{block_end}")
        orf_key = (row["chrom"], str(row["chromStart"]), str(row["chromEnd"]),
                   row["strand"], "|".join(blocks))
        
        label = {}
        for col in label_columns:
            label[col] = row[col]
        
        if orf_key in labels:
            print(orf_key)
        else:
            labels[orf_key] = label


def value_post_process(value, col_name):
    if col_name == "price":
        return max(17, value)
    elif col_name.startswith("dist_"):
        return -math.log10(value)
    else:
        return value


def read_features(data_path, features, labels, label_name):
    feature_set = pd.read_csv(data_path, sep='\t')
    columns_to_exclude = ["chrom_id", "orf_start", "orf_end", "strand", "exon_blocks", "orf_sequence"]
    feature_columns = [col_name for col_name in feature_set.columns if col_name not in columns_to_exclude]
    
    count_pos, count_neg = 0, 0
    for index, row in feature_set.iterrows():
        orf_key = (row["chrom_id"], str(row["orf_start"]), str(row["orf_end"]),
                   row["strand"], row["exon_blocks"])
        if orf_key not in labels:
            continue
            
        label = labels[orf_key][label_name]
        if label >= 2:
            binary_label = 1.0
            count_pos += 1
        else:
            binary_label = 0.0
            count_neg += 1

        current_features = []
        for col_name in feature_columns:
            current_features.append(float(row[col_name]))
            
        features.append((orf_key, current_features, binary_label))     
    print(count_pos, count_neg)
    

def read_features_no_labels(data_path, features):
    feature_set = pd.read_csv(data_path, sep='\t')
    columns_to_exclude = ["chrom_id", "orf_start", "orf_end", "strand", "exon_blocks", "orf_sequence"]
    feature_columns = [col_name for col_name in feature_set.columns if col_name not in columns_to_exclude]
    
    for index, row in feature_set.iterrows():
        orf_key = (row["chrom_id"], str(row["orf_start"]), str(row["orf_end"]),
                   row["strand"], row["exon_blocks"])
        current_features = []
        for col_name in feature_columns:
            current_features.append(float(row[col_name])) 
        features.append((orf_key, current_features, 0))     
    
    
def get_dataset(features, validation_chroms, test_chroms):
    training_data_x, training_data_y = [], []
    validation_data_x, validation_data_y = [], []
    test_data_x, test_data_y = [], []
    for orf_key, current_features, label in features:
        if orf_key[0] in validation_chroms:
            validation_data_x.append(current_features)
            validation_data_y.append(label)
        elif orf_key[0] in test_chroms:
            test_data_x.append(current_features)
            test_data_y.append(label)
        else:
            training_data_x.append(current_features)
            training_data_y.append(label)

    return (training_data_x, training_data_y,
            validation_data_x, validation_data_y,
            test_data_x, test_data_y)


def generate_orf_id(row):
    """
    """
    
    block_sizes = [int(val) for val in row["blockSizes"].split(",")[:-1]]
    block_starts = [int(val) for val in row["blockStarts"].split(",")[:-1]]
    blocks = []
    for block_idx in range(row["blockCount"]):
        block_start = row["chromStart"] + block_starts[block_idx]
        block_end = row["chromStart"] + block_starts[block_idx] + block_sizes[block_idx]
        blocks.append(f"{block_start}-{block_end}")
    
    return f'{row.chrom}_{row.chromStart}_{row.chromEnd}_{row.strand}_{"|".join(blocks)}'


def load_features(data_dir, datasets):
    """
    """

    feature_dfs = []
    for dataset in datasets:
        tmp_df = pd.read_csv(data_dir.joinpath(f'{dataset}'), sep='\t')
        tmp_df['dataset'] = '-'.join(dataset.split('_')[:-2])
        feature_dfs.append(tmp_df)
    
    feature_df = pd.concat(feature_dfs)
    feature_df['orf_idx_str'] =  feature_df.apply(lambda x: f'{x.chrom_id}_{x.orf_start}_{x.orf_end}_{x.strand}_{x.exon_blocks}', axis=1)

    return feature_df


def load_truth_datasets(truth_df, data_dir, overwrite, dataset_names=['iPSC', 'MB1', 'Gaertner']):
    """
    """
    datasets = {}

    X_dfs = []
    y_arrays = []
    orf_ids = []

    if not overwrite and data_dir.joinpath('datasets.pkl').exists():
        with open(data_dir.joinpath('datasets.pkl'), 'rb') as file:
            datasets = pickle.load(file)
            return datasets

    for dataset_name in dataset_names:
        dataset_file_names = []
        with open(data_dir.joinpath(f'{dataset_name}.txt'), 'r') as infile:
            for line in infile.readlines():
                exp_name = line.rstrip('\n')
                dataset_file_names.append(f"{exp_name}_orf_features.csv")

        truth_label = f'score.({dataset_name})'
        feature_df = load_features(data_dir, dataset_file_names)

        #drop_cols = ['chrom_id', 'orf_start', 'orf_end', 'strand', 'exon_blocks', 'dataset',
        #             'orf_sequence']
        drop_cols = ['orf_start', 'orf_end', 'strand', 'exon_blocks', 'dataset', 'orf_sequence']
        feature_df = feature_df.drop(columns=drop_cols)
        #feature_df = feature_df.groupby('orf_idx_str').mean()
        
        merged_df = feature_df.merge(truth_df[['orf_id', truth_label]], left_on='orf_idx_str', right_on='orf_id')
        merged_df.reset_index(inplace=True)

        y = merged_df[truth_label].copy()
        y[y > 0] = 1
        y = y.values

        drop_cols = ['orf_id', 'orf_idx_str', truth_label, 'index']
        numeric_feat_df = merged_df.drop(columns=drop_cols)
        X = numeric_feat_df

        X_dfs.append(X)
        y_arrays.append(y)
        orf_ids.append(merged_df['orf_id'])

        datasets[dataset_name] = Dataset(X, y, dataset_name, merged_df['orf_id'])

    datasets['all'] = Dataset(pd.concat(X_dfs), np.concatenate(y_arrays), 'all',  pd.concat(orf_ids))
    datasets['all'].X.reset_index(inplace=True, drop=True)

    with open(data_dir.joinpath('datasets.pkl'), 'wb') as file:
        pickle.dump(datasets, file)

    return datasets


def download_feature_set(experiment_name, local_path):
    """
    """
    orf_path = f"s3://velia-piperuns-dev/{experiment_name}/output/{experiment_name}_found_by_any_caller.csv"
    local_path = f"{local_folder}/{experiment_name}_found_by_any_caller.csv"
    cmd = f"aws s3 cp {orf_path} {local_path}"
    print(cmd)
    os.system(cmd)
    return local_path


def load_orfs(experiment_name, local_folder, merged_orfs_info, orf_callers):
    """
    """
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


def merge_orfs(local_folder="./data", orf_callers=["price", "ribotish", "ribocode"]):
    """
    """

    merged_orfs_info = {}
    count = 0
    for experiment_name in open(f"{local_folder}/experiments.txt"):
        experiment_name = experiment_name.split("PRE")[1].strip()[:-1]
        count += load_orfs(experiment_name, local_folder, merged_orfs_info, orf_callers)

    pd.DataFrame(merged_orfs_info.values()).reset_index().to_csv(
        f"{local_folder}/merged_orfs_found_by_any_caller.csv",
        sep="\t",
        index=False) 
