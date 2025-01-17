import numpy as np
import pandas as pd
import os


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
        feature_dfs.append(pd.read_csv(data_dir.joinpath(f'{dataset}'), sep='\t'))
    
    feature_df = pd.concat(feature_dfs)
    feature_df['orf_id'] = feature_df.apply(lambda x: f'{x.chrom_id}_{x.orf_start}_{x.orf_end}_{x.strand}_{x.exon_blocks}', axis=1)
    feature_df.set_index('orf_id', inplace=True)
    feature_df = feature_df.select_dtypes(include='number')
    feature_df = feature_df.groupby('orf_id').mean()
    
    return feature_df


class Dataset:
    def __init__(self, X, y, name, orf_ids, model=None, feature_importance_df=None):
        self.X = X
        self.y = y
        self.name = name
        self.orf_ids = orf_ids
        self.model = model
        self.feature_importance_df = feature_importance_df