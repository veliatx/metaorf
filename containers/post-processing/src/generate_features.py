"""
This program computes the by eye features for the given ORFs and (optionally) plots
the distributions of their values.
"""

import os
import pandas as pd
import pyBigWig
import math
import numpy as np
from veliadb.base import Session, Assembly

import compute_features_utils
import orf_utils


def load_big_wig_file(bigwig_pos_file, bigwig_neg_file):
    bw_pos = pyBigWig.open(bigwig_pos_file)
    bw_neg = pyBigWig.open(bigwig_neg_file)
    return {"+": bw_pos, "-": bw_neg}


def load_from_merged_orf_file(merged_orf_calls_file):
    header = []
    orfs = {}
    for line in open(merged_orf_calls_file):
        elements = line.strip("\n").split("\t")
        if not header:
            header = elements
            continue
        orf = {}
        for col_idx, col_name in enumerate(header):
            orf[col_name] = elements[col_idx]
        orf["orf_start"], orf["orf_end"] = int(orf["orf_start"]), int(orf["orf_end"])
        orf_key = (orf["chrom_id"], str(orf["orf_start"]), str(orf["orf_end"]), orf["strand"], orf["exon_blocks"])
        orfs[orf_key] = orf
    return orfs


def add_tis_transformer_feature(orfs, tis_transformer_file, callers):
    if not os.path.exists(tis_transformer_file):
        for orf in orfs.values():
            orf["tis_transformer_score"] = ""
        return

    score_map = {}
    for line in open(tis_transformer_file):
        transcript_id, chrom_id, pos, tis_transformer_score = line.strip().split()
        pos = int(pos)
        tis_transformer_score = float(tis_transformer_score)
        if transcript_id not in score_map:
            score_map[transcript_id] = {}
        score_map[transcript_id][(chrom_id, pos)] = tis_transformer_score
        
    for orf in orfs.values():
        chrom_id = orf["chrom_id"]
        if orf["strand"] == "+":
            pos = orf["orf_start"]
        elif orf["strand"] == "-":
            pos = orf["orf_end"] - 1
            
        tis_transformer_score = 0
        for caller in callers:
            transcript_id = orf[f"transcript_id_{caller}"]
            if transcript_id in score_map and (chrom_id, pos) in score_map[transcript_id]:
                tis_transformer_score = max(tis_transformer_score, score_map[transcript_id][(chrom_id, pos)])
        orf["tis_transformer_score"] = tis_transformer_score


def exon_block_to_veliadb_style(exon_blocks, exon_delimiter = '|'):
    block_sizes = []
    chrom_starts = []
    for e in exon_blocks.split(exon_delimiter):
        e_start, e_end = e.split('-')
        block_sizes.append(str(int(e_end)-int(e_start)))
        chrom_starts.append(str(int(e_start)+1)) # Add 1 to starts to convert to 1-indexed
    return ';'.join(chrom_starts), ';'.join(block_sizes)


def metaorf_call_to_veliadb_hash_string(row, chrom_id_to_assembly_id, bp_df_ids):
    if chrom_id_to_assembly_id:
        assembly_id = chrom_id_to_assembly_id[row['chrom_id']]
    else:
        assembly_id = row['chrom_id']
    chrom_starts, block_sizes = exon_block_to_veliadb_style(row['exon_blocks'])
    chrom_start = int(row['orf_start'])+1 # Add 1 to start to convert to 1-indexed
    hash_string = f"{assembly_id}_{chrom_start}_{row['orf_end']}_{row['strand']}_{chrom_starts}_{block_sizes}"
    
    if hash_string in bp_df_ids:
        return hash_string
    else:
        return ""


def add_big_prot_ids(orfs, organism):
    if organism == "human":
        ucsc_style2assembly_id = {}
        for a in Session().query(Assembly).all():
            if a.ucsc_style_name != 'na':
                if "_" in a.ucsc_style_name:
                    chrom_name = a.genbank_accession
                else:
                    chrom_name = a.ucsc_style_name        
                ucsc_style2assembly_id[chrom_name] = a.id

        bp_df = pd.read_csv('s3://velia-data-dev/VDC_004_annotation/big_prot/v0.9.3/orfset_v0.9.3_full_minlen_15_maxlen_999999999_orfs.csv.gz')
    elif organism == "mouse":
        ucsc_style2assembly_id = {}
        bp_df = pd.read_csv('s3://velia-data-dev/VDC_004_annotation/mouse_prot/orfset_mouse_v0.7.1_minlen_15_maxlen_200_orfs.csv.gz')
    bp_df_ids = set(bp_df["orf.orf_idx_str"])

    for orf in orfs.values():
        orf["bigprot_id"] = metaorf_call_to_veliadb_hash_string(orf, ucsc_style2assembly_id, bp_df_ids)


def add_qc_features(orfs, qc_file):
    with open(qc_file) as tsv_file:
        tsv = pd.read_table(tsv_file)
        for orf in orfs.values():
            orf["size_peak_frac"] = tsv['size_peak_frac'][0]
            orf["size_gini"] = tsv['size_gini'][0]
            orf["periodicity_score"] = tsv['periodicity_score'][0]


def save_into_file(orfs, coverages, drop_values, callers, feature_file_path):
    coverage_features, drop_features = None, None
    with open(feature_file_path, "w") as ofile:
        for orf_key, orf in orfs.items():
            # add header to the first line
            if coverage_features is None:
                header = "chrom_id\torf_start\torf_end\tstrand\texon_blocks\torf_sequence"
                coverage_features = [feature_name for feature_name in coverages[orf_key]]
                header = header + "\t" + "\t".join(coverage_features)
                drop_features = [feature_name for feature_name in drop_values[orf_key]]
                header = header + "\t" + "\t".join(drop_features)
                header = header + "\t" + "\t".join(callers)
                header = header + "\t" + "tis_transformer_score"
                header = header + "\t" + "size_peak_frac"
                header = header + "\t" + "size_gini"
                header = header + "\t" + "periodicity_score"
                header = header + "\t" + "bigprot_id" + "\n"
                ofile.write(header)

            line_to_write = "\t".join([str(val) for val in orf_key])
            line_to_write = line_to_write + f"\t{orf['orf_sequence']}"
            line_to_write = line_to_write + "\t" + "\t".join([str(val) for val in coverages[orf_key].values()])
            line_to_write = line_to_write + "\t" + "\t".join([str(val) for val in drop_values[orf_key].values()])

            scores = [orf[f"orf_score_{caller}"]
                        if orf[f"orf_score_{caller}"].strip() else "0" for caller in callers]
            line_to_write = line_to_write + "\t" + "\t".join(scores)
            line_to_write = line_to_write + "\t" + str(orf["tis_transformer_score"])
            line_to_write = line_to_write + "\t" + str(orf["size_peak_frac"])
            line_to_write = line_to_write + "\t" + str(orf["size_gini"])
            line_to_write = line_to_write + "\t" + str(orf["periodicity_score"])
            line_to_write = line_to_write + "\t" + orf["bigprot_id"] + "\n"
            ofile.write(line_to_write)


def get_by_eye_features(orfs_of_interest_file, experiment_name, data_dir, callers,
                        transcript_coordinates):
    bigwig = load_big_wig_file(
        bigwig_pos_file=f"{data_dir}/aligned/{experiment_name}_all.psite.bed.sorted.bam.pos.norm.bw",
        bigwig_neg_file=f"{data_dir}/aligned/{experiment_name}_all.psite.bed.sorted.bam.neg.norm.bw")
    orfs = load_from_merged_orf_file(merged_orf_calls_file=orfs_of_interest_file)

    coverages, min_cov_value = compute_features_utils.get_coverages(orfs, bigwig)
    drop_values = compute_features_utils.get_signal_drop_values(orfs, callers, bigwig, transcript_coordinates, min_cov_value)
    return orfs, coverages, drop_values


def generate_features_main(experiment_name, data_dir, transcript_bed_file,
                            tis_transformer_file, callers, organism, qc_file):
    transcript_coordinates = orf_utils.read_transcript_bed_file(transcript_bed_file)
    orfs, coverages, drop_values = get_by_eye_features(
        orfs_of_interest_file=f"{data_dir}/{experiment_name}_found_by_any_caller.csv",
        experiment_name=experiment_name,
        data_dir=data_dir,
        callers=callers,
        transcript_coordinates=transcript_coordinates)
    
    add_tis_transformer_feature(orfs, tis_transformer_file, callers)
    add_big_prot_ids(orfs, organism)
    add_qc_features(orfs, qc_file)
    save_into_file(orfs, coverages, drop_values, callers,
                   feature_file_path=f"{data_dir}/{experiment_name}_orf_features.csv")
    
