"""
This program computes the by eye features for the given ORFs and (optionally) plots
the distributions of their values.
"""

import os
import pandas as pd
import pyBigWig
import math
import numpy as np

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
    score_map = {}
    for line in open(tis_transfromer_file):
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
        elif orf["strand"] == "-"
            pos = orf["orf_end"] - 1
            
        tis_transformer_score = 0
        for caller in callers:
            transcript_id = orf[f"transcript_id_{caller}"]
            if transcript_id in score_map and (chrom_id, pos) in score_map[transcript_id]:
                tis_transformer_score = max(tis_transformer_score, score_map[transcript_id][(chrom_id, pos)])
        orf["tis_transformer_score"] = tis_transformer_score


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
                header = header + "\t" + "tis_transformer_score" + "\n"
                ofile.write(header)

            line_to_write = "\t".join([str(val) for val in orf_key])
            line_to_write = line_to_write + f"\t{orf['orf_sequence']}"
            line_to_write = line_to_write + "\t" + "\t".join([str(val) for val in coverages[orf_key].values()])
            line_to_write = line_to_write + "\t" + "\t".join([str(val) for val in drop_values[orf_key].values()])

            scores = [orf[f"orf_score_{caller}"]
                        if orf[f"orf_score_{caller}"].strip() else "0" for caller in callers]
            line_to_write = line_to_write + "\t" + "\t".join(scores)
            line_to_write = line_to_write + "\t" + orf["tis_transformer_score"] + "\n"
            ofile.write(line_to_write)


def get_by_eye_features(orfs_of_interest_file, experiment_name, data_dir, callers,
                        transcript_coordinates):
    bigwig = load_big_wig_file(
        bigwig_pos_file=f"{data_dir}/aligned/{experiment_name}_all.psite.bed.sorted.bam.pos.norm.bw",
        bigwig_neg_file=f"{data_dir}/aligned/{experiment_name}_all.psite.bed.sorted.bam.neg.norm.bw")
    orfs = load_from_merged_orf_file(merged_orf_calls_file=orfs_of_interest_file)

    coverages = compute_features_utils.get_coverages(orfs, bigwig)
    drop_values = compute_features_utils.get_signal_drop_values(orfs, callers, bigwig, transcript_coordinates)
    return orfs, coverages, drop_values


def generate_features_main(experiment_name, data_dir, transcript_bed_file, tis_transfromer_file, callers):
    transcript_coordinates = orf_utils.read_transcript_bed_file(transcript_bed_file)
    orfs, coverages, drop_values = get_by_eye_features(
        orfs_of_interest_file=f"{data_dir}/{experiment_name}_found_by_any_caller.csv",
        experiment_name=experiment_name,
        data_dir=data_dir,
        callers=callers,
        transcript_coordinates=transcript_coordinates)
    add_tis_transformer_feature(orfs, tis_transformer_file, callers)
    save_into_file(orfs, coverages, drop_values, callers,
                   feature_file_path=f"{data_dir}/{experiment_name}_orf_features.csv")
    
