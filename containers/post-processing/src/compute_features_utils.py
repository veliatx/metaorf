"""This script defines functions processing coverage information for each ORF."""

import math
import numpy as np

import compute_drop_utils


def compute_frame_difference(frames, idx_1, idx_2, min_cov_value=0, max_limit=2):
    if frames[idx_1] <= min_cov_value and frames[idx_2] <= min_cov_value:
        return 1
    elif frames[idx_1] <= min_cov_value:
        return max_limit
    else:
        return min(max_limit, frames[idx_2] / frames[idx_1])


def compute_empty_codon_number(read_coverage, min_cov_value):
    n_empty = 0
    for codon_idx in range(len(read_coverage)//3):
        if sum(read_coverage[3*codon_idx:3*(codon_idx+1)]) == 0:
            n_empty += 1
    return -math.log10(max(n_empty, 1))


def get_empty_length_fixed_position(read_coverage, min_cov_value, start_idx):
    current_violations = 0
    current_end = start_idx
    for idx in range(start_idx, len(read_coverage)):
        if read_coverage[idx] > min_cov_value:
            current_violations += 1
            if current_violations > 1 and (current_violations / float(current_end - start_idx)) > 0.05:
                break
        current_end += 1
    return current_end - start_idx


def compute_empty_length(read_coverage, min_cov_value):
    longest_length = -np.inf
    for idx in range(len(read_coverage)):
        current_longest = get_empty_length_fixed_position(read_coverage, min_cov_value, idx)
        longest_length = max(longest_length, current_longest)
    return -math.log10(max(longest_length, 1))


def get_coverages(orfs, bigwig):
    """Computes the coverage values for the entire ORF as well as all reading frames."""

    min_cov_value = np.inf
    for orf in orfs.values():
        for block in orf["exon_blocks"].split("|"):
            start, end = [int(elem) for elem in block.split("-")]
            read_coverage = [cov_val if (not math.isnan(cov_val)) else 0 
                            for cov_val in bigwig[orf["strand"]].values(orf["chrom_id"], start, end)]

            for cov_val in read_coverage:
                if cov_val > 0 and cov_val < min_cov_value:
                    min_cov_value = cov_val

    coverages = {}
    for orf_key, orf in orfs.items():
        sum_genome_coverage = 0
        for cov_val in bigwig[orf["strand"]].values(orf["chrom_id"], orf["orf_start"], orf["orf_end"]):
            if math.isnan(cov_val):
                cov_val = 0
            sum_genome_coverage += cov_val

        # get the ribo-seq signals for the ORF
        read_coverage_whole_orf = []
        for block in orf["exon_blocks"].split("|"):
            start, end = [int(elem) for elem in block.split("-")]
            read_coverage = [cov_val if (not math.isnan(cov_val)) else 0 
                            for cov_val in bigwig[orf["strand"]].values(orf["chrom_id"], start, end)]
            read_coverage_whole_orf += read_coverage

        # coverage for the entire ORF
        coverage_sum = sum(read_coverage_whole_orf)
        coverage_mean = coverage_sum / float(len(read_coverage_whole_orf))
        coverage_std = np.std(read_coverage_whole_orf)

        if orf["strand"] == "-":
            read_coverage_whole_orf = read_coverage_whole_orf[::-1]
        
        # coverage for different frame reading frames
        frames = [0 for _ in range(3)]
        for pos_idx, cov_val in enumerate(read_coverage_whole_orf):
            frames[pos_idx % 3] += cov_val
        frames = [frames[idx]*3.0/float(len(read_coverage_whole_orf)) for idx in range(3)]

        frames_first_60 = [0 for _ in range(3)]
        for pos_idx, cov_val in enumerate(read_coverage_whole_orf[:60]):
            frames_first_60[pos_idx % 3] += cov_val
        frames_first_60 = [frames_first_60[idx]*3.0/float(len(read_coverage_whole_orf[:60])) for idx in range(3)]

        frames_last_60 = [0 for _ in range(3)]
        for pos_idx, cov_val in enumerate(read_coverage_whole_orf[-63:-3]):
            frames_last_60[pos_idx % 3] += cov_val
        frames_last_60 = [frames_last_60[idx]*3.0/float(len(read_coverage_whole_orf[-63:-3])) for idx in range(3)]

        # save into map
        coverages[orf_key] = {
            "mean": -math.log10(max(coverage_mean, 1e-10)),
            "sum": -math.log10(max(coverage_sum, 1e-10)),
            "std": coverage_std,
            "n_reads_orf_vs_genome": coverage_sum / float(sum_genome_coverage) if sum_genome_coverage != 0 else 0,
            "pos_1_vs_0": compute_frame_difference(read_coverage_whole_orf[:3], 0, 1, max_limit=10),
            "pos_2_vs_0": compute_frame_difference(read_coverage_whole_orf[:3], 0, 2, max_limit=10),
            "frames_1_vs_0": compute_frame_difference(frames, 0, 1),
            "frames_2_vs_0": compute_frame_difference(frames, 0, 2),
            "periodicity_first_60_1_vs_0": compute_frame_difference(frames_first_60, 0, 1, min_cov_value/6.0),
            "periodicity_first_60_2_vs_0": compute_frame_difference(frames_first_60, 0, 2, min_cov_value/6.0),
            "periodicity_last_60_1_vs_0": compute_frame_difference(frames_last_60, 0, 1, min_cov_value/6.0),
            "periodicity_last_60_2_vs_0": compute_frame_difference(frames_last_60, 0, 2, min_cov_value/6.0),
            "n_empty_codons": compute_empty_codon_number(read_coverage_whole_orf, min_cov_value),
            "longest_empty_length_whole": compute_empty_length(read_coverage_whole_orf, min_cov_value),
            "longest_empty_length_first_30": compute_empty_length(read_coverage_whole_orf[:30], min_cov_value),
            "longest_empty_length_last_30": compute_empty_length(read_coverage_whole_orf[-33:-3], min_cov_value),}
    return coverages, min_cov_value


def get_signal_drop_values(orfs, callers, bigwig, transcript_coordinates, min_cov_value):
    drop = {}
    for orf_key, orf in orfs.items():
        # get the ribo-seq signals for the ORF
        read_coverage_whole_orf = []
        for block in orf["exon_blocks"].split("|"):
            start, end = [int(elem) for elem in block.split("-")]
            read_coverage = [cov_val if (not math.isnan(cov_val)) else 0 
                            for cov_val in bigwig[orf["strand"]].values(orf["chrom_id"], start, end)]
            read_coverage_whole_orf += read_coverage
        
        # get the transcript_id of this ORF - first non empty transcript_id annotation
        transcript_id = [orf[f"transcript_id_{caller}"] for caller in callers if orf[f"transcript_id_{caller}"].strip()][0]
        
                
        # get the ribo-seq signals around the ORF
        (read_coverage_start_codon_context,
        read_coverage_start_codon_immediate_context,
        read_coverage_stop_codon_context,
        read_coverage_orf_qc) = compute_drop_utils.get_context_signals(transcript_coordinates[transcript_id], orf, bigwig)
        
        # sanity check
        if read_coverage_whole_orf != read_coverage_orf_qc:
            exit((read_coverage_whole_orf, read_coverage_orf_qc))
        
        # compute drop values
        if orf["strand"] == "+":
            drop_features = compute_drop_utils.compute_drop_features(
                read_coverage_start_codon_context, read_coverage_start_codon_immediate_context,
                read_coverage_whole_orf, read_coverage_stop_codon_context, min_cov_value)
        elif orf["strand"] == "-":
            drop_features = compute_drop_utils.compute_drop_features(
                read_coverage_start_codon_context, read_coverage_start_codon_immediate_context,
                read_coverage_whole_orf[::-1], read_coverage_stop_codon_context, min_cov_value)
        else:
            exit("Incorrect strand symbol.")

        drop[orf_key] = drop_features
    return drop
