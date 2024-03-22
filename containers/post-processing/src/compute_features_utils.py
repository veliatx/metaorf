"""This script defines functions processing coverage information for each ORF."""

import math
import numpy as np

import compute_drop_utils


def compute_frame_difference(frames, idx_1, idx_2, max_limit=2):
    if frames[idx_1] == 0:
        return max_limit
    return min(max_limit, frames[idx_2] / frames[idx_1])


def get_coverages(orfs, bigwig):
    """Computes the coverage values for the entire ORF as well as all reading frames."""

    coverages = {}
    for orf_key, orf in orfs.items():
        # get the ribo-seq signals for the ORF
        read_coverage_whole_orf = []
        for block in orf["exon_blocks"].split("|"):
            start, end = [int(elem) for elem in block.split("-")]
            read_coverage = [cov_val if (not math.isnan(cov_val)) else 0 
                            for cov_val in bigwig[orf["strand"]].values(orf["chrom_id"], start, end)]
            read_coverage_whole_orf += read_coverage

        # coverage for the entire ORF
        coverage_sum = sum(read_coverage_whole_orf)
        transcript_coverage_sum = sum([
            cov_val if (not math.isnan(cov_val)) else 0 
            for cov_val in bigwig[orf["strand"]].values(orf["chrom_id"], orf["orf_start"], orf["orf_end"])])
        coverage_mean = coverage_sum / float(len(read_coverage_whole_orf))
        coverage_std = np.std(read_coverage_whole_orf)
        
        # coverage for different frame reading frames
        frames = [0 for _ in range(3)]
        for pos_idx, cov_val in enumerate(read_coverage_whole_orf):
            frames[pos_idx % 3] += cov_val
        if orf["strand"] == "-":
            frames = frames[::-1]

        # save into map
        coverages[orf_key] = {
            "mean": coverage_mean,
            "sum": coverage_sum,
            "std": coverage_std,
            "n_reads_orf_vs_transcript": coverage_sum / float(transcript_coverage_sum) if transcript_coverage_sum != 0 else 0,
            "pos_1_vs_0": compute_frame_difference(read_coverage_whole_orf[:3], 0, 1, max_limit=10),
            "pos_2_vs_0": compute_frame_difference(read_coverage_whole_orf[:3], 0, 2, max_limit=10),
            "frames_1_vs_0": compute_frame_difference(frames, 0, 1),
            "frames_2_vs_0": compute_frame_difference(frames, 0, 2)}
    return coverages


def get_signal_drop_values(orfs, callers, bigwig, transcript_coordinates):
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
                read_coverage_whole_orf, read_coverage_stop_codon_context)
        elif orf["strand"] == "-":
            drop_features = compute_drop_utils.compute_drop_features(
                read_coverage_start_codon_context, read_coverage_start_codon_immediate_context,
                read_coverage_whole_orf[::-1], read_coverage_stop_codon_context)
        else:
            exit("Incorrect strand symbol.")

        drop[orf_key] = drop_features
    return drop