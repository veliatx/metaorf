"""This script defines functions compute the signal drop values for each ORF."""
import math


def get_context_signals(transcript, orf, bigwig):
    """Gets the converage values for the context positions for the given ORF."""

    # left context
    read_coverage_left_context = []
    for block_idx in range(transcript["block_count"]):
        block_start = transcript["start"] + transcript["block_starts"][block_idx]
        block_end = transcript["start"] + transcript["block_starts"][block_idx] + transcript["block_sizes"][block_idx]
        
        if block_end < orf["orf_start"]:
            read_coverage_left_context += bigwig[orf["strand"]].values(orf["chrom_id"], block_start, block_end)
        elif block_end == orf["orf_start"]:
            exit("block_end cannot be eqaul to orf_start")
        else:
            if block_start != orf["orf_start"]:
                read_coverage_left_context += bigwig[orf["strand"]].values(orf["chrom_id"], block_start, orf["orf_start"])
            break
    read_coverage_left_context = [cov_val if (not math.isnan(cov_val)) else 0 for cov_val in read_coverage_left_context]
    
    # ORF
    read_coverage_orf = []
    for block_idx in range(transcript["block_count"]):
        block_start = transcript["start"] + transcript["block_starts"][block_idx]
        block_end = transcript["start"] + transcript["block_starts"][block_idx] + transcript["block_sizes"][block_idx]
        
        if block_end <= orf["orf_start"]:
            continue
        else:
            if block_end < orf["orf_end"]:
                read_coverage_orf += bigwig[orf["strand"]].values(orf["chrom_id"],
                                                                  max(orf["orf_start"], block_start),
                                                                  block_end)
            else:
                read_coverage_orf += bigwig[orf["strand"]].values(orf["chrom_id"],
                                                                  max(orf["orf_start"], block_start),
                                                                  orf["orf_end"])
                break
    read_coverage_orf = [cov_val if (not math.isnan(cov_val)) else 0 for cov_val in read_coverage_orf]
    
    # right context
    read_coverage_right_context = []
    for block_idx in range(transcript["block_count"]):
        block_start = transcript["start"] + transcript["block_starts"][block_idx]
        block_end = transcript["start"] + transcript["block_starts"][block_idx] + transcript["block_sizes"][block_idx]
        
        if block_end <= orf["orf_end"]:
            continue
        else:
            if block_start < orf["orf_end"]:
                read_coverage_right_context += bigwig[orf["strand"]].values(orf["chrom_id"], orf["orf_end"], block_end)
            elif block_start == orf["orf_end"]:
                exit("block_start cannot be eqaul to orf_end")
            else:
                read_coverage_right_context += bigwig[orf["strand"]].values(orf["chrom_id"], block_start, block_end)
    read_coverage_right_context = [cov_val if (not math.isnan(cov_val)) else 0 for cov_val in read_coverage_right_context]
        
    if orf["strand"] == "+":
        start_codon_context = read_coverage_left_context[-60:-3]
        start_codon_immediate_context = read_coverage_left_context[-3:]
        stop_codon_context = read_coverage_right_context[:60]
    elif orf["strand"] == "-":
        start_codon_context = read_coverage_right_context[3:60][::-1]
        start_codon_immediate_context = read_coverage_right_context[:3][::-1]
        stop_codon_context = read_coverage_left_context[-60:][::-1]
    else:
        exit("Unexpected strand symbol.")
    return start_codon_context, start_codon_immediate_context, stop_codon_context, read_coverage_orf

    
def compute_mean_signal(signals):
    if not signals:
        return 0
    else:
        return sum(signals) / float(len(signals))


def compute_max_signal(signals):
    if not signals:
        return 0
    else:
        return max(signals)


def compute_frame_difference(frames, min_cov_value, max_limit=2):
    if max(frames) < min_cov_value*5:
        return 0
    elif frames[0] == 0:
        return max_limit
    else:
        return min(max_limit, (frames[1] + frames[2]) / frames[0])


def compute_periodicity(read_coverage, min_cov_value, start_context=False):
    n_units = len(read_coverage) // 3
    if start_context:
        read_coverage = read_coverage[-3*n_units:]
    else:
        read_coverage = read_coverage[:3*n_units]

    frame_coverage = [0, 0, 0]
    in_frame_bins = 0
    n_bins = 0
    for idx in range(n_units):
        n_bins += 1
        for offset in range(3):
            frame_coverage[offset] += read_coverage[idx*3+offset]
        if (read_coverage[idx*3] > read_coverage[idx*3+1]) and (read_coverage[idx*3] > read_coverage[idx*3+2]):
            in_frame_bins += 1

    if n_bins == 0:
        in_frame_bins_ratio = 0
    else:
        in_frame_bins_ratio = in_frame_bins / float(n_bins)
    return compute_frame_difference(frame_coverage, min_cov_value), in_frame_bins_ratio


def compute_drop_mean(devidend, divisor, min_cov_value):
    devidend_mean = compute_mean_signal(devidend)
    divisor_mean = compute_mean_signal(divisor)

    if divisor_mean < min_cov_value / 2.0:
        return 0
    elif devidend_mean == 0:
        return 3
    else:
        return min(3, max(-3, -math.log10(devidend_mean/divisor_mean)))


def compute_drop_max(devidend, divisor, min_cov_value):
    devidend_max = compute_max_signal(devidend)
    divisor_max = compute_max_signal(divisor)
    
    if divisor_max < min_cov_value:
        return 0
    elif devidend_max == 0:
        return 3
    else:
        return min(3, max(-3, -math.log10(devidend_max/divisor_max)))
    

def compute_drop_values(devidend, divisor, min_cov_value):
    return compute_drop_mean(devidend, divisor, min_cov_value), compute_drop_max(devidend, divisor, min_cov_value)


def get_distance_of_nearest_peak(read_coverage, cmp_signal, default_max_limit=60):
    distance = default_max_limit
    for pos, signal in enumerate(read_coverage):
        if signal >= cmp_signal:
            distance = pos+1
            break
    distance = -math.log10(min(default_max_limit, distance))
    return distance


def compute_nearest_peak(read_coverage_start_codon_context, read_coverage_whole_orf):    
    start_codon_peak = read_coverage_whole_orf[0]
    dist_neg_100 = get_distance_of_nearest_peak(read_coverage_start_codon_context[::-1], start_codon_peak)
    dist_neg_150 = get_distance_of_nearest_peak(read_coverage_start_codon_context[::-1], start_codon_peak*1.5)
    dist_pos_100 = get_distance_of_nearest_peak(read_coverage_whole_orf[1:60], start_codon_peak)
    dist_pos_150 = get_distance_of_nearest_peak(read_coverage_whole_orf[1:60], start_codon_peak*1.5)
    return dist_neg_100, dist_neg_150, dist_pos_100, dist_pos_150


def compute_drop_features(read_coverage_start_codon_context, read_coverage_start_codon_immediate_context,
                          read_coverage_whole_orf, read_coverage_stop_codon_context, min_cov_value):
    # start codon
    five_utr_vs_cds_mean, five_utr_vs_cds_max = compute_drop_values(
        read_coverage_start_codon_context,
        read_coverage_start_codon_immediate_context+read_coverage_whole_orf[:60],
        min_cov_value)
    five_utr_vs_start_codon_mean, five_utr_vs_start_codon_max = compute_drop_values(
        read_coverage_start_codon_context,
        read_coverage_start_codon_immediate_context+read_coverage_whole_orf[:6],
        min_cov_value)
    cds_utr_vs_start_codon_mean, cds_utr_vs_start_codon_max = compute_drop_values(
        read_coverage_whole_orf[6:60],
        read_coverage_start_codon_immediate_context+read_coverage_whole_orf[:6],
        min_cov_value)

    five_utr_periodicity, five_utr_in_frame_bins = compute_periodicity(
        read_coverage_start_codon_context[30:]+read_coverage_start_codon_immediate_context,
        min_cov_value)
        
    # stop codon
    three_utr_vs_cds_mean, three_utr_vs_cds_max = compute_drop_values(
        read_coverage_stop_codon_context+read_coverage_whole_orf[-3:],
        read_coverage_whole_orf[-60:-3],
        min_cov_value)
    three_utr_vs_stop_codon_mean, three_utr_vs_stop_codon_max = compute_drop_values(
        read_coverage_stop_codon_context+read_coverage_whole_orf[-3:],
        read_coverage_whole_orf[-6:-3],
        min_cov_value)
    cds_utr_vs_stop_codon_mean, cds_utr_vs_stop_codon_max = compute_drop_values(
        read_coverage_whole_orf[-60:-6],
        read_coverage_whole_orf[-6:-3],
        min_cov_value)
    three_utr_periodicity, three_utr_in_frame_bins = compute_periodicity(
        read_coverage_whole_orf[-3:] + read_coverage_stop_codon_context[:27],
        min_cov_value)
    
    # nearest peak
    dist_neg_100, dist_neg_150, dist_pos_100, dist_pos_150 = compute_nearest_peak(
        read_coverage_start_codon_context+read_coverage_start_codon_immediate_context,
        read_coverage_whole_orf)
    
    return {"five_utr_vs_cds_mean": five_utr_vs_cds_mean,
            "five_utr_vs_cds_max": five_utr_vs_cds_max,
            "five_utr_vs_start_codon_mean": five_utr_vs_start_codon_mean,
            "five_utr_vs_start_codon_max": five_utr_vs_start_codon_max,
            "cds_utr_vs_start_codon_mean": cds_utr_vs_start_codon_mean,
            "cds_utr_vs_start_codon_max": cds_utr_vs_start_codon_max,
            "three_utr_vs_cds_mean": three_utr_vs_cds_mean,
            "three_utr_vs_cds_max": three_utr_vs_cds_max,
            "three_utr_vs_stop_codon_mean": three_utr_vs_stop_codon_mean,
            "three_utr_vs_stop_codon_max": three_utr_vs_stop_codon_max,
            "cds_utr_vs_stop_codon_mean": cds_utr_vs_stop_codon_mean,
            "cds_utr_vs_stop_codon_max": cds_utr_vs_stop_codon_max,
            "dist_neg_100": dist_neg_100,
            "dist_neg_150": dist_neg_150,
            "dist_pos_100": dist_pos_100,
            "dist_pos_150": dist_pos_150,
            "five_utr_periodicity": five_utr_periodicity,
            "five_utr_in_frame_bins": five_utr_in_frame_bins,
            "three_utr_periodicity": three_utr_periodicity,
            "three_utr_in_frame_bins": three_utr_in_frame_bins,
           }
