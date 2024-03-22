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


def compute_drop_mean(devidend, divisor):
    devidend_mean = compute_mean_signal(devidend)
    divisor_mean = compute_mean_signal(divisor)

    if divisor_mean == 0:
        return -3
    elif devidend_mean == 0:
        return 3
    else:
        return min(3, max(-3, -math.log10(devidend_mean/divisor_mean)))


def compute_drop_max(devidend, divisor):
    devidend_max = compute_max_signal(devidend)
    divisor_max = compute_max_signal(divisor)
    
    if divisor_max == 0:
        return -3
    elif devidend_max == 0:
        return 3
    else:
        return min(3, max(-3, -math.log10(devidend_max/divisor_max)))
    

def compute_drop_values(devidend, divisor):
    return compute_drop_mean(devidend, divisor), compute_drop_max(devidend, divisor)


def get_distance_of_nearest_peak(read_coverage, cmp_signal):
    distance = len(read_coverage)
    for pos, signal in enumerate(read_coverage):
        if signal >= cmp_signal:
            distance = pos+1
            break
    return distance


def compute_nearest_peak(read_coverage_start_codon_context, read_coverage_whole_orf):    
    start_codon_peak = read_coverage_whole_orf[0]
    dist_neg_100 = get_distance_of_nearest_peak(read_coverage_start_codon_context[::-1], start_codon_peak)
    dist_neg_150 = get_distance_of_nearest_peak(read_coverage_start_codon_context[::-1], start_codon_peak*1.5)
    dist_pos_100 = get_distance_of_nearest_peak(read_coverage_whole_orf[1:], start_codon_peak)
    dist_pos_150 = get_distance_of_nearest_peak(read_coverage_whole_orf[1:], start_codon_peak*1.5)
    return dist_neg_100, dist_neg_150, dist_pos_100, dist_pos_150


def compute_drop_features(read_coverage_start_codon_context, read_coverage_start_codon_immediate_context,
                          read_coverage_whole_orf, read_coverage_stop_codon_context):
    # start codon
    five_utr_vs_cds_mean, five_utr_vs_cds_max = compute_drop_values(
        read_coverage_start_codon_context,
        read_coverage_start_codon_immediate_context+read_coverage_whole_orf[:60])
    five_utr_vs_start_codon_mean, five_utr_vs_start_codon_max = compute_drop_values(
        read_coverage_start_codon_context,
        read_coverage_start_codon_immediate_context+read_coverage_whole_orf[:6])
    cds_utr_vs_start_codon_mean, cds_utr_vs_start_codon_max = compute_drop_values(
        read_coverage_whole_orf[6:60],
        read_coverage_start_codon_immediate_context+read_coverage_whole_orf[:6])
        
    # stop codon
    three_utr_vs_cds_mean, three_utr_vs_cds_max = compute_drop_values(
        read_coverage_stop_codon_context+read_coverage_whole_orf[-3:],
        read_coverage_whole_orf[-60:-3])
    three_utr_vs_stop_codon_mean, three_utr_vs_stop_codon_max = compute_drop_values(
        read_coverage_stop_codon_context+read_coverage_whole_orf[-3:],
        read_coverage_whole_orf[-6:-3])
    cds_utr_vs_stop_codon_mean, cds_utr_vs_stop_codon_max = compute_drop_values(
        read_coverage_whole_orf[-60:-6],
        read_coverage_whole_orf[-6:-3])
    
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
           }