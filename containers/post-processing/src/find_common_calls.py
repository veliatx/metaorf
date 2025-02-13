"""This program identifies the overlaps among the following orf-callers:
    RiboCode
    RiboTish
    PRICE
"""

import os
import pandas as pd
import orf_utils


def combine_info(info_list, idx):
    return ",".join([str(info[idx]) for info in info_list])

def find_overlaps(all_orf_calls, orf_callers, transcript_list):
    """Identifies the overlaps among calls from different ORF callers."""

    merged_orf_calls = {}
    merged_orf_calls_full_list = {}
    for orf_caller in orf_callers:
        orf_calls = all_orf_calls[orf_caller]
        for key in orf_calls:
            orf_score_list = sorted(orf_calls[key], key=lambda x:x[3], reverse=True)
            if transcript_list is not None:
                orf_score_list_filtered = [orf_info for orf_info in orf_score_list if orf_info[1] in transcript_list]
            else:
                orf_score_list_filtered = orf_score_list

            if orf_score_list_filtered:
                if key not in merged_orf_calls:
                    merged_orf_calls[key] = {}
                merged_orf_calls[key][orf_caller] = {
                    "gene_id": orf_score_list_filtered[0][0],
                    "transcript_id": orf_score_list_filtered[0][1],
                    "orf_type": orf_score_list_filtered[0][2],
                    "score": orf_score_list_filtered[0][3]}
            
            if key not in merged_orf_calls_full_list:
                merged_orf_calls_full_list[key] = {}
            merged_orf_calls_full_list[key][orf_caller] = {
                "gene_id": combine_info(orf_score_list_filtered, 0),
                "transcript_id": combine_info(orf_score_list_filtered, 1),
                "orf_type": combine_info(orf_score_list_filtered, 2),
                "score": combine_info(orf_score_list_filtered, 3)}
    return merged_orf_calls, merged_orf_calls_full_list


def create_bed_line(chrom_id, orf_start, orf_end, strand, exon_blocks,
                   gene_id, transcript_id):
    bed_line = (f"{chrom_id}\t{orf_start}\t{orf_end}\t{gene_id+';'+transcript_id}\t0\t{strand}"
                f"\t{orf_start}\t{orf_end}\t0\t{len(exon_blocks)}")
    starts, sizes = [], []
    for block in exon_blocks.split("|"):
        start, end = [int(elem) for elem in block.split("-")]
        sizes.append(str(end-start))
        starts.append(str(start-orf_start))
    bed_line += f"\t{','.join(sizes)}\t{','.join(starts)}\n"
    return bed_line


def write_overlaps_into_file(merged_orf_calls, overlap_csv_file, orf_callers):
    """For each cell_type/replicate/experiment, writes the merged ORF calls into file."""

    with open(overlap_csv_file, "w") as overlap_output:
        # add the table header
        header = ("chrom_id\torf_start\torf_end\tstrand\torf_sequence\texon_blocks")
        for orf_caller in orf_callers:
            header += (f"\torf_score_{orf_caller}\tORF_type_{orf_caller}\t"
                       f"gene_id_{orf_caller}\ttranscript_id_{orf_caller}")
        overlap_output.write(header+"\n")

        for key in merged_orf_calls:
            (chrom_id, orf_start, orf_end, strand, orf_sequence, exon_blocks) = key
            if len(orf_sequence) > 600:
                continue

            # write into file
            line_to_write = (f"{chrom_id}\t{orf_start}\t{orf_end}\t{strand}\t{orf_sequence}\t{exon_blocks}")
            for orf_caller in orf_callers:
                if orf_caller in merged_orf_calls[key]:
                    orf_call_info = merged_orf_calls[key][orf_caller]
                    line_to_write += (f"\t{orf_call_info['score']}\t{orf_call_info['orf_type']}"
                                      f"\t{orf_call_info['gene_id']}\t{orf_call_info['transcript_id']}")
                else:
                    line_to_write += (f"\t \t \t \t ")
            overlap_output.write(line_to_write+"\n")
    return


def write_overlaps_into_bed(merged_orf_calls, overlap_bed_file, orf_callers):
    with open(overlap_bed_file, "w") as overlap_bed_output:
        for key in merged_orf_calls:
            (chrom_id, orf_start, orf_end, strand, orf_sequence, exon_blocks) = key
            if len(orf_sequence) > 600:
                continue

            # write into file
            caller_count = sum([1 for orf_caller in orf_callers if orf_caller in merged_orf_calls[key]])
            if caller_count == len(orf_callers):
                bed_file_line = create_bed_line(
                    chrom_id, orf_start, orf_end, strand, exon_blocks,
                    merged_orf_calls[key][orf_callers[0]]['gene_id'], merged_orf_calls[key][orf_callers[0]]['transcript_id'])
                overlap_bed_output.write(bed_file_line)
    return


def read_transcript_list(transcript_list_file, rna_seq_name, threshold):
    if rna_seq_name is None or rna_seq_name == "":
        return None
    
    transcript_table = pd.read_csv(transcript_list_file, sep="\t")
    transcript_table.set_index("tx", inplace=True)

    return set(transcript_table[transcript_table[rna_seq_name] >= threshold].index)


def find_common_calls_main(data_path, experiment_name, transcript_bed_file, reference_fasta, orf_callers,
                           transcript_list_file, rna_seq_name, transcript_tpm_threshold):
    transcript_coordinates = orf_utils.read_transcript_bed_file(transcript_bed_file)
    ref_sequences = orf_utils.read_genome(reference_fasta)

    if "ribocode" in orf_callers:
        ribocode_result_quick_check = orf_utils.ribocode_result_processor(
            orf_result_file=f"{data_path}/ribocode_results/ribocode_results.txt",
            ).format_results(
                ref_sequences, transcript_coordinates)
    else:
        ribocode_result_quick_check = None

    if "ribotish" in orf_callers:
        ribotish_result_quick_check = orf_utils.ribotish_result_processor(
            orf_result_file=f"{data_path}/ribotish_results/pred.txt",
            ).format_results()
    else:
        ribotish_result_quick_check = None

    if "price" in orf_callers:
        price_result_quick_check = orf_utils.price_result_processor(
            orf_result_file=f"{data_path}/price_results/{experiment_name}.orfs.tsv",
            ).format_results(ref_sequences)
    else:
        price_result_quick_check = None

    if "orfquant" in orf_callers:
        orfquant_result_quick_check = orf_utils.orfquant_result_processor(
            orf_result_file=f"{data_path}/orfquant_results/{experiment_name}_final_ORFquant_results.csv",
            ).format_results(ref_sequences, transcript_coordinates)
    else:
        orfquant_result_quick_check = None

    all_orf_calls = {"price": price_result_quick_check,
                     "ribotish": ribotish_result_quick_check,
                     "ribocode": ribocode_result_quick_check,
                     "orfquant": orfquant_result_quick_check}
    overlap_csv_file = f"{data_path}/{experiment_name}_found_by_any_caller.csv"
    full_transcript_info_csv_file = f"{data_path}/{experiment_name}_found_by_any_caller_all_transcripts.csv"
    overlap_bed_file = f"{data_path}/{experiment_name}_found_by_all_caller.bed"
    transcript_list = read_transcript_list(transcript_list_file, rna_seq_name, transcript_tpm_threshold)
    merged_orf_calls, merged_orf_calls_full_list = find_overlaps(all_orf_calls, orf_callers, transcript_list)
    write_overlaps_into_file(merged_orf_calls, overlap_csv_file, orf_callers)
    write_overlaps_into_bed(merged_orf_calls, overlap_bed_file, orf_callers)
    write_overlaps_into_file(merged_orf_calls_full_list, full_transcript_info_csv_file, orf_callers)
