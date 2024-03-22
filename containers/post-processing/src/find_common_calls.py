"""This program identifies the overlaps among the following orf-callers:
    RiboCode
    RiboTish
    PRICE
"""

import os
import copy

import orf_utils


def find_overlaps(all_orf_calls, orf_callers):
    """Identifies the overlaps among calls from different ORF callers."""

    merged_orf_calls = {}
    for orf_caller in orf_callers:
        orf_calls = all_orf_calls[orf_caller]
        for key in orf_calls:
            if key not in merged_orf_calls:
                merged_orf_calls[key] = {}

            # construct orf info into a dictionary
            merged_orf_calls[key][orf_caller] = {
                    "gene_id": orf_calls[key][0],
                    "transcript_id": orf_calls[key][1],
                    "orf_type": orf_calls[key][2],
                    "score": orf_calls[key][3]}
    return merged_orf_calls


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


def write_overlaps_into_file(merged_orf_calls, overlap_csv_file, overlap_bed_file, orf_callers):
    """For each cell_type/replicate/experiment, writes the merged ORF calls into file."""

    with open(overlap_csv_file, "w") as overlap_output, open(overlap_bed_file, "w") as overlap_bed_output:
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
            caller_count = 0
            line_to_write = (f"{chrom_id}\t{orf_start}\t{orf_end}\t{strand}\t{orf_sequence}\t{exon_blocks}")
            for orf_caller in orf_callers:
                if orf_caller in merged_orf_calls[key]:
                    orf_call_info = merged_orf_calls[key][orf_caller]
                    line_to_write += (f"\t{orf_call_info['score']}\t{orf_call_info['orf_type']}"
                                      f"\t{orf_call_info['gene_id']}\t{orf_call_info['transcript_id']}")
                    caller_count += 1
                else:
                    line_to_write += (f"\t \t \t \t ")
            overlap_output.write(line_to_write+"\n")
            
            if caller_count == len(orf_callers):
                bed_file_line = create_bed_line(
                    chrom_id, orf_start, orf_end, strand, exon_blocks,
                    merged_orf_calls[key][orf_callers[0]]['gene_id'], merged_orf_calls[key][orf_callers[0]]['transcript_id'])
                overlap_bed_output.write(bed_file_line)
    return


def find_common_calls_main(data_path, experiment_name, transcript_bed_file, reference_fasta, orf_callers):
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
    overlap_bed_file = f"{data_path}/{experiment_name}_found_by_all_caller.bed"
    merged_orf_calls = find_overlaps(all_orf_calls, orf_callers)
    write_overlaps_into_file(merged_orf_calls, overlap_csv_file, overlap_bed_file, orf_callers)
