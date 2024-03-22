"""This script defines functions that read and format orf-calling results
from various orf-callers, including:
* ribocode
* ribotish
* price
* orfquant
"""

import math
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd 


def read_genome(ref_path_list):
    ref_dict = {}
    for ref_path in ref_path_list:
        for seq in SeqIO.parse(ref_path, "fasta"):
            chromID = seq.id
            chromSeq = (str(seq.seq)).upper()
            ref_dict[chromID] = chromSeq
    return ref_dict


def read_transcript_bed_file(transcript_bed_file):
    """Returns a dictionary that maps transcript_ids to their genomic locations."""

    transcripts = {}
    for line in open(transcript_bed_file):
        elements = line.strip().split()
        # print(elements)
        transcripts[elements[3]] = {
            "chrom_id": elements[0],
            "start": int(elements[1]),
            "end": int(elements[2]),
            "strand": elements[5],
            "block_count": int(elements[9]),
            "block_sizes": [int(size) for size in elements[10].split(",") if size],
            "block_starts": [int(start_pos) for start_pos in elements[11].split(",") if start_pos],
        }
    return transcripts


class orf_result_processor:
    def __init__(self, orf_result_file):
        self.orf_result_file = orf_result_file


    def format_results(self, *args, **kwargs):
        """Read ORF-calling results and format the results into a dictionary object.

        The dictionary object is constructed as:
            key -> (chrom_id, orf_start, orf_end, strand, nucleotide_seq, transcript_id) 
            content -> (gene_id, ORF_type, orf_score, exon_blocks)
        """
        pass
    

    def _get_block_info(self, transcript, tstart, tstop, ref_sequences):
        """Gets the block coordinates and sequences of the given ORF.

        Given the information of which transcript the ORF of interest maps to
        and the start/stop on the transcript, this function can compute the 
        corresponding genomic coordinates and sequences. Blocks will be presented
        in the format as start_1-end_1(exclusive, 0-based)|start_2-end_2|.....
        """

        if transcript["strand"] == "-":
            transcript_length = sum(transcript["block_sizes"])
            tstart, tstop = (transcript_length-tstop), (transcript_length-tstart)

        blocks = []
        current_transcript_length = 0
        orf_sequence = ""
        orf_start, orf_end = None, None
        for block_index in range(transcript["block_count"]):
            if tstart - current_transcript_length >= transcript["block_sizes"][block_index]:
                current_transcript_length += transcript["block_sizes"][block_index]
                continue
            current_start = max(0, tstart - current_transcript_length)
            current_end = min(tstop - current_transcript_length, transcript["block_sizes"][block_index])
            current_abs_start = current_start + transcript["block_starts"][block_index] + transcript["start"]
            current_abs_end = current_end + transcript["block_starts"][block_index] + transcript["start"]
            if orf_start is None:
                orf_start = current_abs_start
            orf_end = current_abs_end
            blocks.append(f"{current_abs_start}-{current_abs_end}")
            orf_sequence += ref_sequences[transcript["chrom_id"]][current_abs_start:current_abs_end]
            if (tstop - current_transcript_length) <= transcript["block_sizes"][block_index]:
                break
            else:
                current_transcript_length += transcript["block_sizes"][block_index]

        if transcript["strand"] == "-":
            orf_sequence = str((Seq(orf_sequence)).reverse_complement())
        return "|".join(blocks), orf_sequence, orf_start, orf_end


class ribocode_result_processor(orf_result_processor):
    def __init__(self, orf_result_file):
        super().__init__(orf_result_file)


    def format_results(self, ref_sequences, transcript_coordinates):
        """
        format of the returns
            map:
                key: chrom_id orf_start orf_end strand nucleotide_seq transcript_id
                content: gene_id ORF_type pval_combined exon_blocks
        """
        
        ribocode_result_quick_check = {}
        column_names = []
        for line in open(self.orf_result_file):
            elements = line.strip().split("\t")
            
            if not column_names:
                column_names = [elem for elem in elements]
            else:    
                orf = {}
                for idx, elem in enumerate(elements):
                    orf[column_names[idx]] = elem
                    
                if orf["strand"] == "+":
                    orf_start = int(orf["ORF_gstart"]) - 1
                    orf_end = int(orf["ORF_gstop"])
                else:
                    orf_start = int(orf["ORF_gstop"]) - 1
                    orf_end = int(orf["ORF_gstart"])
                
                orf_tstart = int(orf["ORF_tstart"]) - 1
                orf_tstop = int(orf["ORF_tstop"])
                orf_blocks, orf_sequence, _, _ = self._get_block_info(
                    transcript_coordinates[orf["transcript_id"]], orf_tstart, orf_tstop, ref_sequences)
                
                log_score = -math.log10(float(orf["pval_combined"])) if float(orf["pval_combined"]) != 0 else 999
                ribocode_result_quick_check[
                    (orf["chrom"], orf_start, orf_end, orf["strand"], orf_sequence, orf_blocks)] = [
                        orf["gene_id"],
                        orf["transcript_id"],
                        orf["ORF_type"],
                        log_score]
        return ribocode_result_quick_check


class ribotish_result_processor(orf_result_processor):
    def __init__(self, orf_result_file):
        super().__init__(orf_result_file)


    def format_results(self):
        """
        format of the returns
            map:
                key: chrom_id orf_start orf_end strand nucleotide_seq transcript_id
                content: gene_id ORF_type pval_combined exon_blocks
        """
    
        ribotish_result_quick_check = {}
        column_names = []
        for line in open(self.orf_result_file):
            elements = line.strip().split("\t")
            
            if not column_names:
                column_names = [elem for elem in elements]
            else:
                orf = {}
                for idx, elem in enumerate(elements):
                    orf[column_names[idx]] = elem
                    
                chrom_id, start_end, strand = orf["GenomePos"].split(":")
                orf_start, orf_end = start_end.split("-")
                log_score = -math.log10(float(orf["RiboPvalue"])) if float(orf["RiboPvalue"]) != 0 else 999
                orf_blocks = "|".join(orf["Blocks"].split(","))
                ribotish_result_quick_check[
                    (chrom_id, int(orf_start), int(orf_end), strand, orf["Seq"], orf_blocks)] = [
                    orf["Gid"],
                    orf["Tid"],
                    orf["TisType"],
                    log_score]
        return ribotish_result_quick_check
    

class price_result_processor(orf_result_processor):
    def __init__(self, orf_result_file):
        super().__init__(orf_result_file)


    def format_results(self, ref_sequences):
        """
        format of the returns
            map:
                key: chrom_id orf_start orf_end strand nucleotide_seq transcript_id
                content: gene_id ORF_type pval_combined exon_blocks
        """

        price_result_quick_check = {}
        column_names = []
        for line in open(self.orf_result_file):
            elements = line.strip().split("\t")
            
            if not column_names:
                column_names = [elem for elem in elements]
            else:    
                orf = {}
                for idx, elem in enumerate(elements):
                    orf[column_names[idx]] = elem

                if float(orf["p value"]) > 0.1:
                    continue
                    
                id_elements = orf["Id"].split("_")
                tid = "_".join(id_elements[:-2])
                orf_type = id_elements[-2]
                    
                chrom_id_strand, start_end_list = orf["Location"].split(":")
                chrom_id = chrom_id_strand[:-1]
                if chrom_id in [str(chr_num) for chr_num in range(1, 23)]+["X", "Y"]:
                    chrom_id = "chr"+chrom_id
                
                if chrom_id not in ref_sequences:
                    print(chrom_id)
                    continue
                
                strand = chrom_id_strand[-1]
                orf_start, orf_end = None, None
                orf_sequence = ""
                for start_end in start_end_list.split("|"):
                    start, end = start_end.split("-")
                    if orf_start is None:
                        orf_start = start
                    orf_end = end
                    orf_sequence += ref_sequences[chrom_id][int(start):int(end)]
                
                if strand == "-":
                    orf_sequence = str((Seq(orf_sequence)).reverse_complement())
                
                log_score = -math.log10(float(orf["p value"])) if float(orf["p value"]) != 0 else 999
                price_result_quick_check[
                    (chrom_id, int(orf_start), int(orf_end), strand, orf_sequence, start_end_list)] = [
                        orf["Gene"],
                        tid,
                        orf_type,
                        log_score]
        return price_result_quick_check


class orfquant_result_processor(orf_result_processor):
    def __init__(self, orf_result_file):
        super().__init__(orf_result_file)


    def format_results(self, ref_sequences, transcript_coordinates):
        """
        format of the returns
            map:
                key: chrom_id orf_start orf_end strand nucleotide_seq transcript_id
                content: gene_id ORF_type pval_combined exon_blocks
        """
        
        orfquant_result_quick_check = {}        
        for _, orf in pd.read_csv(self.orf_result_file).iterrows():            
            orf_tstart = orf["start"] - 1
            orf_tstop = orf["end"]
            orf_blocks, orf_sequence, orf_start, orf_end = self._get_block_info(
                transcript_coordinates[orf["transcript_id"]], orf_tstart, orf_tstop, ref_sequences)
            
            log_score = -math.log10(orf["pval"]) if orf["pval"] != 0 else 999
            orfquant_result_quick_check[
                (orf["region.seqnames"], orf_start, orf_end, orf["region.strand"], orf_sequence, orf_blocks)] = [
                    orf["gene_id"],
                    orf["transcript_id"],
                    orf["ORF_category_Tx"],
                    log_score]
        return orfquant_result_quick_check

