"""Re-formats TIS transformer results into the BED format.
"""



if __name__ == "__main__":
    tis_result_file = ""
    samples = []
    for group_idx in range(1, 7):
        load_tis_data(tis_result_file %group_idx, samples)

    output_bed_file = ""
    write_into_bed_file(samples, output_bed_file)
