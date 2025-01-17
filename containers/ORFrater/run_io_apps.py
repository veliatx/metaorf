"""This program implements an application for the ORF-RATER I/O.
"""

import os
import random
import tempfile
import shutil
import math
import streamlit as st
import pandas as pd
from datetime import datetime
from io import StringIO
import numpy as np
import boto3
import json
from smart_open import open



box_content = ("\"input_dir\": \"/mount/efs/riboseq_callers/data/ORFrater/input_batch_tmp_%d_%d\",\n"
               "\"output_dir\": \"/mount/efs/riboseq_callers/data/ORFrater/output_batch_tmp_%d_%d\",\n"
               "\"annotation_dir\": \"/mount/efs/riboseq_callers/data/ORFrater/human_GRCh38_p14\",\n"
               "\"contaminant_genomes\": \"hg38_rRNA.fa,hg38_tRNA.fa\",\n"
               "\"reference_genomes\": \"GRCh38.p14.genome.fa\",\n"
               "\"genome_annotation_prefix\": \"veliadb_v1\",\n"
               "\"contaminant_genome_index\": \"star\",\n"
               "\"reference_genome_index\": \"star\",\n"
               "\"transcriptome_bed_file\": \"veliadb_v1.bed\",\n"
               "\"pseudogenes\": \"transcript_ids_on_pseudogenes.txt\",\n"
               "\"gene_names\": \"genename_mapping.txt\",\n"
               "\"multimap\": \"3\",\n")

def increase_experiment():
    st.session_state["n_experiments"] += 1


def decrease_experiment():
    st.session_state["n_experiments"] = max(0, st.session_state["n_experiments"]-1)


def increase_experiment_to_compare():
    st.session_state["n_experiments_to_compare"] += 1


def decrease_experiment_to_compare():
    st.session_state["n_experiments_to_compare"] = max(0, st.session_state["n_experiments_to_compare"]-1)


def click_submit():
    st.session_state["submitted"] = True
    

def check_result():
    st.session_state["check_result"] = True


def check_qc_result():
    st.session_state["check_qc_result"] = True


def check_comparison():
    st.session_state["check_comparison"] = True


def get_description(experiment):
    with open(f"s3://velia-piperuns-dev/{experiment}/{experiment}.json") as json_file:
        parameters = json.load(json_file)
    return parameters["note"]


def process_data_frame(df):
    df = df.dropna(axis=0, how="all", ignore_index=True)
    df = df.T
    header = df.iloc[0]
    df = df[1:]
    df.columns = header
    return df


def process_csv(uploaded_csv, state):
    """Reads parameters from a CSV file."""

    experiment_flags = []
    if uploaded_csv is None:
        return experiment_flags

    sample_df = pd.DataFrame()
    job_df = pd.DataFrame()
    reading_status = {"sample": False, "job": False}
    for index, row in pd.read_csv(uploaded_csv, header=None).iterrows():
        if np.all([str(val).strip() in ["", "nan"] for _, val in row.items()]):
            continue

        if str(row.iloc[0]).strip().startswith("["):
            for content_type in reading_status:
                reading_status[content_type] = False
            if str(row.iloc[0]).strip() == "[Samples]":
                reading_status["sample"] = True
            if str(row.iloc[0]).strip() == "[Riboseq_jobs]":
                reading_status["job"] = True
            continue

        if reading_status["sample"]:
            sample_df = pd.concat([sample_df, row], ignore_index = True, axis=1)
        if reading_status["job"]:
            job_df = pd.concat([job_df, row], ignore_index = True, axis=1)

    sample_df = process_data_frame(sample_df)
    st.subheader("Samples")
    st.write(sample_df)
    st.subheader("ORF calling jobs")
    job_df = process_data_frame(job_df)
    st.write(job_df)

    samples = {} 
    for index, row in sample_df.iterrows():
        sample_name = row["sample_id"]

        if sample_name not in samples:
            samples[sample_name] = []

        if row['containing_folder'].startswith("s3"):
            sample_path_list = []
            for fastq_file in row['R1_fastq_file'].split(","):
                sample_path_list.append(os.path.join(row['containing_folder'], fastq_file))
            sample_path = ",".join(sample_path_list)
        else:
            sample_path_list = []
            for fastq_file in row['R1_fastq_file'].split(","):
                sample_path_list.append(os.path.join(
                    "s3://velia-data-dev/VDC_003_ngs/primary/raw_data/",
                    row['containing_folder'],
                    "Analysis/1/Data/fastq",
                    fastq_file))
            sample_path = ",".join(sample_path_list)
        if row["organism"] == "human":
            flags = box_content
            if state["skip_orfcalling"]:
                flags += "\"skip_orfcalling\": \"%s\",\n" %state["skip_orfcalling"]
        else:
            flags = None
            st.write(f"The organism {row['organism']} is not supported at the moment.")

        if (row["umi"] is None or 
            (isinstance(row["umi"], float) and math.isnan(row["umi"])) or
            (isinstance(row["umi"], str) and not row["umi"].strip())):
            row["umi"] = ""
        else:
            row["umi"] = row["umi"].strip()
        samples[sample_name].append((sample_path, flags, row["adaptor_sequence"].strip(), row["umi"]))

    experiment_flags = []
    for index, row in job_df.iterrows():
        job_name = "VPR_orfcalling_%s" %(datetime.today().strftime("%Y%m%d%H%M%S"))
        if state["skip_orfcalling"]:
            job_name += "_mapping_only"
        sample_paths_combined = []
        tis_paths_combined = []
        flags = None
        chx_adapter_sequence, chx_umi = None, None
        tis_adapter_sequence, tis_umi = None, None
        note = ""

        if str(row["CHX"]).strip() not in ["", "nan", "NA"]:
            for sample_name in row["CHX"].split(","):
                sample_name = sample_name.strip()
                job_name += f"_{sample_name}"
                for (sample_path, sample_flags, adapter_sequence, umi) in samples[sample_name]:
                    sample_paths_combined.append(sample_path)
                    flags = sample_flags %(state["tmp_file_random_idx"], index, state["tmp_file_random_idx"], index)
                    chx_adapter_sequence = adapter_sequence
                    chx_umi = umi
            flags += (f"\"adapter_sequence\": \"{chx_adapter_sequence}\",\n"
                      f"\"umi\": \"{chx_umi}\",\n")

        if str(row["TIS"]).strip() not in ["", "nan", "NA"]:
            for sample_name in str(row["TIS"]).split(","):
                sample_name = sample_name.strip()
                job_name += f"_tis_{sample_name}"
                for (sample_path, sample_flags, adapter_sequence, umi) in samples[sample_name]:
                    tis_paths_combined.append(sample_path)
                    if flags is None:
                        flags = sample_flags %(state["tmp_file_random_idx"], index, state["tmp_file_random_idx"], index)
                    tis_adapter_sequence = adapter_sequence
                    tis_umi = umi
            flags += (f"\"adapter_sequence_tis\": \"{tis_adapter_sequence}\",\n"
                      f"\"umi_tis\": \"{tis_umi}\",\n")
        
        if str(row["Note"]).strip() not in ["", "nan", "NA"]:
            note = str(row["Note"]).strip()

        experiment_flags.append((job_name, note, ",".join(sample_paths_combined),
                                 ",".join(tis_paths_combined), flags))

    st.subheader("Jobs to submit")
    data = {"job_name":[], "note":[], "CHX_samples":[], "TIS_samples":[], "flags":[]}
    for (job_name, note, sample_paths, tis_paths, flags) in experiment_flags:
        data["job_name"].append(job_name)
        data["note"].append(note)
        data["CHX_samples"].append(sample_paths)
        data["TIS_samples"].append(tis_paths)
        data["flags"].append(flags)
    parameters_df = pd.DataFrame(data=data)
    st.write(parameters_df)
    st.download_button(
            label="Download job parameters as CSV",
            data=parameters_df.to_csv().encode('utf-8'),
            file_name=f"parameters_{state['tmp_file_random_idx']}.csv",
            mime="text/csv")
    return experiment_flags


def add_text_areas(state):
    """Add text areas and collect user-input data."""

    experiment_flags = []
    for count in range(state['n_experiments']):
        name = st.text_input(
                label=f"Experiment {count} name",
                value="BL11_PBMCL1_2_S7_R1_001")
        note = st.text_input(label=f"Experiment {count} note", value="")
        sample = st.text_input(
                label=f"Experiment {count} sample path",
                value="s3://velia-data-dev/VDC_003_ngs/primary/raw_data/230222_VH01198_15_AAAY5CMHV/Analysis/1/Data/fastq/BL11_PBMCL1_2_S7_R1_001.fastq.gz")
        adapter_sequence = st.text_input(label=f"Experiment {count} Adapter sequence - CHX", value="auto")
        umi = st.text_input(label=f"Experiment {count} UMI - CHX", value="")
        tis = st.text_input(
                label=f"Experiment {count} TIS path (Leave this empty if TIS reads are not available)", value="")
        adapter_sequence_tis = st.text_input(label=f"Experiment {count} Adapter sequence - TIS", value="auto")
        umi_tis = st.text_input(label=f"Experiment {count} UMI - TIS", value="")
        flags = st.text_area(
                label=f"Experiment {count} additional parameters",
                value=box_content %(state["tmp_file_random_idx"], count, state["tmp_file_random_idx"], count),
                height=300)

        if sample:
            flags += (f"\"adapter_sequence\": \"{adapter_sequence}\",\n"
                      f"\"umi\": \"{umi}\",\n")
        if tis:
            flags += (f"\"adapter_sequence_tis\": \"{adapter_sequence_tis}\",\n"
                      f"\"umi_tis\": \"{umi_tis}\",\n")
        if state["skip_orfcalling"]:
            flags += "\"skip_orfcalling\": \"%s\",\n" %state["skip_orfcalling"]

        experiment_flags.append((name, note, sample, tis, flags))
        st.divider()
    return experiment_flags


def add_experiments(state):
    """Add text areas and collect user-input data."""

    experiment_list = []
    for count in range(state['n_experiments_to_compare']):
        option_name = st.selectbox(
            f'Experiment {count} name',
            [""]+st.session_state["job_names"],
            disabled=state["check_comparison"])
        if option_name:
            experiment_list.append(option_name)
    return experiment_list


def write_into_json(temp_dir, name, note, sample, tis, flags):
    """Writes parameters into json."""

    with open(os.path.join(temp_dir, name+".json"), "w") as ofile:
        ofile.write("{\n\t")
        ofile.write(f"\"experiment_name\": \"{name}\",\n\t")
        ofile.write(f"\"note\": \"{note}\",\n\t")
        ofile.write(f"\"sample_paths_s3\": \"{sample}\",\n\t")
        if tis.strip():
            ofile.write(f"\"tis_paths\": \"{tis}\",\n\t")
        flags = flags.replace(",\n", ",\n\t")
        if flags.endswith(",\n\t"):
            flags = flags[:-3] + "\n"
        ofile.write(flags)
        ofile.write("}\n")


def get_job_names_on_s3():
    """Retrieves top-level folder names."""

    folder_names = []

    client = boto3.client('s3')
    paginator = client.get_paginator('list_objects')
    result = paginator.paginate(Bucket='velia-piperuns-dev', Delimiter='/')
    for prefix in result.search('CommonPrefixes'):
        folder_name = prefix.get('Prefix')
        if folder_name.endswith("/"):
            folder_name = folder_name[:-1]

        if folder_name.startswith("VPR_orfcalling"):
            folder_names.append(folder_name)
    return folder_names


def check_tis_exists(state):
    s3 = boto3.resource('s3')
    obj = s3.Object('velia-piperuns-dev', f'{state["qc_job_name_to_check"]}/output/dup_stats_tis.log')
    try:
        obj.load()
        return True
    except:
        return False
    

def read_features(experiment_list, feature_name):
    feature_list = []
    for experiment in experiment_list:
        with open(f"s3://velia-piperuns-dev/{experiment}/output/aligned/{experiment}_qc.tsv") as tsv_file:
            tsv = pd.read_table(tsv_file)
            feature_list.append(tsv[feature_name].values)
    return feature_list


def get_log_data(experiment_list, log_name):
    features = {}
    for experiment in experiment_list:
        for log_line in open(f"s3://velia-piperuns-dev/{experiment}/output/{log_name}"):
            key, value = [dt.strip() for dt in log_line.strip().split("|")]
            if key not in features:
                features[key] = []
            value = float(value) if "." in value else int(value)
            features[key].append([value])
    return features


if __name__ == "__main__":
    # states
    if "n_experiments" not in st.session_state:
        st.session_state["n_experiments"] = 0
    if "tmp_file_random_idx" not in st.session_state:
        st.session_state["tmp_file_random_idx"] = random.randint(0, 1e10)
    if "submitted" not in st.session_state:
        st.session_state["submitted"] = False
    if "skip_orfcalling" not in st.session_state:
        st.session_state["skip_orfcalling"] = False
    if "job_uploaded" not in st.session_state:
        st.session_state["job_uploaded"] = False
    if "check_result" not in st.session_state:
        st.session_state["check_result"] = False
    if "job_name_to_check" not in st.session_state:
        st.session_state["job_name_to_check"] = None
    if "check_qc_result" not in st.session_state:
        st.session_state["check_qc_result"] = False
    if "qc_job_name_to_check" not in st.session_state:
        st.session_state["qc_job_name_to_check"] = None
    if "job_names" not in st.session_state:
        st.session_state["job_names"] = get_job_names_on_s3()
    if "n_experiments_to_compare" not in st.session_state:
        st.session_state['n_experiments_to_compare'] = 0
    if "check_comparison" not in st.session_state:
        st.session_state['check_comparison'] = False

    tab1, tab2, tab3, tab4 = st.tabs(["Submit", "Result", "QC", "Comparison"])
    with tab1:
        # input widgets
        st.title(":dna: Submit ORF-calling jobs")
        st.divider()

        st.session_state["skip_orfcalling"] = st.toggle('Skip orf-calling')
        st.header('Load jobs from a sample sheet')
        uploaded_csv = st.file_uploader("Choose a sample sheet")
        experiment_flags_uploaded = process_csv(uploaded_csv, st.session_state)
        st.divider()

        st.header('Manually add jobs')
        experiment_flags_manual = add_text_areas(st.session_state)
        st.button("Add a new job", on_click=increase_experiment)
        st.button("Remove the last job", on_click=decrease_experiment)
        st.divider()

        st.header('Submit jobs to AWS :sunglasses:')

        if not st.session_state["submitted"]:
            st.button(':red[**Submit**]', on_click=click_submit)
        elif not st.session_state["job_uploaded"]:
            temp_dir = "parameter_%d" %st.session_state["tmp_file_random_idx"]
            os.mkdir(temp_dir)
            if experiment_flags_manual:
                for i in range(st.session_state['n_experiments']):
                    name, note, sample, tis, flags = experiment_flags_manual[i]
                    write_into_json(temp_dir, name, note, sample, tis, flags)
                st.write(f'Submitted! Refresh this page to submit new jobs.')
            elif experiment_flags_uploaded:
                for i in range(len(experiment_flags_uploaded)):
                    name, note, sample, tis, flags = experiment_flags_uploaded[i]
                    write_into_json(temp_dir, name, note, sample, tis, flags)
                st.write(f'Submitted! Refresh this page to submit new jobs.')
            else:
                st.write("No job is specified. Refresh this page to try again.")

            cmd = ("python launch_batch_job_auto.py "
                   "--experiment_parameter_folder %s " %temp_dir)
            os.system(cmd)
            shutil.rmtree(temp_dir)
            st.session_state["job_uploaded"] = True
        else:
            st.write(f'Submitted! Refresh this page to submit new jobs.')

    with tab2:
        st.title(":mag: Check ORF-calling results")
        st.divider()

        # s3 = boto3.resource('s3')
        if not st.session_state["check_result"]:
            st.session_state["job_name_to_check"] = st.text_input("Type a job name to check", "")
            st.write("OR")
            option = st.selectbox('Select a job name to check', [""]+st.session_state["job_names"])
            if option:
                st.session_state["job_name_to_check"] = option
            st.button('Check results', on_click=check_result)
        else:
            try:
                st.write(f'**Experiment name:** {st.session_state["job_name_to_check"]}')
                st.write(f'**Experiment description:** {get_description(st.session_state["job_name_to_check"])}')
                st.divider()
                st.header('Download ORF-RATER results')

                # ORFRATER h5
                with open(f"s3://velia-piperuns-dev/{st.session_state['job_name_to_check']}"
                          "/output/orfrater_results/orfratings.h5", "rb") as h5_file:
                    st.download_button(
                            label="Download ORF-RATER result (h5 format)",
                            data=h5_file,
                            file_name="orfratings.h5",
                            mime='application/x-hdf5')
                    st.write(f"File location: s3://velia-piperuns-dev/{st.session_state['job_name_to_check']}"
                             f"/output/orfrater_results/orfratings.h5")
                    st.divider()

                # ORFRATER bed file
                with open(f"s3://velia-piperuns-dev/{st.session_state['job_name_to_check']}"
                          "/output/orfrater_results/orfratings.bed") as bed_file:
                    st.download_button(
                            label="Download ORF-RATER result (bed format)",
                            data=bed_file,
                            file_name="orfratings.bed")
                    st.write(f"File location: s3://velia-piperuns-dev/{st.session_state['job_name_to_check']}"
                             f"/output/orfrater_results/orfratings.bed")
                    st.divider()

                # table of quantified translation values
                with open(f"s3://velia-piperuns-dev/{st.session_state['job_name_to_check']}"
                          "/output/orfrater_results/quant.h5", "rb") as h5_file:
                    st.download_button(
                            label="Download quantified translation values (h5 format)",
                            data=h5_file,
                            file_name="quant.h5",
                            mime='application/x-hdf5')
                    st.write(f"File location: s3://velia-piperuns-dev/{st.session_state['job_name_to_check']}"
                             f"/output/orfrater_results/quant.h5")
            except Exception as error:
                print(error)
                st.write(f"Some (or all) of the ORF calling results for {st.session_state['job_name_to_check']} "
                         f"cannot be found. Please refresh this page to submit a new query.")
    with tab3:
        st.title(":wrench: Check quality control results")
        st.divider()

        s3 = boto3.resource('s3')
        if not st.session_state["check_qc_result"]:
            st.session_state["qc_job_name_to_check"] = st.text_input("Type a job name to check its QC results", "")
            st.write("OR")
            option_qc = st.selectbox('Select a job name to check its QC results', [""]+st.session_state["job_names"])
            if option_qc:
                st.session_state["qc_job_name_to_check"] = option_qc
            st.button('Check QC results', on_click=check_qc_result)
        else:
            try:
                st.write(f'**Experiment name:** {st.session_state["qc_job_name_to_check"]}')
                st.write(f'**Experiment description:** {get_description(st.session_state["qc_job_name_to_check"])}')
                st.divider()
                st.header("Read mapping")
                st.subheader("Number of reads at different stages")
                log_list = []
                for log_line in open(f"s3://velia-piperuns-dev/{st.session_state['qc_job_name_to_check']}"
                                     "/output/pipeline_stats.log"):
                    log_list.append([dt.strip() for dt in log_line.strip().split("|")])
                st.write(pd.DataFrame(np.array(log_list)))
                st.divider()

                st.header("CHX samples")
                st.subheader("Duplication")
                log_list = []
                for log_line in open(f"s3://velia-piperuns-dev/{st.session_state['qc_job_name_to_check']}"
                                     "/output/dup_stats.log"):
                    log_list.append([dt.strip() for dt in log_line.strip().split("|")])
                st.write(pd.DataFrame(np.array(log_list)))

                st.subheader("Footprints, coverage, and periodicity")
                with open(f"s3://velia-piperuns-dev/{st.session_state['qc_job_name_to_check']}"
                          f"/output/aligned/{st.session_state['qc_job_name_to_check']}_qc.pdf", "rb") as pdf_file:
                    st.download_button(
                            label="Download QC figures",
                            data=pdf_file,
                            mime='application/octet-stream',
                            file_name=f"{st.session_state['qc_job_name_to_check']}_qc.pdf")
                    
                with open(f"s3://velia-piperuns-dev/{st.session_state['qc_job_name_to_check']}"
                          f"/output/aligned/{st.session_state['qc_job_name_to_check']}_qc.tsv") as tsv_file:
                    st.download_button(
                            label="Download QC table",
                            data=tsv_file,
                            file_name=f"{st.session_state['qc_job_name_to_check']}_qc.tsv")
                st.divider()
                
                if check_tis_exists(st.session_state):
                    st.header("TIS samples")
                    st.subheader("Duplication")
                    log_list = []
                    for log_line in open(f"s3://velia-piperuns-dev/{st.session_state['qc_job_name_to_check']}"
                                        "/output/dup_stats_tis.log"):
                        log_list.append([dt.strip() for dt in log_line.strip().split("|")])
                    st.write(pd.DataFrame(np.array(log_list)))

                    st.subheader("Footprints, coverage, and periodicity")
                    with open(f"s3://velia-piperuns-dev/{st.session_state['qc_job_name_to_check']}"
                            f"/output/aligned_tis/{st.session_state['qc_job_name_to_check']}_tis_qc.pdf", "rb") as pdf_file:
                        st.download_button(
                                label="Download QC figures",
                                data=pdf_file,
                                mime='application/octet-stream',
                                file_name=f"{st.session_state['qc_job_name_to_check']}_tis_qc.pdf")
                        
                    with open(f"s3://velia-piperuns-dev/{st.session_state['qc_job_name_to_check']}"
                            f"/output/aligned_tis/{st.session_state['qc_job_name_to_check']}_tis_qc.tsv") as tsv_file:
                        st.download_button(
                                label="Download QC table",
                                data=tsv_file,
                                file_name=f"{st.session_state['qc_job_name_to_check']}_tis_qc.tsv")
            except Exception as error:
                print(error)
                st.write(f"Some (or all) of the QC results for {st.session_state['qc_job_name_to_check']} "
                         f"cannot be found. Please refresh this page to submit a new query.")
    with tab4:
        st.title(":scales: Compare results across experiments")
        st.divider()

        st.header('Add experiments to compare')
        experiment_list = add_experiments(st.session_state)
        st.button("Add a new experiment", on_click=increase_experiment_to_compare)
        st.button("Remove the last experiment", on_click=decrease_experiment_to_compare)
        st.divider()

        if not st.session_state["check_comparison"]:
            st.button('Compare', on_click=check_comparison)
        else:
            st.header("Experiment index")
            st.write(pd.DataFrame([[experiment] for idx, experiment in enumerate(experiment_list)], columns=["Experiment name"]))
            x_arr = [idx for idx, _ in enumerate(experiment_list)]
            st.divider()

            st.header("Riboseq data quality")
            for feature_name in ["periodicity_score", "size_peak", "size_gini",
                                "cds_v_5utr", "cds_v_3utr", "cds_v_utr"]:
                st.subheader(feature_name)
                data = pd.DataFrame(read_features(experiment_list, feature_name), index=x_arr, columns=[feature_name])
                st.bar_chart(data)
                st.divider()

            st.header("Read mapping")
            for key, value in get_log_data(experiment_list, "pipeline_stats.log").items():
                st.subheader(key)
                if len(value) < len(x_arr):
                    continue
                data = pd.DataFrame(value, index=x_arr, columns=[key])
                st.bar_chart(data)
                st.divider()

            st.header("UMI duplicate")
            for key, value in get_log_data(experiment_list, "dup_stats.log").items():
                st.subheader(key)
                if len(value) < len(x_arr):
                    st.write("Not Available.")
                    continue
                data = pd.DataFrame(value, index=x_arr, columns=[key])
                st.bar_chart(data)
                st.divider()
            #print(read_features(experiment_list))
            #data = pd.DataFrame([10, 20], index=x_arr)
            #st.bar_chart(data)


##############################
