# MetaORF

MetaORF is an end-to-end pipeline designed to identify small open reading frames (sORFs) using ribo-seq data.

## Installation

After downloading this repository to your local machine, you can install MetaORF using PIP

`pip install .`

Next, you need to build the Docker images listed in the `contains` folder. These images set up the ribo-seq callers and the data processing pipeline.

## Run
### Configure data paths and experiment parameters
MetaORF can simultaneously process multiple ribo-seq FASTQ files from different cell types or tissues. Each cell type or tissue is treated as an individual experiment. To set up these experiments, create a sample sheet and, for each experiment, specify:
* Name of the experiment
* Data path to the Ribo-seq FASTQ file(s)
* Adaptor sequence
* Organism (current, only human and mouse are supported)
* UMI (leaving it empty if UMI is not used)
Each row of the sample sheet specifies one individual experiment.

### Prepare data
Use the following command to prepare the data:

`metaorf sample_sheets/<your sample sheet>.csv`

MetaORF uses AWS Batch to ensure computational parallelism for efficiency. It processes the FASTQ files, extracts features, and stores them in formatted output files. Definitions of the extracted features can be found in `metaorf_modeling/feature_definition_0416.txt`.

### Run model
Once the datasets are ready, use the scripts in the `metaorf_modeling` folder to identify and prioritize ORFs for each experiment:
* 00_merge_orfs.ipynb: Merges ORF candidates from all experiments into a single table in CSV format.
* 01_truthset_eval.ipynb: Loads the truth set and trains a gradient boosting model for ORF prioritization.
* 02_predict_ensemble.ipynb: Makes inferences and evaluates the confidence of ORF calls.

## License
MetaORF is licensed under the MIT License. See the LICENSE file for more details.
