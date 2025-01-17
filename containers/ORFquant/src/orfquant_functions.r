## required libraries:
# conflicted, R.utils, tidyverse, multidplyr, ineq, IRanges, GenomicRanges, GenomicFeatures, ORFik
# install.packages(c("conflicted","R.utils","tidyverse","multidplyr","ineq","BiocManager"))
# BiocManager::install(c("IRanges", "GenomicRanges", "GenomicFeatures", "ORFik"))

##########################################################################################
# Setup
help_text <- "
# Canonical transcript annotation
-g[tf]        The path to a `gtf` file for canonical transcript definitions.  This or 
              `-in_rdata` must be provided.

-f[asta]      The path to the reference genome fasta file.

# Bam file data
-b[am]        A comma-separated list of bam files to process.
-n[ame]       Name of the transcript annotation GTF file.

# Output file
-o[out]       Output files prefix
-x[experiment_name]       Name of the current experiment.

-h[elp]       This help message
"

# quietly load pacman
suppressMessages(suppressWarnings(library(pacman, warn.conflicts=FALSE, quietly=TRUE, verbose=FALSE)))
## stop and quit
exit <- function() { invokeRestart("abort") }
full_stop <- function(..., status=0) {stop(paste0(...), call.=FALSE); exit()}
## turn off legends in ggplot
nolegend <- function() theme(legend.position="none")

## Get command line arguments
p_load(R.utils, update=FALSE, install=TRUE)
config <- R.utils::cmdArgs(unique=TRUE)
configured <- names(config)
if (length(configured) == 0 ||
    any(configured == "h") ||
    any(configured == "help")) {
    cat(help_text)
    exit()
}

### make short arguments long
names(config)[ names(config) == "g"] <- "gtf"
names(config)[ names(config) == "f"] <- "fasta"
names(config)[ names(config) == "b"] <- "bam"
names(config)[ names(config) == "n"] <- "name"
names(config)[ names(config) == "o"] <- "out"
names(config)[ names(config) == "x"] <- "experiment_name"
configured <- names(config)


suppressPackageStartupMessages(library("RiboseQC"))
suppressPackageStartupMessages(library("ORFquant"))

cat("Prepare annotation files...\n")
annotation_directory = paste(config[['out']], "orfquant_annot", sep="/")
prepare_annotation_files(
    annotation_directory = annotation_directory,
    twobit_file = config[['fasta']],
    gtf_file = config[['gtf']],
    annotation_name = config[['name']],
    export_bed_tables_TxDb = F,
    forge_BSgenome = T,
    create_TxDb = T)

cat("Run RiboseQC analysis...\n")
RiboseQC_analysis(
    annotation_file = paste(annotation_directory, paste(config[['name']], "gtf_Rannot", sep="."), sep="/"),
    bam_files = config[['bam']],
    report_file = paste(config[['out']], paste(config[['name']], "html", sep="."), sep="/"),
    write_tmp_files = F)


###install.packages(
###    paste(annotation_directory, "BSgenome.Homo.sapiens.veliadbv11", sep="/"),
###    repos=NULL, type="source")
###library("BSgenome.Homo.sapiens.veliadbv11")

cat("Run ORFquant...\n")
run_ORFquant(
    for_ORFquant_file = paste(config[['bam']], "for_SaTAnn", sep="_"),
    annotation_file = paste(annotation_directory, paste(config[['name']], "gtf_Rannot", sep="."), sep="/"),
    n_cores = 16,
    prefix = paste(config[['out']], config[['experiment_name']], sep="/"),
    canonical_start_only = F)

cat("Save results into CSV...\n")
orfquant_result_file <- paste(config[['out']], paste(config[['experiment_name']], "final_ORFquant_results", sep="_"), sep="/")
load(orfquant_result_file)
orfs <- ORFquant_results$ORFs_tx
columns_to_drop <- c(
    "unique_features_reads",
    "adj_unique_features_reads",
    "scaling_factors",
    "compatible_with",
    "compatible_with_longest",
    "NMD_candidate_compatible_txs",
    "Distance_to_lastExEx_compatible_txs",
    "longest_ORF")
mcols(orfs) <- mcols(orfs)[, !(names(mcols(orfs)) %in% columns_to_drop)]
write.csv(orfs, file=paste(orfquant_result_file, "csv", sep="."))
cat("Finished.\n")
