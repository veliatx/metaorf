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
-m[ap]        A tab-separated text file containing two named columns, `name` and `file`, 
              that contain the sample names and file paths to the `bam` files to be 
              processed.  This or `-bam` must be provided.
-b[am]        A comma-separated list of bam files to process.  This or `-map` must be 
              provided.
-n[ame]       An optional comma-separated list of sample names in the same order and 
              number as `-bam` files.  Ignored if `-map` is provided.

# Output file
-o[out]       Output files prefix

# Parallel procesing
-p[rocesses]  The number of parallel processes. Default is total cores - 1.

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
### do we want help?
if (length(configured) == 0 ||
    any(configured == "h") ||
    any(configured == "help")) {
    cat(help_text)
    exit()
}
### make short arguments long
names(config)[ names(config) == "g"] <- "gtf"
names(config)[ names(config) == "f"] <- "fasta"
names(config)[ names(config) == "m"] <- "map"
names(config)[ names(config) == "b"] <- "bam"
names(config)[ names(config) == "n"] <- "name"
names(config)[ names(config) == "o"] <- "out"
names(config)[ names(config) == "p"] <- "n_cpu"

configured <- names(config)

cat("Loading R libraries...\n")
## make sure prefered versions of functions are used
p_load(conflicted, update=FALSE, install=TRUE)

## data table libraries
p_load(tidyverse, update=FALSE, install=TRUE)
p_load(data.table, update=FALSE, install=TRUE)
p_load(gridExtra, update=FALSE, install=TRUE)

## multi processing
p_load(parallel, update=FALSE, install=TRUE)


## genome-related libraries
p_load(GenomicRanges, update=FALSE, install=TRUE)
p_load(GenomicFeatures, update=FALSE, install=TRUE)
p_load(ORFik, update=FALSE, install=TRUE)
#p_load(BSgenome.Hsapiens.UCSC.hg38, update=FALSE, install=TRUE)

## statistics
p_load(ineq, update=FALSE, install=TRUE) # Gini


conflict_prefer("collapse", "IRanges", quiet=TRUE)
conflict_prefer("select", "dplyr", quiet=TRUE)
conflict_prefer("slice", "dplyr", quiet=TRUE)
conflict_prefer("filter", "dplyr", quiet=TRUE)
conflict_prefer("rename", "dplyr", quiet=TRUE)
conflict_prefer("last", "dplyr", quiet=TRUE)


##########################################################################################
# Check file names and arguments
cat("Checking script arguments...\n")
for (x in c("gtf","map","in_rdata")) {
    if (any(configured == x) && ! file.exists(config[[x]]))
        full_stop("The ",x," file '",config[[x]],"' could not be read")
}
if (!any(configured == "gtf") && !any(configured == "in_rdata")) full_stop("One of 'gtf' or 'in_rdata' must be specified")
if (any(configured == "gtf") && any(configured == "in_rdata")) warn("Both 'gtf' and 'in_rdata' are specified -- 'gtf' will be ignored")
if (!any(configured == "bam") && !any(configured == "map")) full_stop("One of 'map' or 'bam' must be specified")
if (any(configured == "bam") && any(configured == "map")) warn("Both 'bam' and 'map' are specified -- 'bam' will be ignored")
if (!any(configured == "out")) full_stop("Output file prefix must be specified")

# check output file locations
if (file.access(base::dirname(config[["out"]]), mode=2) == -1) full_stop("Directory for output is not writeable")
if (file.exists(paste0(config[["out"]],".tsv"))) full_stop("TSV output file already exists -- will not overwrite")
if (file.exists(paste0(config[["out"]],".pdf"))) full_stop("PDF output file already exists -- will not overwrite")
#pdf(file=paste0(config[["out"]],".pdf"), onefile=TRUE, paper="letter")

# set number of processors for multiprocessing
max_proc <- parallel::detectCores()
requested_proc <- as.numeric(config[["n_cpu"]])
if (!is.na(requested_proc) && length(requested_proc)>0 && requested_proc > 0 ) {
    n_proc <- min(max_proc, requested_proc)
} else {
    n_proc <- max(1, max_proc - 1)  
}
cat("Running on", n_proc,"processors\n")
options(mc.cores=n_proc)



if (any(configured == "map")) {
    bam_map <- read_tsv(config[["map"]], show_col_types = FALSE)
    if (!nrow(bam_map)>0 || !any(colnames(bam_map)=="name") || !any(colnames(bam_map)=="file")) {
        full_stop("The bam mapping file is misconfigured.  See `help`.")
    }
} else if (any(configured == "bam")) {
    bam_map <- tibble(file=unlist(stringr::str_split(config[['bam']], ",")))
    if (any(configured == "name")) {
        # use the provided sample names
        names <- unlist(stringr::str_split(config[['name']], ","))
        if (length(names) != nrow(bam_map))
            full_stop("The number of sample names does not match the number of bam files")
        bam_map <- bam_map %>% mutate(name=names)
    } else {
        # or generate sample names from bam file names
        bam_map <- bam_map %>% mutate(name=str_match(file,"(\\w[^/.]+)\\.bam")[,2])
        if (any(is.na(bam_map$name)))
            full_stop("Sample names could not be parsed from the bam file names")
    }
}
if (length(unique(bam_map$name)) < length(bam_map$name)) full_stop("The sample names must be unique")
# check that all the bam files are good
bam_map <- bam_map %>% mutate(file_ok = map_lgl(file, file.exists))
if (any(bam_map$file_ok == FALSE)) full_stop("One or more bam files listed in the bam mapping file can not be read")


##########################################################################################
# Read inputs
## Read canonical transcript information
cat("Reading canonical transcript information...\n")
if (!str_detect(config[['gtf']], stringr::regex("\\.gtf$", ignore_case=TRUE))) {
    full_stop("The specified gtf file does not end in '.gtf'")
}
cat(" → Parsing gtf file\n")
txdb <- suppressWarnings(loadTxdb(paste(config[['gtf']], ".db", sep="")))
# txdb <- suppressWarnings(loadTxdb(config[['gtf']]))
# txdb <- makeTxdbFromGenome(config[['gtf']], config[['fasta']], organism = "Homo sapiens", optimize = TRUE, return = TRUE)
tx_Names <- filterTranscripts(txdb, minFiveUTR=100, minCDS=100, minThreeUTR=100, longestPerGene=TRUE) # valid transcripts
utr_5 <- loadRegion(txdb, part="leaders", names.keep=tx_Names)
utr_3 <- loadRegion(txdb, part="trailers", names.keep=tx_Names)
cds   <- loadRegion(txdb, part="cds", names.keep=tx_Names)
tx  <- loadRegion(txdb, part="tx", names.keep=tx_Names)
windows_Start <- startRegion(cds, tx, TRUE, upstream = 50, 49)


##########################################################################################
# Process data
## Read bam files
cat("Loading bam files...\n")
data <- bam_map %>% mutate(
    bam = mclapply(file, readBam, mc.cores=n_proc), 
    widths = map(bam, width),
    good_reads = map_dbl(widths, ~sum(.x >= 26 & .x <= 34))
)
bad_bams <- data %>% filter(good_reads < 1000)
if (nrow(bad_bams) > 0 ) {
  cat(nrow(bad_bams),"bam files dropped - not RiboSeq or too few reads\n")
  data <- data %>% filter(good_reads >= 1000) %>% select(-widths)
}

cat("Calculating offsets...\n")
data$shifts <- lapply(data$bam, function(.x) detectRibosomeShifts(.x, txdb, firstN=60)) # does not work parallelized

cat("Measuring ribosome footprints...\n")
data <- data %>% 
  mutate(
    shifted = map2(bam, shifts, shiftFootprints),
    footprints = map(bam, ~{.x %>% readWidths() %>% table() %>% as.data.frame() %>% mutate(Size=as.numeric(as.character(`.`)))}),
    size_peak = map_dbl(footprints, ~{.x %>% slice_max(Freq, n=1, with_ties=FALSE) %>% pull(Size)}),
    # fraction reads +/- 1 from peak
    size_peak_frac=map2_dbl(footprints, size_peak, ~{sum(.x$Freq[ .x$Size >= .y-1 & .x$Size <= .y+1  ]) / sum(.x$Freq)}),
    # gini value of size distribution
    size_gini=map_dbl(footprints, ~Gini(.x$Freq[.x$Freq>10]))
    )

p1 <- data %>% 
  select(name, footprints) %>% 
  unnest(footprints) %>% 
  # collapse values with footprint > 40
  mutate(Freq=if_else(Size==40, sum(Freq[Size>40]),Freq)) %>% 
  filter(Size <= 40) %>% 
  rename(Footprint=Size, Frequency=Freq, Sample=name) %>% 
  ggplot(aes(x=Footprint, y=Frequency, fill=Sample)) + geom_col(alpha=0.75) + 
  facet_grid(Sample~., scales="free_y") +
  ggtitle("Footprint Size Distribution") +
  theme_bw() + nolegend()
#print(p1)

## transcript coverage
cat("Calculating transcript location coverage...\n")
data <- data %>% mutate(coverage_table = map(shifted, function(.x) {
    coverage_5utr <- metaWindow(.x, utr_5, scoring = NULL, feature = "utr_5", scaleTo=100, withFrames=FALSE)
    coverage_cds <-  metaWindow(.x, cds,   scoring = NULL, feature = "cds",   scaleTo=100, withFrames=FALSE)
    coverage_3utr <- metaWindow(.x, utr_3, scoring = NULL, feature = "utr_3", scaleTo=100, withFrames=FALSE)
    coverage <- rbindlist(list(coverage_5utr, coverage_cds, coverage_3utr))
    coverage_by_gene <- coverage[, `:=`(total = sum(score), positions=sum(score>0)), by=genes]
    good_genes <- coverage_by_gene$positions >= 5 & coverage_by_gene$total >= 5 & coverage_by_gene$total < 1000
    if (sum(good_genes) < 500) warning("Too few genes for accurate coverage analysis!", call.=FALSE)
    coverage <- coverage %>% as_tibble() %>% filter( genes %in% coverage_by_gene$genes[good_genes])
    # drop last nucleotide of utr_5
    coverage <- coverage %>% filter(!(feature=="utr_5" & position == 100))
    summarized_coverage <- coverage %>% 
        group_by(genes) %>% mutate(total=sum(score)) %>% ungroup() %>%
        mutate(norm_score=score/total*100, feature=factor(feature, levels=c("utr_5","cds","utr_3"))) %>% 
        group_by(feature, position) %>%
        summarize(summed_score=sum(score), summed_norm=sum(norm_score), .groups="drop")
    summarized_coverage <- summarized_coverage %>% arrange(feature,position)
    summarized_coverage
}))
rm(utr_5, cds, utr_3)

bad <- data %>% mutate(coverage_rows=map_dbl(coverage_table, nrow)) %>% filter(coverage_rows==0)
if (nrow(bad) > 0) {
  warning(paste0("!!!!! ",nrow(bad)," samples dropped from analysis"))
  data <- data %>% mutate(coverage_rows=map_dbl(coverage_table, nrow)) %>% filter(coverage_rows>0)
}

data <- data %>% mutate(
    summed_coverage = map(coverage_table, ~{.x %>% filter(position >= 3 & position <= 98) %>% group_by(feature) %>% summarize(sum=sum(summed_score), .groups="drop")}),
    cds_v_5utr = map_dbl(summed_coverage, ~{.x$sum[.x$feature=="cds"]/(.x$sum[.x$feature=="cds"]+.x$sum[.x$feature=="utr_5"])}),
    cds_v_3utr = map_dbl(summed_coverage, ~{.x$sum[.x$feature=="cds"]/(.x$sum[.x$feature=="cds"]+.x$sum[.x$feature=="utr_3"])}),
    cds_v_utr =  map_dbl(summed_coverage, ~{.x$sum[.x$feature=="cds"]/(.x$sum[.x$feature=="cds"]+.x$sum[.x$feature=="utr_3"]+.x$sum[.x$feature=="utr_5"])})
) %>% 
    select(-summed_coverage)

p2 <- data %>% 
    collect() %>%
    select(Sample=name, coverage_table) %>% 
    unnest(coverage_table) %>% 
    rename(Coverage=summed_norm) %>% 
    ggplot(aes(x=position, fill=feature)) + 
    geom_area(aes(y=Coverage)) + 
    geom_line(aes(y=Coverage), color="black", size=0.5) + 
    facet_grid(Sample~feature, scales="free_y") + 
    theme_bw() + nolegend() +
    ggtitle("Transcript Coverage")
#print(p2)


## periodicity
cat("Calculating periodicity...\n")
data <- data %>% mutate(
    periods = map(shifted, ~metaWindow(.x, windows_Start, withFrames=TRUE)),
    periodicity_score = map_dbl(periods, ~{
        .x %>% 
            filter(position >= 0) %>% 
            mutate(
                codon_0 = if_else(frame == 0, 1, 0),
                codon = cumsum(codon_0)
            ) %>% 
            group_by(codon) %>% 
            summarize( fraction_0 = score[frame==0]/sum(score), .groups="drop") %>% 
            summarize( fraction_0 = mean(fraction_0) ) %>% 
            pull(fraction_0)
    })
)

p3 <- data %>% 
    select(Sample=name, periods) %>% 
    unnest(periods) %>% 
    rename(Coverage=score, Frame=frame, Position=position) %>% 
    filter(Position %inrange% c(-6,29)) %>% 
    mutate(Frame=factor(Frame)) %>% 
    ggplot(aes(x=Position, fill=Frame)) + 
    geom_col(aes(y=Coverage)) + 
    facet_grid(Sample~., scale="free_y") + 
    ggtitle("Codon Periodicity") +
    theme_bw() +
    theme(legend.position="bottom")
#print(p3)

##########################################################################################
# Save results
cat("Saving results...\n")

# extract offset data into wide format
wide_shifts <- data %>% 
  select(name, shifts) %>% 
  unnest_wider(shifts) %>% 
  unnest(c(fraction, offsets_start)) %>% 
  filter(fraction >= 26 & fraction <= 33) %>% 
  pivot_wider(names_from=fraction, values_from=offsets_start, names_prefix="Offset.")

## save tsv
data %>% left_join(wide_shifts, by="name") %>% 
    mutate(reads=map_dbl(bam, length)) %>% 
    select(file, name, reads, size_peak, size_peak_frac, size_gini, cds_v_5utr, cds_v_3utr, cds_v_utr, periodicity_score, starts_with("Offset")) %>%
    readr::write_tsv(file=paste0(config[["out"]],".tsv"))

## finish pdf
combined_plots <- do.call(marrangeGrob, args = list(grobs = list(p1, p2, p3), ncol=1, nrow=1))
ggsave(paste0(config[["out"]],".pdf"), combined_plots, width=8.5, height=11)


cat("Finished.\n")
