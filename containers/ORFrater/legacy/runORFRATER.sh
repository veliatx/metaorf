#!/bin/bash

echo_prefix="" # set to echo for troubleshooting

# conda
conda_env=orfrater
path_to_orfrater="conda run -n $conda_env python /home/ec2-user/ORF-RATER"
conda_source=/home/ec2-user/anaconda3/etc/profile.d/conda.sh
source $conda_source
conda run -n $conda_env python --version

exclude_genes=(AC006548.28 ADGRV1 AHNAK AHNAK2 AK4 ATP8B1 CADM2 CHS.10076 CHS.23448 DST EPPK1 FAT FGF14 HELLPAR HMCN1 HMCN2 HYDIN KCNQ1OT1 KIAA1109 KMT2D LINC-PINT LOC101927919 LOC105377700 LOC440434 LRIF1 MACF1 MDN1 MUC16 MUC4 MYRIP NEB NM_001164462.1 OBSCN PCLO PLSCR1 PPM1A PPP2R2A PRH1-TAS2R14 RNF213 RNF219 RP3-394A18.1 RRN3P3 RUNX1 RYR1 RYR2 SCN7A SMURF1 SYNE1 SYNE2 TMEM101 TRPS1 TTN UBR4 ZNF253)

if [ $# -eq 0 ]; then
        echo "ORF-RATER pipeline. If there are multiple initiation site (TIS) bam files, should use samtools to combine them first."
        echo "Usage: run_orfrater.sh [options] Riboseq.bam [Riboseq2.bam]" 
        echo ""
        echo -e "\tOptions:"
        echo -e "\t-----------------"
        echo -e "\t-b)\tPath to BED file. Be sure it corresponds to the GTF file used for alignment."
        echo -e "\t-s)\tpseudogenes list. Text file list of transcript IDs annotated as pseudogenes, one per line."
        echo -e "\t-f)\tPath to genome fasta file. Be sure it corresponds to the GTF file used for alignment."
        echo -e "\t-g)\tPath to file with gene name lookup table. Formatted as a two-column tab-delimited file: the first column with the transcript ID and second column with the gene names."
        echo -e "\t-m)\tMax length of footprint to consider"
        echo -e "\t-n)\tMin length of footprint to consider"
        echo -e "\t-t)\tWhether to include TIS (e.g. Harringtonine or LTM) and the path to the transcriptome-aligned TIS BAM file."
        echo -e "\t-d)\tWorking directory. Uses directory of fastq file if not suppled."
        echo -e "\t-a)\tAdapter sequence (original:AGATCGGAAGAGCACACGTCTGAAC), illumina adapter: AGATCGGAAGAGCACACGTCT (default)"
        echo -e "\t-x)\tAdapter clipper program ([fastx]/atria/skewer)"
        echo -e "\t-G)\tGenome to use (one of [hg19], hg38, mm10, mm9, rn5, rn6 or full path to STAR index directory"
        echo -e "\t-C)\tContaminant genome to counterscreen against, determined from `Genome` short code if not provided"
        exit 0
else
        bedfile="/genomes/hg38/hg38_annotation/hg38.ncbiRefSeq_no-altChr.bed"
        pseudogenes=""
        fasta="/genomes/hg38/hg38_annotation/hg38_no-altChr.fa"
        genenames=""
        tis=0
        maxLength=34
        minLength=24
        bamfile_tis=""
        # output="orfrater"
        PARAMS=""
        clipSeq="AGATCGGAAGAGCACACGTCT"
        clipper="FASTX"
        genome="hg19"
        contaminant_path=""
        declare -a MAYBE_FASTQS

        while (( "$#" )); do
                case "$1" in
                -b)
                        bedfile=$2
                        shift 2
                ;;
                -s)
                        pseudogenes=$2
                        shift 2
                ;;
                -f)
                        fasta=$2
                        shift 2
                ;;
                -g)
                        genenames=$2
                        shift 2
                ;;
                -m)
                maxLength=$2
                shift 2
                ;;
                -n)
                minLength=$2
                shift 2
                ;;
                -t)
                        tis=1
                        echo "Running with initiation site enriched dataset (i.e. Harr or LTM)"
                        bamfile_tis=$2
                        shift 2
                ;;
                ###########################
                # Gene Cutler customization 2021.10.11
                -d)
                        custom_DIR=$2
                        shift 2
                        ;;
                -x)
                        if [[ $2 == "Atria" || $2 == "atria" || $2 == "a" || $2 == "A" ]]; then
                                clipper="ATRIA"
                        elif [[ $2 == "Skewer" || $2 == "skewer" || $2 == "S" || $2 == "s" ]]; then
                                clipper="SKEWER"
                        else
                                clipper="FASTX"
                        fi
                        shift 2
                        ;;
                -a)
                        clipSeq=$2
                        shift 2
                        ;;
                -G)
                        genome=$2
                        shift 2
                        ;;
                -C)
                        contaminant_path=$2
                        shift 2
                        ;;
                --) # end argument parsing
                        shift
                        break
                ;;
                -*|--*=) # unsupported flags
                        echo "Error: Unsupported flag $1" >&2
                        exit 1
                ;;
                *) # preserve positional arguments
                        MAYBE_FASTQS+=($1)
                        shift
                ;;
                esac
        done
fi



die_if_empty() {
        if [ ! -e $1 ]; then
                echo "Error: File $1 doesn't exist!  Exiting."
                exit 1
        elif [ ! -s $1 ]; then
                echo "Error: File $1 is empty!  Exiting." 
                exit 1
        else
                return 1
        fi
}

# Function to trim, filter, and align fastq to bam
# input = full filename of fastq file
# output = final bam file saved to $BAM_FROM_FASTQ
make_bam_for_fastq() {
        local FULL_FASTQ=$1
        echo "============  Generating aligned bam for $FULL_FASTQ"
        local FASTQ_FILENAME=${FULL_FASTQ##/*/}
        local GZIPPED=0

        if [[ $FULL_FASTQ =~ .fastq.gz$ ]]; then 
                GZIPPED=1; 
                FASTQ_BASENAME=${FASTQ_FILENAME%%.fastq.gz}
        else
                FASTQ_BASENAME=${FASTQ_FILENAME%%.fastq}
        fi

        local FINAL_BAM="${OUTDIR}/${FASTQ_BASENAME}.filtered.aligned.unique.bam"
        if [ -s $FINAL_BAM ]; then
                echo "Fastq supplied, but BAM already exists.  Using $FINAL_BAM"
                BAM_FROM_FASTQ=$FINAL_BAM
                return 0
        fi


        # set the path to the STAR genome indices
        if [[ $genome == "hg19" || $genome == "hg38" || $genome == "mm10" || $genome == "mm9" || $genome == "rn5" ||  $genome == "rn6" ]]; then
                genome_path=/genomes/${genome}/${genome}.star
                contaminant_path=/genomes/${genome}/${genome}cont.star
        elif [ -d "$genome" ] && [ -e "$genome/SAindex" ]; then
                genome_path=$genome
        else
                echo "Error: No valid STAR aligner genome reference found"
                exit 1
        fi

        if [ ! -d $contaminant_path ] || [ ! -e $contaminant_path/SAindex ]; then
                echo "Error: No valid STAR aligner contaminant reference found"
                exit 1
        fi



        if [ -e "${OUTDIR}/${FASTQ_BASENAME}.trimmed.fastq.gz" ]; then
                echo "Trimmed fastq already exists... skipping"
        else
                echo "Trimming $FULL_FASTQ"
                start_timestamp=$(date +%s)
                if [ $clipper == "ATRIA" ]; then
                        echo "Trimming with Atria"
                        WORK_DIR=`mktemp -d -p "$OUTDIR"`
                        # Atria produces a lot of messages, so dump those
                        /usr/local/bin/atria -q 33 --polyG -a $clipSeq --length-range 20:500 -n 10 -t 12 -r $FULL_FASTQ -o $WORK_DIR 2> /dev/null 
                        /usr/local/bin/atria -q 20 -a $clipSeq --length-range 20:500 -n 15 -t 12 -r $FULL_FASTQ -o "$WORK_DIR" 2> /dev/null 
                        if [ -e "$WORK_DIR" ] && [ -e "$WORK_DIR"/*.atria.log ]; then
                                mv "$WORK_DIR/$FASTQ_BASENAME.atria.fastq.gz" "${OUTDIR}/${FASTQ_BASENAME}.trimmed.fastq.gz"
                                /bin/rm -rf "$WORK_DIR"
                        else
                                echo "Error: No output directory for Atria found"
                                exit 1
                        fi
                elif [ $clipper == "SKEWER" ]; then
                        echo "Trimming with Skewer"
                        /usr/local/bin/skewer -x $clipSeq --mode any -r 0.05 --min 20 -n --quiet --stdout --threads 8 $FULL_FASTQ | pigz -c --fast > "${OUTDIR}/${FASTQ_BASENAME}.trimmed.fastq.gz"
                else 
                        echo "Trimming with Fastx"
                        if [ $GZIPPED -eq 0 ]; then
                                cat $FULL_FASTQ | fastx_clipper -Q33 -l 20 -n -v -c -a $clipSeq | fastx_trimmer -Q33 -f 1 2>${OUTDIR}/trim.log > "${OUTDIR}/${FASTQ_BASENAME}.trimmed.fastq"
                        else
                                pigz -dc $FULL_FASTQ | fastx_clipper -Q33 -l 20 -n -v -c -a $clipSeq | fastx_trimmer -Q33 -f 1 2>${OUTDIR}/trim.log > "${OUTDIR}/${FASTQ_BASENAME}.trimmed.fastq"
                        fi
                        pigz "${OUTDIR}/${FASTQ_BASENAME}.trimmed.fastq"
                fi
                ###########################
                curr_timestamp=$(date +%s)
                elapsed_time=$(expr $curr_timestamp - $start_timestamp)
                echo "Done trimming... trimming took $(expr $elapsed_time / 60) minutes"
        fi

# fastq alignment
        die_if_empty "${OUTDIR}/${FASTQ_BASENAME}.trimmed.fastq.gz"
        if [ -s "${OUTDIR}/${FASTQ_BASENAME}.filtered.aligned.unique.bam" ]; then
                echo "Looks like you already have alignments for $FULL_FASTQ. Skipping..."
        else
                echo "Aligning to contaminant database (tRNA/rRNA)"
                ###########################
                # Gene Cutler customization 2021.10.11
# align to contaminant database
                STAR --readFilesCommand zcat --outSAMstrandField intronMotif --outReadsUnmapped Fastx \
                                --genomeDir $contaminant_path --runThreadN 18 \
                                --readFilesIn "${OUTDIR}/${FASTQ_BASENAME}.trimmed.fastq.gz" --outFileNamePrefix "${OUTDIR}/${FASTQ_BASENAME}.contaminant."
                echo "Aligning unaligned to $genome with refGene annotations"
# align to full genome
                STAR --outSAMstrandField intronMotif --genomeDir $genome_path --runThreadN 18 \
                                --readFilesIn "${OUTDIR}/${FASTQ_BASENAME}.contaminant.Unmapped.out.mate1" \
                                --outFileNamePrefix "${OUTDIR}/${FASTQ_BASENAME}.filtered." \
                                --outFilterMismatchNmax 2 --outFilterMultimapNmax 5 --outSAMprimaryFlag AllBestScore \
                                --chimScoreSeparation 10 --chimScoreMin 20 --chimSegmentMin 15 --outSAMattributes All
                echo "Done aligning ${OUTDIR}/${FASTQ_BASENAME}.contaminant.Unmapped.out.mate1"
                pigz --force "${OUTDIR}/${FASTQ_BASENAME}.contaminant.Unmapped.out.mate1"
                /bin/rm -f "${OUTDIR}/${FASTQ_BASENAME}."*.Log.*
        fi

# filter and sort alignment results into bam file
        if [ -s "${OUTDIR}/${FASTQ_BASENAME}.filtered.aligned.unique.bam" ]; then
                echo "Seems you already have the filtered sorted bams. Skipping... "
        else
                samtools view -@ 6 -bS "${OUTDIR}/${FASTQ_BASENAME}.filtered.Aligned.out.sam" | samtools sort -@ 6 -l 0 -o "${OUTDIR}/${FASTQ_BASENAME}.sorted.bam"
                samtools view -@ 6 -bS -F 0X100 "${OUTDIR}/${FASTQ_BASENAME}.sorted.bam" > "${OUTDIR}/${FASTQ_BASENAME}.noSecond.bam"
                samtools view -@ 6 -q 10 -b "${OUTDIR}/${FASTQ_BASENAME}.noSecond.bam" > $FINAL_BAM
                /bin/rm "${OUTDIR}/${FASTQ_BASENAME}.filtered.Aligned.out.sam" \
                        "${OUTDIR}/${FASTQ_BASENAME}.sorted.bam" \
                        "${OUTDIR}/${FASTQ_BASENAME}.contaminant.Aligned.out.sam" \
                        "${OUTDIR}/${FASTQ_BASENAME}.contaminant.SJ.out.tab" \
                        "${OUTDIR}/${FASTQ_BASENAME}.filtered.Chimeric.out.sam" \
                        "${OUTDIR}/${FASTQ_BASENAME}.filtered.Chimeric.out.junction" \
                        "${OUTDIR}/${FASTQ_BASENAME}.filtered.SJ.out.tab"
        fi

        echo "============  Completed generating aligned bam for $FULL_FASTQ: $FINAL_BAM"
        BAM_FROM_FASTQ=$FINAL_BAM
        return 0
}



echo "Running ORF-RATER"
echo "Using BED file: $bedfile"
echo "Using pseudogene file: $pseudogenes"
echo "Using fasta file: $fasta"
echo "Using gene names file: $genenames"
echo "Min footprint length: $minLength"
echo "Max footprint length: $maxLength"

echo "Input files:"
for BAM in ${MAYBE_FASTQS[@]}; do echo ">>>> $BAM"; done
echo "-----------------"
echo 

if [ ! -z $custom_DIR ]; then 
        OUTDIR=${custom_DIR%%/}
        mkdir -p $OUTDIR
else 
        OUTDIR="."
fi


#       if input files are fastqs, they will be aligned here
#       resulting filenames of new or input bams saved to $BAMS array
declare -a BAMS
for MAYBE_FASTQ in ${MAYBE_FASTQS[@]}; do
#               if "BAM" files are actually fastqs, run the alignment pipeline from RunRiboSeq
        if [[ $MAYBE_FASTQ =~ .fastq(.gz)?$ ]]; then
                BAM_FROM_FASTQ=""
                make_bam_for_fastq $MAYBE_FASTQ
                die_if_empty $BAM_FROM_FASTQ
                BAMS+=($BAM_FROM_FASTQ)
        elif [[ $MAYBE_FASTQ =~ .bam$ ]]; then
                # the file was already a bam file and not a fastq file
                BAMS+=($MAYBE_FASTQ)
        else
                echo "Format of $MAYBE_FASTQ is not recognized.  Must be fastq, fastq.gz, or bam.  Exiting!"
                exit 1
        fi
done


# index bams
for BAM in ${BAMS[@]}; do
        BAI="$BAM.bai"
        if [ ! -e $BAI ]; then
                echo "Indexing $BAM"
                samtools index -b -@ 6 $BAM $BAI
        fi
done


# Find most common P-site offset. Should run for each dataset separately. 
# 1. accumulate arrays of input and output files
declare -a BAM_ARRAY
declare -a OFFSETFILES
declare -a FILTERED_OFFSETFILES
declare -a TALLYFILES
ALL_BAMS=${BAMS[@]}
if [ -n $bamfile_tis ]; then ALL_BAMS+=($bamfile_tis); fi
for BAM in ${ALL_BAMS[@]}; do
        BASENAME=${BAM##*/}
        BASENAME=${BASENAME%%.*}
        BAM_ARRAY+=($BAM)
        TALLYFILES+=("$BASENAME.tallies.txt")
        OFFSETFILES+=("$BASENAME.offsets_auto.txt")
done

# 2. run in parallel on these inputs/outputs
echo "Running psite_trimmed.py on ${BAM_ARRAY[*]}"
parallel -P 20% "source $conda_source  ~/anaconda3/etc/profile.d/conda.sh && \
                $echo_prefix $path_to_orfrater/psite_trimmed.py --force \
                --subdir $OUTDIR --offsetfile '{1}' \
                        --tallyfile '{2}' -p 3 --cdsbed $bedfile --minrdlen 18 --maxrdlen 40 '{3}'" \
                        :::  ${OFFSETFILES[@]} :::+ ${TALLYFILES[@]} :::+ ${BAM_ARRAY[@]}
# 3. filter resulting offset files
for OFFSET in ${OFFSETFILES[@]}; do
        OFFSETFILE_FILTERED=${OFFSET%%.offsets_auto.txt}.offsets.txt
        FILTERED_OFFSETFILES+=($OFFSETFILE_FILTERED)
        awk '$1 >= a && $1 <= b { print }' a="$minLength" b="$maxLength" < $OUTDIR/$OFFSET > $OUTDIR/$OFFSETFILE_FILTERED
done

# Prune transcripts, remove pseudogenes if file provided
echo "Running prune_transcripts.py"
if [ -n "$pseudogenes" ]; then PSEUDOGENES_ARG=("--pseudogenes" "$pseudogenes"); fi
$echo_prefix $path_to_orfrater/prune_transcripts.py $fasta $BAMS --force --peakfrac 1 --minreads 20 \
    --minlen 28 --maxlen 31 --inbed $bedfile --outbed "$OUTDIR/transcripts.bed" ${PSEUDOGENES_ARG[@]} \
    --summarytable  "$OUTDIR/tidsummary.txt" -vvp 8
if [ ! -e "$OUTDIR/transcripts.bed"  ]; then echo "prune_transcripts.py failed to make a bed file.  Exiting!"; exit 1; fi
# Create transcript families
echo "Running make_tfams.py"
if [ -n "$genenames" ]; then GENENAMES_ARG=("-g" "$genenames"); fi
$echo_prefix $path_to_orfrater/make_tfams.py -v ${GENENAMES_ARG[@]} --force --tfamstem "$OUTDIR/tfams" --inbed "$OUTDIR/transcripts.bed" 

echo "Running find_orfs_and_types.py"
# lots of python errors to STDERR
$echo_prefix $path_to_orfrater/find_orfs_and_types.py $fasta --tfamstem "$OUTDIR/tfams" --inbed "$OUTDIR/transcripts.bed" --orfstore "$OUTDIR/orf.h5" --codons NTG -v -p 24 2> $OUTDIR/find_orfs_and_types.errors.txt
die_if_empty "$OUTDIR/tfams.txt"

if [ $tis -eq 1 ]; then
        mkdir $OUTDIR/TIS
        BASENAME=${bamfile_tis##*/}
        BASENAME=${bamfile_tis%%.*}
        OFFSETFILE=$BASENAME.offsets.txt
        cp $OUTDIR/$OFFSETFILE $OUTDIR/TIS/offsets.txt
        echo "Running regress_orfs.py with initiation site enriched dataset"
        $echo_prefix $path_to_orfrater/regress_orfs.py $bamfile_tis --startonly --startcount 1 -vvp 16 \
                --orfstore              $OUTDIR/orf.h5 \
                --inbed                 $OUTDIR/transcripts.bed \
                --offsetfile    offsets.txt \
                --exclude $exclude_genes \
                 2>/dev/null

        echo "Running regress_orfs.py on regular dataset"
        mkdir $OUTDIR/ND
        cp $OUTDIR/${FILTERED_OFFSETFILES[0]} $OUTDIR/ND/offsets.txt
        $echo_prefix $path_to_orfrater/regress_orfs.py $BAMS --subdir $OUTDIR/ND --restrictbystarts $OUTDIR/TIS -vvp 16 \
                --orfstore              $OUTDIR/orf.h5 \
                --inbed                 $OUTDIR/transcripts.bed \
                --offsetfile    offsets.txt \
                --exclude ${exclude_genes[@]} \
                 2>/dev/null
        die_if_empty "$OUTDIR/ND/regression.h5"

        echo "Running rate_regression_output.py"
        $echo_prefix $path_to_orfrater/rate_regression_output.py $OUTDIR/ND $OUTDIR/TIS --minperleaf 1 --minforestscore 0.5 -vp 16 \
                --ratingsfile $OUTDIR/orfratings.h5 \
                2> $OUTDIR/rate_regression_output.errors.txt

else
        echo "No TIS provided. Running regress_orfs.py on regular dataset"
        mkdir $OUTDIR/ND
        cp $OUTDIR/${FILTERED_OFFSETFILES[0]} $OUTDIR/ND/offsets.txt
          # use fewer cores at this step
        $echo_prefix $path_to_orfrater/regress_orfs.py $BAMS --subdir $OUTDIR/ND -vvp 8 \
                --orfstore              $OUTDIR/orf.h5 \
                --inbed                 $OUTDIR/transcripts.bed \
                --offsetfile    offsets.txt \
                --exclude ${exclude_genes[@]} \
                 2> $OUTDIR/regress_orfs.errors.txt
        die_if_empty "$OUTDIR/ND/regression.h5"

        echo "Running rate_regression_output.py"
        $echo_prefix $path_to_orfrater/rate_regression_output.py $OUTDIR/ND --minperleaf 1 --minforestscore 0.5 -vp 16 \
                --ratingsfile $OUTDIR/orfratings.h5 \
                 2> $OUTDIR/rate_regression_output.errors.txt
        die_if_empty "$OUTDIR/orfratings.h5"
fi

echo "Running make_orf_bed.py"
$echo_prefix $path_to_orfrater/make_orf_bed.py --ratingsfile $OUTDIR/orfratings.h5 --minrating 0.1 -c --inbed $OUTDIR/transcripts.bed --outbed $OUTDIR/orfratings.bed

echo "Running quantify_orfs.py"
# full BAM names will generate errors, convert to basenames (no path, no extension) and supply with the --names flag
BASENAMES=""
for BAM in ${BAMS[@]}; do BASENAME=${BAM##*/}; BASENAME=${BASENAME%%.*}; BASENAMES="$BASENAMES $BASENAME"; done
BASENAMES=${BASENAMES##*( )}
$echo_prefix $path_to_orfrater/quantify_orfs.py ${BAMS[@]} --names $BASENAMES \
        --force \
        -vp 16 \
        --subdir $OUTDIR \
        --offsetfile ND/offsets.txt \
        --metagenefile ND/metagene.txt \
        --quantfile quant.h5 \
        --ratingsfile $OUTDIR/orfratings.h5 \
        --inbed  $OUTDIR/transcripts.bed \
        --minrating 0.4