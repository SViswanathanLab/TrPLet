#!/usr/bin/env bash

# activate environment
activate="${snakemake_params[conda_env]}"
source $activate TrPLet

exec 2> "${snakemake_log[0]}"  # send all stderr from this script to the log file

count_path="${snakemake_input[counts]}"
ref_path="rsem_ref/gencode_v38"
count_file=$(basename "$count_path")
sample_name="${count_file%_Aligned.toTranscriptome.out.bam}"
THREADS=${snakemake_params[threads]}

OUTDIR=$(dirname "${snakemake_output[0]}")

mkdir -p $OUTDIR # create output directory if it does not exist


rsem-calculate-expression -p $THREADS \
--bam \
--paired-end \
--strandedness none \
--append-names \
"${count_path}" "${ref_path}" "${OUTDIR}/${sample_name}"
