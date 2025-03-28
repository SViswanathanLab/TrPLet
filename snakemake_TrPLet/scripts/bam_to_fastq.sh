#!/usr/bin/env bash

# activate environment
activate="${snakemake_params[conda_env]}"
source $activate TrPLet

exec 2> "${snakemake_log[0]}"  # send all stderr from this script to the log file

read="${snakemake_input[bam_files]}"  # don't double-quote this - we want word splitting

THREADS=${snakemake_params[threads]}
Sorted_bam="${snakemake_output[0]}"
r1="${snakemake_output[1]}"
r2="${snakemake_output[2]}"

OUTDIR=$(dirname "${snakemake_output[0]}")
mkdir -p $OUTDIR # Create directory if it doesn't exist

samtools sort -@ $THREADS -n $read -o $Sorted_bam
samtools fastq -@ $THREADS $Sorted_bam -1 $r1 -2 $r2 -0 /dev/null -s /dev/null -n
