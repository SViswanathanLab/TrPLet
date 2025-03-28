#!/usr/bin/env bash

# activate environment
activate="${snakemake_params[conda_env]}"
source $activate TrPLet

exec 2> "${snakemake_log[0]}" 

reads=(${snakemake_input[reads]})  
r1="${reads[0]}"
r2="${reads[1]}"

GENOMEDIR="${snakemake_input[index]}"
GTFFILE="${snakemake_input[gtf]}"

OUTDIR=$(dirname "${snakemake_output[0]}")
sample_name=$(basename "$OUTDIR")
THREADS=${snakemake_params[threads]}

# Create output directory if it does not exist
mkdir -p $OUTDIR

if [[ "$r1" == *".gz"* ]]; then
  STAR --runThreadN $THREADS \
  --genomeDir $GENOMEDIR \
  --readFilesIn "${r1}" "${r2}" \
  --outFileNamePrefix "${OUTDIR}/${sample_name}_" \
  --readFilesCommand zcat \
  --outSAMattributes All \
  --outSAMtype BAM Unsorted \
  --quantMode TranscriptomeSAM \
  --sjdbGTFfile $GTFFILE \
  --outReadsUnmapped Fastx \
  --outMultimapperOrder Random
else
  STAR --runThreadN $THREADS \
  --genomeDir $GENOMEDIR \
  --readFilesIn "${r1}" "${r2}" \
  --outFileNamePrefix "${OUTDIR}/${sample_name}_" \
  --outSAMattributes All \
  --outSAMtype BAM Unsorted \
  --quantMode TranscriptomeSAM \
  --sjdbGTFfile $GTFFILE \
  --outReadsUnmapped Fastx \
  --outMultimapperOrder Random
fi
