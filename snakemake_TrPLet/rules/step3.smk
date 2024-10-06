rule batchCorrect:
    input:
        "results/step2/7samples_TCGA11373_merged_counts.tsv.gz",
    output:
        "results/step3/batch_corrected_Wang_TCGA_counts_with_lineageCOV.tsv.gz",
    log:
        "logs/batchCorrect.log"
    params:
        batchCorrect=batchCorrect_script,
        threads=10,
    threads: 10
    resources:
        mem_mb=256000
    shell:
        """
        mkdir -p results/step3

        export PATH=$HOME/miniforge3/envs/TrPLet/bin:$PATH

        Rscript {params.batchCorrect} --merged_count {input} --output {output} 2> {log}

        """
