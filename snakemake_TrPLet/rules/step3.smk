rule batchCorrect:
    input:
        mergedCount="results/step2/merged_with_TCGA11373_counts.tsv.gz",
        lineageInfo="data/sample_lineageInfo.txt"
    output:
        "results/step3/batch_corrected_TCGA_counts_with_lineageCOV.tsv.gz",
    log:
        "logs/batchCorrect.log"
    params:
        batchCorrect=batchCorrect_script,
        threads=10,
        exportPath=exportPath,
    threads: 10
    resources:
        mem_mb=256000
    shell:
        """
        mkdir -p results/step3

        # export PATH=$HOME/miniforge3/envs/TrPLet/bin:$PATH
        export PATH={params.exportPath}:$PATH

        Rscript {params.batchCorrect} --merged_count {input.mergedCount} --output {output} --lineages {input.lineageInfo} 2> {log}

        """
