rule rsem_ref:
    input:
        star_path=STAR,
        rsem_path=rsem,
    output:
        "data/gencode.v38.annotation.gtf",
        "data/GRCh38.primary_assembly.genome.fa",
        directory("rsem_ref"),
        multiext("rsem_ref/gencode_v38", ".chrlist", ".n2g.idx.fa", ".transcripts.fa", ".grp", ".seq", ".idx.fa", ".ti", "Log.out"),
    log:
        "logs/rsem_ref.log"
    threads: 30
    params:
        threads=30
    shell:
        """
        source /etc/profile.d/modules.sh

        exec 2> {log}

        module load {input.star_path}
        module load {input.rsem_path}

        wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz
        gunzip -c gencode.v38.annotation.gtf.gz > data/gencode.v38.annotation.gtf


        wget -c --tries=80 --timeout=60 https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh38.primary_assembly.genome.fa.gz
        gunzip -c GRCh38.primary_assembly.genome.fa.gz > data/GRCh38.primary_assembly.genome.fa

        mkdir -p rsem_ref
        rsem-prepare-reference -p {threads} \
        --gtf data/gencode.v38.annotation.gtf \
        --star \
        data/GRCh38.primary_assembly.genome.fa \
        rsem_ref/gencode_v38

        """



rule Align:
    input:
        star_path=STAR,
        index="rsem_ref",
        gtf="data/gencode.v38.annotation.gtf",
        reads=["data/{Sample}"+fq1_suffix, "data/{Sample}"+fq2_suffix],
    output:
        "results/step1/STAR_results/{Sample}/{Sample}_Aligned.toTranscriptome.out.bam",
        "results/step1/STAR_results/{Sample}/{Sample}_ReadsPerGene.out.tab",
        "results/step1/STAR_results/{Sample}/{Sample}_Log.final.out",
    log:
        "logs/star/{Sample}.log"
    threads: 30
    params:
        threads=30
    script:
        "../scripts/starAlign.sh"



rule rsem_quant:
    input:
        rsem_ref=multiext("rsem_ref/gencode_v38", ".chrlist", ".n2g.idx.fa", ".transcripts.fa", ".grp", ".seq", ".idx.fa", ".ti"),
        rsem_path=rsem,
        counts="results/step1/STAR_results/{Sample}/{Sample}_Aligned.toTranscriptome.out.bam",
    output:
        "results/step1/rsem_results/{Sample}/{Sample}.isoforms.results",
    log:
        "logs/rsem_quant/{Sample}.log"
    threads: 30
    params:
        threads=30
    script:
        "../scripts/rsemQuant.sh"



rule isoform_count:
    input:
        isoform_quants=expand("results/step1/rsem_results/{Sample}/{Sample}.isoforms.results", Sample=Samples),
    output:
        "results/step1/isoform_count_matrix.txt"
    log:
        "logs/isoform_count.log"
    script:
        "../scripts/rsem_isoform_countMatrix.py"
