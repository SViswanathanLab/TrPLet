rule bam_to_fastq:
    input:
        bam_files="data/{Sample}.bam",
    output:
        "results/bam_to_fastq/{Sample}/{Sample}_sorted.bam",
        "results/bam_to_fastq/{Sample}/{Sample}_R1.fastq.gz",
        "results/bam_to_fastq/{Sample}/{Sample}_R2.fastq.gz",
    log:
        "logs/bam_to_fastq/{Sample}.log"
    threads: 30
    params:
        threads=30,
        conda_env=config['conda_env'],
    script:
        "../scripts/bam_to_fastq.sh"

rule rsem_ref:
    output:
        "data/gencode.v38.annotation.gtf",
        "data/GRCh38.primary_assembly.genome.fa",
        directory("rsem_ref"),
        multiext("rsem_ref/gencode_v38", ".chrlist", ".n2g.idx.fa", ".transcripts.fa", ".grp", ".seq", ".idx.fa", ".ti"),
    log:
        "logs/rsem_ref.log"
    threads: 30
    params:
        threads=30,
        conda_env=config['conda_env'],
    shell:
        """
        source {params.conda_env} TrPLet
        #source /etc/profile.d/modules.sh

        exec 2> {log}

        wget -c http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz
        gunzip -c gencode.v38.annotation.gtf.gz > data/gencode.v38.annotation.gtf


        wget -c --tries=80 --timeout=60 http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh38.primary_assembly.genome.fa.gz
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
        index="rsem_ref",
        gtf="data/gencode.v38.annotation.gtf",
        reads=["results/bam_to_fastq/{Sample}/{Sample}"+fq1_suffix, "results/bam_to_fastq/{Sample}/{Sample}"+fq2_suffix],
    output:
        "results/step1/STAR_results/{Sample}/{Sample}_Aligned.toTranscriptome.out.bam",
        "results/step1/STAR_results/{Sample}/{Sample}_Log.final.out",
    log:
        "logs/star/{Sample}.log"
    threads: 30
    params:
        threads=30,
        conda_env=config['conda_env'],
    script:
        "../scripts/starAlign.sh"



rule rsem_quant:
    input:
        rsem_ref=multiext("rsem_ref/gencode_v38", ".chrlist", ".n2g.idx.fa", ".transcripts.fa", ".grp", ".seq", ".idx.fa", ".ti"),
        counts="results/step1/STAR_results/{Sample}/{Sample}_Aligned.toTranscriptome.out.bam",
    output:
        "results/step1/rsem_results/{Sample}/{Sample}.isoforms.results",
        "results/step1/rsem_results/{Sample}/{Sample}.genes.results",
    log:
        "logs/rsem_quant/{Sample}.log"
    threads: 30
    params:
        threads=30,
        conda_env=config['conda_env'],
    script:
        "../scripts/rsemQuant.sh"
