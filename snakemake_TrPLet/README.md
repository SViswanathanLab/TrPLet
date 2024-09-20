## Directory structure
```
.
├── config            # Contains sample sheet (sample_sheet.tsv) and config file (config.yaml)
├── rules             # Snakemake rules
├── scripts           # Scripts to run each step 
├── README.md
├── Snakefile         # Snakemake workflow
└── TrPLet.yml        # env file 

```

## Usage
* Step1: Download the snakemake_TrPLet package
* Step2: Construct env from TrPLet.yml
  ```
  conda env create -f TrPLet.yml -n TrPLet
  conda activate TrPLet
  ```
* Step3: Prepare inputs
  * Create a folder named ```data``` within ```snakemake_TrPLet```
  * Copy fq1&fq2 files into folder ```data```
  * Prepare ```config/sample_sheet.tsv```
    * Paired-end data is assumed.
    * 3 types of RNAseq data formats are accommodated: **.fastq.gz, .fq.gz, .fastq**
    * The ```config/sample_sheet.tsv``` file is an example sample sheet.
    * Modify ```config/sample_sheet.tsv``` to have the first column consisting of sample names, the second column consisting of fq1 file names, and the third column consisting of fq2 file names. Each column is separated by **one space**. 
    * The fq1 & fq2 file names must contain the full sample names.
* Step4: Run snakemake pipeline
  ```
  cd snakemake_TrPLet
  snakemake --unlock
  snakemake --executor cluster-generic --jobs 20 --latency-wait 60 --cluster-generic-submit-cmd "qsub -l h_vmem=256G, -pe pvm 32 -o $HOME/snakemake_TrPLet/joblogs/ -e $HOME/snakemake_TrPLet/joblogs/"
  ```
  * The running takes long, so it is usual the command execution is interrupted, need to rerun the above commands to generate all results expected. 
  
  
