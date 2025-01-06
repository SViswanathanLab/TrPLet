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
  cd snakemake_TrPLet
  conda env create -f TrPLet.yml -n TrPLet
  ```
* Step3: Prepare inputs
  * Create a folder named ```data``` within ```snakemake_TrPLet```
  * Copy fq1&fq2 files into folder ```data```
  * Prepare ```config/sample_sheet.tsv```
    * Paired-end data is assumed.
    * gzipped and non-gzipped RNAseq data formats are accommodated: **.fastq.gz, .fastq**
    * The ```config/sample_sheet.tsv``` file is an example sample sheet.
    * Modify ```config/sample_sheet.tsv``` to have the first column consisting of sample names, the second column consisting of fq1 file names, and the third column consisting of fq2 file names. Each column is separated by **one space**. 
    * The fq1 & fq2 file names must contain the full sample names.
      
* Step4: Run snakemake pipeline in the background without requiring terminal connection all the time (since the analysis takes long)
  ```
  screen -S TrPLet
  conda activate TrPLet # activate TrPLet environment
  snakemake --unlock
  snakemake --executor cluster-generic --jobs 20 --latency-wait 60 --cluster-generic-submit-cmd "qsub -l h_vmem=256G, -pe pvm 32 -o $HOME/snakemake_TrPLet/joblogs/ -e $HOME/snakemake_TrPLet/joblogs/"
  ```
  * Press ```Ctrl+A``` and then ```D``` to detach the screen as needed. Use ```screen -ls``` to see the list of screens. Now closing the terminal or losing connection to the cluster should not interrupt the snakemake job. 
    
  
  
