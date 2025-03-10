## Directory structure
```
.
├── config            # Contains sample sheet (sample_sheet.csv) and config file (config.yaml)
├── rules             # Snakemake rules
├── scripts           # Scripts to run each step 
├── README.md
├── Snakefile         # Snakemake workflow
└── TrPLet_env.yml    # env file 

```

## Usage
### Step1: Download the snakemake_TrPLet package
### Step2: Construct env from TrPLet.yml
  ```
  cd snakemake_TrPLet
  conda env create -f TrPLet_env.yml -n TrPLet
  ```
### Step3: Prepare inputs
  * Create a folder named ```data``` within ```snakemake_TrPLet```
  * Copy RNA-seq files into folder ```data```

### Step4: Prepare ```config/sample_sheet.csv```
  * Paired-end data is assumed.
  * The following data file formats are accommodated: **.fastq.gz, .fq.gz, .fastq, .fq, .bam**
  * The ```config/sample_sheet.csv``` file is an example sample sheet, modify it as needed so that 
      * 1st column: sample names
      * 2nd column: fq1 file names (if using .bam as input, put the name of .fastq.gz files generated from .bam)
      * 3rd column: fq2 file names (if using .bam as input, put the name of .fastq.gz files generated from .bam)
      * 4th column: RNA-seq file format - choose between **bam** and **fq** (all inputs need to be in the same file format)
      * 5th column: type of inputs - choose between **Non-TCGA** or **Cellline** (all inputs need to be in the same type)
      * 6th column: if the inputs belong to **Non-TCGA**, put the **integer code for lineage** based on the table below; if the inputs belong to **Cellline**, put **NA**
        
        ![Image](https://github.com/user-attachments/assets/affe87db-9482-4333-bc9e-1259a9061b7d)
    
    **Each column is separated by one comma.** 
    **The fq1 & fq2 file names must contain the full sample names.**
  
  Example 1 - .bam, None-TCGA: 
  ```
  SRR13786145,SRR13786145_1.fastq.gz,SRR13786145_2.fastq.gz,bam,None-TCGA,11
  SRR13786146,SRR13786146_1.fastq.gz,SRR13786146_2.fastq.gz,bam,None-TCGA,11
  SRR13786147,SRR13786147_1.fastq.gz,SRR13786147_2.fastq.gz,bam,None-TCGA,11
  ```

  Example 2 - .fq.gz, Cellline: 
  ```
  FUUR-1_rep1,FUUR-1_rep1_1.fq.gz,FUUR-1_rep1_2.fq.gz,fq,Cellline,NA
  FUUR-1_rep2,FUUR-1_rep2_1.fq.gz,FUUR-1_rep2_2.fq.gz,fq,Cellline,NA
  FUUR-1_rep3,FUUR-1_rep3_1.fq.gz,FUUR-1_rep3_2.fq.gz,fq,Cellline,NA
  ```
      
### Step5: On the job submission node (argos-qsub1), run snakemake pipeline in the background without requiring terminal connection all the time (since the analysis takes long)
  ```
  screen -S TrPLet
  conda activate TrPLet # activate TrPLet environment
  
  snakemake --unlock
  snakemake --executor cluster-generic --jobs 20 --latency-wait 60 --cluster-generic-submit-cmd "qsub -l h_vmem=256G, -pe pvm 16 -o $HOME/snakemake_TrPLet/joblogs/ -e $HOME/snakemake_TrPLet/joblogs/"

  ```
  * Press ```Ctrl+A``` and then ```D``` to detach the screen as needed. Use ```screen -ls``` to see the list of screens. Now closing the terminal or losing connection to the cluster should not interrupt the snakemake job.
    

### Output
* The results are saved in the folder ```snakemake_TrPLet/results```.
    * For Non-TCGA workflow, the final result is stored as ```snakemake_TrPLet/results/step5_to_7/dep_predictions_Zscoredlog2TPMp1_bestMODEL_all_genes_geq0.2.csv```.
    * For Cellline workflow, the final result is stored as ```snakemake_TrPLet/results/step3/dep_predictions_Zscoredlog2TPMp1_bestMODEL_all_genes_geq0.2.csv```.
    
**The above commands are designed for the DFCI argos cluster. Modifications may be required to adapt the pipeline for use on other clusters.**
