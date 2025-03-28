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
  * Copy RNA-seq files (.bam or .fq or gene-level matrix) into folder ```data```
  * Download files [here](https://www.dropbox.com/scl/fo/wmfhqzzspjfqiezqhajk5/ADAUJUQgCYhnCuI3aOLxHYs?rlkey=v0fmkxyn0cqwm1o9xsqpfdwxk&st=7sc4rj0l&dl=0) into folder ```data```

### Step4: Prepare ```config/sample_sheet.csv```
  * Paired-end data is assumed.
  * The following data file formats are accommodated: **.fastq.gz, .fq.gz, .fastq, .fq, .bam, gene_count_matrix.txt**
      * If the input is gene-level count matrix from RNA-seq, we need to have the file named gene_count_matrix.txt and in the following format:
        ```
        gene_id  effective_length  Sample1  Sample2	
        ENSG00000249352.3_7SK  #  #  #
        ENSG00000249352.3_7SK  #  #  #
        ENSG00000249352.3_7SK  #  #  #
        ENSG00000249352.3_7SK  #  #  #
        ```
        **The names of the columns and file must strictly follow the required format.**
  * The ```config/sample_sheet.csv``` file is an example sample sheet, modify it as needed so that 
      * 1st column: sample names
      * 2nd column: fq1 file names (if using .bam as input, put the name of .fastq.gz files generated from .bam)
      * 3rd column: fq2 file names (if using .bam as input, put the name of .fastq.gz files generated from .bam)
      * 4th column: RNA-seq file format - choose from **bam, fq, or matrix** (all inputs need to be in the same file format)
      * 5th column: type of inputs - choose between **Tumor** or **Cellline** (all inputs need to be in the same type)
      * 6th column: if the inputs belong to **Tumor**, put the **integer code for lineage** based on the table below; if the inputs belong to **Cellline**, put **NA**
        
        ![Image](https://github.com/user-attachments/assets/affe87db-9482-4333-bc9e-1259a9061b7d)
    
    **Each column is separated by one comma.** 
    **The fq1 & fq2 file names must contain the full sample names.**
    **Ensure there is no empty row at the bottom.**
  
  Example 1 - .bam, Tumor: 
  ```
  SRR13786145,SRR13786145_1.fastq.gz,SRR13786145_2.fastq.gz,bam,Tumor,11
  SRR13786146,SRR13786146_1.fastq.gz,SRR13786146_2.fastq.gz,bam,Tumor,11
  SRR13786147,SRR13786147_1.fastq.gz,SRR13786147_2.fastq.gz,bam,Tumor,11
  ```

  Example 2 - .fq.gz, Cellline: 
  ```
  FUUR-1_rep1,FUUR-1_rep1_1.fq.gz,FUUR-1_rep1_2.fq.gz,fq,Cellline,NA
  FUUR-1_rep2,FUUR-1_rep2_1.fq.gz,FUUR-1_rep2_2.fq.gz,fq,Cellline,NA
  FUUR-1_rep3,FUUR-1_rep3_1.fq.gz,FUUR-1_rep3_2.fq.gz,fq,Cellline,NA
  ```

  Example 3 - matrix, Tumor: 
  ```
  LS-GD-0030,LS-GD-0030_1.fastq.gz,LS-GD-0030_2.fastq.gz,matrix,Tumor,6
  LS-GD-0100,LS-GD-0100_1.fastq.gz,LS-GD-0100_2.fastq.gz,matrix,Tumor,6
  LS-GD-0031,LS-GD-0031_1.fastq.gz,LS-GD-0031_2.fastq.gz,matrix,Tumor,6
  LS-GD-0101,LS-GD-0101_1.fastq.gz,LS-GD-0101_2.fastq.gz,matrix,Tumor,6
  ```
      
### Step5: On the job submission node (argos-qsub1), run snakemake pipeline in the background without requiring terminal connection all the time (since the analysis takes long)
  ```
  screen -S TrPLet
  conda activate TrPLet # activate TrPLet environment
  
  snakemake --unlock
  snakemake --executor cluster-generic --jobs 20 --latency-wait 60 --cluster-generic-submit-cmd "qsub -l h_vmem=256G, -pe pvm 16 -o $HOME/snakemake_TrPLet/joblogs/ -e $HOME/snakemake_TrPLet/joblogs/"

  ```
  * Press ```Ctrl+A``` and then ```D``` to detach the screen as needed. Use ```screen -ls``` to see the list of screens. Now closing the terminal or losing connection to the cluster should not interrupt the snakemake job.
  * Adjust the memory or number of nodes to request as needed. May also need to adjust the chunk size in step3_batchcorrect.R when processing large number of samples.
  * If the pipeline run is interrupted, need to rerun the following to resume from the point where it stopped.
   ```
   screen -r SCREENID
   snakemake --unlock
   snakemake --executor cluster-generic --jobs 20 --latency-wait 60 --cluster-generic-submit-cmd "qsub -l h_vmem=256G, -pe pvm 16 -o $HOME/snakemake_TrPLet/joblogs/ -e $HOME/snakemake_TrPLet/joblogs/"

   ```   

### Output
* The results are saved in the folder ```snakemake_TrPLet/results```.
    * For Tumor workflow, the final result is stored as ```snakemake_TrPLet/results/step5_to_7/dep_predictions_Zscoredlog2TPMp1_bestMODEL_all_genes_geq0.2.csv```.
    * For Cellline workflow, the final result is stored as ```snakemake_TrPLet/results/step3/dep_predictions_Zscoredlog2TPMp1_bestMODEL_all_genes_geq0.2.csv```.
    
**The above commands are designed for the DFCI argos cluster. Modifications may be required to adapt the pipeline for use on other clusters.**
