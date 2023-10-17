# Running Workflow on OSCAR

## Directory Structure

 * **0_data:** is an empty directory in which to download sequneces and metadata from GISAID for analyses.
 * **1_scripts:** contains shell scripts to run the pipeline as reflected in ```/covid19_analysis/1_scripts``` the singularity image can be pulled directly to oscar or your local machine using ```singularity pull covid19.sif docker://ericsalomaki/covid_new_pango:05092023``` from the `1_scripts` directory.
 * **2_metadata:** contains the ```Dockerfile``` that was used to create the container for running the pipeline, a GFF file, QC rules file, and the reference fasta file and genbank file.
 * **3_results** will be created while the pipeline is running and results will be written to ```/covid19_analysis/3_results/${YYYYMMDD}```


## Running Pipeline via Oscar Slurm Batch Submission  
  
To run the covid pipeline, navigate to ```/PATH/TO/CLONED/REPO/covid19_analysis/1_scripts/``` and run:   
```
sbatch run_slurm.sh /ABSOLUTE/PATH/TO/SEQUENCE/DATA/covid_sequences.fasta
```  
Results will be produced in ```/covid19_analysis/3_results/${YYYYMMDD}```

A run with ~20,000 input sequences takes roughly 30 minutes to complete the primary pangolin analyses and produce figures on Oscar with 24 threads and 128G RAM allocated, however the IQ-tree analysis will run for several days. If incomplete, IQ-tree uses checkpoints and therefore the analysis can be continued beyond the allocated time, if necessary.

  
## Running Pipeline via Oscar Interactive Session

To run thie pipeline in an interact session, first enter a screen `screen -S JOBNAME` and then initiate an interact session with enough resources (`interact -t 24:00:00 -n 24 -m 128G`)
  
Navigate to the `1_scripts` directory:  
```
cd /PATH/TO/CLONED/REPO/covid19_analysis/1_scripts
```
  
Enter the singularity container and mount the parent directory:

```
singularity exec -B /ABSOLUTE/PATH/TO/CLONED/REPO/covid19_analysis/ /PATH/TO/CLONED/REPO/covid19_analysis/1_scripts/covid19.sif bash 
```  

Once inside the container, run:

``` 
bash run.sh /ABSOLUTE/PATH/TO/SEQUENCE/DATA/covid_sequences.fasta
```

To leave the screen use `ctl + a + d` and to return use `screen -r JOBNAME`  
  
Results will be produced in `/PATH/TO/CLONED/REPO/covid19_analysis/3_results/${YYYYMMDD}`

## Example Usage for Oscar
```
sbatch /PATH/TO/CLONED/REPO/covid19_analysis/1_scripts/run_slurm.sh /PATH/TO/CLONED/REPO/covid19_analysis/0_data/sequenceData.fasta
```