## Running Pipeline via Oscar Slurm Batch Submission  
  
To run the covid pipeline, navigate to ```/gpfs/data/ris3/1_scripts/``` and run:   
```
sbatch /gpfs/data/ris3/0_data/gisaid_20220926/run_slurm.sh /gpfs/data/ris3/PATH/TO/SEQUENCE/DATA
```  
Results will be produced in ```/gpfs/data/ris3/3_results/${YYYYMMDD}```

A run with ~20,000 input sequences takes roughly 8 hours to complete

  
## Running Pipeline via Oscar Interactive Session

To run thie pipeline in an interact session, first enter a screen `screen -S JOBNAME` and then initiate an interact session with enough resources (`interact -t 24:00:00 -n 24 -m 128G`)
  
Navigate to the `1_scripts` directory:  
```
cd /gpfs/data/ris3/1_scripts
```
  
Enter the singularity container and mount the parent directory:

```
singularity exec -B /gpfs/data/ris3/ /gpfs/data/ris3/1_scripts/covid12162022_latest.sif bash 
```  

Once inside the container, run:

``` 
bash run.sh /gpfs/data/ris3/PATH/TO/SEQUENCE/DATA
```

To leave the screen use `ctl + a + d` and to return use `screen -r JOBNAME`  
  
Results will be produced in `/gpfs/data/ris3/3_results/${YYYYMMDD}`


The `run.sh` script and usage instructions are modified from that on https://github.com/kantorlab/covid-pipeline 
