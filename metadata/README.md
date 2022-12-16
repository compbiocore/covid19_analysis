To run thie pipeline on Oscar, first enter an interact session with enough resources (```interact -t 1:00:00 -n 4``` is plenty for test_run.sh; for actual data with run.sh I have been using ```interact -t 12:00:00 -n 24```)


To enter the singularity container and mount the data directory use:

```
singularity exec -B /gpfs/data/ris3/0_data/ covid12162022_latest.sif bash
```

Once inside the container, simply run:

Test Run (will produce test_results directory):

```
bash _test_run.sh
```
will produce results in test_results directory


Real Data (will produce results directory)
```
bash run.sh
```

As new data is downloaded, the path to the data will need to be updated in the run.sh script

