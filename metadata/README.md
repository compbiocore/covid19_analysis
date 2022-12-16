To enter the singularity container and mount the data directory use:

```
singularity exec -B /gpfs/data/ris3/0_data/ covid12162022_latest.sif bash
```

Once inside the container, simply run:

```
bash run.sh
```

As new data is downloaded, the path to the data will need to be updated in the run.sh script
