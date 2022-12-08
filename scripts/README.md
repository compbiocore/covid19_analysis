# Running

The `run.sh` script and usage instructions are modified from the one included in https://github.com/kantorlab/covid-pipeline to account for updates to Nextclade and Nextalign commands

The input file of sequences in FASTA format should be named `ri_sequence.fa` (edit the `run.sh` script to use a different input filename).

The pipeline is designed to be submitted as a batch job to a SLURM cluster:

    sbatch run.sh

If you are running locally instead of on a cluster, simply execute the script in bash:

    bash run.sh

Results are written to the `results` subdirectory.