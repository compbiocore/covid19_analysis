#!/bin/env bash
#SBATCH -n 24
#SBATCH -t 1-00:00:00
#SBATCH -J covid_pipeline
#SBATCH --mem=128G

singularity exec -B /gpfs/data/ris3/ /gpfs/data/ris3/1_scripts/covid12162022_latest.sif bash run.sh $1
