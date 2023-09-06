#!/bin/env bash
#SBATCH -n 24
#SBATCH -N 1
#SBATCH -t 10-00:00:00
#SBATCH -J covid_pipeline
#SBATCH --mem=128G

pth=$(cd ../ && pwd)

singularity exec -B /oscar/data/ris3/ ${pth}/1_scripts/covid_new_pango.sif bash run.sh $1
