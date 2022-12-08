#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=20G
#SBATCH -t 12:00:00

set -e
BIN=bin

mkdir -p results

# Run pangolin
pangolin ri_sequences.fa -o results/pangolin --alignment --no-temp

# Run nextclade
$BIN/nextclade dataset get -n sars-cov-2 -o nextclade_dataset
$BIN/nextclade run \
	--input-dataset 'nextclade_dataset' \
	--output-json 'results/nextclade.json' \
	--output-csv 'results/nextclade.csv' \
	--output-tsv 'results/nextclade.tsv' \
	--output-tree 'results/nextclade.auspice.json' \
	--input-qc-config 'src/qcRulesConfig.json' \
  ri_sequences.fa \
>results/nextclade.log

# Run nextalign
$BIN/nextalign run \
	--genemap=src/genemap.gff \
	--genes=E,M,N,ORF10,ORF14,ORF1a,ORF1b,ORF3a,ORF6,ORF7a,ORF7b,ORF8,ORF9b,S \
	--output-all=results/nextalign \
	--input-ref=src/reference.fasta \
	ri_sequences.fa