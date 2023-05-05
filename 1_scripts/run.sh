set -e

day=$(date "+%Y%m%d")

mkdir -p ../3_results/${day}
cd ../3_results/${day}

# Run pangolin
pangolin $1 -o ../3_results/${day}/pangolin --alignment --no-temp

# Run nextclade
nextclade dataset get -n sars-cov-2 -o nextclade_dataset
nextclade run \
	--input-dataset nextclade_dataset \
	--output-json ../3_results/${day}/nextclade.json \
	--output-csv ../3_results/${day}/nextclade.csv \
	--output-tsv ../3_results/${day}/nextclade.tsv \
	--output-tree ../3_results/${day}/nextclade.auspice.json \
	--input-qc-config ../2_metadata/qcRulesConfig.json \
	$1
	>../3_results/${day}/nextclade.log

# Run nextalign
nextalign run \
	--genemap=../2_metadata/genemap.gff \
	--genes=E,M,N,ORF10,ORF14,ORF1a,ORF1b,ORF3a,ORF6,ORF7a,ORF7b,ORF8,ORF9b,S \
	--output-all=../3_results/${day}/nextalign \
	--input-ref=../2_metadata/reference.fasta \
	$1

python ../1_scripts/qc.py
Rscript ../1_scripts/figures.R
Rscript ../1_scripts/num-voc-voi.R
Rscript ../1_scripts/new_perc_mutations.R
