set -e

day=$(date "+%Y%m%d")

mkdir -p /gpfs/data/ris3/3_results/${day}
cd /gpfs/data/ris3/3_results/${day}

# Run pangolin
pangolin $1 -o /gpfs/data/ris3/3_results/${day}/pangolin --alignment --no-temp

# Run nextclade
nextclade dataset get -n sars-cov-2 -o nextclade_dataset
nextclade run \
	--input-dataset nextclade_dataset \
	--output-json /gpfs/data/ris3/3_results/${day}/nextclade.json \
	--output-csv /gpfs/data/ris3/3_results/${day}/nextclade.csv \
	--output-tsv /gpfs/data/ris3/3_results/${day}/nextclade.tsv \
	--output-tree /gpfs/data/ris3/3_results/${day}/nextclade.auspice.json \
	--input-qc-config /gpfs/data/ris3/2_metadata/qcRulesConfig.json \
	$1
	>/gpfs/data/ris3/3_results/${day}/nextclade.log

# Run nextalign
nextalign run \
	--genemap=/gpfs/data/ris3/2_metadata/genemap.gff \
	--genes=E,M,N,ORF10,ORF14,ORF1a,ORF1b,ORF3a,ORF6,ORF7a,ORF7b,ORF8,ORF9b,S \
	--output-all=/gpfs/data/ris3/3_results/${day}/nextalign \
	--input-ref=/gpfs/data/ris3/2_metadata/reference.fasta \
	$1
