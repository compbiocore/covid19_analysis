set -e


mkdir -p results

# Run pangolin
pangolin /gpfs/data/ris3/0_data/gisaid_20220926/sequenceData.fasta -o results/pangolin --alignment --no-temp

# Run nextclade
nextclade dataset get -n sars-cov-2 -o nextclade_dataset
nextclade run \
	--input-dataset 'nextclade_dataset' \
	--output-json 'results/nextclade.json' \
	--output-csv 'results/nextclade.csv' \
	--output-tsv 'results/nextclade.tsv' \
	--output-tree 'results/nextclade.auspice.json' \
	--input-qc-config 'src/qcRulesConfig.json' \
	/gpfs/data/ris3/0_data/gisaid_20220926/sequenceData.fasta \
	>results/nextclade.log

# Run nextalign
nextalign run \
	--genemap=src/genemap.gff \
	--genes=E,M,N,ORF10,ORF14,ORF1a,ORF1b,ORF3a,ORF6,ORF7a,ORF7b,ORF8,ORF9b,S \
	--output-all=results/nextalign \
	--input-ref=src/reference.fasta \
	/gpfs/data/ris3/0_data/gisaid_20220926/sequenceData.fasta
