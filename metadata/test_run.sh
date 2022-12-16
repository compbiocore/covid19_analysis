set -e

mkdir -p test_results

# Run pangolin
pangolin /gpfs/data/ris3/0_data/gisaid_20220926/singTest/test_sequences.fa -o test_results/pangolin --alignment --no-temp

# Run nextclade
nextclade dataset get -n sars-cov-2 -o nextclade_dataset
nextclade run \
	--input-dataset 'nextclade_dataset' \
	--output-json 'test_results/nextclade.json' \
	--output-csv 'test_results/nextclade.csv' \
	--output-tsv 'test_results/nextclade.tsv' \
	--output-tree 'test_results/nextclade.auspice.json' \
	--input-qc-config 'src/qcRulesConfig.json' \
	/gpfs/data/ris3/0_data/gisaid_20220926/singTest/test_sequences.fa \
	> test_results/nextclade.log

# Run nextalign
nextalign run \
	--genemap=src/genemap.gff \
	--genes=E,M,N,ORF10,ORF14,ORF1a,ORF1b,ORF3a,ORF6,ORF7a,ORF7b,ORF8,ORF9b,S \
	--output-all=test_results/nextalign \
	--input-ref=src/reference.fasta \
	/gpfs/data/ris3/0_data/gisaid_20220926/singTest/test_sequences.fa
