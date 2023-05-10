set -e

day=$(date "+%Y%m%d")
pth=$(cd ../ && pwd)

mkdir -p ${pth}/3_results/${day}
cd ${pth}/3_results/${day}

# Run pangolin
pangolin $1 -o ${pth}/3_results/${day}/pangolin --alignment --no-temp

# Run nextclade
nextclade dataset get -n sars-cov-2 -o nextclade_dataset
nextclade run \
	--input-dataset nextclade_dataset \
	--output-json ${pth}/3_results/${day}/nextclade.json \
	--output-csv ${pth}/3_results/${day}/nextclade.csv \
	--output-tsv ${pth}/3_results/${day}/nextclade.tsv \
	--output-tree ${pth}/3_results/${day}/nextclade.auspice.json \
	--input-qc-config ${pth}/2_metadata/qcRulesConfig.json \
	$1
	>${pth}/3_results/${day}/nextclade.log

# Run nextalign
nextalign run \
	--genemap=${pth}/2_metadata/genemap.gff \
	--genes=E,M,N,ORF10,ORF14,ORF1a,ORF1b,ORF3a,ORF6,ORF7a,ORF7b,ORF8,ORF9b,S \
	--output-all=${pth}/3_results/${day}/nextalign \
	--input-ref=${pth}/2_metadata/reference.fasta \
	$1

python3 ${pth}/1_scripts/qc.py ${day} ${pth}
Rscript ${pth}/1_scripts/figures.R ${day} ${pth}
Rscript ${pth}/1_scripts/RISHL_only_plots.R ${day} ${pth}
Rscript ${pth}/1_scripts/new_perc_mutations.R ${day} ${pth}

echo "Initial Analyses Complete, Running IQtree"

awk '/^>/{f=!d[$1];d[$1]=1}f' nextalign.aligned.fasta > nextalign.aligned.no_dups.fasta
iqtree -s ${pth}/3_results/${day}/nextalign/nextalign.aligned.no_dups.fasta --prefix ${pth}/3_results/${day}/iqtree2 -st DNA -m GTR+F -t 24 --mem 128G

