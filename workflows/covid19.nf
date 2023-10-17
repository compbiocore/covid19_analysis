#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.download_reads = false
params.run_analysis = true

process downloadGISAID {
  container 'cowmoo/covid_pipeline:latest'

  publishDir "$params.out_dir/gisaid/", mode: 'copy', overwrite: false

  output:
    path "gisaid.fasta", emit: fasta
    path "gisaid.csv", emit: csv
    path "sra_run.txt", emit: runs

  script:
    """
    export GISAIDR_USERNAME='${params.username}'
    export GISAIDR_PASSWORD='${params.password}'
    Rscript /data/gisaid_download.R
    """
}

process runAnalysisPipeline {
    container 'ericsalomaki/covid_new_pango:05092023'

    publishDir "$params.out_dir/analysis/", mode: 'copy', overwrite: false

    input:
      file(gisaid_fasta)

    output:
        path "*.*"

    script:
     """
        set -e

        day=$(date "+%Y%m%d")
        pth=$(pwd)

        mkdir -p ${pth}/3_results/${day}
        cd ${pth}/3_results/${day}

        # Run pangolin
        pangolin ${gisaid_fasta} -o ${pth}/3_results/${day}/pangolin --alignment --no-temp

        # Run nextclade
        nextclade dataset get -n sars-cov-2 -o nextclade_dataset
        nextclade run \
            --input-dataset nextclade_dataset \
            --output-json ${pth}/3_results/${day}/nextclade.json \
            --output-csv ${pth}/3_results/${day}/nextclade.csv \
            --output-tsv ${pth}/3_results/${day}/nextclade.tsv \
            --output-tree ${pth}/3_results/${day}/nextclade.auspice.json \
            --input-qc-config ${pth}/2_metadata/qcRulesConfig.json \
            ${gisaid_fasta}
            >${pth}/3_results/${day}/nextclade.log

        # Run nextalign
        nextalign run \
            --genemap=${pth}/2_metadata/genemap.gff \
            --genes=E,M,N,ORF10,ORF14,ORF1a,ORF1b,ORF3a,ORF6,ORF7a,ORF7b,ORF8,ORF9b,S \
          --output-all=${pth}/3_results/${day}/nextalign \
            --input-ref=${pth}/2_metadata/reference.fasta \
            ${gisaid_fasta}

        python3 ${pth}/1_scripts/qc.py ${day} ${pth}
        Rscript ${pth}/1_scripts/figures.R ${day} ${pth}
        Rscript ${pth}/1_scripts/RISHL_only_plots.R ${day} ${pth}
        Rscript ${pth}/1_scripts/new_perc_mutations.R ${day} ${pth}
        rm ${pth}/3_results/${day}/Rplots.pdf
        echo "Initial Analyses Complete, Running IQtree"
     """
}

process downloadSRA {
  container 'cowmoo/covid_pipeline:latest'

  publishDir "$params.out_dir/gisaid/", mode: 'copy', overwrite: false

  input:
    file(sra_run)

  output:
    path "*.fastq"

  script:
    """
    #cat ${sra_run} | xargs fastq-dump --split-files
    sed -n 1,5p ${sra_run} | xargs fastq-dump --split-files
    """
}


workflow {
  gisaid_result = downloadGISAID()

  if (params.download_reads) {
    downloadSRA(gisaid_result.runs)
  }

  if (params.run_analysis) {
    runAnalysisPipeline(gisaid_result.fasta)
  }
}