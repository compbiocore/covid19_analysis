#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process downloadGISAID {
  container 'cowmoo/covid_pipeline:latest'

  publishDir "$params.out_dir/gisaid/", mode: 'copy', overwrite: false

  output:
    path "gisaid.fasta"
    path "gisaid.csv"

  script:
    """
    export GISAIDR_USERNAME='${params.username}'
    export GISAIDR_PASSWORD='${params.password}'
    Rscript /data/gisaid_download.R
    """
}

workflow {
  downloadGISAID()
}