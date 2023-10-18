# Running Workflow via Nextflow

The following documentation details on how to run the Covid19 analysis pipeline using Nextflow on any computing environment.

## Installation

### 1. Check out Github repo
First, check out the Github repo:

```commandline
git clone https://github.com/compbiocore/covid19_analysis.git
```

### 2. Install Nextflow and Singularity

#### Option A: On Any Computing Environment

If you do not have Singularity already; you can install it by referring to the [Singularity installation guide](https://docs.sylabs.io/guides/3.0/user-guide/installation.html) here.

If you do not have Nextflow already; you can install it by referring to the [Nextflow installation guide](https://www.nextflow.io/docs/latest/getstarted.html#installation) here.

After installing Singularity, ensure that in your Nextflow configuration file, you have enabled Singularity in Nextflow. You can refer to the [Singularity configuration guide](https://www.nextflow.io/docs/edge/container.html#id24) here; or in another words, add the following block in the `nextflow.config` file that Nextflow is sourcing:
```commandline
...
singularity {
    enabled = true
}
```

#### Option B: On Brown OSCAR Computing Environment

If you are on Brown OSCAR computing environment, you can simply install Nextflow and Singularity computing environment by following the [set up instructions here](https://github.com/compbiocore/workflows_on_OSCAR). And then to initialize the Nextflow environment, simply type in:
```commandline
nextflow_start
```


## Running the Nextflow Workflow

Once you have finished installing (or already have the requisites satisfied), you can run the Nextflow pipeline with the following command:

```
cd $PROJECT_REPO
nextflow run $PROJECT_REPO/workflows/covid19.nf \
--output_dir $OUTPUT_DIR --username $GISAID_USER --password='$GISAID_PASSWORD' \
--project_github $PROJECT_REPO
```  

## Output Directory

Below is a brief walk-through and explaination of all the workflow workproducts: 

#### Output 1: GISAID Sequence Files and Metadata

In `$OUTPUT_DIR/gisaid`:
 - `gisaid.fasta`, the sequence containing for all sequences downloaded from GISAID given a certain geolocation (e.g., USA/Rhode Island). 
 - `gisaid.csv`, the GISAID metadata file for all the sequences given the certain geolocation
 - `sra_run.txt`, all of the SRA id's linked to the GISAID sequences in this workflow. 

#### Output 2: Analysis Files

