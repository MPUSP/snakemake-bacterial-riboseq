# snakemake-bacterial-riboseq

![Platform](https://img.shields.io/badge/platform-all-green)
[![Snakemake](https://img.shields.io/badge/snakemake-≥7.0.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/MPUSP/snakemake-bacterial-riboseq/workflows/Tests/badge.svg?branch=main)](https://github.com/MPUSP/snakemake-bacterial-riboseq/actions?query=branch%3Amain+workflow%3ATests)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)

---

A Snakemake workflow for the analysis of bacterial riboseq data.

- [snakemake-bacterial-riboseq](#snakemake-bacterial-riboseq)
  - [Usage](#usage)
  - [Workflow overview](#workflow-overview)
  - [Installation](#installation)
    - [Additional tools](#additional-tools)
  - [Running the workflow](#running-the-workflow)
    - [Input data](#input-data)
      - [Reference genome](#reference-genome)
      - [Read data](#read-data)
    - [Execution](#execution)
    - [Parameters](#parameters)
  - [Authors](#authors)
  - [References](#references)

## Usage

The usage of this workflow is described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/?usage=MPUSP%2Fsnakemake-bacterial-riboseq).

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) repository and its DOI (see above).

## Workflow overview

TODO: include first part of the figure here.

## Installation

**Step 1: Clone this repository**

```bash
git clone https://github.com/MPUSP/snakemake-bacterial-riboseq.git
cd snakemake-bacterial-riboseq
```

**Step 2: Install dependencies**

It is recommended to install snakemake and run the workflow with `conda`, `mamba` or `micromamba`.

```bash
# download Miniconda3 installer
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
# install Conda (respond by 'yes')
bash miniconda.sh
# update Conda
conda update -y conda
# install Mamba
conda install -n base -c conda-forge -y mamba
```

**Step 3: Create snakemake environment**

This step creates a new conda environment called `snakemake-bacterial-riboseq`.

```bash
# create new environment with dependencies & activate it
mamba create -c conda-forge -c bioconda -n snakemake-bacterial-riboseq snakemake pandas
conda activate snakemake-bacterial-riboseq
```

### Additional tools

**Important note:**

All other dependencies for the workflow are **automatically pulled as `conda` environments** by snakemake, when running the workflow with the `--use-conda` parameter (recommended).


## Running the workflow

### Input data

#### Reference genome

An NCBI Refseq ID, e.g. `GCF_000006945.2`. Find your genome assembly and corresponding ID on [NCBI genomes](https://www.ncbi.nlm.nih.gov/data-hub/genome/). Alternatively use a custom pair of `*.fasta` file and `*.gff` file that describe the genome of choice.

Important requirements when using custom `*.fasta` and `*.gff` files:

- `*.gff` genome annotation must have the same chromosome/region name as the `*.fasta` file (example: `NC_003197.2`)
- `*.gff` genome annotation must have `gene` and `CDS` type annotation that is automatically parsed to extract transcripts
- all chromosomes/regions in the `*.gff` genome annotation must be present in the `*.fasta` sequence
- but not all sequences in the `*.fasta` file need to have annotated genes in the `*.gff` file

#### Read data

Ribosome footprint sequencing data in `*.fastq.gz` format. The currently supported input data are **single-end, strand-specific reads**. Input data files are supplied via a mandatory table, whose location is indicated in the `config.yml` file (default: `samples.tsv`). The sample sheet has the following layout:

| sample   | condition | replicate | lib_prep | data_folder | fq1                      |
| -------- | --------- | --------- | -------- | ----------- | ------------------------ |
| RPF-RTP1 | RPF-RTP   | 1         | mpusp    | data        | RPF-RTP1_R1_001.fastq.gz |
| RPF-RTP2 | RPF-RTP   | 2         | mpusp    | data        | RPF-RTP2_R1_001.fastq.gz |

Some configuration parameters of the pipeline may be specific for your data and library preparation protocol. The options should be adjusted in the `config.yml` file. For example:

- Minimum and maximum read length after adapter removal (see option `cutadapt: default`). Here, the test data has a minimum read length of 15 + 7 = 22 (2 nt on 5'end + 5 nt on 3'end), and a maximum of 45 + 7 = 52.
- Unique molecular identifiers (UMIs). For example, the protocol by [McGlincy & Ingolia, 2017](https://doi.org/10.1016/J.YMETH.2017.05.028) creates a UMI that is located on both the 5'-end (2 nt) and the 3'-end (5 nt). These UMIs are extracted with `umi_tools` (see options `umi_extraction: method` and `pattern`).

### Execution

To run the workflow from command line, change the working directory.

```bash
cd path/to/snakemake-bacterial-riboseq
```

Adjust the global and module-specific options in the default config file `config/config.yml`.
Before running the entire workflow, you can perform a dry run using:

```bash
snakemake --dry-run
```

To run the complete workflow with test files using **`conda`**, execute the following command. The definition of the number of compute cores is mandatory.

```bash
snakemake --cores 10 --use-conda --directory .test
```

### Parameters


## Authors

- Dr. Rina Ahmed-Begrich
  - Affiliation: [Max-Planck-Unit for the Science of Pathogens](https://www.mpusp.mpg.de/) (MPUSP), Berlin, Germany
  - ORCID profile: https://orcid.org/0000-0002-0656-1795
- Dr. Michael Jahn
  - Affiliation: [Max-Planck-Unit for the Science of Pathogens](https://www.mpusp.mpg.de/) (MPUSP), Berlin, Germany
  - ORCID profile: https://orcid.org/0000-0002-3913-153X
  - github page: https://github.com/m-jahn


Visit the MPUSP github page at https://github.com/MPUSP for more info on this workflow and other projects.

## References

- Essential tools are linked in the top section of this document
- The sequencing library preparation is based on the publication:
> McGlincy, N. J., & Ingolia, N. T. _Transcriptome-wide measurement of translation by ribosome profiling_. Methods, 126, 112–129, **2017**. https://doi.org/10.1016/J.YMETH.2017.05.028.
