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
mamba env create -c conda-forge -c bioconda -n snakemake-bacterial-riboseq snakemake pandas
conda activate snakemake-bacterial-riboseq
```

### Additional tools

**Important note:**

All other dependencies for the workflow are **automatically pulled as `conda` environments** by snakemake, when running the workflow with the `--use-conda` parameter (recommended).


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
