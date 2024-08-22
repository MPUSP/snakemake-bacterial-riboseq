# snakemake-bacterial-riboseq

## Usage

The usage of this workflow is described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/?usage=MPUSP%2Fsnakemake-bacterial-riboseq).

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) repository and its DOI (see above).

## Workflow overview

This workflow is a best-practice workflow for the analysis of ribosome footprint sequencing (Ribo-Seq) data.
The workflow is built using [snakemake](https://snakemake.readthedocs.io/en/stable/) and consists of the following steps:

 1. Obtain genome database in `fasta` and `gff` format (`python`, [NCBI Datasets](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/))
    1. Using automatic download from NCBI with a `RefSeq` ID
    2. Using user-supplied files
 2. Check quality of input sequencing data (`FastQC`)
 3. Cut adapters and filter by length and/or sequencing quality score (`cutadapt`)
 4. Deduplicate reads by unique molecular identifier (UMI, `umi_tools`)
 5. Map reads to the reference genome (`STAR aligner`)
 6. Sort and index for aligned seq data (`samtools`)
 7. Filter reads by feature type (`bedtools`)
 8. Generate summary report for all processing steps (`MultiQC`)
 9. Shift ribo-seq reads according to the ribosome's P-site alignment (`R`, `ORFik`)
 10. Calculate basic gene-wise statistics such as RPKM (`R`, `ORFik`)
 11. Return report as HTML and PDF files (`R markdown`, `weasyprint`)

If you want to contribute, report issues, or suggest features, please get in touch on [github](https://github.com/MPUSP/snakemake-bacterial-riboseq).

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

**Note:**

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
| RPF-RTP1 | RPF-RTP   | 1         | McGlincy | data        | RPF-RTP1_R1_001.fastq.gz |
| RPF-RTP2 | RPF-RTP   | 2         | McGlincy | data        | RPF-RTP2_R1_001.fastq.gz |

Some configuration parameters of the pipeline may be specific for your data and library preparation protocol. The options should be adjusted in the `config.yml` file. For example:

- Minimum and maximum read length after adapter removal (see option `cutadapt: default`). Here, the test data has a minimum read length of 15 + 7 = 22 (2 nt on 5'end + 5 nt on 3'end), and a maximum of 45 + 7 = 52.
- Unique molecular identifiers (UMIs). For example, the protocol by [McGlincy & Ingolia, 2017](https://doi.org/10.1016/J.YMETH.2017.05.028) creates a UMI that is located on both the 5'-end (2 nt) and the 3'-end (5 nt). These UMIs are extracted with `umi_tools` (see options `umi_extraction: method` and `pattern`).

### Execution

To run the workflow from command line, change the working directory.

```bash
cd path/to/snakemake-bacterial-riboseq
```

Adjust options in the default config file `config/config.yml`.
Before running the entire workflow, you can perform a dry run using:

```bash
snakemake --dry-run
```

To run the complete workflow with test files using **`conda`**, execute the following command. The definition of the number of compute cores is mandatory.

```bash
snakemake --cores 10 --use-conda --directory .test
```

### Parameters

This table lists all parameters that can be used to run the workflow.

| parameter              | type | details                                     | default                                      |
| ---------------------- | ---- | ------------------------------------------- | -------------------------------------------- |
| **samplesheet**        |      |                                             |                                              |
| path                   | str  | path to samplesheet, mandatory              | "config/samples.tsv"                         |
| **get_genome**         |      |                                             |                                              |
| database               | str  | one of `manual`, `ncbi`                     | `ncbi`                                       |
| assembly               | str  | RefSeq ID                                   | `GCF_000006785.2`                            |
| fasta                  | str  | optional path to fasta file                 | Null                                         |
| gff                    | str  | optional path to gff file                   | Null                                         |
| gff_source_type        | str  | list of name/value pairs for GFF source     | see config file                              |
| **cutadapt**           |      |                                             |                                              |
| fivep_adapter          | str  | sequence of the 5' adapter                  | Null                                         |
| threep_adapter         | str  | sequence of the 3' adapter                  | `ATCGTAGATCGGAAGAGCACACGTCTGAA`              |
| default                | str  | additional options passed to `cutadapt`     | [`-q 10 `, `-m 22 `, `-M 52`, `--overlap=3`] |
| **umi_extraction**     |      |                                             |                                              |
| method                 | str  | one of `string` or `regex`, see manual      | `regex`                                      |
| pattern                | str  | string or regular expression                | `^(?P<umi_0>.{5}).*(?P<umi_1>.{2})$`         |
| **umi_dedup**          |      |                                             |                                              |
| options                | str  | default options for deduplication           | see config file                              |
| **star**               |      |                                             |                                              |
| index                  | str  | location of genome index; if Null, is made  | Null                                         |
| genomeSAindexNbases    | num  | length of pre-indexing string, see STAR man | 9                                            |
| multi                  | num  | max number of loci read is allowed to map   | 10                                           |
| sam_multi              | num  | max number of alignments reported for read  | 1                                            |
| intron_max             | num  | max length of intron; 0 = automatic choice  | 1                                            |
| default                | str  | default options for STAR aligner            | see config file                              |
| **extract_features**   |      |                                             |                                              |
| biotypes               | str  | biotypes to exclude from mapping            | [`rRNA`, `tRNA`]                             |
| CDS                    | str  | CDS type to include for mapping             | [`protein_coding`]                           |
| **bedtools_intersect** |      |                                             |                                              |
| defaults               | str  | remove hits, sense strand, min overlap 20%  | [`-v `, `-s `, `-f 0.2`]                     |
| **annotate_orfs**      |      |                                             |                                              |
| window_size            | num  | size of 5'-UTR added to CDS                 | 30                                           |
| **shift_reads**        |      |                                             |                                              |
| window_size            | num  | start codon window to determine shift       | 30                                           |
| read_length            | num  | size range of reads to use for shifting     | [27, 45]                                     |
| end_alignment          | str  | end used for alignment of RiboSeq reads     | `3prime`                                     |
| shift_table            | str  | optional table with offsets per read length | Null                                         |
| export_bigwig          | str  | export shifted reads as bam file            | True                                         |
| export_ofst            | str  | export shifted reads as ofst file           | False                                        |
| skip_shifting          | str  | skip read shifting entirely                 | False                                        |
| skip_length_filter     | str  | skip filtering reads by length              | False                                        |
| **multiqc**            |      |                                             |                                              |
| config                 | str  | path to multiqc config                      | `config/multiqc_config.yml`                  |
| **report**             |      |                                             |                                              |
| export_figures         | bool | export figures as `.svg` and `.png`         | True                                         |
| export_dir             | str  | sub-directory for figure export             | `figures/`                                   |
| figure_width           | num  | standard figure width in px                 | 875                                          |
| figure_height          | num  | standard figure height in px                | 500                                          |
| figure_resolution      | num  | standard figure resolution in dpi           | 125                                          |
