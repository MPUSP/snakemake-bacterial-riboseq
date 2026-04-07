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

| sample   | condition | replicate | fq1                           |
| -------- | --------- | --------- | ----------------------------- |
| RPF-RTP1 | RPF-RTP   | 1         | data/RPF-RTP1_R1_001.fastq.gz |
| RPF-RTP2 | RPF-RTP   | 2         | data/RPF-RTP2_R1_001.fastq.gz |

Some configuration parameters of the pipeline may be specific for your data and library preparation protocol. The options should be adjusted in the `config.yml` file. For example:

- Minimum and maximum read length after adapter removal (see option `cutadapt: default`). Here, the test data has a minimum read length of 15 + 7 = 22 (2 nt on 5'end + 5 nt on 3'end), and a maximum of 45 + 7 = 52.
- Unique molecular identifiers (UMIs). For example, the protocol by [McGlincy &amp; Ingolia, 2017](https://doi.org/10.1016/J.YMETH.2017.05.028) creates a UMI that is located on both the 5'-end (2 nt) and the 3'-end (5 nt). These UMIs are extracted with `umi_tools` (see options `umi_extraction: method` and `pattern`).

Example configuration files for different sequencing protocols can be found in `resources/protocols/`.
