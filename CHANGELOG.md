# Changelog

## 1.0.0 (2024-08-02)


### Features

* adapted read shifting to pipeline, make bam export default ([f29761d](https://github.com/MPUSP/snakemake-bacterial-riboseq/commit/f29761d024216ec3b546900b846372326043ab79))
* added description to run examples ([4b31486](https://github.com/MPUSP/snakemake-bacterial-riboseq/commit/4b31486ee45b98c2d79ea67c81929b54e7dd5118))
* added report for workflow results ([31bb75a](https://github.com/MPUSP/snakemake-bacterial-riboseq/commit/31bb75a66540a4aaab503322a2f6335059d4e1be))
* added rule to calculate basic gene-wise statistics ([18e316c](https://github.com/MPUSP/snakemake-bacterial-riboseq/commit/18e316c7692ee68318abe6c39c5a0804398d17f3))
* added wf overview figure, closes [#4](https://github.com/MPUSP/snakemake-bacterial-riboseq/issues/4) ([dcc8a1a](https://github.com/MPUSP/snakemake-bacterial-riboseq/commit/dcc8a1a21736149ad97e7e382cbf91f0ac58ca12))
* bring genome work on par with current state of small-orfs pipe ([8bd026c](https://github.com/MPUSP/snakemake-bacterial-riboseq/commit/8bd026c6e2492f4c997e6bcdf612a26e74144eb4))
* prepared wf to run test data; removed unused conf options; snakefmt + linting ([d351622](https://github.com/MPUSP/snakemake-bacterial-riboseq/commit/d35162222fc4e3008a366d7d1b8df16cdc35c168))


### Bug Fixes

* added dev branch to gh actions workflow ([051de2c](https://github.com/MPUSP/snakemake-bacterial-riboseq/commit/051de2c9e9d9b4b0062b427b6e17bca3224353d8))
* added missing channels for env definitions ([5184b38](https://github.com/MPUSP/snakemake-bacterial-riboseq/commit/5184b3819a58be808308635357d0da63d29f9da9))
* added missing options for gff source type ([2ee01d5](https://github.com/MPUSP/snakemake-bacterial-riboseq/commit/2ee01d523a884b94e377780c5aa4c295786432e6))
* bug in PCA plotting vars ([bfa3cd6](https://github.com/MPUSP/snakemake-bacterial-riboseq/commit/bfa3cd69598d247f5d34ce9d4e10b9703eb04d82))
* changed protocol name ([6bc4dc6](https://github.com/MPUSP/snakemake-bacterial-riboseq/commit/6bc4dc6a989a29af9cfb79cf64198a0ffaa69c26))
* changed protocol name in README ([1702026](https://github.com/MPUSP/snakemake-bacterial-riboseq/commit/17020268796bd6561bed705a80382454f7901d04))
* cleaned plot_mapping_length script and added logging output, closes [#3](https://github.com/MPUSP/snakemake-bacterial-riboseq/issues/3) ([5fb31af](https://github.com/MPUSP/snakemake-bacterial-riboseq/commit/5fb31afc1f339779ad6cded58308e09f1855a028))
* make fastqc run quietly (no progress printed to terminal) ([930b55c](https://github.com/MPUSP/snakemake-bacterial-riboseq/commit/930b55c0d68c22c11125d7ea4ff07a116e8e3f6a))
* remove sORF related steps ([1c2ff00](https://github.com/MPUSP/snakemake-bacterial-riboseq/commit/1c2ff00c16ba8239e63d08b5d70b325bd73f218b))
* removed run_workflow.sh script. ([0742005](https://github.com/MPUSP/snakemake-bacterial-riboseq/commit/074200535f0acb93f0ba496530a70eb4f9baff03))
* removed unnecessary params ([b357ab9](https://github.com/MPUSP/snakemake-bacterial-riboseq/commit/b357ab97fff99744fa9bd1049703d9c83f52ab12))
* snakefmt error ([a4dd408](https://github.com/MPUSP/snakemake-bacterial-riboseq/commit/a4dd4085d2761689955917ada0a6f0c4baf19e28))
* updated release please workflow ([b71b2e4](https://github.com/MPUSP/snakemake-bacterial-riboseq/commit/b71b2e49941988fb83f4a73167974b1b1b21b3dc))
