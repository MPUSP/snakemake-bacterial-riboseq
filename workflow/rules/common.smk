# import basic packages
import pandas as pd
from snakemake.utils import validate
from os import path


# read sample sheet
samples = (
    pd.read_csv(config["samplesheet"], sep="\t", dtype={"sample": str})
    .set_index("sample", drop=False)
    .sort_index()
)


# validate sample sheet and config file
validate(samples, schema="../../config/schemas/samples.schema.yml")
validate(config, schema="../../config/schemas/config.schema.yml")


# define final targets of the workflow
def get_final_output():
    targets = []
    targets.append("results/multiqc/multiqc_report.html")
    targets.append("results/report/report.pdf")
    targets.append(
        expand(
            "results/{mapping_status}/length_dist/{sample}_length_dist.tsv",
            mapping_status=["mapped", "deduplicated", "filtered_bam"],
            sample=samples.index,
        )
    )
    return targets


# get fastq files
def get_fastq(wildcards):
    if hasattr(wildcards, "status"):
        if wildcards.status == "raw":
            return samples.loc[wildcards.sample]["fq1"]
        if wildcards.status == "clipped":
            return f"results/clipped/{wildcards.sample}.fastq.gz"
    else:
        return samples.loc[wildcards.sample]["fq1"]


# get bam files
def get_bam(wildcards):
    if wildcards.mapping_status == "mapped":
        return expand(
            os.path.join("results", "mapped", "{sample}.bam"), sample=wildcards.sample
        )
    if wildcards.mapping_status == "deduplicated":
        return expand(
            os.path.join("results", "deduplicated", "{sample}.bam"),
            sample=wildcards.sample,
        )
    if wildcards.mapping_status == "filtered_bam":
        return expand(
            os.path.join("results", "filtered_bam", "{sample}.bam"),
            sample=wildcards.sample,
        )
