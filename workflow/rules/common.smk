# import basic packages
import pandas as pd
from snakemake.utils import validate
from os import path


samples = (
    pd.read_csv(config["samplesheet"], sep="\t", dtype={"sample": str})
    .set_index("sample", drop=False)
    .sort_index()
)

# TODO: write validation schema
# validate(SAMPLES, schema="../config/schemas/samples.schema.yml")


def get_final_output():
    targets = []
    targets.append("results/multiqc/multiqc_report.html")
    targets.append("results/feature_stats/feature_stats.csv")
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
    if wildcards.status == "raw":
        return expand(
            "{input_dir}/{sample}",
            input_dir=samples.loc[wildcards.sample]["data_folder"],
            sample=samples.loc[wildcards.sample]["fq1"],
        )
    if wildcards.status == "clipped":
        return expand("results/clipped/{sample}.fastq.gz", sample=samples.index)


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
