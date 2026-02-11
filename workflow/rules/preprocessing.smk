# ----------------------------------------------------- #
# Riboseq preprocessing:                                #
# ----------------------------------------------------- #


# module to make QC report
# -----------------------------------------------------
rule fastqc:
    input:
        get_fastq,
    output:
        html="results/fastqc_{status}/{sample}_fastqc.html",
        zip="results/fastqc_{status}/{sample}_fastqc.zip",
    message:
        """--- Checking fastq files with FastQC."""
    log:
        "results/fastqc_{status}/log/{sample}.log",
    threads: max(1, int(workflow.cores * 0.25))
    resources:
        mem_mb=1024,
    wrapper:
        "v7.6.0/bio/fastqc"


# module to trim adapters from reads
# -----------------------------------------------------
rule cutadapt:
    input:
        get_fastq,
    output:
        fastq="results/clipped/{sample}.fastq.gz",
        qc="results/clipped/{sample}.qc.txt",
    params:
        adapters=config["cutadapt"]["adapters"],
        extra=config["cutadapt"]["default"],
    message:
        """--- Trim adapters from reads."""
    threads: max(1, int(workflow.cores * 0.25))
    log:
        "results/clipped/log/{sample}.log",
    wrapper:
        "v7.9.0/bio/cutadapt/se"


# module to extract UMIs and attach to read name
# -----------------------------------------------------
rule umi_extraction:
    input:
        fastq="results/clipped/{sample}.fastq.gz",
    output:
        fastq="results/umi_extraction/{sample}.fastq.gz",
    conda:
        "../envs/umitools.yml"
    message:
        """--- Extracting UMIs."""
    params:
        method=config["umi_extraction"]["method"],
        pattern=lambda wc: config["umi_extraction"]["pattern"],
    log:
        path="results/umi_extraction/log/{sample}.log",
    shell:
        "umi_tools extract "
        "--extract-method={params.method} "
        "--bc-pattern='{params.pattern}' "
        "--stdin {input.fastq} "
        "--stdout {output.fastq} > {log.path}"


# module to fetch genome from NCBI or Ensemble
# -----------------------------------------------------
rule get_genome:
    input:
        fasta=lambda wildcards: (
            config["get_genome"]["fasta"]
            if config["get_genome"]["database"] == "manual"
            else []
        ),
        gff=lambda wildcards: (
            config["get_genome"]["gff"]
            if config["get_genome"]["database"] == "manual"
            else []
        ),
    output:
        fasta="results/get_genome/genome.fasta",
        gff="results/get_genome/genome.gff",
        fai="results/get_genome/genome.fasta.fai",
    params:
        database=config["get_genome"]["database"],
        assembly=config["get_genome"]["assembly"],
        gff_source_types=config["get_genome"]["gff_source_type"],
    message:
        "--- Parsing genome GFF and FASTA files"
    log:
        path="results/get_genome/log/get_genome.log",
    wrapper:
        "https://raw.githubusercontent.com/MPUSP/mpusp-snakemake-wrappers/refs/heads/main/get_genome"


# module to map reads to ref genome using STAR aligner
# -----------------------------------------------------
rule create_star_index:
    input:
        genome="results/get_genome/genome.fasta",
    output:
        path=directory("results/get_genome/index"),
    conda:
        "../envs/star.yml"
    message:
        """--- STAR index creation."""
    params:
        index=config["star"]["index"],
        indexNbases=config["star"]["genomeSAindexNbases"],
    log:
        path="results/get_genome/log/star_index.log",
    shell:
        "if [ {params.index} == None ]; then "
        "mkdir {output.path};"
        "STAR --runMode genomeGenerate "
        "--genomeDir {output.path} "
        "--genomeFastaFiles {input.genome} "
        "--genomeSAindexNbases {params.indexNbases} > {log.path}; "
        "rm -f ./Log.out; "
        "else "
        "ln -s {params.index} {output.path}; "
        "echo 'made symbolic link from {params.index} to {output.path}' > {log.path}; "
        "fi;"


# module to map reads to ref genome using STAR aligner
# -----------------------------------------------------
rule star_mapping:
    input:
        fastq=rules.umi_extraction.output,
        genome=rules.create_star_index.output,
    output:
        bam="results/mapped/unsorted/{sample}.bam",
    conda:
        "../envs/star.yml"
    message:
        """--- STAR mapping."""
    params:
        default=config["star"]["default"],
        multi=config["star"]["multi"],
        sam_multi=config["star"]["sam_multi"],
        intron_max=config["star"]["intron_max"],
        outprefix=lambda w, output: f"{os.path.splitext(output.bam)[0]}_",
    log:
        path="results/mapped/log/{sample}.log",
    threads: max(1, int(workflow.cores * 0.25))
    shell:
        "STAR "
        "--runThreadN {threads} "
        "--genomeDir {input.genome} "
        "--readFilesIn {input.fastq} "
        "{params.default} "
        "--outFilterMultimapNmax {params.multi} "
        "--alignIntronMax {params.intron_max} "
        "--outSAMmultNmax {params.sam_multi} "
        "--outFileNamePrefix {params.outprefix} "
        "> {output.bam} 2> {log.path}"


# module to sort and index bam file using samtools
# -----------------------------------------------------
rule mapping_sorted_bam:
    input:
        rules.star_mapping.output.bam,
    output:
        bam="results/mapped/{sample}.bam",
        bai="results/mapped/{sample}.bam.bai",
    conda:
        "../envs/samtools.yml"
    log:
        "results/mapped/log/samtools_{sample}.log",
    message:
        """--- Samtools sort and index bam files."""
    params:
        tmp="results/mapped/sort_{sample}_tmp",
    threads: max(1, int(workflow.cores * 0.25))
    shell:
        "samtools sort -@ {threads} -O bam -T {params.tmp} -o {output.bam} {input} 2> {log}; "
        "samtools index -@ {threads} {output.bam} 2>> {log}"


# module to deduplicate reads
# -----------------------------------------------------
rule umi_dedup:
    input:
        bam=rules.mapping_sorted_bam.output.bam,
        bai=rules.mapping_sorted_bam.output.bai,
    output:
        bam="results/deduplicated/{sample}.bam",
        bai="results/deduplicated/{sample}.bam.bai",
    conda:
        "../envs/umitools.yml"
    message:
        """--- UMI tools deduplication."""
    params:
        tmp="results/deduplicated/sort_{sample}_tmp",
        default=config["umi_dedup"],
    threads: max(1, int(workflow.cores * 0.25))
    log:
        path="results/deduplicated/log/{sample}.log",
        stderr="results/deduplicated/log/{sample}.stderr",
        stats="results/deduplicated/log/{sample}_umi_stats.txt",
    shell:
        "umi_tools dedup {params.default} --stdin={input.bam} --output-stats={log.stats} --log={log.path} 2> {log.stderr} | "
        "samtools sort -@ {threads} -O bam -T {params.tmp} -o {output.bam}; "
        "samtools index {output.bam}"


# module to extract selected biotypes from gff file
# -----------------------------------------------------
rule extract_features:
    input:
        gff=rules.get_genome.output.gff,
    output:
        gff="results/get_genome/selected_features.gff",
    log:
        path="results/get_genome/log/extract_features.log",
    message:
        """--- Removing selected biotype features from genome annotation."""
    params:
        features=config["extract_features"]["biotypes"],
    conda:
        "../envs/extract_features.yml"
    script:
        "../scripts/extract_features.py"


# module to filter bam file
# -----------------------------------------------------
rule filter_bam:
    input:
        bam=rules.umi_dedup.output.bam,
        bai=rules.umi_dedup.output.bai,
        gff=rules.extract_features.output.gff,
    output:
        bam="results/filtered_bam/{sample}.bam",
        bai="results/filtered_bam/{sample}.bam.bai",
        stats="results/filtered_bam/{sample}_stats.txt",
    log:
        path="results/filtered_bam/log/{sample}.log",
    message:
        """--- Removing reads mapping to selected biotype regions."""
    conda:
        "../envs/filter_bam.yml"
    params:
        defaults=config["bedtools_intersect"]["defaults"],
    threads: max(1, int(workflow.cores * 0.25))
    shell:
        "intersectBed -abam {input.bam} -b {input.gff} {params.defaults} | "
        "samtools sort -@ {threads} > {output.bam} 2> {log.path}; "
        "samtools index -@ {threads} {output.bam} 2>> {log.path}; "
        "samtools stats -@ {threads} {output.bam} > {output.stats}"


# module to extract mapping length of bam file
# -----------------------------------------------------
rule extract_mapping_length:
    input:
        bam=get_bam,
    output:
        tsv="results/{mapping_status}/length_dist/{sample}_length_dist.tsv",
        pdf="results/{mapping_status}/length_dist/{sample}_length_dist.pdf",
    message:
        """--- Extract mapping length of BAM input: {wildcards.mapping_status}/{wildcards.sample}.bam"""
    conda:
        "../envs/plot_mapping_length.yml"
    log:
        path="results/{mapping_status}/length_dist/log/{sample}_length_dist.log",
    script:
        "../scripts/plot_mapping_length.py"


# module to run multiQC on input + processed files
# -----------------------------------------------------
rule multiqc:
    input:
        expand(
            "results/fastqc_{status}/{sample}_fastqc.html",
            sample=samples.index,
            status=config["multiqc"]["fastqc_stage"],
        ),
        expand("results/clipped/{sample}.fastq.gz", sample=samples.index),
        expand(
            "results/umi_extraction/{sample}.fastq.gz",
            sample=samples.index,
        ),
        expand(
            "results/mapped/unsorted/{sample}.bam",
            sample=samples.index,
        ),
        expand(
            "results/{mapping_status}/{sample}.bam",
            mapping_status=["deduplicated", "filtered_bam"],
            sample=samples.index,
        ),
    output:
        report="results/multiqc/multiqc_report.html",
    conda:
        "../envs/multiqc.yml"
    message:
        """--- Generating MultiQC report for seq data."""
    params:
        config=config["multiqc"]["config"],
    log:
        path="results/multiqc/log/multiqc.log",
    shell:
        "outdir=`echo {output.report} | cut -f 1-2 -d /`; "
        "multiqc -c {params.config} --force --verbose --dirs --outdir ${{outdir}} results &> {log.path}"
