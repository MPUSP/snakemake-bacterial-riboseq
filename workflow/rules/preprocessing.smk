# ----------------------------------------------------- #
# Riboseq preprocessing:                                #
# ----------------------------------------------------- #

# module to make QC report
# -----------------------------------------------------
rule fastqc:
    input:
        get_fastq,
    output:
        report=directory("results/fastqc_{status}/{sample}"),
    conda:
        "../envs/fastqc.yml",
    message:
        """--- Checking fastq files with FastQC."""
    log:
        "results/fastqc_{status}/log/{sample}.log",
    threads: int(workflow.cores * 0.2), # assign 20% of max cores
    shell:
        "mkdir -p {output.report};"
        "fastqc --nogroup -extract --threads {threads} -o {output.report} {input} > {log}"


# module to trim adapters from reads
# -----------------------------------------------------
rule cutadapt:
    input:
        fastq=lambda wc: expand("{input_dir}/{sample}",
            input_dir=samples.loc[wc.sample]["data_folder"],
            sample=samples.loc[wc.sample]["fq1"]),
    output:
        fastq="results/clipped/{sample}.fastq.gz",
    conda:
        "../envs/cutadapt.yml"
    message:
        """--- Trim adapters from reads."""
    params:
        fivep_adapter=config["cutadapt"]["fivep_adapter"],
        threep_adapter=config["cutadapt"]["threep_adapter"],
        default=config["cutadapt"]["default"],
    log:
        stdout="results/clipped/log/{sample}.log",
        stderr="results/clipped/log/{sample}.stderr",
    threads: int(workflow.cores * 0.4), # assign 40% of max cores
    shell:
        "if [ {params.fivep_adapter} != None ]; then "
        "fivep=`echo -g {params.fivep_adapter}`; "
        "else fivep=''; "
        "fi; "
        "if [ {params.threep_adapter} != None ]; then "
        "threep=`echo -a {params.threep_adapter}`; "
        "else threep=''; "
        "fi; "
        "cutadapt ${{fivep}} ${{threep}} "
        "{params.default} --cores {threads} "
        "-o {output.fastq} {input.fastq} > {log.stdout} 2> {log.stderr}"


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
    output:
        path=directory("results/get_genome"),
        fasta="results/get_genome/genome.fasta",
        gff="results/get_genome/genome.gff",
    conda:
        "../envs/get_genome.yml"
    message:
        """--- Parsing genome GFF and FASTA files."""
    params:
        database=config["get_genome"]["database"],
        assembly=config["get_genome"]["assembly"],
        fasta=config["get_genome"]["fasta"],
        gff=config["get_genome"]["gff"],
    log:
        path="results/get_genome/log/get_genome.log",
    script:
        "../scripts/get_genome.py"


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
    threads: int(workflow.cores * 0.2), # assign 20% of max cores
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
        "results/mapped/log/samtools_{sample}.log"
        # os.path.join(output_dir, "mapping/log/samtools_sort_{sample}.stderr")
    message: """--- Samtools sort and index bam files."""
    params:
        tmp="results/mapped/sort_{sample}_tmp",
    threads: int(workflow.cores * 0.2), # assign 20% of max cores
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
    threads: int(workflow.cores * 0.2), # assign 20% of max cores
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
    threads: int(workflow.cores * 0.2), # assign 20% of max cores
    shell:
        "intersectBed -abam {input.bam} -b {input.gff} {params.defaults} | "
        "samtools sort -@ {threads} > {output.bam} 2> {log.path}; "
        "samtools index -@ {threads} {output.bam} 2>> {log.path}; "
        "samtools stats -@ {threads} {output.bam} > {output.stats}"


# module to extract mapping length of bam file
# -----------------------------------------------------
rule extract_mapping_length:
    input:
        get_bam,
    output:
        "results/{mapping_status}/length_dist/{sample}_length_dist.tsv",
    message: """--- Extract mapping length of BAM input: {wildcards.mapping_status}/{wildcards.sample}.bam"""
    conda:
        "../envs/plot_mapping_length.yml",
    log:
        os.path.join("results", "{mapping_status}", "length_dist", "log", "{sample}.log"),
    params:
        outdir=os.path.join("results","{mapping_status}", "length_dist"),
    shell:
        "workflow/scripts/plot_mapping_length.py -b {input} -o {params.outdir} -s {wildcards.sample} 2> {log}"


# module to run multiQC on input + processed files
# -----------------------------------------------------
rule multiqc:
    input:
        expand(
            "results/fastqc_clipped/{sample}_fastqc.html",
            sample=samples.index),
        expand(
            "results/clipped/{sample}.fastq.gz",
            sample=samples.index),
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
