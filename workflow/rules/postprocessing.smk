# ----------------------------------------------------- #
# Riboseq postprocessing:                               #
# ----------------------------------------------------- #


# module to extract selected biotypes from gff file
# -----------------------------------------------------
rule extract_mRNA_features:
    input:
        gff=rules.get_genome.output.gff,
    output:
        gff="results/get_genome/mRNA_features.gff",
    log:
        path="results/get_genome/log/mRNA_features.log",
    message:
        """--- Removing selected biotype features from genome annotation."""
    params:
        features=config["extract_features"]["CDS"],
    conda:
        "../envs/extract_features.yml"
    script:
        "../scripts/extract_features.py"


# module to process annotated ORFs
# -----------------------------------------------------
rule annotate_orfs:
    params:
        config["annotate_orfs"],
        max_cores=workflow.cores,
    input:
        fasta=rules.get_genome.output.fasta,
        gff=rules.get_genome.output.gff,
    output:
        df_annotated_orfs="results/annotate_orfs/df_annotated_orfs.csv",
        granges="results/annotate_orfs/granges.RData",
        txdb="results/get_genome/genome.gff.db",
    conda:
        "../envs/r_orfik.yml"
    message:
        """--- Preparing annotated ORFs."""
    log:
        path="results/annotate_orfs/log/log.txt",
    script:
        "../scripts/annotate_orfs.R"


# module to filter and shift reads
# -----------------------------------------------------
rule shift_reads:
    input:
        bam="results/filtered_bam/{sample}.bam",
        bai="results/filtered_bam/{sample}.bam.bai",
        granges=rules.annotate_orfs.output.granges,
    output:
        bam="results/shift_reads/{sample}.bam",
        png_preshift="results/shift_reads/{sample}_preshift.png",
        png_postshift="results/shift_reads/{sample}_postshift.png",
        shift="results/shift_reads/{sample}_shift.csv",
    conda:
        "../envs/r_orfik.yml"
    message:
        """--- Shifting Ribo-Seq reads."""
    params:
        config["shift_reads"],
    threads: workflow.cores
    log:
        path="results/shift_reads/log/{sample}.log",
    script:
        "../scripts/shift_reads.R"


# module to calculate feature-wise statistics
# -----------------------------------------------------
rule feature_stats:
    input:
        samples=config["samplesheet"],
        fasta=rules.get_genome.output.fasta,
        gff_rna=rules.extract_features.output.gff,
        gff_cds=rules.extract_mRNA_features.output.gff,
        granges=rules.annotate_orfs.output.granges,
        orfs=rules.annotate_orfs.output.df_annotated_orfs,
        txdb=rules.annotate_orfs.output.txdb,
        bam_filtered=expand("results/filtered_bam/{sample}.bam", sample=samples.index),
        bam_shifted=expand("results/shift_reads/{sample}.bam", sample=samples.index),
    output:
        csv="results/feature_stats/feature_stats.csv",
    conda:
        "../envs/r_orfik.yml"
    message:
        """--- Calculating feature-wise statistics."""
    threads: workflow.cores
    log:
        path="results/feature_stats/log/stats.log",
    script:
        "../scripts/feature_stats.R"


# module to report results as HTML notebook
# -----------------------------------------------------
rule report_html:
    input:
        orfs_annotated=rules.annotate_orfs.output.df_annotated_orfs,
        orfs_features=rules.feature_stats.output.csv,
        granges=rules.annotate_orfs.output.granges,
        fasta=rules.get_genome.output.fasta,
        bam_filtered=expand("results/filtered_bam/{sample}.bam", sample=samples.index),
        bam_shifted=expand("results/shift_reads/{sample}.bam", sample=samples.index),
    output:
        html="results/report/report.html",
    conda:
        "../envs/r_orfik.yml"
    message:
        """--- Writing HTML report with workflow results."""
    params:
        config["report"],
    log:
        path="results/report/log/report_html.log",
    script:
        "../notebooks/report.Rmd"


# module to convert HTML to PDF output
# -----------------------------------------------------
rule report_pdf:
    input:
        html=rules.report_html.output.html,
    output:
        pdf="results/report/report.pdf",
    conda:
        "../envs/report_pdf.yml"
    message:
        """--- Converting HTML report to PDF."""
    log:
        path="results/report/log/report_pdf.log",
    shell:
        "weasyprint -v {input.html} {output.pdf} &> {log.path}"
