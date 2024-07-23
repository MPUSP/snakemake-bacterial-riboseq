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


# module to process annotated and predicted ORFs
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
    conda:
        "../envs/r_orfik.yml"
    message:
        """--- Listing annotated and potential new ORFs."""
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
        ofst_orig="results/shift_reads/{sample}_orig.ofst",
        ofst_filt="results/shift_reads/{sample}_filt.ofst",
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
