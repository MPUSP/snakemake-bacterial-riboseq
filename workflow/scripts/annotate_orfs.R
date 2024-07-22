# LOAD PACKAGES
# ------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(Biostrings)
  library(GenomicFeatures)
  library(ORFik)
})


# CONFIGURATION
# ------------------------------
sm_params <- snakemake@params[[1]]
max_cores <- snakemake@params[[2]]
for (param in names(sm_params)) {
  assign(param, sm_params[[param]], envir = .GlobalEnv)
}


# STAGE 1 : FILE PREPARATION
# ------------------------------
# import genome sequence
messages <- c("importing genome sequence and annotation")
genome_fasta <- snakemake@input[["fasta"]]
genome_gff <- snakemake@input[["gff"]]
genome_dna <- readDNAStringSet(genome_fasta)
indexFa(genome_fasta)
organism <- str_remove_all(names(genome_dna)[1], "^[a-zA-Z0-9_]+.[0-9]+ |, .*")

# import genome annotation with chromosome metadata
txdb <- makeTxdbFromGenome(
  genome_gff,
  genome_fasta,
  organism = organism,
  optimize = TRUE,
  pseudo_5UTRS_if_needed = window_size,
  return = TRUE
)

# harmonize metadata
seqinfo(genome_dna) <- seqinfo(txdb)

# extract all transcripts and cds
list_cds <- loadRegion(txdb, "cds", by = "tx")
list_tx <- loadRegion(txdb, "mrna", by = "tx")

# make leader and start codon regions
list_leader <- startCodons(list_cds) %>%
  extendLeaders(extension = window_size)
list_start <- extendTrailers(list_leader, window_size)

# extract sequences for all CDSs
list_cds_seq <- extractTranscriptSeqs(genome_dna, list_cds)

# STAGE 1 : FIND ALL POSSIBLE ORFS
# ------------------------------

# list of all possible ORFs
list_orfs <- ORFik::findORFsFasta(
  genome_dna,
  startCodon = startDefinition(orf_start_codon_table),
  stopCodon = paste(orf_stop_codon, collapse = "|"),
  minimumLength = floor(sorf_min_length / 3),
  longestORF = orf_longest_only,
  is.circular = FALSE
)

# sort and rename predicted ORFs
list_orfs <- sort(list_orfs)
names(list_orfs) <- paste0(
  "ORF",
  formatC(seq_along(list_orfs),
    flag = "0",
    width = nchar(length(list_orfs))
  )
)

# coerce ORF GRanges to data frame
df_annotated_orfs <- list_cds %>%
  as_tibble() %>%
  dplyr::filter(exon_rank == 1) %>%
  mutate(
    sequence = as.character(list_cds_seq),
    start_codon = str_sub(sequence, 1, 3),
    stop_codon = str_sub(sequence, -3, -1),
    intergenic = FALSE,
    intragenic = FALSE,
    partial_overlap = FALSE
  )

# export results
write_csv(df_annotated_orfs, file = snakemake@output[["df_annotated_orfs"]])
save(
  list = c("list_cds", "list_tx", "list_orfs", "list_leader", "list_start"),
  file = snakemake@output[["granges"]]
)

# export log
write_lines(
  file = snakemake@log[["path"]],
  x = paste0("ANNOTATE_ORFS: ", messages)
)
