# LOAD PACKAGES
# ------------------------------

suppressWarnings({
  suppressPackageStartupMessages({
    library(tidyverse)
    library(Biostrings)
    library(GenomicFeatures)
    library(GenomeInfoDb)
    library(ORFik)
  })
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
genome_fa_index <- indexFa(genome_fasta)

# try to look up organism and taxonomy ID
df_taxonomy <- loadTaxonomyDb()
header_gff <- read_lines(genome_gff, n_max = 10)
header_species <- na.omit(str_match(header_gff, ";species=.*"))[1]

if (!is.na(header_species)) {
  messages <- c(messages, "trying to guess taxonomy ID from GFF file")
  taxon_id <- header_species %>%
    str_extract("wwwtax.cgi(%3F|\\?)id(%3D|=)[0-9]+") %>%
    str_extract("[0-9]+$") %>%
    as.numeric()
  messages <- c(messages, paste("extracted taxon ID:", taxon_id))
  tax_entry <- filter(df_taxonomy, tax_id == taxon_id)
  if (!is.na(tax_entry[1, "genus"])) {
    genome_name <- paste(tax_entry[1, "genus"], tax_entry[1, "species"])
    messages <- c(messages, paste("found genome name:", genome_name)) 
  } else {
    genome_name <- NA
    messages <- c(messages, paste0("the extracted Taxonomy ID '", taxon_id, "' is not valid"))
  }
} else {
  genome_name <- NA
}

if (is.na(genome_name)) {
  messages <- c(messages, "trying to guess genome name from fasta file")
  genome_name <- str_remove_all(
    names(genome_dna)[1],
    paste0(genome_seqlevels[1], " |chromosome|\\,.*")
  )
  messages <- c(messages, paste("extracted genome name:", genome_name))
  df_taxonomy$name <- paste(df_taxonomy$genus, df_taxonomy$species)
  tax_entry <- filter(df_taxonomy, name == genome_name)
  taxon_id <- as.numeric(tax_entry[1, "tax_id"])
  if (!is.na(taxon_id)) {
    messages <- c(messages, paste("found taxon ID:", taxon_id))
  } else {
    genome_name <- NA
  }
}

if (is.na(genome_name)) {
  messages <- c(messages, "taxonomy guessing failed: falling back to arbitrary taxon ID")
  genome_name <- "root"
  taxon_id <- 1
}

# import genome annotation with chromosome metadata
quiet_txdb <- quietly(makeTxdbFromGenome)
txdb <- quiet_txdb(
  genome_gff,
  genome_fasta,
  organism = genome_name,
  optimize = TRUE,
  pseudo_5UTRS_if_needed = window_size,
  return = TRUE
)
messages <- append(messages, paste0("warning: ", txdb$warnings))
messages <- append(messages, txdb$messages)

# harmonize metadata
seqinfo(genome_dna) <- seqinfo(txdb$result)

# extract all transcripts and cds
list_cds <- loadRegion(txdb$result, "cds", by = "tx")
list_tx <- loadRegion(txdb$result, "mrna", by = "tx")

# parse genome gff file GFF
df_gff <- genome_gff %>%
  read_tsv(
    comment = "#",
    col_names = c(
      "seqnames", "source", "feature", "start", "end",
      "score", "strand", "frame", "attribute"
    ),
    show_col_types = FALSE
  ) %>%
  mutate(
    width = end - start,
    name = str_extract(attribute, ";Name=[a-zA-Z0-9_-]+;") %>%
      str_remove_all("Name=|;"),
    id = str_extract(attribute, "ID=gene-[a-zA-Z0-9_-]+;") %>%
      str_remove_all("ID=gene-|;")
  ) %>%
  filter(name %in% names(list_cds))

# make leader and start codon regions
list_leader <- startCodons(list_cds) %>%
  extendLeaders(extension = window_size)
list_start <- extendTrailers(list_leader, window_size)

# extract sequences for all CDSs
list_cds_seq <- extractTranscriptSeqs(genome_dna, list_cds)

# coerce ORF GRanges to data frame
df_annotated_orfs <- list_cds %>%
  as_tibble() %>%
  dplyr::filter(exon_rank == 1) %>%
  mutate(
    sequence = as.character(list_cds_seq),
    start_codon = str_sub(sequence, 1, 3),
    stop_codon = str_sub(sequence, -3, -1)
  )

# export results
write_csv(df_annotated_orfs, file = snakemake@output[["df_annotated_orfs"]])
save(
  list = c("list_cds", "list_tx", "list_leader", "list_start"),
  file = snakemake@output[["granges"]]
)

# export log
write_lines(
  file = snakemake@log[["path"]],
  x = paste0("ANNOTATE_ORFS: ", messages)
)
