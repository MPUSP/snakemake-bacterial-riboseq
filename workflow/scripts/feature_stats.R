# LOAD PACKAGES
# ------------------------------

suppressWarnings({
  suppressPackageStartupMessages({
    library(tidyverse)
    library(Biostrings)
    library(GenomicFeatures)
    library(GenomeInfoDb)
    library(ORFik)
    library(parallel)
    library(doParallel)
    library(AnnotationDbi)
  })
})


# IMPORT FILES
# ------------------------------

# import genome sequence
messages <- c("importing genome sequence and annotation")
samplesheet <- snakemake@input[["samples"]]
genome_fasta <- snakemake@input[["fasta"]]
genome_gff_rna <- snakemake@input[["gff_rna"]]
genome_gff_cds <- snakemake@input[["gff_cds"]]
genome_dna <- readDNAStringSet(genome_fasta)
genome_fa_index <- indexFa(genome_fasta)

# import annotated ORF table + granges + txdb
messages <- append(messages, "importing annotated ORFs, table and Granges")
df_annotated <- read_csv(snakemake@input[["orfs"]], show_col_types = FALSE)
load(snakemake@input[["granges"]])
txdb <- loadDb(snakemake@input[["txdb"]])

# import original and processed read data
messages <- append(messages, "reading unshifted BAM files")
bam <- lapply(snakemake@input[["bam_filtered"]], ORFik::readBam)
names(bam) <- str_split_i(snakemake@input[["bam_filtered"]], "\\/", -1)

messages <- append(messages, "reading shifted and end-aligned BAM files")
bam_shift <- lapply(snakemake@input[["bam_shifted"]], ORFik::readBam)
names(bam_shift) <- str_split_i(snakemake@input[["bam_shifted"]], "\\/", -1)


# STATISTICAL ANALYSIS
# ------------------------------

# harmonize metadata
seqinfo(genome_dna) <- seqinfo(txdb)

# generate new instance of DB file on the fly
output_temp_db <- str_remove(snakemake@output[["csv"]], "[a-zA-Z_]*.csv$")
messages <- append(messages, paste0(
  "temporary TxDb files stored in ", output_temp_db
))

make_temp_txdb <- function(txdb, file) {
  db_file <- paste0(
    output_temp_db, "txdb_",
    format(Sys.time(), "%Y%m%d_%X"), "_",
    formatC(file, flag = "0", width = 3)
  )
  saveDb(txdb, file = db_file)
  loadDb(db_file)
}

# compute features from riboseq data using customized ORFik function
compute_features <- function(grl, RFP, RNA = NULL, Gtf) {
  computeFeaturesQuietly <- purrr::quietly(ORFik::computeFeatures)
  cl <- makePSOCKcluster(min(c(length(RFP), snakemake@threads)))
  registerDoParallel(cl)
  result <- foreach(file = seq_along(RFP), .combine = bind_rows) %do% {
    computeFeaturesQuietly(
      grl = grl,
      RFP = RFP[[file]],
      RNA = RNA[[file]],
      Gtf = make_temp_txdb(Gtf, file),
      riboStart = snakemake@config$shift_reads$read_length[1],
      riboStop = snakemake@config$shift_reads$read_length[2],
      uorfFeatures = TRUE,
      sequenceFeatures = FALSE
    )$result %>%
      mutate(
        sample_name = names(RFP)[file],
        seqnames = unlist(grl) %>%
          {
            .[!duplicated(names(.))]
          } %>%
          seqnames() %>% as.character(),
        group_name = names(grl),
        cpmRFP = countRFP / sum(countRFP) * 10^6
      ) %>%
      dplyr::select(14, 15, 13, 1, 16, 2:12)
  }
  stopCluster(cl)
  return(result)
}

messages <- append(messages, "calculating gene-wise statistics")
df_stats <- compute_features(
  grl = list_cds,
  RFP = bam_shift,
  RNA = NULL,
  Gtf = txdb
) %>%
  as_tibble()

# add sample and replicate information to feature table
df_stats <- mutate(df_stats, sample_name = str_remove(sample_name, "\\.bam$"))
df_samples <- read_tsv(samplesheet, show_col_types = FALSE) %>%
  dplyr::rename(sample_name = sample) %>%
  dplyr::select(sample_name, condition, replicate)
df_stats <- left_join(
  df_stats, df_samples,
  by = "sample_name"
)
df_stats <- df_stats %>%
  relocate(seqnames, group_name, sample_name, condition, replicate)


# EXPORT RESULTS
# ------------------------------

# export main table
write_csv(df_stats, file = snakemake@output[["csv"]])

# remove temporary txdb files
messages <- append(messages, "removing temporary txdb files")
temp_files <- list.files(path = output_temp_db, pattern = "txdb_", full.names = TRUE)
output_dump <- file.remove(temp_files)

# export log
write_lines(
  file = snakemake@log[["path"]],
  x = paste0("FEATURE STATS: ", messages)
)
