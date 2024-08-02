# LOAD PACKAGES
# ------------------------------

suppressWarnings({
  suppressPackageStartupMessages({
    library(tidyverse)
    library(Biostrings)
    library(GenomicFeatures)
    library(ORFik)
    library(rtracklayer)
  })
})


# CONFIGURATION
# ------------------------------
#
sm_params <- snakemake@params[[1]]
window_size <- sm_params$window_size
export_ofst <- sm_params$export_ofst
export_bigwig <- sm_params$export_bigwig
max_cores <- snakemake@threads
input_file <- snakemake@input$bam
input_shift <- sm_params$shift_table
output_bam_filt <- snakemake@output$bam
output_ofst_filt <- str_replace(output_bam_filt, ".bam$", ".ofst")
output_bigwig_filt <- str_remove(output_bam_filt, ".bam$")
output_png_preshift <- snakemake@output$png_preshift
output_png_postshift <- snakemake@output$png_postshift
output_shift <- snakemake@output$shift
read_length <- sm_params$read_length
end_alignment <- sm_params$end_alignment
skip_shifting <- sm_params$skip_shifting
skip_length_filter <- sm_params$skip_length_filter


# STAGE 1 : IMPORT SEQ DATA
# ------------------------------
# import *.bam files
messages <- c()
messages <- append(messages, paste0("importing sequencing data: ", input_file))

# import orf annotation
messages <- append(messages, "import ORF annotation in GRanges format")
load(snakemake@input$granges)

# read bam file
messages <- append(messages, "read unshifted BAM files")
bam <- ORFik::readBam(input_file)

# filter bam file by read lengths
if (!skip_length_filter) {
  messages <- append(
    messages,
    paste0("filtering reads by length (", read_length[1], " to ", read_length[2], " nt)")
  )
  bam_filtered <- bam[readWidths(bam) %in% seq(read_length[1], read_length[2])]
} else {
  bam_filtered <- bam
}

# collapse reads to single base position
bam_onebase <- convertToOneBasedRanges(bam_filtered, method = end_alignment, addSizeColumn = TRUE)


# STAGE 2 : SHIFT READS TO P-SITE
# ------------------------------
# function to determine read density in window around start codon
get_density <- function(footprints, cds, transcript, window_size) {
  windowPerReadLength(
    grl = cds,
    tx = transcript,
    reads = footprints,
    pShifted = FALSE,
    drop.zero.dt = TRUE,
    scoring = "zscore",
    upstream = window_size,
    downstream = window_size
  )
}

# function to shift reads aligned to the 3' end instead of 5' (default)
shiftFootprints_3prime <- function(footprints, shifts, sort = TRUE) {
  if (!is(shifts, "data.frame"))
    stop("shifts must be data.frame/data.table")
  if (nrow(shifts) == 0)
    stop("No shifts found in data.frame")
  if (is.null(shifts$fraction) | is.null(shifts$offsets_start))
    stop("Either fraction (read lengths) or offsets_start (shifts by nt) column in shifts is not set!")
  selected_lengths <- shifts$fraction
  selected_shifts <- shifts$offsets_start
  lengths_all <- readWidths(footprints, along.reference = TRUE)
  shifts_all <- selected_shifts[match(lengths_all, selected_lengths)]
  direction <- ifelse(strand(footprints) == "+", -1, 1)
  shifts_all <- shifts_all * direction
  shifted <- shift(footprints, shifts_all)
  shifted <- sortSeqlevels(shifted)
  if (sort) {
    messages <- append(messages, "Sorting shifted footprints...")
    shifted <- sort(shifted)
  }
  return(shifted)
}

# function to fall back to generic shift table when everything fails
fall_back <- function(df) {
  shift_fail <- c(
    "Automatic determination of shift offsets failed. ",
    "Using a generic shift table as fall back. ",
    "Shift table is exported to ./results/shift_reads/<sample>_shift.csv"
  )
  warning(paste0(shift_fail, collapse = "\n"))
  df$shift_model <- NA
  if (end_alignment == "5prime") {
    df$offsets_start <- (df$fraction - 16) * -1
  } else {
    df$offsets_start <- 16
  }
  return(df)
}

# function to determine read shift offsets from density in window
get_offset <- function(footprints, cds, transcript, window_size) {
  df_density <- get_density(footprints, cds, transcript, window_size) %>%
    as_tibble() %>%
    arrange(fraction, position)
  df_shift <- df_density %>%
    complete(position, fraction) %>%
    group_by(fraction) %>%
    mutate(score = replace(score, is.na(score), 0)) %>%
    mutate(
      reference = filter(
        ., fraction == round(sum(read_length) / 2)
      )$score[1:n()],
      offset = filter(
        ., fraction == round(sum(read_length) / 2),
        score == max(score)
      )$position
    )
  if (sum(df_shift$score == 0) / nrow(df_shift) > 0.5) {
    df_shift <- df_shift %>%
      summarize(
        .groups = "drop",
        offset = NA,
        shift = NA
      )
  } else {
    df_shift <- df_shift %>%
      summarize(
        .groups = "drop",
        offset = offset[1],
        shift = position[which.max(
          ccf(reference, score, lag.max = window_size, plot = FALSE)$acf
        )]
      )
    if (end_alignment == "5prime") {
      df_shift <- df_shift %>%
        mutate(shift = ifelse(c(0, diff(shift) == 1), shift, NA))
    }
  }
  if (sum(!is.na(df_shift$shift)) < 5) {
    df_result <- fall_back(df_shift)
  } else {
    df_result <- df_shift %>%
      mutate(
        shift_model = round(predict(lm(shift ~ fraction, .), .)),
        offsets_start = offset - shift_model
      )
    if (end_alignment == "5prime" & any(df_result$offsets_start > 0)) {
      df_result <- fall_back(df_shift)
    }
  }
  df_result
}

# function to calculate and plot heatmap of start site coverage
plot_hitmap <- function(
    footprints, cds, transcript, title = NULL,
    window_size = 30, plot = TRUE) {
  if (any(width(footprints[1:100]) > 1)) {
    footprints <- convertToOneBasedRanges(
      footprints,
      method = end_alignment,
      addSizeColumn = TRUE
    )
  }
  # check TSS window size and remove too short/long windows
  windows = startRegion(cds, transcript, TRUE, window_size, window_size)
  windowSize <- (2 * window_size) + 1
  matching_window <- which(widthPerGroup(windows) == windowSize)

  hitmap <- windowPerReadLength(
    grl = cds[matching_window],
    tx = transcript[matching_window],
    reads = footprints,
    pShifted = FALSE,
    drop.zero.dt = TRUE,
    upstream = window_size,
    downstream = window_size
  )
  if (plot) {
    coverageHeatMap(
      hitmap,
      title = title,
      colors = c("white", "#B3B3B3", "#E6AB02", "#E7298A", "#7570B3"),
      scoring = "zscore", addFracPlot = TRUE
    )
  } else {
    hitmap
  }
}

# shift bam files
if (skip_shifting) {
  messages <- append(messages, "shifting not performed, reads are only end-aligned")
  bam_filtered_shifted <- bam_onebase
} else {
  # read optional shift table
  if (!is.null(input_shift)) {
    messages <- append(messages, paste0("importing shift table: ", input_shift))
    df_shift <- read_csv(input_shift, show_col_types = FALSE)
  } else {
    messages <- append(messages, "calculating shift table")
    df_shift <- get_offset(bam_onebase, list_cds, list_tx, window_size)
  }
  if (end_alignment == "5prime") {
    messages <- append(messages, "shifting footprints and reducing to 5' end")
    bam_filtered_shifted <- shiftFootprints(bam_filtered, df_shift)
  } else if (end_alignment == "3prime") {
    messages <- append(messages, "shifting footprints and reducing to 3' end")
    bam_filtered_shifted <- shiftFootprints_3prime(bam_onebase, df_shift)
  }
}


# STAGE 3 : EXPORT RESULTS
# ------------------------
messages <- append(messages, "exporting all read output files")

# export shifted/filtered bam file in conventional bam format (optional)
rtracklayer::export(object = bam_filtered_shifted, con = output_bam_filt, format = "bam")

# export shifted/filtered bam file in space-efficient osft format
if (export_ofst) {
  export.ofst(bam_filtered_shifted, file = output_ofst_filt)
}

# export shifted/filtered bam file in bigwig track format (optional)
if (export_bigwig) {
  export.bigWig(bam_filtered_shifted, file = output_bigwig_filt)
}

# export shift table
messages <- append(messages, "exporting shift table")
write_csv(df_shift, file = output_shift)

# plot heatmap of start site coverage before and after shifting
messages <- append(messages, "exporting heatmaps of start codon coverage before and after shifting")
png(filename = output_png_preshift, width = 875, height = 563, res = 120)
plot_hitmap(
  bam_onebase,
  list_cds, list_tx,
  title = paste0("File: ", input_file, " -- before shifting")
)
invisible(capture.output(dev.off()))
png(filename = output_png_postshift, width = 875, height = 563, res = 120)
plot_hitmap(
  bam_filtered_shifted,
  list_cds, list_tx,
  title = paste0("File: ", input_file, " -- after shifting")
)
invisible(capture.output(dev.off()))

# export log
messages <- append(messages, "export of shifted read files complete")
write_lines(
  file = snakemake@log[["path"]],
  x = paste0("SHIFT_READS: ", messages)
)
