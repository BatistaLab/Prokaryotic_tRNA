setwd("")
# working directory containing reference.fa, reference.fai,
# and subfolders for each condition containing replicate BAMs + indexes
library(data.table)
library(plyr)
library(tidyverse)
library(Rsamtools)

# taken from
# https://seqqc.wordpress.com/2015/03/10/calculate-nucelotide-frequency-with-rsamtools-pileup/
pileupFreq <- function(pileupres) {
  nucleotides <- levels(pileupres$nucleotide)
  res <- split(pileupres, pileupres$seqnames)
  res <- lapply(res, function(x) {
    split(x, x$pos)
  })
  res <- lapply(res, function(positionsplit) {
    nuctab <- lapply(positionsplit, function(each) {
      chr <- as.character(unique(each$seqnames))
      pos <- as.character(unique(each$pos))
      tablecounts <- sapply(
        nucleotides,
        function(n) {
          sum(each$count[each$nucleotide == n])
        }
      )
      c(chr, pos, tablecounts)
    })
    nuctab <-
      data.frame(do.call("rbind", nuctab), stringsAsFactors = FALSE)
    rownames(nuctab) <- NULL
    nuctab
  })
  res <- data.frame(do.call("rbind", res), stringsAsFactors = FALSE)
  # modified to not crash if no reads are found for the reference
  if (nrow(res) == 0) {
    return(0)
  } else {

  }
  rownames(res) <- NULL
  colnames(res) <-
    c("seqnames", "start", levels(pileupres$nucleotide))
  res[3:ncol(res)] <- apply(res[3:ncol(res)], 2, as.numeric)
  res
}

# list of tRNAs entered manually
# because reading from txt file gives weird results
list_of_experiments <- as.list(list.dirs(recursive = FALSE))
list_of_tRNAs <- list()
p_param <-
  PileupParam(
    max_depth = 1000000,
    distinguish_strands = TRUE,
    min_nucleotide_depth = 0
  )
tRNA_seqs <-
  open(
    FaFile(
      "reference.fa", # replace with specific reference
      "reference.fa.fai"
      # replace with index generated from specific reference using samtools-faidx
    )
  )
tRNA_lens <- seqlengths(tRNA_seqs)

# performs pileup for each tRNA
# for each experimental condition across replicates
# indexes experiment folders
# creates summary dataframe
# iterates through set list of tRNAs
# for every tRNA, gets a reference of positions and bases from index files
# indexes replicates
# performs pileup per experiment per tRNA per replicate
# assembles replicates horizontally, joined on tRNA index to reference
# final metrics are calculated per tRNA
# each tRNA is assembled vertically into the final format
for (experiment in list_of_experiments) {
  all_tRNAs_list <- list()
  pileup_output <- data.frame(matrix(
    ncol = 34,
    nrow = 0,
    dimnames = list(
      NULL,
      c(
        "tRNA_index",
        "A_rep1",
        "C_rep1",
        "G_rep1",
        "T_rep1",
        "N_rep1",
        "depth_rep1",
        "A_freq_rep1",
        "C_freq_rep1",
        "G_freq_rep1",
        "T_freq_rep1",
        "N_freq_rep1",
        "A_rep2",
        "C_rep2",
        "G_rep2",
        "T_rep2",
        "N_rep2",
        "depth_rep2",
        "A_freq_rep2",
        "C_freq_rep2",
        "G_freq_rep2",
        "T_freq_rep2",
        "N_freq_rep2",
        "A_rep3",
        "C_rep3",
        "G_rep3",
        "T_rep3",
        "N_rep3",
        "depth_rep3",
        "A_freq_rep3",
        "C_freq_rep3",
        "G_freq_rep3",
        "T_freq_rep3",
        "N_freq_rep3"
      )
    )
  ))

  for (tRNA in list_of_tRNAs) {
    tRNA_length <- tRNA_lens[[tRNA]]
    tRNA_reference <-
      getSeq(
        tRNA_seqs, GRanges(
          seqnames = tRNA, IRanges(start = 1, end = tRNA_length)
        )
      )
    reference_dataframe <-
      data.frame(ref_base = unlist(strsplit(
        as.character(tRNA_reference, use.names = FALSE),
        split = ""
      )[[1]]))
    reference_dataframe$tRNA_index <- 1:nrow(reference_dataframe)
    reference_dataframe$tRNA_index <-
      lapply(reference_dataframe$tRNA_index, function(x) {
        paste(names(tRNA_reference), x, sep = "_")
      })
    reference_dataframe <- reference_dataframe[, c("tRNA_index", "ref_base")]
    reference_dataframe <-
      reference_dataframe %>% mutate(
        tRNA_index = as.character(tRNA_index)
      )
    exp_name <- substring(experiment, 3)
    replicate_files <- as.list(list.files(
      path = experiment,
      no.. = TRUE,
      pattern = "(?:.bam$)"
    ))
    param <- ScanBamParam(which = GRanges(
      tRNA,
      IRanges(
        start = 0,
        end = 1000
      )
    ))
    output_list <- list()
    output_list <- append(list(reference_dataframe), output_list)
    for (replicate in replicate_files) {
      rep_file <- paste(experiment, replicate, sep = "/")
      rep_name <- str_extract(rep_file, pattern = "(?:rep.)")
      bf <- BamFile(rep_file)
      rep_pileup <- pileup(
        bf,
        include_deletions = 1,
        include_insertions = 1,
        ignore_query_Ns = FALSE,
        distinguish_nucleotides = TRUE,
        scanBamParam = param,
        pileupParam = p_param
      )
      freqs <- pileupFreq(rep_pileup)
      if (is.data.frame(freqs) == FALSE) {
        next
      }
      freqs <- select(freqs, -c("-", "=", "+"))
      freqs$depth <- rowSums(freqs[, c("A", "C", "G", "T", "N")])
      freqs$A_freq <- freqs$A / freqs$depth
      freqs$C_freq <- freqs$C / freqs$depth
      freqs$G_freq <- freqs$G / freqs$depth
      freqs$T_freq <- freqs$T / freqs$depth
      freqs$N_freq <- freqs$N / freqs$depth
      colnames(freqs) <- paste(colnames(freqs), rep_name, sep = "_")
      ref <- paste("seqnames", rep_name, sep = "_")
      pos <- paste("start", rep_name, sep = "_")
      freqs <- unite(freqs, "tRNA_index", c(ref, pos))
      output_list <- append(list(freqs), output_list)
    }
    # https://stackoverflow.com/a/49346473
    no_of_reps_with_reads <- length(output_list) - 1
    output_df <-
      output_list %>% purrr::reduce(dplyr::full_join, by = "tRNA_index")
    output_df <- output_df %>% replace(is.na(.), 0)
    output_df <- output_df %>%
      mutate(
        A_avg_freq = rowSums(
          select(
            output_df, starts_with("A_freq")
          )
        ) / no_of_reps_with_reads
      )
    output_df <- output_df %>%
      mutate(
        C_avg_freq = rowSums(
          select(output_df, starts_with("C_freq"))
        ) / no_of_reps_with_reads
      )
    output_df <- output_df %>%
      mutate(
        G_avg_freq = rowSums(
          select(output_df, starts_with("G_freq"))
        ) / no_of_reps_with_reads
      )
    output_df <- output_df %>%
      mutate(
        T_avg_freq = rowSums(
          select(output_df, starts_with("T_freq"))
        ) / no_of_reps_with_reads
      )
    output_df <- output_df %>%
      mutate(
        N_avg_freq = rowSums(
          select(output_df, starts_with("N_freq"))
        ) / no_of_reps_with_reads
      )
    output_df <- output_df %>%
      mutate(
        misincorporation_freq = case_when(
          ref_base == "A" ~ 1 - A_avg_freq,
          ref_base == "C" ~ 1 - C_avg_freq,
          ref_base == "T" ~ 1 - T_avg_freq,
          ref_base == "G" ~ 1 - G_avg_freq,
          TRUE ~ 0
        )
      )
    output_df <- output_df %>%
      mutate(
        avg_depth = rowSums(
          select(output_df, starts_with("depth_"))
        ) / no_of_reps_with_reads
      )
    output_df <- output_df %>%
      mutate(
        total_depth = rowSums(
          select(output_df, starts_with("depth"))
        )
      )
    output_df <- output_df %>%
      mutate(
        misincorporation_freq = case_when(
          total_depth == 0 ~ 0,
          TRUE ~ misincorporation_freq
        )
      )
    pileup_output <- rbind.fill(pileup_output, output_df)
    print(paste(tRNA, "finished"))
  }
  print(paste(exp_name, "finished"))
  raw_output_name <-
    paste(exp_name, "collated_freqs.csv", sep = "_")
  summary_output_name <-
    paste(exp_name, "summary_freqs.csv", sep = "_")
  filtered_output_name <-
    paste(exp_name, "filtered_freqs.csv", sep = "_")
  pileup_output <- pileup_output %>% replace(is.na(.), 0)
  raw_pileup_output <-
    pileup_output[, c(1, 41, 35, 29, 18, 7, 42, 43, 2:6, 8:17, 19:28, 30:34, 36:40)]
  summary_pileup_output <-
    pileup_output[, c(
      "tRNA_index",
      "ref_base",
      "misincorporation_freq",
      "avg_depth",
      "A_avg_freq",
      "C_avg_freq",
      "G_avg_freq",
      "T_avg_freq",
      "N_avg_freq",
      "total_depth"
    )]
  summary_pileup_output <- summary_pileup_output %>%
    rowwise() %>%
    mutate(
      tRNA_name = sub("_.*", "", tRNA_index)
    )
  summary_pileup_output <- summary_pileup_output %>%
    rowwise() %>%
    mutate(
      tRNA_position = as.numeric(sub(".*_", "", tRNA_index))
    )
  colnames(raw_pileup_output) <-
    paste(colnames(raw_pileup_output),
      as.character(experiment),
      sep = "_"
    )
  filtered_pileup_output <- summary_pileup_output %>%
    rowwise() %>%
    mutate(
      misincorporation_freq = case_when(
        avg_depth < 100 ~ NA,
        TRUE ~ misincorporation_freq
      )
    )
  write.csv(raw_pileup_output, raw_output_name)
  write.csv(summary_pileup_output, summary_output_name)
  write.csv(filtered_pileup_output, filtered_output_name)
}