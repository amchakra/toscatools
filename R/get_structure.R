# ==========
# Functions to analyse structures
# ==========

#' Get sequence
#'
#' Get sequences for left and right hybrid arms
#'
#' @param hybrids.dt Hybrids data.table
#' @param genome.dt Transcriptome FASTA converted to data.table
#' @return hybrids.dt with \code{L_sequence} and \code{R_sequence} columns
#' @import data.table
#' @export

get_sequence <- function(hybrids.dt, genome.dt) {

  # Get L
  hybrids.dt[, gene_id := L_seqnames]
  setkey(genome.dt, gene_id)
  setkey(hybrids.dt, gene_id)
  seq.dt <- genome.dt[hybrids.dt]
  seq.dt[, `:=`(L_sequence, stringr::str_sub(sequence, start = L_start, end = L_end))]
  seq.dt[, `:=`(c("gene_id", "sequence"), NULL)]

  # Get R
  # Do separately in case on different genes
  seq.dt[, gene_id := R_seqnames]
  setkey(seq.dt, gene_id)
  seq.dt <- genome.dt[seq.dt]
  seq.dt[, `:=`(R_sequence, stringr::str_sub(sequence, start = R_start, end = R_end))]
  seq.dt[, `:=`(c("gene_id", "sequence"), NULL)]

  return(seq.dt)
}

#' Shuffle sequence
#'
#' Shuffles sequences using \code{uShuffle}
#'
#' @param sequence Sequence as character vector
#' @param number Number of shuffled sequences to return
#' @param klet Preserve klet nucleotide frequencies
#' @param seed Seed
#' @return Shuffled sequences as character vector
#' @export

shuffle_sequence <- function(sequence, number = 1, klet = 2, seed = 42) {
  system(paste0("uShuffle -seed ", seed, " -k ", klet, " -n ", number, " -s ", sequence), intern = TRUE)
}

#' Analyses hybrid structure
#'
#' Analyses hybrid structure using \code{RNAduplex} from \code{ViennaRNA}
#'
#' @param name Hybrid name
#' @param L_sequence Left sequence
#' @param R_sequence Right sequence
#' @return data.table with binding energy and dot-bracket structure
#' @import data.table
#' @export

analyse_structure <- function(name, L_sequence, R_sequence) {
  input <- paste0(L_sequence, "\n", R_sequence)
  rnaduplex <- system("RNAduplex --noLP", input = input, intern = TRUE)

  # Get MFE
  rnaduplex <- gsub("\\s+", "_", rnaduplex)
  mfe <- sapply(strsplit(rnaduplex, "_"), "[", 5)
  if (mfe == "(") mfe <- sapply(strsplit(rnaduplex, "_"), "[", 6) # Positives < 10 have an extra space
  mfe <- as.numeric(gsub("\\(|\\)", "", mfe))

  # Get dot.bracket
  db <- sapply(strsplit(rnaduplex, "_"), "[", 1)

  # Arms widths
  l_width <- nchar(L_sequence)
  r_width <- nchar(R_sequence)

  # Get from/to positions
  l <- sapply(strsplit(rnaduplex, "_"), "[", 2)
  r <- sapply(strsplit(rnaduplex, "_"), "[", 4)
  l_from <- as.numeric(sapply(strsplit(l, ","), "[", 1))
  l_to <- as.numeric(sapply(strsplit(l, ","), "[", 2))
  r_from <- as.numeric(sapply(strsplit(r, ","), "[", 1))
  r_to <- as.numeric(sapply(strsplit(r, ","), "[", 2))

  # If no structure found e.g.
  # L_sequence <- "ACACAACACACACACACACACA"
  # R_sequence <- "ACACACACACACACACACACAAAA"
  # then l = "1,1" and r = "0,0"

  if (r == "0,0") {
    db <- paste0(c(rep(".", l_width)), "&", c(rep(".", r_width)))
  } else {

    # Pad db with . according to from/to
    l_db <- sapply(strsplit(db, "&"), "[", 1)
    r_db <- sapply(strsplit(db, "&"), "[", 2)
    l_db <- unlist(strsplit(l_db, "")) # convert string to vector of characters
    r_db <- unlist(strsplit(r_db, ""))
    l_db <- c(rep(".", l_from - 1L), l_db, rep(".", l_width - l_to))
    r_db <- c(rep(".", r_from - 1L), r_db, rep(".", r_width - r_to))

    stopifnot(all(length(l_db == l_width), length(r_db) == r_width))
    db <- paste0(c(l_db, "&", r_db), collapse = "")
  }

  return(data.table(name = name, mfe = mfe, structure = db))
}

#' Get binding energy
#'
#' Calculates binding energy using \code{RNAduplex} from \code{ViennaRNA}
#'
#' @param name Hybrid name
#' @param L_sequence Left sequence
#' @param R_sequence Right sequence
#' @return data.table with binding energy
#' @import data.table
#' @export

get_mfe <- function(name, L_sequence, R_sequence) {
  input <- paste0(L_sequence, "\n", R_sequence)
  rnaduplex <- system("RNAduplex --noLP", input = input, intern = TRUE)

  # Get MFE
  rnaduplex <- gsub("\\s+", "_", rnaduplex)
  mfe <- sapply(strsplit(rnaduplex, "_"), "[", 5)
  if (mfe == "(") mfe <- sapply(strsplit(rnaduplex, "_"), "[", 6) # Positives < 10 have an extra space
  mfe <- as.numeric(gsub("\\(|\\)", "", mfe))

  return(data.table(name = name, mfe = mfe))
}

#' Get shuffles binding energies
#'
#' Calculates mean and standard deviation of 100 shuffled binding energies using \code{RNAduplex} from \code{ViennaRNA}
#'
#' @param name Hybrid name
#' @param L_sequence Left sequence
#' @param R_sequence Right sequence
#' @return data.table with mean and standard deviation of shuffles binding energies
#' @import data.table
#' @export

get_shuffled_mfe <- function(name, L_sequence, R_sequence) {
  L <- shuffle_sequence(L_sequence, number = 100, klet = 2)
  R <- shuffle_sequence(R_sequence, number = 100, klet = 2)

  shuffled_mfe.dt <- rbindlist(lapply(seq_along(L), function(i) get_mfe(name, L[i], R[i])))
  shuffled_mfe.dt[, `:=`(
    mean_shuffled_mfe = mean(mfe, na.rm = TRUE),
    sd_shuffled_mfe = sd(mfe, na.rm = TRUE)
  ),
  by = name
  ]

  return(unique(shuffled_mfe.dt[, .(name, mean_shuffled_mfe, sd_shuffled_mfe)]))
}
