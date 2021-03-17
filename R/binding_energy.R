

#' Title
#'
#' @param hybrids.dt
#' @param genome.dt
#'
#' @return
#'
#' @import data.table
#' @export

GetSequence <- function(hybrids.dt, genome.dt) {

  # Get L
  hybrids.dt[, gene_id := L_seqnames]
  setkey(genome.dt, gene_id)
  setkey(hybrids.dt, gene_id)
  seq.dt <- genome.dt[hybrids.dt]
  seq.dt[, `:=`(L_sequence, stringr::str_sub(sequence, start = L_start, end = L_end))]
  seq.dt[, `:=` (c("gene_id", "sequence"), NULL)]

  # Get R
  # Do separately in case on different genes
  seq.dt[, gene_id := R_seqnames]
  setkey(seq.dt, gene_id)
  seq.dt <- genome.dt[seq.dt]
  seq.dt[, `:=`(R_sequence, stringr::str_sub(sequence, start = R_start, end = R_end))]
  seq.dt[, `:=` (c("gene_id", "sequence"), NULL)]

  return(seq.dt)

}



#' Title
#'
#' @param sequence
#' @param klet
#' @param seed
#'
#' @return
#' @export

ShuffleSequence <- function(sequence, number = 1, klet = 2, seed = 42) {

  system(paste0("uShuffle -seed ", seed, " -k ", klet, " -n ", number, " -s ", sequence), intern = TRUE)

}

#' Title
#'
#' @param sequence1
#' @param sequence2
#'
#' @return
#' @export

GetMFE <- function(L_sequence, R_sequence) {

  input <- paste0(L_sequence, "\n", R_sequence)
  rnaduplex <- system("RNAduplex --noLP", input = input, intern = TRUE)

  rnaduplex <- gsub("\\s+", "_", rnaduplex) 
  mfe <- sapply(strsplit(rnaduplex, "_"), "[", 5) # Need to fix for positives < 10 ? as they have an extra space
  mfe <- as.numeric(gsub("\\(|\\)", "", mfe))
  
  # Get dot.bracket
  db <- sapply(strsplit(rnaduplex, "_"), "[", 1)
  
  # Get from/to positions
  l <- sapply(strsplit(rnaduplex, "_"), "[", 2)
  r <- sapply(strsplit(rnaduplex, "_"), "[", 4)
  l_from <- as.numeric(sapply(strsplit(l, ","), "[", 1))
  l_to <- as.numeric(sapply(strsplit(l, ","), "[", 2))
  r_from <- as.numeric(sapply(strsplit(r, ","), "[", 1))
  r_to <- as.numeric(sapply(strsplit(r, ","), "[", 2))
  l_width <- nchar(L_sequence)
  r_width <- nchar(R_sequence)
  
  # Pad db with . according to from/to
  l_db <- sapply(strsplit(db, "&"), "[", 1)
  r_db <- sapply(strsplit(db, "&"), "[", 2)
  l_db <- unlist(strsplit(l_db, "")) # convert string to vector of characters
  r_db <- unlist(strsplit(r_db, ""))
  l_db <- c(rep(".", l_from - 1L), l_db, rep(".", l_width - l_to))
  r_db <- c(rep(".", r_from - 1L), r_db, rep(".", r_width - r_to))
  stopifnot(all(length(l_db == l_width), length(r_db) == r_width))
  db <- paste0(c(l_db, "&", r_db), collapse = "")
  # if(mfe > 0) mfe <- 0

  return(list(mfe = mfe, db = db))

}
