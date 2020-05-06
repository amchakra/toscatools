

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

ShuffleSequence <- function(sequence, klet = 2, seed = 42) {

  system(paste0("uShuffle -seed ", seed, " -k ", klet, " -n 1 -s ", sequence), intern = TRUE)

}

#' Title
#'
#' @param sequence1
#' @param sequence2
#'
#' @return
#' @export

GetMFE <- function(sequence1, sequence2) {

  input <- paste0(sequence1, "\n", sequence2)
  rnaduplex <- system("RNAduplex --noLP", input = input, intern = TRUE)

  rnaduplex <- gsub("\\s+", "_", rnaduplex)
  mfe <- sapply(strsplit(rnaduplex, "_"), "[", 5)

  mfe <- as.numeric(gsub("\\(|\\)", "", mfe))
  # if(mfe > 0) mfe <- 0

  return(mfe)

}
