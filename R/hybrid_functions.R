# ==========
# Hybrid utility functions
# ==========

#' Reorients hybrids
#'
#' Reorients hybrids, so left arm always before right arm for intragenic hybrids and left arm seqname alphabetically before right arm seqname for intergenic
#'
#' @param hybrids.dt Hybrids data.table
#' @return Reoriented hybrids.dt
#' @import data.table
#' @export

reorient_hybrids <- function(hybrids.dt) {

  # First do starts
  correct.dt <- hybrids.dt[L_start <= R_start]
  incorrect.dt <- hybrids.dt[L_start > R_start]

  renamed <- gsub("^L_", "X_", names(incorrect.dt))
  renamed <- gsub("^R_", "L_", renamed)
  renamed <- gsub("^X_", "R_", renamed)

  setnames(incorrect.dt, renamed)

  reoriented.dt <- rbindlist(list(correct.dt, incorrect.dt), use.names = TRUE)

  stopifnot(all(reoriented.dt$L_start <= reoriented.dt$R_start))

  # Then do subject (to make sure intergenics in same order)
  correct.dt <- reoriented.dt[L_seqnames <= R_seqnames]
  incorrect.dt <- reoriented.dt[L_seqnames > R_seqnames]

  renamed <- gsub("^L_", "X_", names(incorrect.dt))
  renamed <- gsub("^R_", "L_", renamed)
  renamed <- gsub("^X_", "R_", renamed)

  setnames(incorrect.dt, renamed)

  reoriented.dt <- rbindlist(list(correct.dt, incorrect.dt), use.names = TRUE)
  stopifnot(all(reoriented.dt$L_subject <= reoriented.dt$R_subject))
  stopifnot(nrow(reoriented.dt) == nrow(hybrids.dt))

  return(reoriented.dt)
}


#' Converts GRangesList to hybrid data.table
#'
#' @param hybrids.grl GRangesList of length 2 with L and R GRanges
#' @return hybrids data.table
#' @import data.table
#' @export

convert_to_datatable <- function(hybrids.grl) {
  L.dt <- as.data.table(hybrids.grl$L)
  R.dt <- as.data.table(hybrids.grl$R)
  setkey(L.dt, name)
  setkey(R.dt, name)

  hybrids.dt <- merge(L.dt, R.dt, by = "name")

  # Rename L
  L.names <- grep(".x", names(hybrids.dt), value = TRUE)
  L.renamed <- paste0("L_", gsub(".x", "", L.names))
  setnames(hybrids.dt, L.names, L.renamed)

  # Rename R
  R.names <- grep(".y", names(hybrids.dt), value = TRUE)
  R.renamed <- paste0("R_", gsub(".y", "", R.names))
  setnames(hybrids.dt, R.names, R.renamed)

  return(hybrids.dt)
}


#' Converts hybrids data.table to GRanges
#'
#' Converts hybrids data.table to left or right arm GRanges, either using transcriptomic or genomic coordinates
#'
#' @param hybrids.dt hybrids data.table
#' @param arm Convert left or right arm
#' @param genomic Use genomic coordinates
#' @return GRanges of left or right hybrid arm
#' @import data.table
#' @export


convert_to_granges <- function(hybrids.dt, arm = c("L", "R"), genomic = FALSE) {
  if (!arm %in% c("L", "R")) stop("[ERROR]: arm should be L or R")

  if (!genomic) {
    if (arm == "L") {
      L.dt <- hybrids.dt[, grep("^R_", names(hybrids.dt), invert = TRUE), with = FALSE]
      setnames(L.dt, names(L.dt), gsub("^L_", "", names(L.dt)))
      L.gr <- GenomicRanges::GRanges(L.dt)

      return(L.gr)
    } else if (arm == "R") {
      R.dt <- hybrids.dt[, grep("^L_", names(hybrids.dt), invert = TRUE), with = FALSE]
      setnames(R.dt, names(R.dt), gsub("^R_", "", names(R.dt)))
      R.gr <- GenomicRanges::GRanges(R.dt)

      return(R.gr)
    }
  } else if (genomic) {
    if (!any(grepl("_genomic_", names(hybrids.dt)))) stop("Genomic coordinates should have been calculated")
    if (arm == "L") {
      L.dt <- hybrids.dt[, grep("^R_", names(hybrids.dt), invert = TRUE), with = FALSE]
      setnames(L.dt, names(L.dt), gsub("^L_genomic_", "", names(L.dt)))
      L.gr <- GenomicRanges::GRanges(L.dt)

      return(L.gr)
    } else if (arm == "R") {
      R.dt <- hybrids.dt[, grep("^L_", names(hybrids.dt), invert = TRUE), with = FALSE]
      setnames(R.dt, names(R.dt), gsub("^R_genomic_", "", names(R.dt)))
      R.gr <- GenomicRanges::GRanges(R.dt)

      return(R.gr)
    }
  }
}
