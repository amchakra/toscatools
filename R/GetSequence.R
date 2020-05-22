

#' GetSequence
#'
#' @param gr
#' @param fasta
#'
#' @return
#' @export

GetSequenceDS <- function(gr, fasta) {

  # Load as DNA stringset if filename supplied
  # Should add some checks

  if(is.character(fasta)) {

    fasta <- Biostrings::readDNAStringSet(fasta)

  }

  seq <- fasta[gr]
  gr$sequence <- seq # currently not outputting as character

  return(gr)

}



#' Title
#'
#' @param hybrids.dt
#' @param genome.dt
#'
#' @import data.table
#'
#' @return
#' @export

# Should only do this for intragenic?

GetSJMotifs <- function(hybrids.dt, genome.dt) {

  genomic.dt <- hybrids.dt[orientation == "genomic"]
  reverse.dt <- hybrids.dt[orientation == "reverse"]

  genomic.flank.dt <- genomic.dt[, `:=` (L_start = L_end + 1,
                                         L_end = L_end + 2,
                                         R_start = R_start - 2,
                                         R_end = R_start - 1)]

  reverse.flank.dt <- reverse.dt[, `:=` (L_start = L_start - 2,
                                         L_end = L_start - 1,
                                         R_start = R_end + 1,
                                         R_end = R_end + 2)]

  flank.dt <- rbindlist(list(genomic.flank.dt, reverse.flank.dt))

  flank.seq.dt <- GetSequence(flank.dt, genome.dt)
  flank.seq.dt$sj <- paste0(flank.seq.dt$L_sequence, flank.seq.dt$R_sequence)

  stopifnot(all(nchar(flank.seq.dt[orientation == "genomic"]$sj) == 4))
  stopifnot(all(nchar(flank.seq.dt[orientation == "reverse"]$sj) %in% 0:4)) # Some chimeric include e.g reverse orientation at the start of a gene (e.g. rRNA), so no splice motif

  setkey(hybrids.dt, name)
  setkey(flank.seq.dt, name)

  hybrids.dt <- merge(hybrids.dt, flank.seq.dt[, .(name, sj)], by = "name", all.x = TRUE)

  return(hybrids.dt)

}
