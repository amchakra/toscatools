

#' Title
#'
#' @param hybrids.x.dt
#' @param hybrids.y.dt
#'
#' @return
#' @export
#' @import data.table
#'
#' @examples

find_hybrid_overlaps <- function(hybrids.x.dt, hybrids.y.dt) {

  # Make sure uniformly oriented
  hybrids.x.dt <- reorient_hybrids(hybrids.x.dt)
  hybrids.y.dt <- reorient_hybrids(hybrids.y.dt)

  # Convert to GRanges
  L.x.gr <- ConvertToGRanges(hybrids.x.dt, "L")
  R.x.gr <- ConvertToGRanges(hybrids.x.dt, "R")

  L.y.gr <- ConvertToGRanges(hybrids.y.dt, "L")
  R.y.gr <- ConvertToGRanges(hybrids.y.dt, "R")

  all.seqlevels <- unique(c(hybrids.x.dt$L_seqnames, hybrids.x.dt$R_seqnames, hybrids.y.dt$L_seqnames, hybrids.y.dt$R_seqnames))
  GenomeInfoDb::seqlevels(L.x.gr) <- all.seqlevels
  GenomeInfoDb::seqlevels(R.x.gr) <- all.seqlevels
  GenomeInfoDb::seqlevels(L.y.gr) <- all.seqlevels
  GenomeInfoDb::seqlevels(R.y.gr) <- all.seqlevels

  # Get overlaps
  L.ol <- as.data.table(GenomicRanges::findOverlaps(L.x.gr, L.y.gr))
  R.ol <- as.data.table(GenomicRanges::findOverlaps(R.x.gr, R.y.gr))

  L.ol[, q_s := paste0(queryHits, "_", subjectHits)]
  R.ol[, q_s := paste0(queryHits, "_", subjectHits)]

  # Both arms overlap
  ol <- L.ol[L.ol$q_s %in% R.ol$q_s]
  ol[, q_s := NULL]

  return(ol)

}
