
#' Title
#'
#' @param aligned.bam
#'
#' @return
#'
#' @export

ExtractGenomicReads <- function(aligned.bam) {

  ga <- GenomicAlignments::readGAlignments(aligned.bam, use.names = TRUE, param = Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isMinusStrand = FALSE)))
  ga <- ga[GenomicAlignments::njunc(ga) == 1]
  grl <- GenomicAlignments::grglist(ga)
  gr <- unlist(grl)
  gr$name <- names(gr)
  names(gr) <- NULL

  L.gr <- gr[c(TRUE, FALSE)]
  R.gr <- gr[c(FALSE, TRUE)]

  L.gr$q_side <- "L"
  R.gr$q_side <- "R"

  return(GenomicRanges::GRangesList(L = L.gr, R = R.gr))

}
