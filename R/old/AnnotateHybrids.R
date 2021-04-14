
#' Title
#'
#' @param hybrids.dt
#' @param ref.gr
#'
#' @return
#'
#' @import data.table

AnnotateHybrids <- function(hybrids.dt, ref.gr) {

  setkey(hybrids.dt, name, id)

  L.gr <- ConvertToGRanges(hybrids.dt, arm = "L")
  R.gr <- ConvertToGRanges(hybrids.dt, arm = "R")

  L.gr <- AnnotateHybridArm(arm.gr = L.gr, ref.gr = ref.gr)
  R.gr <- AnnotateHybridArm(arm.gr = R.gr, ref.gr = ref.gr)

  L.annot.dt <- as.data.table(L.gr)[, .(name, annotation)]
  L.annot <- stringr::str_split_fixed(L.annot.dt$annotation, pattern = "\\|", 5)
  L.annot.dt <- cbind(L.annot.dt, L.annot)
  setnames(L.annot.dt, c("annotation", "V1", "V2", "V3", "V4", "V5"), c("L_annotation", "L_transcriptid", "L_geneid", "L_genename", "L_biotype", "L_region"))
  setkey(L.annot.dt, name)
  annotated.hybrids.dt <- merge(hybrids.dt, L.annot.dt, by = "name")

  R.annot.dt <- as.data.table(R.gr)[, .(name, annotation)]
  R.annot <- stringr::str_split_fixed(R.annot.dt$annotation, pattern = "\\|", 5)
  R.annot.dt <- cbind(R.annot.dt, R.annot)
  setnames(R.annot.dt, c("annotation", "V1", "V2", "V3", "V4", "V5"), c("R_annotation", "R_transcriptid", "R_geneid", "R_genename", "R_biotype", "R_region"))
  setkey(R.annot.dt, name)
  annotated.hybrids.dt <- merge(annotated.hybrids.dt, R.annot.dt, by = "name")

  return(annotated.hybrids.dt)

}

#' Title
#'
#' @param arm.gr
#' @param ref.gr
#'
#' @return

AnnotateHybridArm <- function(arm.gr, ref.gr) {

  ol <- GenomicRanges::findOverlaps(GenomicRanges::resize(arm.gr, width = 1, fix = "start"), ref.gr)
  stopifnot(all(seq_along(arm.gr) %in% S4Vectors::queryHits(ol)))
  stopifnot(anyDuplicated(S4Vectors::queryHits(ol)) == 0)

  arm.gr$annotation <- as.character(NA)
  arm.gr[S4Vectors::queryHits(ol)]$annotation <- ref.gr[S4Vectors::subjectHits(ol)]$annotation
  stopifnot(all(!is.na(arm.gr$annotation)))

  return(arm.gr)

}
