

#' Title
#'
#' @param hybrids.dt
#' @param regions.gr
#'
#' @return
#' @export
#' @import data.table
#'
#' @examples

annotate_hybrids <- function(hybrids.dt, regions.gr) {

  # ==========
  # Left
  # ==========

  L.gr <- primavera::convert_to_granges(hybrids.dt, arm = "L", genomic = TRUE)

  # Sync seqlevels
  common.seqlevels <- unique(c(GenomeInfoDb::seqlevels(L.gr),
                               GenomeInfoDb::seqlevels(regions.gr)))
  GenomeInfoDb::seqlevels(L.gr) <- common.seqlevels
  GenomeInfoDb::seqlevels(regions.gr) <- common.seqlevels

  ol <- GenomicRanges::findOverlaps(GenomicRanges::resize(L.gr, width = 1, fix = "center"), regions.gr)
  stopifnot(all(!duplicated(S4Vectors::queryHits(ol))))

  match.gr <- regions.gr[S4Vectors::subjectHits(ol)]
  hybrids.dt[S4Vectors::queryHits(ol), `:=` (L_region = match.gr$type,
                                  L_gene_id = match.gr$gene_id,
                                  L_gene_name = match.gr$gene_name,
                                  L_biotype = match.gr$biotype)]

  # bug fix for iCount missing a couple of intergenic regions
  hybrids.dt[is.na(L_region), `:=` (L_region = "intergenic",
                                    L_gene_id = ".",
                                    L_gene_name = "None",
                                    L_biotype = "")]

  # ==========
  # Right
  # ==========

  R.gr <- primavera::convert_to_granges(hybrids.dt, arm = "R", genomic = TRUE)

  # Sync seqlevels
  common.seqlevels <- unique(c(GenomeInfoDb::seqlevels(R.gr),
                               GenomeInfoDb::seqlevels(regions.gr)))
  GenomeInfoDb::seqlevels(R.gr) <- common.seqlevels
  GenomeInfoDb::seqlevels(regions.gr) <- common.seqlevels

  ol <- GenomicRanges::findOverlaps(GenomicRanges::resize(R.gr, width = 1, fix = "center"), regions.gr)
  stopifnot(all(!duplicated(S4Vectors::queryHits(ol))))

  match.gr <- regions.gr[S4Vectors::subjectHits(ol)]
  hybrids.dt[S4Vectors::queryHits(ol), `:=` (R_region = match.gr$type,
                                  R_gene_id = match.gr$gene_id,
                                  R_gene_name = match.gr$gene_name,
                                  R_biotype = match.gr$biotype)]

  # bug fix for iCount missing a couple of intergenic regions
  hybrids.dt[is.na(R_region), `:=` (R_region = "intergenic",
                                    R_gene_id = ".",
                                    R_gene_name = "None",
                                    R_biotype = "")]

  # rRNA
  hybrids.dt[grep("rRNA|rDNA", L_seqnames), `:=` (L_gene_id = L_seqnames,
                                                  L_region = "rRNA",
                                                  L_gene_name = L_seqnames,
                                                  L_biotype = "rRNA")]
  hybrids.dt[grep("rRNA|rDNA", R_seqnames), `:=` (R_gene_id = R_seqnames,
                                                  R_region = "rRNA",
                                                  R_gene_name = R_seqnames,
                                                  R_biotype = "rRNA")]

  stopifnot(all(!is.na(c(hybrids.dt$L_region, hybrids.dt$R_region))))
  return(hybrids.dt)

}
