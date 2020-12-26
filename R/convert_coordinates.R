#' Title
#'
#' @param hybrids.dt
#' @param genes.gr
#' @param filename
#' @param sam_tag
#'
#' @return
#' @export
#' @import data.table
#'
#' @examples

convert_coordinates <- function(hybrids.dt, genes.gr, filename, sam_tag = TRUE) {

  # Get genes.dt
  genes.dt <- as.data.table(genes.gr)[, .(fasta_id, seqnames, start, end, strand)]
  setkey(genes.dt, fasta_id)
  setkey(hybrids.dt, L_seqnames)

  coord.dt <- merge(hybrids.dt, genes.dt, by.x = "L_seqnames", by.y = "fasta_id")

  # Do positive strand genes
  coord.dt[strand == "+", `:=` (L.start = L_start + start - 1,
                                L.end = L_end + start - 1,
                                R.start = R_start + start - 1,
                                R.end = R_end + start - 1),
           by = name]

  # Do negative strand genes
  # Flip R and L for BEDgraph
  coord.dt[strand == "-", `:=` (R.end = end - L_start,
                                R.start = end - L_end,
                                L.end = end - R_start,
                                L.start = end - R_end),
           by = name]

  # Calculate bedgraph blocks
  # NB BED is 0 based
  coord.dt[, `:=` (chrom = seqnames,
                   chromStart = L.start - 1,
                   chromEnd = R.end,
                   name = name,
                   score = 0,
                   strand = strand,
                   thickStart = L.start - 1,
                   thickEnd = R.end,
                   itemRgb = 0,
                   blockCount = 2,
                   blockSizes = paste0(L.end - L.start + 1, ",", R.end - R.start + 1),
                   blockStarts = paste0("0,", R.start - L.start)),
           by = name]

  bed.dt <- coord.dt[, .(chrom,	chromStart, chromEnd, name, score, strand, thickStart, thickEnd, itemRgb, blockCount, blockSizes, blockStarts)]

  if(sam_tag == TRUE) {

    if("cluster" %in% names(coord.dt)) bed.dt$name <- paste0(bed.dt$name, "_", coord.dt$cluster)
    if("orientation" %in% names(coord.dt)) bed.dt$name <- paste0(bed.dt$name, "_", coord.dt$orientation)
    if("mfe" %in% names(coord.dt)) bed.dt$name <- paste0(bed.dt$name, "_", coord.dt$mfe)

  }

  fwrite(bed.dt, file = filename, sep = "\t", quote = FALSE, col.names = FALSE, scipen = 999) # Needed to add scipen as otherwise some coordinates end up scientific

}
