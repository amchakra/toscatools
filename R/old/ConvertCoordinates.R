

#' Title
#'
#' @param seq.dt
#' @param genes.gr
#' @param cores
#'
#' @return
#' @import data.table

ConvertCoordinates <- function(seq.dt, genes.gr, cores) {

  cl <- parallel::makeForkCluster(cores)
  g.grl <- parallel::parLapply(cl = cl, 1:nrow(seq.dt), function(i) {

    message(i)

    gene <- seq.dt[i]$L_seqnames

    gene.gr <- genes.gr[genes.gr$fasta_id == gene]

    if(as.character(GenomicRanges::strand(gene.gr)) == "+") {

      gene.start <- as.integer(GenomicRanges::start(gene.gr))

      L.start = seq.dt[i]$L_start + gene.start - 1
      L.end = seq.dt[i]$L_end + gene.start - 1
      R.start = seq.dt[i]$R_start + gene.start - 1
      R.end = seq.dt[i]$R_end + gene.start - 1

      L.gr <- GenomicRanges::GRanges(seqnames = as.character(GenomicRanges::seqnames(gene.gr)),
                      ranges = IRanges::IRanges(start = L.start,
                                       end = L.end),
                      strand = "+")

      R.gr <- GenomicRanges::GRanges(seqnames = as.character(GenomicRanges::seqnames(gene.gr)),
                      ranges = IRanges::IRanges(start = R.start,
                                       end = R.end),
                      strand = "+")

    } else if(as.character(GenomicRanges::strand(gene.gr)) == "-") {

      gene.start <- as.integer(GenomicRanges::end(gene.gr))
      L.start <- gene.start - seq.dt[i]$L_start
      L.end <- gene.start - seq.dt[i]$L_end
      R.start <- gene.start - seq.dt[i]$R_start
      R.end <- gene.start - seq.dt[i]$R_end

      L.gr <- GenomicRanges::GRanges(seqnames = as.character(GenomicRanges::seqnames(gene.gr)),
                      ranges = IRanges::IRanges(start = L.end,
                                       end = L.start),
                      strand = "-")

      R.gr <- GenomicRanges::GRanges(seqnames = as.character(GenomicRanges::seqnames(gene.gr)),
                      ranges = IRanges::IRanges(start = R.end,
                                       end = R.start),
                      strand = "-")

    }

    L.gr$name <- seq.dt[i]$name
    R.gr$name <- seq.dt[i]$name

    return(c(L.gr, R.gr))

  })
  parallel::stopCluster(cl)

  g.grl <- sort(GenomicRanges::GRangesList(g.grl))
  return(g.grl)

}

#' Title
#'
#' @param seq.dt
#' @param genes.gr
#' @param filename
#' @param sam_tag
#'
#' @return
#' @import data.table

ExportGenomicBED <- function(seq.dt, genes.gr, filename, sam_tag = TRUE) {

  # Get genes.dt
  genes.dt <- as.data.table(genes.gr)[, .(fasta_id, seqnames, start, end, strand)]
  setkey(genes.dt, fasta_id)
  setkey(seq.dt, L_seqnames)

  coord.dt <- merge(seq.dt, genes.dt, by.x = "L_seqnames", by.y = "fasta_id")

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


#' Title
#'
#' @param g.grl
#' @param hybrids.dt
#' @param filename
#' @param sam_tag
#'
#' @return
#'
#' @import data.table

ExportBED <- function(g.grl, hybrids.dt, filename, sam_tag = TRUE) {

  rtracklayer::export.bed(g.grl, filename)

  bed.dt <- fread(filename,
                  col.names = c("chrom",	"chromStart", "chromEnd", "name", "score", "strand", "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts"))

  # TODO need to generalise. Currently needs MFE
  bed.dt$name <- hybrids.dt$name
  if(sam_tag == TRUE) {

    if("cluster" %in% names(hybrids.dt)) bed.dt$name <- paste0(bed.dt$name, "_", hybrids.dt$cluster)
    if("orientation" %in% names(hybrids.dt)) bed.dt$name <- paste0(bed.dt$name, "_", hybrids.dt$orientation)
    if("mfe" %in% names(hybrids.dt)) bed.dt$name <- paste0(bed.dt$name, "_", hybrids.dt$mfe)

  }
    # bed.dt$name <- paste0(bed.dt$name, "_", hybrids.dt$cluster, "_", hybrids.dt$orientation, "_", hybrids.dt$mfe)

  # bed.dt[, `:=` (name = hybrids.dt$name,
  #                # score = hybrids.dt$cluster,
  #                cluster = as.character(hybrids.dt$cluster))]
  # bed.dt[is.na(cluster), cluster := "None"]
  #
  # colour.dt <- data.table(cluster = 1:max(hybrids.dt$cluster, na.rm = TRUE),
  #                         colour = RColorBrewer::brewer.pal(9, "Set1"))
  # black.dt <- data.table(cluster = "None",
  #                        colour = "#000001")
  # colour.dt <- rbind(colour.dt, black.dt)
  # colour.dt <- cbind(colour.dt, t(col2rgb(colour.dt$colour)))
  # colour.dt[, rgb := paste0(red, ",", green, ",", blue)]
  # # colour.dt[cluster == "None", rgb := NA]
  #
  # bed.dt <- merge(bed.dt, colour.dt[, .(cluster, rgb)], by = "cluster")
  # bed.dt[, itemRgb := rgb]
  # bed.dt[, `:=` (cluster = NULL, rgb = NULL)]

  fwrite(bed.dt, file = filename, sep = "\t", quote = FALSE, col.names = FALSE)

}




