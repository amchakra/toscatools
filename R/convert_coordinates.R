#' Convert coordinates
#'
#' Convert from transcriptomic to genomic coordinates
#'
#' @param hybrids.dt Hybrids data.table
#' @param genes.gr Transcript gtf as GRanges
#' @return \code{hybrids.dt} with genomic coordinates in \code{[L|R]_genomic_*} columns
#' @export
#' @import data.table

convert_coordinates <- function(hybrids.dt, genes.gr) {

  # Get genes.dt
  genes.dt <- as.data.table(genes.gr)[, .(fasta_id, seqnames, start, end, strand)]
  setkey(genes.dt, fasta_id)
  setkey(hybrids.dt, L_seqnames)

  coord.dt <- merge(hybrids.dt, genes.dt, by.x = "L_seqnames", by.y = "fasta_id")

  # Do L
  setkey(hybrids.dt, L_seqnames)
  coord.dt <- merge(hybrids.dt, genes.dt, by.x = "L_seqnames", by.y = "fasta_id")

  coord.dt[strand == "+", `:=`(
    L_genomic_seqnames = seqnames,
    L_genomic_start = L_start + start - 1,
    L_genomic_end = L_end + start - 1,
    L_genomic_strand = "+"
  ),
  by = name
  ]
  coord.dt[strand == "-", `:=`(
    L_genomic_seqnames = seqnames,
    L_genomic_start = end - L_end,
    L_genomic_end = end - L_start,
    L_genomic_strand = "-"
  ),
  by = name
  ]

  # Do R
  coord.dt[, `:=`(seqnames = NULL, start = NULL, end = NULL, strand = NULL)]
  setkey(coord.dt, R_seqnames)
  coord.dt <- merge(coord.dt, genes.dt, by.x = "R_seqnames", by.y = "fasta_id")

  coord.dt[strand == "+", `:=`(
    R_genomic_seqnames = seqnames,
    R_genomic_start = R_start + start - 1,
    R_genomic_end = R_end + start - 1,
    R_genomic_strand = "+"
  ),
  by = name
  ]
  coord.dt[strand == "-", `:=`(
    R_genomic_seqnames = seqnames,
    R_genomic_start = end - R_end,
    R_genomic_end = end - R_start,
    R_genomic_strand = "-"
  ),
  by = name
  ]

  # Tidy up
  coord.dt[, `:=`(seqnames = NULL, start = NULL, end = NULL, strand = NULL)]

  col.names <- names(coord.dt)
  ordered.col.names <- c(
    grep("^L|^R", col.names, invert = TRUE, value = TRUE),
    grep("^L", col.names, value = TRUE),
    grep("^R", col.names, value = TRUE)
  )
  stopifnot(all(ordered.col.names %in% col.names))
  setcolorder(coord.dt, ordered.col.names)

  return(coord.dt)
}

#' Export BED
#'
#' Calculate and export BED12 file with genomic coordinates
#'
#' @param hybrids.dt Hybrids data.table with genomic coordinates
#' @param filename BED filename
#' @param sam_tag Add SAM tags to \code{name} column in BED12
#' @return BED12 file exported
#' @export
#' @import data.table

export_genomic_bed <- function(hybrids.dt, filename, sam_tag = TRUE) {
  if (!all(hybrids.dt$L_seqnames == hybrids.dt$R_seqnames)) stop("Hybrids should all be intragenic")
  if (any(grepl("rRNA|rDNA", c(hybrids.dt$L_seqnames, hybrids.dt$R_seqnames)))) stop("Hybrids should not include rRNA")
  if (!"_genomic_" %in% names(hybrids.dt)) stop("Genomic coordinates should have been calculated")

  coord.dt <- hybrids.dt[, `:=`(
    L.start = L_genomic_start,
    L.end = L_genomic_end,
    R.start = R_genomic_start,
    R.end = R_genomic_end,
    seqnames = L_genomic_seqnames,
    strand = L_genomic_strand
  ),
  by = name
  ]

  # # Flip R and L for BEDgraph for negative strand
  coord.dt[strand == "-", `:=`(
    R.start = L_genomic_start,
    R.end = L_genomic_end,
    L.start = R_genomic_start,
    L.end = R_genomic_end
  ),
  by = name
  ]

  # Calculate bedgraph blocks
  # NB BED is 0 based
  coord.dt[, `:=`(
    chrom = seqnames,
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
    blockStarts = paste0("0,", R.start - L.start)
  ),
  by = name
  ]

  bed.dt <- coord.dt[, .(chrom, chromStart, chromEnd, name, score, strand, thickStart, thickEnd, itemRgb, blockCount, blockSizes, blockStarts)]

  if (sam_tag == TRUE) {
    if ("sample" %in% names(coord.dt)) {
      bed.dt$name <- paste0(bed.dt$name, "_", coord.dt$sample)
    } else {
      bed.dt$name <- paste0(bed.dt$name, "_")
    }
    if ("cluster" %in% names(coord.dt)) {
      bed.dt$name <- paste0(bed.dt$name, "_", coord.dt$cluster)
    } else {
      bed.dt$name <- paste0(bed.dt$name, "_")
    }
    if ("orientation" %in% names(coord.dt)) {
      bed.dt$name <- paste0(bed.dt$name, "_", coord.dt$orientation)
    } else {
      bed.dt$name <- paste0(bed.dt$name, "_")
    }
    if ("mfe" %in% names(coord.dt)) {
      bed.dt$name <- paste0(bed.dt$name, "_", coord.dt$mfe)
    } else {
      bed.dt$name <- paste0(bed.dt$name, "_")
    }
  }

  fwrite(bed.dt, file = filename, sep = "\t", quote = FALSE, col.names = FALSE, scipen = 999) # Needed to add scipen as otherwise some coordinates end up scientific
}
