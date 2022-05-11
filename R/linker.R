#' Title
#'
#' @param bam BAM file
#' @param threads Threads
#' @param export_fasta Whether to export a FASTA/Q file
#' @param fasta_filename FASTA/Q filename
#'
#' @return
#' @export
#'
#' @examples

convert_linkerbam_to_hybrids <- function(bam, threads = 4, export_fasta = FALSE, fasta_filename) {

  ga <- GenomicAlignments::readGAlignments(bam,
                                           use.names = TRUE,
                                           param = Rsamtools::ScanBamParam(what = c("seq", "qual")))

  ga <- ga[GenomicAlignments::njunc(ga) == 0]

  cig <- GenomicAlignments::cigarRangesAlongQuerySpace(GenomicAlignments::cigar(ga),
                                                       after.soft.clipping = TRUE,
                                                       drop.empty.ranges = TRUE,
                                                       reduce.ranges = TRUE)
  stopifnot(all(S4Vectors::elementNROWS(cig) == 1))
  cig <- unlist(IRanges::IRangesList(cig))

  # Probably a vectorised way of doing this byt seq[cig] doesn't work?
  S4Vectors::mcols(ga)$seq <- Biostrings::DNAStringSet(parallel::mclapply(seq_along(cig), function(i) { S4Vectors::mcols(ga)$seq[[i]][cig[i]] }, mc.cores = threads))
  S4Vectors::mcols(ga)$qual <- Biostrings::BStringSet(parallel::mclapply(seq_along(cig), function(i) { S4Vectors::mcols(ga)$qual[[i]][cig[i]] }, mc.cores = threads))

  dt <- as.data.table(ga, keep.rownames = TRUE)
  dt[, `:=` (name = gsub("^L\\.|^R\\.", "", rn),
             arm = tstrsplit(rn, "\\.")[[1]])][, rn := NULL]

  hybrids.dt <- merge(dt[arm == "L"], dt[arm == "R"], by = "name")
  hybrids.dt[, seq := paste0(seq.x, seq.y)]
  hybrids.dt[, seq_width := nchar(seq)]

  if(export_fasta) {

    hybrids.fa <- Biostrings::DNAStringSet(paste0(hybrids.dt$seq.x, hybrids.dt$seq.y))
    names(hybrids.fa) <- hybrids.dt$name

    hybrids.qual <- Biostrings::BStringSet(paste0(hybrids.dt$qual.x, hybrids.dt$qual.y))

    Biostrings::writeXStringSet(hybrids.fa, fasta_filename,
                                format = "fastq",
                                qualities = hybrids.qual,
                                compress = TRUE)

  }

  return(hybrids.dt)

}


