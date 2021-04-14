

#' Title
#'
#' @param chimeric.junction
#'
#' @return
#'
#' @import data.table

ExtractChimericReads <- function(chimeric.junction) {

  chimeric.dt <- fread(chimeric.junction,
                       col.names = c("chr_donorA", "brkpt_donorA", "strand_donorA", "chr_acceptorB", "brkpt_acceptorB", "strand_acceptorB",
                                     "junction_type", "repeat_left_lenA", "repeat_right_lenB",
                                     "read_name", "start_alnA", "cigar_alnA", "start_alnB", "cigar_alnB"))


  # Some might be triplexes/circular RNAs. Ignore for now.....
  chimeric.dt <- chimeric.dt[!grepl("N", cigar_alnA)]
  chimeric.dt <- chimeric.dt[!grepl("N", cigar_alnB)]

  # Might need to introduce QC to ensure they don't overlap in query space
  L.gr <- with(chimeric.dt, GenomicRanges::GRanges(seqnames = chr_donorA,
                                    ranges = GenomicRanges::shift(unlist(GenomicAlignments::cigarRangesAlongReferenceSpace(chimeric.dt$cigar_alnA, ops = c("M", "D"), reduce.ranges = TRUE)),
                                                   start_alnA - 1), # include D to ignore short deletions
                                    strand = "+",
                                    name = read_name,
                                    q_side = "L"))

  R.gr <- with(chimeric.dt, GenomicRanges::GRanges(seqnames = chr_acceptorB,
                                    ranges = GenomicRanges::shift(unlist(GenomicAlignments::cigarRangesAlongReferenceSpace(chimeric.dt$cigar_alnB, ops = c("M", "D"), reduce.ranges = TRUE)),
                                                   start_alnB - 1),
                                    strand = "+",
                                    name = read_name,
                                    q_side = "R"))

  return(GenomicRanges::GRangesList(L = L.gr, R = R.gr))

}
