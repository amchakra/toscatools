
#' Title
#'
#' @param aligned.bam
#' @param chimeric.junction
#'
#' @return
#'
#' @import data.table
#' @export
#'
#' @examples

ExtractHybrids <- function(aligned.bam, chimeric.junction) {

  # Read in data
  message("Reading in data")
  ga <- GenomicAlignments::readGAlignments(aligned.bam, use.names = TRUE, param = Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isMinusStrand = FALSE)))
  chimeric.dt <- fread(chimeric.junction,
                       col.names = c("chr_donorA", "brkpt_donorA", "strand_donorA", "chr_acceptorB", "brkpt_acceptorB", "strand_acceptorB",
                                     "junction_type", "repeat_left_lenA", "repeat_right_lenB",
                                     "read_name", "start_alnA", "cigar_alnA", "start_alnB", "cigar_alnB"))
  chimeric.dt <- chimeric.dt[!read_name %in% names(ga)] # Remove chimeric reads that are also in aligned - this is slow!

  # ==========
  # Extract genomic orientation
  # ==========
  message("Analysing genomic orientation reads")
  ga <- ga[GenomicAlignments::njunc(ga) == 1] # Ignore triplexes etc.
  grl <- GenomicAlignments::grglist(ga)
  gr <- unlist(grl)
  gr$name <- names(gr)
  names(gr) <- NULL

  L.gr <- gr[c(TRUE, FALSE)]
  R.gr <- gr[c(FALSE, TRUE)]

  L.gr$q_side <- "L"
  R.gr$q_side <- "R"

  aligned.grl <- GenomicRanges::GRangesList(L = L.gr, R = R.gr)
  aligned.dt <- ConvertToDataTable(aligned.grl)
  setkey(aligned.dt, name)

  # ==========
  # Extract chimeric orientation
  # ==========
  message("Analysing reverse orientation reads")

  # Some might be triplexes/circular RNAs. Ignore for now.....
  chimeric.dt <- chimeric.dt[!grepl("N", cigar_alnA)]
  chimeric.dt <- chimeric.dt[!grepl("N", cigar_alnB)]

  # Might need to introduce QC to ensure they don't overlap in query space
  L.gr <- with(chimeric.dt, GenomicRanges::GRanges(seqnames = chr_donorA,
                                                   ranges = GenomicRanges::shift(unlist(GenomicAlignments::cigarRangesAlongReferenceSpace(chimeric.dt$cigar_alnA, ops = c("M", "D"),
                                                                                                                                          reduce.ranges = TRUE)), start_alnA - 1), # include D to ignore short deletions
                                                   strand = "+",
                                                   name = read_name,
                                                   q_side = "L"))

  R.gr <- with(chimeric.dt, GenomicRanges::GRanges(seqnames = chr_acceptorB,
                                                   ranges = GenomicRanges::shift(unlist(GenomicAlignments::cigarRangesAlongReferenceSpace(chimeric.dt$cigar_alnB, ops = c("M", "D"),
                                                                                                                                          reduce.ranges = TRUE)), start_alnB - 1),
                                                   strand = "+",
                                                   name = read_name,
                                                   q_side = "R"))

  chimeric.grl <- GenomicRanges::GRangesList(L = L.gr, R = R.gr)
  chimeric.dt <- ConvertToDataTable(chimeric.grl)
  setkey(chimeric.dt, name)

  # Identify chimeric overlapping
  ol <- GenomicRanges::pintersect(chimeric.grl$L, chimeric.grl$R)
  ol.dt <- as.data.table(ol)[, .(name, hit)]
  setnames(ol.dt, "hit", "overlapping_hybrid")
  setkey(ol.dt, name)

  chimeric.dt <- merge(chimeric.dt, ol.dt, by = "name", all.x = TRUE)
  stopifnot(all(!is.na(chimeric.dt$overlapping_hybrid))) # Make sure all assigned

  # ==========
  # Extract chimeric orientation
  # ==========

  # Now combine and output
  aligned.dt[, orientation := "genomic"]
  chimeric.dt[, orientation := "reverse"]
  hybrids.dt <- rbindlist(list(aligned.dt, chimeric.dt), fill = TRUE)

  return(hybrids.dt)

}

#' Title
#'
#' @param hybrids.dt
#'
#' @return
#'
#' @import data.table
#' @export
#'

ReorientHybrids <- function(hybrids.dt) {

  correct.dt <- hybrids.dt[L_start <= R_start]
  incorrect.dt <- hybrids.dt[L_start > R_start]

  renamed <- gsub("^L_", "X_", names(incorrect.dt))
  renamed <- gsub("^R_", "L_", renamed)
  renamed <- gsub("^X_", "R_", renamed)

  setnames(incorrect.dt, renamed)

  reoriented.dt <- rbindlist(list(correct.dt, incorrect.dt), use.names = TRUE)

  stopifnot(all(reoriented.dt$L_start <= reoriented.dt$R_start))

  return(reoriented.dt)

}


#' Title
#'
#' @param hybrids.dt
#'
#' @return
#'
#' @import data.table
#' @export

RemovePCRDuplicates <- function(hybrids.dt) {

  # hybrids.dt[, rbc := sub(".*\\:", "", name)]
  hybrids.dt[, rbc := ifelse(grepl("_", name), sub(".*\\:", "", name), sub(".*\\:", "", name))] # Added "_" option for UMI tools
  unique.hybrids.dt <- unique(hybrids.dt, by = c("L_seqnames", "L_start", "L_end", "R_seqnames", "R_start", "R_end", "rbc", "orientation"))
  # unique.hybrids.dt <- unique(hybrids.dt, by = c("L_seqnames", "L_start", "R_seqnames", "R_start", "rbc", "orientation"))

  message("PCR duplication ratio: ", round(nrow(hybrids.dt)/nrow(unique.hybrids.dt), 2))
  return(unique.hybrids.dt)

}
