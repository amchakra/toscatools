
#' Title
#'
#' @param aligned.bam
#' @param chimeric.junction
#'
#' @return
#'
#' @import data.table


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
#' @param aligned.bam
#'
#' @return
#'
#' @import data.table
#'

ExtractHybridsWithinBAM <- function(aligned.bam, remove.repeats = FALSE) {

  ga <- GenomicAlignments::readGAlignments(aligned.bam, use.names = TRUE, param = Rsamtools::ScanBamParam(what = c("flag", "seq")))

  # 0 = mapped positive strand
  # 16 = mapped negative strand
  # 2046 = chimeric second half positive strand
  # 2064 = chimeric second half negative strand

  # ==========
  # Reverse orientation reads
  # ==========

  R.chimeric.ga <- ga[S4Vectors::mcols(ga)$flag %in% c(2048, 2064)]
  R.chimeric.ga <- R.chimeric.ga[GenomicAlignments::njunc(R.chimeric.ga) == 0] # Remove SJ in one arm
  L.chimeric.ga <- ga[S4Vectors::mcols(ga)$flag %in% c(0, 16)]
  L.chimeric.ga <- L.chimeric.ga[names(L.chimeric.ga) %in% names(R.chimeric.ga)]
  L.chimeric.ga <- L.chimeric.ga[GenomicAlignments::njunc(L.chimeric.ga) == 0] # Remove SJ in other arm

  # UMI tools dedup doesn't collapse the supplementary alignment so need to remove missing ones
  # missing <- names(R.chimeric.gr)[!names(R.chimeric.gr) %in% names(L.chimeric.gr)]
  # Also need to do this to filter SJ in one arm or other
  R.chimeric.ga <- R.chimeric.ga[names(R.chimeric.ga) %in% names(L.chimeric.ga)]

  # Sort and check they match up
  L.chimeric.ga <- L.chimeric.ga[sort(names(L.chimeric.ga))]
  R.chimeric.ga <- R.chimeric.ga[sort(names(R.chimeric.ga))]

  stopifnot(all(names(L.chimeric.ga) == names(R.chimeric.ga)))

  # Convert to GRanges
  L.chimeric.gr <- GenomicAlignments::granges(L.chimeric.ga)
  L.chimeric.gr$name <- names(L.chimeric.gr)
  names(L.chimeric.gr) <- NULL
  L.chimeric.gr$q_side <- "L"

  R.chimeric.gr <- GenomicAlignments::granges(R.chimeric.ga)
  R.chimeric.gr$name <- names(R.chimeric.gr)
  names(R.chimeric.gr) <- NULL
  R.chimeric.gr$q_side <- "R"

  # Merge and convert to data.table
  # chimeric.gr <- c(L.chimeric.gr, R.chimeric.gr)
  # dupl <- duplicated(chimeric.gr$name) # Check there is one of each
  # stopifnot(sum(dupl == FALSE) == sum(dupl == TRUE))

  chimeric.grl <- GenomicRanges::GRangesList(L = L.chimeric.gr, R = R.chimeric.gr)
  chimeric.dt <- ConvertToDataTable(chimeric.grl)
  setkey(chimeric.dt, name)

  # Identify chimeric overlapping
  ol <- GenomicRanges::pintersect(chimeric.grl$L, chimeric.grl$R)
  ol.dt <- as.data.table(ol)[, .(name, hit)]
  setnames(ol.dt, "hit", "overlapping_hybrid")
  setkey(ol.dt, name)

  chimeric.dt <- merge(chimeric.dt, ol.dt, by = "name", all.x = TRUE)

  # ==========
  # Genomic orientation reads
  # ==========

  chimeric.reads <- names(ga)[S4Vectors::mcols(ga)$flag %in% c(2048, 2064)]
  aligned.ga <- ga[!names(ga) %in% chimeric.reads]
  aligned.ga <- ga[GenomicAlignments::njunc(ga) == 1] # Ignore triplexes etc.

  # Add in filter
  if(remove.repeats == TRUE) {

    nuc <- c("A", "G", "C", "T")
    filter <- sapply(paste0(expand.grid(nuc, nuc)$Var1, expand.grid(nuc, nuc)$Var2), function(x) paste0(rep(x, each = 5), collapse = ""))
    repeats <- Biostrings::vcountPDict(Biostrings::PDict(filter), S4Vectors::mcols(aligned.ga)$seq, collapse = 2)
    aligned.ga <- aligned.ga[repeats == 0]

  }

  aligned.grl <- GenomicAlignments::grglist(aligned.ga)
  aligned.gr <- unlist(aligned.grl)
  aligned.gr$name <- names(aligned.gr)
  names(aligned.gr) <- NULL

  L.aligned.gr <- aligned.gr[c(TRUE, FALSE)]
  R.aligned.gr <- aligned.gr[c(FALSE, TRUE)]

  L.aligned.gr$q_side <- "L"
  R.aligned.gr$q_side <- "R"

  aligned.grl <- GenomicRanges::GRangesList(L = L.aligned.gr, R = R.aligned.gr)
  aligned.dt <- ConvertToDataTable(aligned.grl)
  setkey(aligned.dt, name)

  # Now combine and output
  aligned.dt[, orientation := "genomic"]
  chimeric.dt[, orientation := "reverse"]
  hybrids.dt <- rbindlist(list(aligned.dt, chimeric.dt), fill = TRUE)
  hybrids.dt <- hybrids.dt[L_strand == "+" & R_strand == "+"] # Both have to be positive strand

  return(hybrids.dt)

}


#' Title
#'
#' @param hybrids.dt
#'
#' @return
#'
#' @import data.table

RemovePCRDuplicates <- function(hybrids.dt) {

  # hybrids.dt[, rbc := sub(".*\\:", "", name)]
  hybrids.dt[, rbc := ifelse(grepl("_", name), sub(".*\\_", "", name), sub(".*\\:", "", name))] # Added "_" option for UMI tools
  unique.hybrids.dt <- unique(hybrids.dt, by = c("L_seqnames", "L_start", "L_end", "R_seqnames", "R_start", "R_end", "rbc", "orientation"))
  # unique.hybrids.dt <- unique(hybrids.dt, by = c("L_seqnames", "L_start", "R_seqnames", "R_start", "rbc", "orientation"))

  message("PCR duplication ratio: ", round(nrow(hybrids.dt)/nrow(unique.hybrids.dt), 2))
  return(unique.hybrids.dt)

}
