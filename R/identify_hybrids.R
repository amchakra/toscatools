# ==========
# Functions to identify hybrids from Blat output
# ==========

#' Load blast8 file
#'
#' Load blast8 file as data.table
#'
#' @param blast8 blast8 file
#' @return blast8 data.table
#' @export
#' @import data.table

load_blast8 <- function(blast8) {
  dt <- data.table::fread(
    input = blast8,
    col.names = c(
      "query",
      "subject",
      "identity",
      "alignment_length",
      "mismatches",
      "gap_openings",
      "q_start",
      "q_end",
      "s_start",
      "s_end",
      "evalue",
      "bit_score"
    )
  )
  setkey(dt, query)
  return(dt)
}

#' Add read lengths
#'
#' Add FASTA read lengths to blast8 data.table
#'
#' @param blast.dt blast8 data.table
#' @param fasta FASTA file
#' @return blast8 data.table with \code{readlength} column
#' @export
#' @import data.table

add_read_lengths <- function(blast.dt, fasta) {

  # Load fasta for read lengths
  fasta <- Biostrings::readDNAStringSet(fasta)
  fasta.dt <- data.table(
    query = names(fasta),
    readlength = BiocGenerics::width(fasta)
  )
  setkey(fasta.dt, query)

  # Add to blast.dt
  setorder(blast.dt, query)
  setkey(blast.dt, query)
  stopifnot(all(blast.dt$query %in% fasta.dt$query))
  blast.dt <- fasta.dt[blast.dt]

  return(blast.dt)
}

#' Calculate blast8 metrics
#'
#' Calcualte blast8 metrics for filtereing
#'
#' @param blast.dt blast8 data.table
#' @return blast8 data.table with \code{min_evalue}, \code{q_length}, \code{unmapped} and \code{mapped} columns
#' @export
#' @import data.table

calculate_blast8_metrics <- function(blast.dt) {

  # Add calculations
  blast.dt[, min_evalue := min(evalue), by = .(query, q_start, q_end)] # minimum e-value for a given partial match
  blast.dt[, q_length := q_end - q_start + 1] # partial match length
  blast.dt[, `:=`(
    unmapped = readlength - q_length,
    mapped = q_length / readlength
  )]
  return(blast.dt)
}

#' Get valid hybrids
#'
#' Calculate valid hybrids from Blat partial alignments
#'
#' @param blast8.query.dt blast8 data.table with read length and metrics
#' @param min_unmapped_length Minimum unmapped read length
#' @param q_minoverlap Maximum overlap in left and right query arms
#' @param q_maxgap Maximum gap between left and right query arms
#' @param s_minoverlap Minimum overlap in left and right subject arms
#' @param xlink_distance Maximum distance from start of the read (i.e. crosslink + 1 nt) before left arm starts
#' @return Hybrids data.table
#' @export
#' @import data.table

get_valid_hybrids <- function(blast.query.dt, min_unmapped_length = 16, q_minoverlap = 4, q_maxgap = 4, s_minoverlap = 0, xlink_distance = 5, max_read_length = 100) {

  # Keep best match for a given query region
  hybrids.dt <- blast.query.dt[evalue == min_evalue]

  # Match up with fasta read length and remove if enough of a continuous match for any hit
  # if(any(hybrids.dt$unmapped < (min_unmapped_length - q_minoverlap) & hybrids.dt$unmapped == 100)) {
  if (any(hybrids.dt$unmapped < (min_unmapped_length - q_minoverlap))) {
    return(data.table())
  } else {

    # Now get combinations
    hybrids.dt[, id := 1:.N]
    hybrids.dt <- merge(hybrids.dt, hybrids.dt, by = c("query"), allow.cartesian = TRUE)
    hybrids.dt <- hybrids.dt[id.y > id.x] # Remove duplicates
    hybrids.dt[, id := paste0(id.x, "_", id.y)]

    # Remove those with significant overlap in the query mappings and too large a gap between the query mappings
    hybrids.dt[, q_ol := min(q_end.x, q_end.y) - max(q_start.x, q_start.y) + 1, by = id]
    hybrids.dt <- hybrids.dt[q_ol <= q_minoverlap][q_ol >= -q_maxgap]

    # Remove those with significant overlap in the subject mappings, if subjects are the same
    hybrids.dt <- hybrids.dt[, s_ol := ifelse(subject.x == subject.y,
      min(s_end.x, s_end.y) - max(s_start.x, s_start.y) + 1,
      0
    ), by = id]
    hybrids.dt <- hybrids.dt[s_ol <= s_minoverlap]

    # Remove those too far away from xlink position
    hybrids.dt <- hybrids.dt[q_start.x < xlink_distance | q_start.y < xlink_distance]

    # Rename columns
    n <- names(hybrids.dt)
    n[grep("\\.x$", n)] <- paste0("L_", gsub("\\.x$", "", n[grep("\\.x$", n)]))
    n[grep("\\.y$", n)] <- paste0("R_", gsub("\\.y$", "", n[grep("\\.y$", n)]))
    n <- gsub("_s_", "_", n)
    n <- gsub("_subject", "_seqnames", n)
    setnames(hybrids.dt, n)
    hybrids.dt[, `:=`(L_strand = "+", R_strand = "+")]

    return(hybrids.dt)
  }
}

#' Find valid overlaps between hybrids
#'
#' @param hybrids.x.dt Hybrids data.table
#' @param hybrids.y.dt Hybrids data.table
#' @return Overlaps data.table
#' @export
#' @import data.table

find_valid_hybrid_overlaps <- function(hybrids.x.dt, hybrids.y.dt) {

  # Make sure uniformly oriented
  hybrids.x.dt <- reorient_hybrids(hybrids.x.dt)
  hybrids.y.dt <- reorient_hybrids(hybrids.y.dt)

  # Convert to GRanges
  L.x.gr <- convert_to_granges(hybrids.x.dt, "L")
  R.x.gr <- convert_to_granges(hybrids.x.dt, "R")

  L.y.gr <- convert_to_granges(hybrids.y.dt, "L")
  R.y.gr <- convert_to_granges(hybrids.y.dt, "R")

  all.seqlevels <- unique(c(hybrids.x.dt$L_seqnames, hybrids.x.dt$R_seqnames, hybrids.y.dt$L_seqnames, hybrids.y.dt$R_seqnames))
  GenomeInfoDb::seqlevels(L.x.gr) <- all.seqlevels
  GenomeInfoDb::seqlevels(R.x.gr) <- all.seqlevels
  GenomeInfoDb::seqlevels(L.y.gr) <- all.seqlevels
  GenomeInfoDb::seqlevels(R.y.gr) <- all.seqlevels

  # Get overlaps
  L.ol <- as.data.table(GenomicRanges::findOverlaps(L.x.gr, L.y.gr))
  R.ol <- as.data.table(GenomicRanges::findOverlaps(R.x.gr, R.y.gr))

  L.ol[, q_s := paste0(queryHits, "_", subjectHits)]
  R.ol[, q_s := paste0(queryHits, "_", subjectHits)]

  # Both arms overlap
  ol <- L.ol[L.ol$q_s %in% R.ol$q_s]
  ol[, q_s := NULL]

  return(ol)

}

#' Filter valid hybrids
#'
#' Filter valid hybrids for each read as
#'  (i) single/unique,
#'  (ii) multiple hits, but one that overlaps a single/unique hit and is selected
#'  (iii) multiple hits that do not overlap single/unique hits and termed ambiguous
#'
#' @param hybrids.dt Hybrids data.table
#' @return Hybrids data.table with \code{hybrid_selection} column
#' @import data.table
#' @export

filter_valid_hybrids <- function(hybrids.dt) {

  hybrids.dt <- reorient_hybrids(hybrids.dt)

  hybrids.dt[, multi := .N, by = query]
  single.hybrids.dt <- hybrids.dt[multi == 1]
  multi.hybrids.dt <- hybrids.dt[multi > 1]

  ol <- find_valid_hybrid_overlaps(multi.hybrids.dt, single.hybrids.dt)
  multi.hybrids.dt <- multi.hybrids.dt[unique(ol$queryHits)]

  # Some multi match more than one single
  multi.hybrids.dt[, q_length_exc_gap := max(L_q_end, R_q_end) - min(L_q_start, R_q_start) + 1 + q_ol, by = .(query, id)]
  multi.hybrids.dt[, total_alignment_length := L_alignment_length + R_alignment_length]

  multi.hybrids.dt <- multi.hybrids.dt[, max_total_alignment_length := max(total_alignment_length), by = query]
  multi.hybrids.dt <- multi.hybrids.dt[total_alignment_length == max_total_alignment_length] # Keep one with with longest total alignment length
  multi.hybrids.dt <- multi.hybrids.dt[, max_q_length_exc_gap := max(q_length_exc_gap), by = query]
  multi.hybrids.dt <- multi.hybrids.dt[q_length_exc_gap == max_q_length_exc_gap] # Keep one with with longest query length minus the gap
  multi.hybrids.dt <- multi.hybrids.dt[, N := .N, by = query][N == 1] # Keep only those with a unique match

  # Write out
  single.hybrids.dt[, hybrid_selection := "single"]
  multi.hybrids.dt[, hybrid_selection := "multi_overlap"]
  valid.hybrids.dt <- rbindlist(list(single.hybrids.dt, multi.hybrids.dt), use.names = TRUE, fill = TRUE)

  other.hybrids.dt <- hybrids.dt[!query %in% valid.hybrids.dt$query]
  other.hybrids.dt[, hybrid_selection := "ambiguous"]

  all.hybrids.dt <- rbindlist(list(valid.hybrids.dt, other.hybrids.dt), use.names = TRUE, fill = TRUE)
  stopifnot(all(hybrids.dt$query %in% all.hybrids.dt$query))

  # Tidy up
  all.hybrids.dt[, `:=` (L_width = L_end - L_start + 1,
                         L_strand = "+",
                         R_width = R_end - R_start + 1,
                         R_strand = "+")]
  all.hybrids.dt[L_seqnames == R_seqnames, orientation := ifelse(L_q_start <= R_q_start, "genomic", "reverse")]
  all.hybrids.dt[, type := ifelse(L_seqnames == R_seqnames, "intragenic", "intergenic")]
  all.hybrids.dt <- all.hybrids.dt[, .(query, orientation, type, hybrid_selection,
                                       L_seqnames, L_q_start, L_q_end, L_start, L_end, L_width, L_strand,
                                       R_seqnames, R_q_start, R_q_end, R_start, R_end, R_width, R_strand)]
  setnames(all.hybrids.dt, c("query", "L_q_start", "L_q_end", "R_q_start", "R_q_end"), c("name", "L_read_start", "L_read_end", "R_read_start", "R_read_end"))

  return(all.hybrids.dt)

}
