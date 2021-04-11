
#' Title
#'
#' @param blast.dt
#'
#' @return
#' @export
#' @import data.table

#' @examples

calculate_blast8_metrics <- function(blast.dt) {

  # Add calculations
  blast.dt[, min_evalue := min(evalue), by = .(query, q_start, q_end)] # minimum e-value for a given partial match
  blast.dt[, q_length := q_end - q_start + 1] # partial match length
  blast.dt[, `:=` (unmapped = readlength - q_length,
                   mapped = q_length/readlength)]
  return(blast.dt)

}

#' Title
#'
#' @param blast8.query.dt
#' @param min_unmapped_length
#' @param q_minoverlap
#' @param q_maxgap
#' @param s_minoverlap
#' @param xlink_distance
#'
#' @return
#' @export
#' @import data.table
#'
#' @examples

get_valid_hybrids <- function(blast.query.dt, min_unmapped_length = 16, q_minoverlap = 4, q_maxgap = 4, s_minoverlap = 0, xlink_distance = 5, max_read_length = 100) {

  # Keep best match for a given query region
  hybrids.dt <- blast.query.dt[evalue == min_evalue]

  # Match up with fasta read length and remove if enough of a continuous match for any hit
  # if(any(hybrids.dt$unmapped < (min_unmapped_length - q_minoverlap) & hybrids.dt$unmapped == 100)) {
  if(any(hybrids.dt$unmapped < (min_unmapped_length - q_minoverlap))) {

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
                                              0), by = id]
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
    hybrids.dt[, `:=` (L_strand = "+", R_strand = "+")]

    return(hybrids.dt)

  }

}

#' Title
#'
#' @param hybrids.dt
#'
#' @return
#' @import data.table
#' @export
#'
#' @examples

filter_valid_hybrids <- function(hybrids.dt) {

  hybrids.dt <- reorient_hybrids(hybrids.dt)

  hybrids.dt[, multi := .N, by = query]
  single.hybrids.dt <- hybrids.dt[multi == 1]
  multi.hybrids.dt <- hybrids.dt[multi > 1]

  ol <- find_hybrid_overlaps_simple(multi.hybrids.dt, single.hybrids.dt)
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
