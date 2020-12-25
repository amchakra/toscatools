
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

get_valid_hybrids <- function(blast.query.dt, min_unmapped_length = 21, q_minoverlap = 4, q_maxgap = 4, s_minoverlap = 0, xlink_distance = 5, max_read_length = 100) {

  # Keep best match for a given query region
  blast.query.dt[, min_evalue := min(evalue), by = .(q_start, q_end)]
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
