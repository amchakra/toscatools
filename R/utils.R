
#' Title
#'
#' @param blast8
#'
#' @return
#' @export
#'
#' @examples

load_blast8 <- function(blast8) {

  dt <- data.table::fread(input = blast8,
                          col.names = c("query",
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
                                        "bit_score"))
  setkey(dt, query)
  return(dt)

}
