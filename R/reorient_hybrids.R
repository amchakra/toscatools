#' Title
#'
#' @param hybrids.dt
#'
#' @return
#'
#' @import data.table
#' @export
#'

reorient_hybrids <- function(hybrids.dt) {

  # First do starts
  correct.dt <- hybrids.dt[L_start <= R_start]
  incorrect.dt <- hybrids.dt[L_start > R_start]

  renamed <- gsub("^L_", "X_", names(incorrect.dt))
  renamed <- gsub("^R_", "L_", renamed)
  renamed <- gsub("^X_", "R_", renamed)

  setnames(incorrect.dt, renamed)

  reoriented.dt <- rbindlist(list(correct.dt, incorrect.dt), use.names = TRUE)

  stopifnot(all(reoriented.dt$L_start <= reoriented.dt$R_start))

  # Then do subject (to make sure intergenics in same order)
  correct.dt <- reoriented.dt[L_seqnames <= R_seqnames]
  incorrect.dt <- reoriented.dt[L_seqnames > R_seqnames]

  renamed <- gsub("^L_", "X_", names(incorrect.dt))
  renamed <- gsub("^R_", "L_", renamed)
  renamed <- gsub("^X_", "R_", renamed)

  setnames(incorrect.dt, renamed)

  reoriented.dt <- rbindlist(list(correct.dt, incorrect.dt), use.names = TRUE)
  stopifnot(all(reoriented.dt$L_subject <= reoriented.dt$R_subject))
  stopifnot(nrow(reoriented.dt) == nrow(hybrids.dt))

  return(reoriented.dt)

}
