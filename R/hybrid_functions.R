

#' Title
#'
#' @param hybrids.grl
#'
#' @return
#'
#' @import data.table
#' @export

ConvertToDataTable <- function(hybrids.grl) {

  L.dt <- as.data.table(hybrids.grl$L)
  R.dt <- as.data.table(hybrids.grl$R)
  setkey(L.dt, name)
  setkey(R.dt, name)

  hybrids.dt <- merge(L.dt, R.dt, by = "name")

  # Rename L
  L.names <- grep(".x", names(hybrids.dt), value = TRUE)
  L.renamed <- paste0("L_", gsub(".x", "", L.names))
  setnames(hybrids.dt, L.names, L.renamed)

  # Rename R
  R.names <- grep(".y", names(hybrids.dt), value = TRUE)
  R.renamed <- paste0("R_", gsub(".y", "", R.names))
  setnames(hybrids.dt, R.names, R.renamed)

  return(hybrids.dt)

}


#' Title
#'
#' @param hybrids.dt
#' @param arm
#'
#' @return
#' @import data.table
#' @export


ConvertToGRanges <- function(hybrids.dt, arm = c("L", "R")) {

  if(arm == "L") {

    L.dt <- hybrids.dt[, grep("^R_", names(hybrids.dt), invert = TRUE), with = FALSE]
    setnames(L.dt, names(L.dt), gsub("^L_", "", names(L.dt)))
    L.gr <- GenomicRanges::GRanges(L.dt)

    return(L.gr)

  } else if(arm == "R") {

    R.dt <- hybrids.dt[, grep("^L_", names(hybrids.dt), invert = TRUE), with = FALSE]
    setnames(R.dt, names(R.dt), gsub("^R_", "", names(L.dt)))
    R.gr <- GenomicRanges::GRanges(R.dt)

    return(R.gr)

  }

}
