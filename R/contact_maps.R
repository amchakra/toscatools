# ==========
# Functions for contact maps
# ==========

#' Get contact map
#'
#' Gets contact map matrix from hybrids.dt (previously subset to given gene/region)
#'
#' @param hybrid.dt Hybrids data.table
#' @param genome.size Size of genome/gene/region
#' @param verbose Monitor progress
#' @return Contact map matrix
#' @export

get_contact_map <- function(hybrid.dt, genome.size, verbose = TRUE) {
  hybrid.dt <- toscatools::reorient_hybrids(hybrid.dt)

  mat <- matrix(data = 0, nrow = genome.size, ncol = genome.size)

  for (i in 1:nrow(hybrid.dt)) {
    if (verbose) if (i %% 1e5 == 0) message(i)
    mat[
      hybrid.dt[i]$L_start:hybrid.dt[i]$L_end,
      hybrid.dt[i]$R_start:hybrid.dt[i]$R_end
    ] <- mat[
      hybrid.dt[i]$L_start:hybrid.dt[i]$L_end,
      hybrid.dt[i]$R_start:hybrid.dt[i]$R_end
    ] + 1
  }

  return(mat)
}

#' Bin contact map matrix
#'
#' @param mat Contact map matrix
#' @param bin.size Size of bin
#' @return Binnecd contact map matrix
#' @export

bin_matrix <- function(mat, bin.size) {
  mat.size <- nrow(mat) / bin.size
  bin.mat <- matrix(data = 0, nrow = mat.size, ncol = mat.size)

  for (i in seq_len(mat.size)) {
    for (j in seq_len(mat.size)) {
      bin.mat[i, j] <- sum(mat[((i * bin.size) - bin.size + 1):(i * bin.size), ((j * bin.size) - bin.size + 1):(j * bin.size)])
    }
  }

  return(bin.mat)
}
