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

#' Bin a rectangular matrix
#'
#' This function bins a rectangular matrix by summing values in each bin. The bin size is adjustable.
#' It loops through each bin, defines the indices for the current bin, sums the values in the current bin, and stores in the binned matrix.
#'
#' @param matrix_data The rectangular matrix to bin.
#' @param bin_size The size of the bins
#' @return A binned matrix.
#'
#' @export

bin_matrix_rect <- function(matrix_data, bin_size) {
  nr <- nrow(matrix_data)
  nc <- ncol(matrix_data)

  # Calculate the number of bins in each dimension
  nr_bins <- ceiling(nr / bin_size)
  nc_bins <- ceiling(nc / bin_size)

  # Create an empty matrix to store the binned values
  binned_matrix <- matrix(0, nrow = nr_bins, ncol = nc_bins)

  # Loop through each bin
  for (i in 1:nr_bins) {
    for (j in 1:nc_bins) {
      # Define the indices for the current bin
      start_row <- (i - 1) * bin_size + 1
      end_row <- min(i * bin_size, nr)
      start_col <- (j - 1) * bin_size + 1
      end_col <- min(j * bin_size, nc)

      # Sum the values in the current bin and store in the binned matrix
      binned_matrix[i, j] <- sum(matrix_data[start_row:end_row, start_col:end_col])
    }
  }

  return(binned_matrix)
}
