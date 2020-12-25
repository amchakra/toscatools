

#' Title
#'
#' @param hybrid.dt
#' @param genome.size
#' @param verbose
#'
#' @return
#' @export
#'
#' @examples

get_contact_map <- function(hybrid.dt, genome.size, verbose = TRUE) {

  hybrid.dt <- primavera::ReorientHybrids(hybrid.dt)

  mat <- matrix(data = 0, nrow = genome.size, ncol = genome.size)

  for(i in 1:nrow(hybrid.dt)) {

    if(verbose) if(i%%1e5 == 0) message(i)
    mat[hybrid.dt[i]$L_start:hybrid.dt[i]$L_end,
        hybrid.dt[i]$R_start:hybrid.dt[i]$R_end] <- mat[hybrid.dt[i]$L_start:hybrid.dt[i]$L_end,
                                                        hybrid.dt[i]$R_start:hybrid.dt[i]$R_end] + 1

  }

  return(mat)

}


#' Title
#'
#' @param mat
#' @param bin.size
#'
#' @return
#' @export
#'
#' @examples

bin_matrix <- function(mat, bin.size) {

  mat.size <- nrow(mat)/bin.size
  bin.mat <- matrix(data = 0, nrow = mat.size, ncol = mat.size)

  for(i in seq_len(mat.size)) {

    for(j in seq_len(mat.size)) {

      bin.mat[i, j] <- sum(mat[((i*bin.size) - bin.size + 1):(i*bin.size), ((j*bin.size) - bin.size + 1):(j*bin.size)])

    }

  }

  return(bin.mat)

}
