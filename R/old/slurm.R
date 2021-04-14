#' Title
#'
#' @param id
#' @param L_sequence
#' @param R_sequence
#'
#' @return

.slurm_GetMFE <- function(id, L_sequence, R_sequence) {

  mfe <- GetMFE(L_sequence, R_sequence)
  names(mfe) <- id
  return(mfe)

}

#' Title
#'
#' @param id
#' @param L_sequence
#' @param R_sequence
#'
#' @return

.slurm_AnalyseMFE <- function(id, L_sequence, R_sequence) {

  # Get MFE
  mfe <- GetMFE(L_sequence, R_sequence)

  L <- ShuffleSequence(L_sequence, number = 100, klet = 2)
  R <- ShuffleSequence(R_sequence, number = 100, klet = 2)

  shuffledmfe <- sapply(seq_along(L), function(i) GetMFE(L[i], R[i]))

  return(list(mfe = mfe,
              shuffled_mean = mean(shuffledmfe),
              shuffled_sd = sd(shuffledmfe)))

}
