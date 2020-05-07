#' Title
#'
#' @param id
#' @param L_sequence
#' @param R_sequence
#'
#' @return
#' @export

.slurm_GetMFE <- function(id, L_sequence, R_sequence) {

  mfe <- GetMFE(L_sequence, R_sequence)
  names(mfe) <- id
  return(mfe)

}
