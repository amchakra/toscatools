#' Title
#'
#' @param id
#' @param L_sequence
#' @param R_sequence
#'
#' @return
#' @export

.slurm_GetMFE <- function(id, L_sequence, R_sequence) {

  binding <- GetMFE(L_sequence, R_sequence)
  mfe <- binding$mfe
  db <- binding$db
  names(mfe) <- id
  names(db) <- id
  
  #return(list(id = id, mfe = mfe, db = db))
  
  mfe.dt <- as.data.table(list(id = id, mfe = mfe, db = db))
  return(mfe.dt)  # returning as data.table could be more convenient? after rslurm can rbindlist

}

#' Title
#'
#' @param id
#' @param L_sequence
#' @param R_sequence
#'
#' @return
#' @export

.slurm_AnalyseMFE <- function(id, L_sequence, R_sequence) {
  # Get MFE
  binding <- GetMFE(L_sequence, R_sequence)
  mfe <- binding$mfe
  db <- binding$db
  
  L <- sapply(seq_along(id), function(i) ShuffleSequence(L_sequence[i], number = 10, klet = 2))
  R <- sapply(seq_along(id), function(i) ShuffleSequence(R_sequence[i], number = 10, klet = 2))
  
  shuffledmfe <- sapply(1:ncol(L),function(i) GetMFE(L[,i], R[,i])$mfe)  #only need mfe for shuffled
  
  results.ls <- list(id = id, mfe = mfe, db = db,
                     shuffled_mean = colMeans(shuffledmfe),
                     shuffled_sd = apply(shuffledmfe, 2, sd))
  
  # return(list(id = id, mfe = mfe, db = db,
  #             shuffled_mean = colMeans(shuffledmfe),
  #             shuffled_sd = apply(shuffledmfe, 2, sd)))
  mfe.dt <- as.data.table(results.ls)
  return(mfe.dt) # returning as data.table could be more convenient? after rslurm can rbindlist
  

  # # Get MFE 
  # mfe <- GetMFE(L_sequence, R_sequence)
  # 
  # L <- ShuffleSequence(L_sequence, number = 100, klet = 2)
  # R <- ShuffleSequence(R_sequence, number = 100, klet = 2)
  # 
  # shuffledmfe <- sapply(seq_along(L), function(i) GetMFE(L[i], R[i]))
  # 
  # return(list(mfe = mfe,
  #             shuffled_mean = mean(shuffledmfe),
  #             shuffled_sd = sd(shuffledmfe)))

}


