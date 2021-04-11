

#' Title
#'
#' @param L_sequence
#' @param R_sequence
#'
#' @return
#'
#' @import data.table
#' @export
#'
#' @examples

get_structure <- function(L_sequence, R_sequence) {

  input <- paste0(L_sequence, "\n", R_sequence)
  rnaduplex <- system("RNAduplex --noLP", input = input, intern = TRUE)

  # Get MFE
  rnaduplex <- gsub("\\s+", "_", rnaduplex)
  mfe <- sapply(strsplit(rnaduplex, "_"), "[", 5)
  if(mfe == "(") mfe <- sapply(strsplit(rnaduplex, "_"), "[", 6) # Positives < 10 have an extra space
  mfe <- as.numeric(gsub("\\(|\\)", "", mfe))

  # Get dot.bracket
  db <- sapply(strsplit(rnaduplex, "_"), "[", 1)

  # Get from/to positions
  l <- sapply(strsplit(rnaduplex, "_"), "[", 2)
  r <- sapply(strsplit(rnaduplex, "_"), "[", 4)
  l_from <- as.numeric(sapply(strsplit(l, ","), "[", 1))
  l_to <- as.numeric(sapply(strsplit(l, ","), "[", 2))
  r_from <- as.numeric(sapply(strsplit(r, ","), "[", 1))
  r_to <- as.numeric(sapply(strsplit(r, ","), "[", 2))
  l_width <- nchar(L_sequence)
  r_width <- nchar(R_sequence)

  # Pad db with . according to from/to
  l_db <- sapply(strsplit(db, "&"), "[", 1)
  r_db <- sapply(strsplit(db, "&"), "[", 2)
  l_db <- unlist(strsplit(l_db, "")) # convert string to vector of characters
  r_db <- unlist(strsplit(r_db, ""))
  l_db <- c(rep(".", l_from - 1L), l_db, rep(".", l_width - l_to))
  r_db <- c(rep(".", r_from - 1L), r_db, rep(".", r_width - r_to))

  stopifnot(all(length(l_db == l_width), length(r_db) == r_width))
  db <- paste0(c(l_db, "&", r_db), collapse = "")

  return(data.table(mfe = mfe, structure = db))

}
