
#' Title
#'
#' @param blast.dt
#' @param fasta
#'
#' @return
#' @export
#' @import data.table
#'
#' @examples

add_read_lengths <- function(blast.dt, fasta) {

  # Load fasta for read lengths
  fasta <- Biostrings::readDNAStringSet(fasta)
  fasta.dt <- data.table(query = names(fasta),
                         readlength = BiocGenerics::width(fasta))
  setkey(fasta.dt, query)

  # Add to blast.dt
  setorder(blast.dt, query)
  setkey(blast.dt, query)
  stopifnot(all(blast.dt$query %in% fasta.dt$query))
  blast.dt <- fasta.dt[blast.dt]

  return(blast.dt)

}
