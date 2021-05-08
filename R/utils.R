

#' Get sequence from NCBI
#'
#' @param accession NCBI accession
#'
#' @return DNAStringSet of sequence
#' @export

get_ncbi_sequence <- function(accession) {

  fa <- rentrez::entrez_fetch(db="nucleotide", id=accession, rettype="fasta")
  lines <- stringr::str_split(fa, "\n")[[1]]
  id <- gsub(">", "", lines[1])
  seq <- paste0(lines[-1], collapse = "")

  # Convert to DNAStringSet
  dss <- Biostrings::DNAStringSet(seq)
  names(dss) <- id

  return(dss)

}
