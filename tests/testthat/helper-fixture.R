# Synthetic genome + Tosca-style reference for coordinate tests. Each reference sequence is the unspliced
# genomic span of a gene in the gene's sense (as tosca/ref/*_annotation.R builds it with getfasta -s).
make_fixture <- function() {
  set.seed(1)
  chrom <- Biostrings::DNAString(paste(sample(c("A", "C", "G", "T"), 1000, replace = TRUE), collapse = ""))
  genes.gr <- GenomicRanges::GRanges(
    seqnames = "chrT",
    ranges = IRanges::IRanges(start = c(101L, 501L), end = c(300L, 760L)),
    strand = c("+", "-"),
    fasta_id = c("geneP", "geneM")
  )
  tx <- Biostrings::DNAStringSet(list(
    geneP = Biostrings::subseq(chrom, 101L, 300L),
    geneM = Biostrings::reverseComplement(Biostrings::subseq(chrom, 501L, 760L))
  ))
  # Arms include the first and last nt of each reference (geneP is 200 nt, geneM 260 nt).
  hybrids.dt <- data.table::data.table(
    name       = paste0("h", 1:5),
    L_seqnames = c("geneP", "geneP", "geneM", "geneM", "geneP"),
    L_start    = c(1L, 150L, 1L, 200L, 40L),
    L_end      = c(20L, 180L, 25L, 230L, 60L),
    R_seqnames = c("geneP", "geneP", "geneM", "geneM", "geneM"),
    R_start    = c(50L, 181L, 100L, 241L, 10L),
    R_end      = c(80L, 200L, 130L, 260L, 35L)
  )
  list(chrom = chrom, genes.gr = genes.gr, tx = tx, hybrids.dt = hybrids.dt)
}

# Sequence at genomic coordinates, read in the gene's sense.
genomic_seq <- function(chrom, start, end, strand) {
  s <- Biostrings::subseq(chrom, start, end)
  as.character(if (strand == "-") Biostrings::reverseComplement(s) else s)
}
