# Toy hybrids with genomic columns already computed. Non-genomic reference entries have no genome
# coordinates: their "genomic" seqname is the entry itself and coordinates are entry-relative.
make_annotation_fixture <- function() {
  entries <- c("geneA", "rRNA_45S", "rDNA", "rRNA_5S", "tRNA-Gly-GCC",
               "snRNA_U1", "snRNA_U4atac", "7SK", "7SL", "YRNA_Y1", "chrM")
  n <- length(entries)
  genomic.seqnames <- ifelse(entries == "geneA", "chrT", entries)
  genomic.start <- ifelse(entries == "geneA", 120L, 10L)
  hybrids.dt <- data.table::data.table(
    name = paste0("h", seq_len(n)),
    L_seqnames = entries,
    L_start = 10L, L_end = 30L,
    R_seqnames = rev(entries),
    R_start = 10L, R_end = 30L,
    L_genomic_seqnames = genomic.seqnames,
    L_genomic_start = genomic.start,
    L_genomic_end = genomic.start + 20L,
    L_genomic_strand = "+",
    R_genomic_seqnames = rev(genomic.seqnames),
    R_genomic_start = rev(genomic.start),
    R_genomic_end = rev(genomic.start) + 20L,
    R_genomic_strand = "+"
  )
  # One gene on chrT: a 5'UTR then CDS
  regions.gr <- GenomicRanges::GRanges(
    seqnames = "chrT",
    ranges = IRanges::IRanges(start = c(101L, 151L), end = c(150L, 300L)),
    strand = "+",
    type = c("UTR5", "CDS"),
    gene_id = "ENSG_A",
    gene_name = "geneA",
    biotype = "protein_coding"
  )
  list(hybrids.dt = hybrids.dt, regions.gr = regions.gr)
}

expected_class <- c(
  "rRNA_45S" = "rRNA", "rDNA" = "rRNA", "rRNA_5S" = "rRNA", "tRNA-Gly-GCC" = "tRNA",
  "snRNA_U1" = "snRNA", "snRNA_U4atac" = "snRNA", "7SK" = "7SK", "7SL" = "7SL",
  "YRNA_Y1" = "Y_RNA", "chrM" = "Mt"
)

test_that("non-genomic entries are labelled by class with the entry name as gene", {
  f <- make_annotation_fixture()
  out <- annotate_hybrids(data.table::copy(f$hybrids.dt), f$regions.gr)
  expect_equal(nrow(out), nrow(f$hybrids.dt))
  for (arm in c("L", "R")) {
    seqnames <- out[[paste0(arm, "_seqnames")]]
    idx <- match(names(expected_class), seqnames)
    expect_identical(out[[paste0(arm, "_region")]][idx], unname(expected_class), label = paste(arm, "region"))
    expect_identical(out[[paste0(arm, "_biotype")]][idx], unname(expected_class), label = paste(arm, "biotype"))
    expect_identical(out[[paste0(arm, "_gene_name")]][idx], names(expected_class), label = paste(arm, "gene_name"))
    expect_identical(out[[paste0(arm, "_gene_id")]][idx], names(expected_class), label = paste(arm, "gene_id"))
  }
})

test_that("a gene arm keeps its regions annotation", {
  f <- make_annotation_fixture()
  out <- annotate_hybrids(data.table::copy(f$hybrids.dt), f$regions.gr)
  # L arm of h1 is geneA at chrT:120-140 (centre in UTR5); R arm of h11 is the same
  x <- out[L_seqnames == "geneA"]
  expect_identical(c(x$L_region, x$L_gene_id, x$L_gene_name, x$L_biotype),
                   c("UTR5", "ENSG_A", "geneA", "protein_coding"))
  y <- out[R_seqnames == "geneA"]
  expect_identical(c(y$R_region, y$R_gene_id, y$R_gene_name, y$R_biotype),
                   c("UTR5", "ENSG_A", "geneA", "protein_coding"))
})

test_that("no arm is left unannotated and no non-genomic entry is called intergenic", {
  f <- make_annotation_fixture()
  out <- annotate_hybrids(data.table::copy(f$hybrids.dt), f$regions.gr)
  expect_false(anyNA(c(out$L_region, out$R_region)))
  expect_false(any(c(out$L_region, out$R_region) == "intergenic"))
  expect_false(any(c(out$L_gene_name, out$R_gene_name) == "None"))
})

test_that("patterns are anchored: names that merely contain a class word are not relabelled", {
  f <- make_annotation_fixture()
  dt <- f$hybrids.dt[1]
  dt[, `:=`(L_seqnames = "geneA_7SKlike", R_seqnames = "chrM_like")]
  dt[, `:=`(R_genomic_seqnames = "chrT", R_genomic_start = 400L, R_genomic_end = 420L)]
  out <- annotate_hybrids(dt, f$regions.gr)
  expect_identical(out$L_region, "UTR5")
  expect_identical(out$R_region, "intergenic")
})
