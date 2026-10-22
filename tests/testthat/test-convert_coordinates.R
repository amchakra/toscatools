test_that("genomic arms carry the same sequence as their transcript arms on both strands", {
  f <- make_fixture()
  out <- convert_coordinates(data.table::copy(f$hybrids.dt), f$genes.gr)
  expect_equal(nrow(out), nrow(f$hybrids.dt))
  for (i in seq_len(nrow(out))) {
    x <- out[i]
    for (arm in c("L", "R")) {
      g <- genomic_seq(f$chrom, x[[paste0(arm, "_genomic_start")]], x[[paste0(arm, "_genomic_end")]],
                       x[[paste0(arm, "_genomic_strand")]])
      t <- as.character(Biostrings::subseq(f$tx[[x[[paste0(arm, "_seqnames")]]]],
                                           x[[paste0(arm, "_start")]], x[[paste0(arm, "_end")]]))
      expect_identical(g, t, label = paste(x$name, arm, "genomic sequence"))
    }
  }
})

test_that("minus-strand arms map to [end - e + 1, end - s + 1]", {
  f <- make_fixture()
  out <- convert_coordinates(data.table::copy(f$hybrids.dt), f$genes.gr)
  x <- out[name == "h3"]                       # geneM: genomic 501-760, minus strand
  expect_identical(c(x$L_genomic_start, x$L_genomic_end), c(760L - 25L + 1L, 760L - 1L + 1L))
  expect_identical(c(x$R_genomic_start, x$R_genomic_end), c(760L - 130L + 1L, 760L - 100L + 1L))
})

test_that("arm widths are preserved and terminal arms land on the gene boundaries", {
  f <- make_fixture()
  out <- convert_coordinates(data.table::copy(f$hybrids.dt), f$genes.gr)
  expect_identical(out$L_genomic_end - out$L_genomic_start, out$L_end - out$L_start)
  expect_identical(out$R_genomic_end - out$R_genomic_start, out$R_end - out$R_start)
  # geneP (+, 101-300): transcript nt 1 -> 101, nt 200 -> 300
  expect_identical(out[name == "h1", L_genomic_start], 101L)
  expect_identical(out[name == "h2", R_genomic_end], 300L)
  # geneM (-, 501-760): transcript nt 1 -> 760, nt 260 -> 501
  expect_identical(out[name == "h3", L_genomic_end], 760L)
  expect_identical(out[name == "h4", R_genomic_start], 501L)
  # nothing outside its gene
  span <- data.table::as.data.table(f$genes.gr)[, .(fasta_id, gstart = start, gend = end)]
  for (arm in c("L", "R")) {
    m <- merge(out, span, by.x = paste0(arm, "_seqnames"), by.y = "fasta_id")
    expect_true(all(m[[paste0(arm, "_genomic_start")]] >= m$gstart))
    expect_true(all(m[[paste0(arm, "_genomic_end")]] <= m$gend))
  }
})

test_that("genomic columns are integer and carry the gene's seqname and strand", {
  f <- make_fixture()
  out <- convert_coordinates(data.table::copy(f$hybrids.dt), f$genes.gr)
  for (col in c("L_genomic_start", "L_genomic_end", "R_genomic_start", "R_genomic_end"))
    expect_type(out[[col]], "integer")
  expect_identical(out$L_genomic_strand, ifelse(out$L_seqnames == "geneM", "-", "+"))
  expect_identical(out$R_genomic_strand, ifelse(out$R_seqnames == "geneM", "-", "+"))
  expect_true(all(as.character(out$L_genomic_seqnames) == "chrT"))
})

test_that("export_genomic_bed writes 0-based BED12 blocks that match the arms", {
  f <- make_fixture()
  out <- convert_coordinates(data.table::copy(f$hybrids.dt), f$genes.gr)
  intra <- out[L_seqnames == R_seqnames]
  bed <- tempfile(fileext = ".bed")
  on.exit(unlink(bed))
  export_genomic_bed(data.table::copy(intra), bed, sam_tag = FALSE)
  b <- data.table::fread(bed, header = FALSE, colClasses = list(character = 11:12),
                         col.names = c("chrom", "chromStart", "chromEnd", "name", "score", "strand", "thickStart",
                                       "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts"))
  expect_equal(nrow(b), nrow(intra))
  for (i in seq_len(nrow(b))) {
    x <- intra[name == b$name[i]]
    starts <- c(x$L_genomic_start, x$R_genomic_start); ends <- c(x$L_genomic_end, x$R_genomic_end)
    o <- order(starts)
    expect_equal(b$chromStart[i], min(starts) - 1)
    expect_equal(b$chromEnd[i], max(ends))
    expect_equal(as.integer(strsplit(b$blockSizes[i], ",")[[1]]), (ends - starts + 1L)[o])
    expect_equal(as.integer(strsplit(b$blockStarts[i], ",")[[1]]), starts[o] - min(starts))
  }
})
