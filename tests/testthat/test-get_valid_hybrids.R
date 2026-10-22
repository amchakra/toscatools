# A single read whose two BLAT hits make one valid intragenic hybrid. `offset` moves both query
# alignments along the read, so the left arm starts at read position 1 + offset.
make_query <- function(offset = 0L) {
  data.table::data.table(
    query = "read1",
    subject = "geneA",
    q_start = c(1L, 31L) + offset, q_end = c(30L, 60L) + offset,
    s_start = c(100L, 300L), s_end = c(129L, 329L),
    evalue = 1e-5, min_evalue = 1e-5,
    unmapped = 30L
  )
}

test_that("the crosslink-distance filter is off by default", {
  expect_s3_class(get_valid_hybrids(make_query(10L)), "data.table")
  expect_s3_class(get_valid_hybrids(make_query(10L), xlink_distance = NA), "data.table")
  expect_s3_class(get_valid_hybrids(make_query(10L), xlink_distance = Inf), "data.table")
})

test_that("with the filter on, a hybrid with an arm at the read start is kept", {
  h <- get_valid_hybrids(make_query(0L), xlink_distance = 5)
  expect_s3_class(h, "data.table")
  expect_equal(nrow(h), 1L)
  expect_equal(h$L_q_start, 1L)
  # read position 4 is still within the window, position 5 is not
  expect_s3_class(get_valid_hybrids(make_query(3L), xlink_distance = 5), "data.table")
})

test_that("with the filter on, a hybrid with both arms away from the read start is discarded", {
  expect_identical(get_valid_hybrids(make_query(4L), xlink_distance = 5), "hybrids_too_far_from_xlink_site")
  expect_identical(get_valid_hybrids(make_query(10L), xlink_distance = 5), "hybrids_too_far_from_xlink_site")
})
