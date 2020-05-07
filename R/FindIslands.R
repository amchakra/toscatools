
#' Title
#'
#' @param hybrids.dt
#' @param minimum.coverage
#' @param minimum.width
#' @param cores
#'
#' @return
#' @export

FindIslands <- function(hybrids.dt, minimum.coverage = 2, minimum.width = 8, cores) {

  # Create GRanges of L reads
  L.gr <- ConvertToGRanges(hybrids.dt, arm = "left")
  seqlevels(L.gr) <- unique(as.character(seqnames(L.gr))) # to speed up coverage and slice
  L.coverage <- coverage(L.gr) # calculate coverage per seqname as RleList
  L.island <- slice(L.coverage, minimum.coverage) # apply minimum coverage threshold and returns RleViewsList

  # Create GRanges of L islands
  L.island.gr <- GRanges(seqnames = Rle(rep(names(L.island), sapply(L.island, length))),
                         ranges = IRanges(start = unlist(start(L.island)), end = unlist(end(L.island))),
                         strand = Rle(rep("+", length(unlist(start(L.island))))),
                         peakcoverage = unlist(viewMaxs(L.island)))
  L.island.gr <- L.island.gr[width(L.island.gr) > minimum.width] # apply minimum island width threshold
  L.island.gr$islandid <- paste0(1:length(L.island.gr), rep("+", length(L.island.gr))) # Need to specify strand for later merging ### maybe change order...

  # Now get the matching R confident regions for each L confident region as a list with each item corresponding to one L confident region
  R.gr <- ConvertToGRanges(hybrid.dt, arm = "right")
  seqlevels(R.gr) <- unique(as.character(seqnames(R.gr))) # not sure if this will speed up or not

  L.island.grl <- split(L.island.gr, L.island.gr$islandid)

  cluster <- makeForkCluster(cores)
  R.island.list <- pblapply(cl = cluster, L.island.grl, GetMatchingRightIsland)
  stopCluster(cluster)
  R.island.gr <- unlist(GRangesList(R.island.list)) # get from list to GRanges (need to go via GRangesList)

  # Now need to get back to data.table format. NB - not all L islands match a R island, some L islands are sometimes an R island (and vice versa - 1%)
  names(L.island.gr) <- NULL
  L.island.gr <- AnnotateHybridArm(L.island.gr) # annotate islands
  L.island.dt <- data.table(as.data.frame(L.island.gr))
  L.island.dt[, arm := "L"]
  setkey(L.island.dt, islandid)

  names(R.island.gr) <- NULL
  R.island.gr <- AnnotateHybridArm(R.island.gr)
  R.island.dt <- data.table(as.data.frame(R.island.gr))
  R.island.dt[, arm := "R"]
  setkey(R.island.dt, islandid)

  # Now merge L and R arms to get back to hybrid data.table
  islands.dt <- merge(L.island.dt, R.island.dt)
  setnames(islands.dt,
           c("islandid", "seqnames.x", "start.x", "end.x", "width.x", "strand.x", "peakcoverage.x", "biotype.x", "region.x", "gene.x", "arm.x",
             "seqnames.y", "start.y", "end.y", "width.y", "strand.y", "peakcoverage.y", "biotype.y", "region.y", "gene.y", "arm.y"),
           c("islandid", "L_seqnames", "L_start", "L_end", "L_width", "L_strand", "L_peakcoverage", "L_biotype", "L_region", "L_gene", "L_arm",
             "R_seqnames", "R_start", "R_end", "R_width", "R_strand", "R_peakcoverage", "R_biotype", "R_region", "R_gene", "R_arm"))

  return(islands.dt)

}

GetMatchingRightIsland <- function(L.island.gr) {

  matching.L.gr <- subsetByOverlaps(L.gr, L.island.gr) # find the L reads that overlap the L island
  matching.R.gr <- R.gr[R.gr$read %in% matching.L.gr$read] # select the matching R read
  seqlevels(matching.R.gr) <- unique(as.character(seqnames(matching.R.gr))) # reduce R seqlevels to 1 to speed up splice operation in loop

  R.coverage <- coverage(matching.R.gr) # calculate coverage per seqname as RleList
  R.island <- slice(R.coverage, minimum.coverage) # apply minimum coverage threshold and returns RleViewsList

  if(length(unlist(start(R.island))) > 0) {
    R.island.gr <- GRanges(seqnames = Rle(rep(names(R.island), sapply(R.island, length))),
                           ranges = IRanges(start = unlist(start(R.island)), end = unlist(end(R.island))),
                           strand = Rle(rep("+", length(unlist(start(R.island))))),
                           peakcoverage = unlist(viewMaxs(R.island)),
                           islandid = L.island.gr$islandid)
    R.island.gr <- R.island.gr[width(R.island.gr) > minimum.width] # apply minimum width threshold
    if(length(R.island.gr) > 0) {
      R.island.gr <- R.island.gr[R.island.gr$peakcoverage == max(R.island.gr$peakcoverage)] # if more than one matching R island select the one with the highest peak coverage first
      R.island.gr <- sort(R.island.gr)[1] # then seek the left-most one (minimizes loop...)
    }
  } else {
    R.island.gr <- GRanges()
  }
  seqlevels(R.island.gr) <- seqlevels(L.island.gr) # match seqlevels for later merging
  return(R.island.gr)
}
