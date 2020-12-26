



#' FindHybridOverlaps
#'
#' @param hybrids.x.dt
#' @param hybrids.y.dt
#' @param percent_overlap
#' @param reciprocal
#'
#' @return
#'
#' @import data.table
#' @export

FindHybridOverlaps <- function(hybrids.x.dt, hybrids.y.dt, percent_overlap = 0.5, reciprocal = FALSE) {

  hybrids.x.dt[, `:=` (m = 1, idx = 1:.N)]
  hybrids.y.dt[, `:=` (m = 1, idx = 1:.N)]

  cj.hybrids.dt <- merge(hybrids.x.dt, hybrids.y.dt, by = "m", allow.cartesian = TRUE)

  if(reciprocal == FALSE) cj.hybrids.dt <- cj.hybrids.dt[name.x < name.y] # Don't do this for edges
  cj.hybrids.dt <- cj.hybrids.dt[name.x != name.y]

  # Calculate overlaps
  cj.hybrids.dt[, `:=` (L_ol = min(L_end.x, L_end.y) - max(L_start.x, L_start.y) + 1,
                        R_ol = min(R_end.x, R_end.y) - max(R_start.x, R_start.y) + 1),
               by = .(name.x, name.y)]
  cj.hybrids.dt <- cj.hybrids.dt[L_ol > 0 & R_ol > 0]

  # If no overlaps, return empty data.table
  if(nrow(cj.hybrids.dt) == 0) return(cj.hybrids.dt[, `:=` (L_p = NA, R_p = NA)][, .(idx.x, name.x, idx.y, name.y, L_p, R_p)])

  # Calculate span
  cj.hybrids.dt[, `:=` (L_sp = max(L_end.x, L_end.y) - min(L_start.x, L_start.y) + 1,
                        R_sp = max(R_end.x, R_end.y) - min(R_start.x, R_start.y) + 1),
                by = .(name.x, name.y)]

  # Calculate percentage overlaps
  cj.hybrids.dt[, `:=` (L_p = L_ol/L_sp,
                        R_p = R_ol/R_sp),
                by = .(name.x, name.y)]

  # Filter based on percentage overlap
  cj.hybrids.dt <- cj.hybrids.dt[L_p > percent_overlap & R_p > percent_overlap]

  hits.dt <- cj.hybrids.dt[, .(idx.x, name.x, idx.y, name.y, L_p, R_p)]
  return(hits.dt)

}

#' Title
#'
#' @param hybrids.dt
#' @param percent_overlap
#'
#' @return
#'
#' @import data.table
#' @export

ClusterHybrids <- function(hybrids.dt, percent_overlap = 0.5) {

  # Get overlap hits
  hits.dt <- FindHybridOverlaps(hybrids.x.dt = hybrids.dt, hybrids.y.dt = hybrids.dt, percent_overlap = percent_overlap, reciprocal = TRUE)

  # If not hits, return all clusters as NA
  if(nrow(hits.dt) == 0) return(hybrids.dt[, cluster := NA])

  # Convert to graph
  n <- as.character(unique(hits.dt$name.x))
  e <- split(as.character(hits.dt$name.y), as.character(hits.dt$name.x))
  NEL <- graph::graphNEL(nodes = n, edgeL = e, edgemode = "undirected")

  # Get connected nodes (i.e. clusters)
  clusters <- RBGL::connectedComp(NEL) # https://support.bioconductor.org/p/85092/
  # table(elementNROWS(clusters))
  clusters.dt <- data.table(cluster = rep(names(clusters), S4Vectors::elementNROWS(clusters)),
                            name = unlist(clusters))
  setkey(clusters.dt, name)

  stopifnot(any(!duplicated(clusters.dt$name))) # Make sure no hybrid is in more than one cluster

  # Merge back
  setkey(hybrids.dt, name)
  hybrids.dt <- merge(hybrids.dt, clusters.dt, by = "name", all.x = TRUE)

  return(hybrids.dt)

}

#' Title
#'
#' @param hybrids.dt
#'
#' @return
#' @import data.table
#' @export

CollapseClustersOld <- function(hybrids.dt) {

  clusters.dt <- hybrids.dt[!is.na(cluster) & !is.infinite(cluster)]
  clusters.dt[, `:=` (L_cluster_start = median(L_start),
                      L_cluster_end = median(L_end),
                      R_cluster_start = median(R_start),
                      R_cluster_end = median(R_end),
                      cluster_mfe = mean(mfe)),
              by = .(L_seqnames, cluster)]

  clusters.dt <- unique(clusters.dt[, .(L_seqnames, L_cluster_start, L_cluster_end, L_strand,
                                        R_seqnames, R_cluster_start, R_cluster_end, R_strand,
                                        cluster, cluster_mfe)])

  setnames(clusters.dt,
           c("L_cluster_start", "L_cluster_end", "R_cluster_start", "R_cluster_end", "cluster_mfe"),
           c("L_start", "L_end", "R_start", "R_end", "mfe"))

  clusters.dt[, name := paste0("C", 1:.N)]

  return(clusters.dt)

}

#' Title
#'
#' @param hybrids.dt
#'
#' @return
#' @import data.table
#' @export

CollapseClusters <- function(hybrids.dt) {

  clusters.dt <- hybrids.dt[!is.na(cluster) & !is.infinite(cluster)]
  clusters.dt[, `:=` (L_cluster_start = floor(median(L_start)),
                      L_cluster_end = ceiling(median(L_end)),
                      R_cluster_start = floor(median(R_start)),
                      R_cluster_end = ceiling(median(R_end))),
              by = .(L_seqnames, cluster)]

  clusters.dt[, count := .N, by = .(L_seqnames, cluster)]
  clusters.dt <- unique(clusters.dt[, .(L_seqnames, L_cluster_start, L_cluster_end, L_strand,
                                        R_seqnames, R_cluster_start, R_cluster_end, R_strand,
                                        cluster, count)])

  setnames(clusters.dt,
           c("L_cluster_start", "L_cluster_end", "R_cluster_start", "R_cluster_end"),
           c("L_start", "L_end", "R_start", "R_end"))

  clusters.dt[, name := paste0("C", 1:.N)]

  return(clusters.dt)

}

#' Title
#'
#' @param hybrid.dt
#' @param sample_size
#' @param percent_overlap
#' @param cores
#'
#' @return
#' @import data.table
#' @export
#'
#' @examples

cluster_hybrids <- function(hybrids.dt, sample_size = -1, percent_overlap = 0.75, cores = 8, verbose = FALSE) {

  hybrids.dt <- primavera::reorient_hybrids(hybrids.dt)

  # Get arm overlaps
  L.gr <- primavera::convert_to_granges(hybrids.dt, arm = "L")
  R.gr <- primavera::convert_to_granges(hybrids.dt, arm = "R")

  L.bed <- tempfile(tmpdir = getwd(), fileext = ".bed")
  R.bed <- tempfile(tmpdir = getwd(), fileext = ".bed")

  L.ol <- tempfile(tmpdir = getwd(), fileext = ".bed")
  R.ol <- tempfile(tmpdir = getwd(), fileext = ".bed")

  if(sample_size == -1) {
    subsample <- 1:length(L.gr)
  } else {
    subsample <- sample(1:length(L.gr), sample_size)
  }

  rtracklayer::export.bed(sort(L.gr[subsample]), L.bed)
  cmd <- paste("bedtools intersect -sorted -wo -s -f", percent_overlap, "-F", percent_overlap, "-e", "-a", L.bed, "-b", L.bed, ">", L.ol)
  if(verbose) message(cmd)
  system(cmd)

  rtracklayer::export.bed(sort(R.gr[subsample]), R.bed)
  cmd <- paste("bedtools intersect -sorted -wo -s -f", percent_overlap, "-F", percent_overlap, "-e", "-a", R.bed, "-b", R.bed, ">", R.ol)
  if(verbose) message(cmd)
  system(cmd)

  # Get overlap graph
  L.dt <- fread(L.ol, select = c(4, 10), col.names = c("q", "s"))
  R.dt <- fread(R.ol, select = c(4, 10), col.names = c("q", "s"))

  # Delete temporary files
  invisible(file.remove(L.bed))
  invisible(file.remove(R.bed))
  invisible(file.remove(L.ol))
  invisible(file.remove(R.ol))

  # Remove self matches
  L.dt <- L.dt[q != s]
  R.dt <- R.dt[q != s]

  # Create list for graph
  # s <- unique(c(L.dt$s, R.dt$s)) # All that have at least one arm overlapping
  s <- intersect(L.dt$s, R.dt$s) # All that have both arms overlapping
  s.list <- parallel::mclapply(s, function(x) {

    return(c(L.dt[s == x]$q, R.dt[s == x]$q)[duplicated(c(L.dt[s == x]$q, R.dt[s == x]$q))])

  }, mc.cores = cores)
  names(s.list) <- s

  NEL <- graph::graphNEL(nodes = as.character(s), edgeL = s.list, edgemode = "undirected")

  # Get connected nodes (i.e. clusters)
  clusters <- RBGL::connectedComp(NEL) # https://support.bioconductor.org/p/85092/
  clusters.dt <- data.table(cluster = rep(names(clusters), S4Vectors::elementNROWS(clusters)),
                            name = unlist(clusters))
  setkey(clusters.dt, name)
  stopifnot(any(!duplicated(clusters.dt$name))) # Make sure no hybrid is in more than one cluster
  clusters.dt[, cluster := paste0("C", cluster)]

  # Merge back
  setkey(hybrids.dt, name)
  hybrids.dt <- merge(hybrids.dt, clusters.dt, by = "name", all.x = TRUE)

  return(hybrids.dt)

}