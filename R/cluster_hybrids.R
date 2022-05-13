# ==========
# Functions to cluster hybrids
# ==========

#' Find overlapping pairs of hybrids
#'
#' Uses \code{bedtools pairtopair} to get all overlaps and then calculates percentage overlaps based on overall span of overlapping left and right arms.
#' Called by \code{cluster_hybrids}
#'
#' @param hybrids.dt Hybrids data.table
#' @param verbose Print \code{bedtools} command
#' @return bedpe.dt with left arm, right arm and mean percentage overlap calculated
#' @import data.table
#' @export

find_hybrid_overlaps <- function(hybrids.dt, verbose = FALSE) {
  hybrids.dt <- toscatools::reorient_hybrids(hybrids.dt)

  # Create BEDPE and get overlaps
  bedpe.colnames <- c("L_seqnames", "L_start", "L_end", "R_seqnames", "R_start", "R_end", "name", "total_count", "L_strand", "R_strand")
  bedpe.dt <- hybrids.dt[, ..bedpe.colnames]
  bedpe.dt[, `:=`(
    L_start = L_start - 1,
    R_start = R_start - 1
  )]

  bedpe <- tempfile(tmpdir = getwd(), fileext = ".bedpe")
  ol <- tempfile(tmpdir = getwd(), fileext = ".bedpe")

  fwrite(bedpe.dt, file = bedpe, sep = "\t", col.names = FALSE)
  cmd <- paste("bedtools pairtopair -rdn -a", bedpe, "-b", bedpe, ">", ol)
  if (verbose) message(cmd)
  system(cmd)

  # Check if there are no overlaps
  if (file.size(ol) != 0) {
    bedpe.dt <- fread(ol, col.names = c(paste0(bedpe.colnames, ".x"), paste0(bedpe.colnames, ".y")))
    # Delete temporary files
    invisible(file.remove(bedpe))
    invisible(file.remove(ol))
  } else {

    # Delete temporary files
    invisible(file.remove(bedpe))
    invisible(file.remove(ol))
    return(data.table())
  }

  # Get calculations and filter
  bedpe.dt[, `:=`(
    L_ol = min(L_end.x, L_end.y) - max(L_start.x, L_start.y) + 1,
    R_ol = min(R_end.x, R_end.y) - max(R_start.x, R_start.y) + 1
  ),
  by = .(name.x, name.y)
  ]

  bedpe.dt[, `:=`(
    L_sp = max(L_end.x, L_end.y) - min(L_start.x, L_start.y) + 1,
    R_sp = max(R_end.x, R_end.y) - min(R_start.x, R_start.y) + 1
  ),
  by = .(name.x, name.y)
  ]

  bedpe.dt[, `:=`(
    L_p = L_ol / L_sp,
    R_p = R_ol / R_sp
  ),
  by = .(name.x, name.y)
  ]

  bedpe.dt[, mean_p := mean(c(L_p, R_p)), by = .(name.x, name.y)]
  return(bedpe.dt)
}

#' Cluster overlapping hybrids
#'
#' Cluster overlapping hybrids using graph network analysis
#'
#' @param hybrid.dt Hybrid data.table
#' @param percent_overlap Minimum percentage overlap for both left and right arms
#' @return
#' @import data.table
#' @export

cluster_hybrids <- function(hybrids.dt, percent_overlap = 0.75, verbose = FALSE) {
  hybrids.bedpe.dt <- find_hybrid_overlaps(hybrids.dt, verbose = verbose)

  if (nrow(hybrids.bedpe.dt) == 0) {
    return(hybrids.dt[, cluster := as.character(NA)])
  }

  sel.bedpe.dt <- hybrids.bedpe.dt[L_p > percent_overlap & R_p > percent_overlap]

  # igraph
  g <- igraph::graph_from_edgelist(el = as.matrix(sel.bedpe.dt[, .(name.x, name.y)]), directed = FALSE)
  igraph::E(g)$weight <- sel.bedpe.dt$mean_p # weight by percent overlap

  c <- igraph::components(g)
  if (verbose) message(c$no, " clusters")

  clusters.dt <- data.table(
    name = names(c$membership),
    cluster = c$membership
  )
  setorder(clusters.dt, cluster)

  # Merge back
  if (nrow(clusters.dt) == 0) clusters.dt[, name := character()] # In case there are no clusters
  setkey(clusters.dt, name)
  if (nrow(clusters.dt) != 0) stopifnot(any(!duplicated(clusters.dt$name))) # Make sure no hybrid is in more than one cluster, but only if there are clusters
  clusters.dt[, tempcluster := paste0("C", cluster)][, cluster := NULL]

  # Order cluster names by number of hybrids
  clusters.order.dt <- clusters.dt[, .N, by = tempcluster]
  setorder(clusters.order.dt, -N, tempcluster)[, cluster := paste0("C", stringr::str_pad(1:.N, width = 3, pad = 0))]
  setnames(clusters.order.dt, "N", "cluster_hybrid_count")
  clusters.dt <- merge(clusters.dt, clusters.order.dt, by = "tempcluster")
  clusters.dt[, tempcluster := NULL]

  # Merge back
  setkey(hybrids.dt, name)
  if ("cluster" %in% names(hybrids.dt)) hybrids.dt[, cluster := NULL] # if clusters already assigned, remove them
  hybrids.clustered.dt <- merge(hybrids.dt, clusters.dt, by = "name", all.x = TRUE)
  hybrids.clustered.dt[is.na(cluster), cluster := "."]

  return(hybrids.clustered.dt)
}

# #' Title
# #'
# #' @param hybrids.dt
# #'
# #' @return
# #' @export
# #' @import data.table
# #'
# #' @examples

# collapse_clusters <- function(hybrids.dt) {

#   clusters.dt <- hybrids.dt[!is.na(cluster) & !is.infinite(cluster)][cluster != ""][cluster != "."]
#   clusters.dt[, `:=` (L_cluster_start = floor(median(L_start)),
#                       L_cluster_end = ceiling(median(L_end)),
#                       R_cluster_start = floor(median(R_start)),
#                       R_cluster_end = ceiling(median(R_end))),
#               by = .(L_seqnames, R_seqnames, cluster)]

#   clusters.dt[, count := .N, by = .(L_seqnames, R_seqnames, cluster)]
#   clusters.dt <- unique(clusters.dt[, .(L_seqnames, L_cluster_start, L_cluster_end, L_strand,
#                                         R_seqnames, R_cluster_start, R_cluster_end, R_strand,
#                                         cluster, count)])

#   setnames(clusters.dt,
#            c("cluster", "L_cluster_start", "L_cluster_end", "R_cluster_start", "R_cluster_end"),
#            c("name", "L_start", "L_end", "R_start", "R_end"))

#   return(clusters.dt)

# }

#' Collapse clusters
#'
#' Collapse overlapping clusters of hybrids into a representative duplex
#'
#' @param hybrids.dt Hybrid data.table
#' @return Clusters data.table
#' @import data.table
#' @export

collapse_clusters <- function(hybrids.dt, mode = c("median", "wide")) {
  if(!"cluster" %in% names(hybrids.dt)) stop("No clusters in hybrids.dt")

  clusters.dt <- hybrids.dt[!is.na(cluster) & !is.infinite(cluster)][cluster != ""][cluster != "."]

  if(mode == "median") {
  clusters.dt[, `:=`(
    L_cluster_start = floor(median(L_start)),
    L_cluster_end = ceiling(median(L_end)),
    R_cluster_start = floor(median(R_start)),
    R_cluster_end = ceiling(median(R_end))
  ),
  by = .(L_seqnames, R_seqnames, cluster)
  ]
  } else if(mode == "wide") {
    clusters.dt[, `:=`(
      L_cluster_start = floor(min(L_start)),
      L_cluster_end = ceiling(max(L_end)),
      R_cluster_start = floor(min(R_start)),
      R_cluster_end = ceiling(max(R_end))
    ),
    by = .(L_seqnames, R_seqnames, cluster)
    ]
  }
  clusters.dt[, count := .N, by = .(L_seqnames, R_seqnames, cluster)]
  clusters.dt <- unique(clusters.dt[, .(
    L_seqnames, L_cluster_start, L_cluster_end, L_strand,
    R_seqnames, R_cluster_start, R_cluster_end, R_strand,
    cluster, count
  )])

  setnames(
    clusters.dt,
    c("L_cluster_start", "L_cluster_end", "R_cluster_start", "R_cluster_end"),
    c("L_start", "L_end", "R_start", "R_end")
  )

  clusters.dt[, name := paste0(cluster, "-", L_seqnames, "-", R_seqnames)]

  return(clusters.dt)
}
