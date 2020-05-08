



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
