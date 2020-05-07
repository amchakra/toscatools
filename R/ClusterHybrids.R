#
#
# ClusterHybrids <- function(hybrids.dt) {
#
#   L.gr <- ConvertToGRanges(hybrids.dt, arm = "L")
#   R.gr <- ConvertToGRanges(hybrids.dt, arm = "R")
#
#   setorder(hybrids.dt, L_seqnames)
#   L.hybrids.dt <- split(hybrids.dt, by = "L_seqnames")
#
#   # Find self hits
#   L.gr <- ConvertToGRanges(L.hybrids.dt$Mapt_ENSMUSG00000018411.17, arm = "L")
#
#   L.ol <- GenomicRanges::findOverlaps(L.gr, drop.self = TRUE, drop.redundant = TRUE)
#
#   # Keep 75% overlap
#   x <- L.gr[queryHits(L.ol)]
#   y <- L.gr[subjectHits(L.ol)]
#   relative_overlap <- width(pintersect(x, y)) / pmin(width(x), width(y))
#   L.ol <- L.ol[relative_overlap >= 0.75]
#
#   clusters <- extractClustersFromSelfHits(L.ol)
#
# }
#
# # https://support.bioconductor.org/p/85092/
# extractClustersFromSelfHits <- function(hits)
# {
#   stopifnot(is(hits, "Hits"))
#   N <- queryLength(hits)
#   stopifnot(N == subjectLength(hits))
#   h <- union(hits, t(hits))
#   qh <- queryHits(h)
#   sh <- subjectHits(h)
#   cid <- cid0 <- seq_len(N)  # cluster ids
#   while (TRUE) {
#     cid2 <- pmin(cid, selectHits(h, "first"))
#     if (identical(cid2, cid))
#       break
#     cid <- cid2
#     h <- Hits(qh, cid[sh], N, N)
#   }
#   unname(splitAsList(cid0, cid))
# }
#
# x <- c("a", "b", "a", "c", "d")
# table <- c("a", "e", "d", "a", "a", "d")
# hits <- findMatches(x, table)  # sorts the hits by query
# hits
#
# selectHits(hits, select="all")  # no-op
#
# selectHits(hits, select="first")
# selectHits(hits, select="first", nodup=TRUE)
#
# selectHits(hits, select="last")
# selectHits(hits, select="last", nodup=TRUE)
#
# selectHits(hits, select="arbitrary")
# selectHits(hits, select="count")
#
# # ========
#
# L.ol <- GenomicRanges::findOverlaps(L.gr, drop.self = FALSE, drop.redundant = FALSE)
# NEL <- split(subjectHits(L.ol), queryHits(L.ol))
# test <- graph::graphNEL(as.character(seq_len(queryLength(L.ol))), NEL)
#
# comp <- RBGL::connectedComp(test)
# comp.dt <- data.table(cluster = rep(names(comp), elementNROWS(comp)),
#                       id = as.integer(unlist(comp)))
#
# L.gr$id <- 1:length(L.gr)
# L.gr$cluster <- comp.dt$cluster[match(L.gr$id, comp.dt$id)]
#
# comp.dt <- as.data.table(comp)
#
# test@partitioning
