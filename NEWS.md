# toscatools 1.0.0.9001

Development version, to be released as 2.0.0. Changes since 1.0.0:

## Bug fixes

* `convert_coordinates()` placed every arm on a minus-strand gene 1 nt too low on the genome. It computed
  `end - L_end` / `end - L_start` (and the same for R) where 1-based closed coordinates need
  `end - L_end + 1` / `end - L_start + 1`. Plus-strand arms were correct. This affects every
  `[L|R]_genomic_start/end` column on a minus-strand gene and everything built from them:
  hybrid and cluster tables, `export_genomic_bed()` BED12 and the BAM made from it, `.bp` arcs, and the
  region/gene annotation of arms that sit at a feature boundary (`annotate_hybrids()` uses arm centres).
  Output from 1.0.0 or earlier should be regenerated, or shifted +1 on the minus strand.
* The genomic coordinate columns are now integer on both strands (they were double).

## Breaking changes

* Removed `convert_linkerbam_to_hybrids()`.
* `get_valid_hybrids()` returns a character discard reason, not an empty `data.table()`, when a read
  yields no valid hybrid. It stops if given more than one query. `max_read_length` (unused) is removed.
  The `xlink_distance` default is now `1e6`, and the crosslink-distance filter is disabled.
* `cluster_hybrids()` gains `cluster_method` (`"components"`, `"leiden"`, `"louvain"`), `weight` and
  `fraction`. The default `fraction = TRUE` uses `bedtools pairtopair -f`
  (`find_hybrid_overlaps_fraction()`), which changes clustering relative to 1.0.0. Pass
  `fraction = FALSE` for the old behaviour. Existing `cluster*` columns are dropped first.
* `collapse_clusters()` validates `mode` with `match.arg()`.

## New features

* `find_hybrid_overlaps_fraction()`, `bin_matrix_rect()`.

## Package

* Added a testthat suite. `convert_coordinates()` and `export_genomic_bed()` are tested on a synthetic
  genome with plus- and minus-strand genes.
* Dependencies are now declared in `Imports` (`Import:` was a typo, so none were installed or checked).
* No compiled code. The Rcpp contact-map kernel added during development (`rcpp_get_contact_map()`) was
  removed before release; `get_contact_map()` remains, in pure R.
