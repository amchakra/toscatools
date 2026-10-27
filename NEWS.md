# toscatools 1.1.0.9002

## New features

* `annotate_hybrids()` labels more non-genomic reference entries (entries with no genome coordinates) in
  the same way as rRNA and tRNA: `gene_id` and `gene_name` are the entry name, and `region` and
  `biotype` are its class. `^snRNA_` -> `snRNA`, `^7SK$` -> `7SK`, `^7SL$` -> `7SL`, `^YRNA_` -> `Y_RNA` and `^chrM$` (the
  mitochondrial genome as one entry) -> `Mt`. Before, these entries were not special-cased, so their
  arms were usually called `intergenic` with gene `None`. rRNA and tRNA labelling is unchanged.

# toscatools 1.1.0.9001

Development version, to be released as 2.0.0. Changes since 1.0.0 (the `dev` branch had been versioned
1.1.0 before these changes without a release; the development version is numbered above it so that
version comparisons order correctly):

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
* `get_valid_hybrids()`'s crosslink-distance filter is now **off by default** (`xlink_distance = NULL`;
  it was 5 in 1.0.0). Set `xlink_distance = 5` for hiCLIP and other iCLIP-style libraries, whose read 5'
  end is the crosslink: a pair is kept only if one arm starts at read position 1-4. Leave it off for
  SPLASH, PARIS or COMRADES. Pairs it removes are reported as `"hybrids_too_far_from_xlink_site"`.
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
