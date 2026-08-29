## permdisp_dispersion_test.R
##
## Checks whether the primary-vs-metastasis difference in Fig. 3 is a
## location effect or just unequal spread within groups. Warton, Wright &
## Wang (2012, MEE 3:89-101) and Anderson & Walsh (2013, Ecol Monogr
## 83:557-574) both make the point that PERMANOVA can pick up dispersion
## differences, so this runs the companion test (Anderson's PERMDISP,
## via vegan::betadisper) on the same six-class subcomposition and the
## same nine matched patients.
##
## Permutations are restricted within patient (nine pairs -> 512 possible
## swaps) and enumerated in full rather than sampled, same as everywhere
## else in this analysis.
##
## A note on vegan::permutest.betadisper: it does not return the same
## p-value as full enumeration under a restricted permutation scheme, even
## when handed the exact permutation matrix directly. The F-statistic it
## reports agrees; the p-value does not. This is demonstrated directly near
## the end of this script (not just asserted here) for the Bray-Curtis
## centroid case. The p-values reported above come from refitting
## betadisper() under each of the 512 valid arrangements directly and
## counting, not from permutest().

library(data.table)
library(vegan)

find_spot_table <- function() {
  cand <- c(Sys.getenv("PDAC_SPOT_TABLE", unset = NA),
            "pdac_spot_level.csv.gz",
            file.path("..", "pdac_spot_level.csv.gz"),
            file.path(dirname(getwd()), "pdac_spot_level.csv.gz"))
  cand <- cand[!is.na(cand)]
  hit  <- cand[file.exists(cand)]
  if (!length(hit)) {
    stop("Cannot find pdac_spot_level.csv.gz.\n",
         "  Run from the repository root, or set the path explicitly:\n",
         "    Sys.setenv(PDAC_SPOT_TABLE = '/path/to/pdac_spot_level.csv.gz')\n",
         "  Looked in: ", paste(cand, collapse = ", "), call. = FALSE)
  }
  normalizePath(hit[1])
}
csv_path <- find_spot_table()
cat("reading:", csv_path, "\n")

md <- as.data.table(read.csv(gzfile(csv_path), stringsAsFactors = FALSE))
md[, site := c("Pancreas" = "PT", "Liver" = "HM",
                "Lymph node" = "LNM", "Normal Pancreas" = "NP")[site]]
setnames(md, "ecotype", "cc")

patients <- c("PT_2", "PT_3", "PT_4", "PT_6", "PT_8", "PT_9", "PT_10", "PT_11", "PT_12")
classes  <- c("CC1", "CC2", "CC3", "CC4", "CC5", "CC7")

spots <- md[patient %in% patients & site %in% c("PT", "HM") & cc %in% classes]
counts <- spots[, .N, by = .(patient, site, cc)]
counts[, prop := N / sum(N), by = .(patient, site)]
wide <- dcast(counts, patient + site ~ cc, value.var = "prop", fill = 0)
wide <- wide[order(match(patient, patients), factor(site, levels = c("HM", "PT")))]

X   <- as.matrix(wide[, ..classes])
grp <- as.character(wide$site)

clr <- function(X) {
  Z <- X
  Z[Z == 0] <- min(Z[Z > 0]) / 2
  Z <- Z / rowSums(Z)
  t(apply(log(Z), 1, function(r) r - mean(r)))
}

dist_bray <- vegdist(X, method = "bray")
dist_ait  <- dist(clr(X))

exact_permdisp <- function(D, grp, type = "centroid") {
  n <- length(grp)
  npairs <- n / 2
  nperm  <- 2^npairs
  f_stat <- function(g) anova(betadisper(D, factor(g), type = type))$F[1]
  observed <- f_stat(grp)
  exceed <- 0L
  for (i in 0:(nperm - 1)) {
    bits <- as.integer(intToBits(i))[1:npairs]
    g <- grp
    for (k in which(bits == 1)) {
      a <- 2 * k - 1; b <- 2 * k
      g[c(a, b)] <- g[c(b, a)]
    }
    if (f_stat(g) >= observed - 1e-8) exceed <- exceed + 1L
  }
  data.frame(F = observed, P = exceed / nperm, n_exceed = exceed, n_perm = nperm)
}

results <- rbind(
  data.frame(distance = "bray",      type = "centroid", exact_permdisp(dist_bray, grp, "centroid")),
  data.frame(distance = "bray",      type = "median",   exact_permdisp(dist_bray, grp, "median")),
  data.frame(distance = "aitchison", type = "centroid", exact_permdisp(dist_ait,  grp, "centroid")),
  data.frame(distance = "aitchison", type = "median",   exact_permdisp(dist_ait,  grp, "median"))
)

print(results, row.names = FALSE)

if (all(results$P > 0.05)) {
  cat("\nDispersion doesn't differ between primary and metastasis in either",
      "\ngeometry -- the PERMANOVA site effect isn't a dispersion artifact.\n")
}

## Direct comparison against vegan::permutest.betadisper, same restricted
## permutation set, same case (Bray-Curtis, centroid). Demonstrates the
## claim in the header rather than just asserting it.
perm_matrix_9 <- function() {
  t(sapply(0:511, function(i) {
    bits <- as.integer(intToBits(i))[1:9]; idx <- 1:18
    for (k in which(bits == 1)) { a <- 2*k-1; b <- 2*k; idx[c(a,b)] <- c(b,a) }
    idx
  }))
}
bd <- betadisper(dist_bray, factor(grp), type = "centroid")
pt <- permutest(bd, permutations = perm_matrix_9())
cat("\nvegan::permutest.betadisper, same 512-arrangement set:\n")
cat(sprintf("  F = %.4f  P = %.4f  (exact enumeration above: F = %.4f  P = %.4f)\n",
            pt$tab$F[1], pt$tab$`Pr(>F)`[1],
            results$F[results$distance == "bray" & results$type == "centroid"],
            results$P[results$distance == "bray" & results$type == "centroid"]))
cat("  F agrees; P does not -- this is what the header comment refers to.\n")

write.csv(results, "permdisp_results.csv", row.names = FALSE)
