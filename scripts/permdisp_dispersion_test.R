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
## when you hand it the exact permutation matrix yourself. The F-statistic
## it reports is fine, the p-value isn't. Checked this two different ways
## before trusting it, so the p-values below come from refitting betadisper
## under each of the 512 valid arrangements directly and counting, not from
## permutest().

library(data.table)
library(vegan)

csv_path <- "/Users/akhaliq/Desktop/rebuttle_natgen/new/pdac_spot_level.csv.gz"

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

write.csv(results, "permdisp_results.csv", row.names = FALSE)
