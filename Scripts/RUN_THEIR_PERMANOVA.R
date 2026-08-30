## =============================================================================
##  RUN_THEIR_PERMANOVA.R
##
##  Reproduces the PERMANOVA in Li & Ma's 02_patient_level_analysis.R.
##
##  Their pipeline reads an 11.6 GB Seurat object (PDAC_Updated.rds) only to
##  extract spot-level metadata. This script skips that step and reads the
##  exported spot table instead. Everything downstream is their code, copied
##  verbatim from lines 86-93 and 159-161 of their script — the same objects,
##  the same names, the same adonis2 call.
##
##  Run:  Rscript RUN_THEIR_PERMANOVA.R
##  Needs: data.table, vegan
##
##  If you would rather run their script end to end from the .rds, do that —
##  this is the shortcut, not a replacement. Either way, report what it gives.
## =============================================================================

## ---- preflight --------------------------------------------------------------
## If either package is missing, everything below fails with confusing errors
## like `could not find function "."` — so stop here instead.
need <- setdiff(c("data.table", "vegan"), rownames(installed.packages()))
if (length(need)) {
  stop("Missing package(s): ", paste(need, collapse = ", "),
       "\n  Run:  install.packages(c(\"data.table\",\"vegan\"))\n")
}
suppressPackageStartupMessages({ library(data.table); library(vegan) })

## Look for the CSV in the working directory, then in ~/Downloads.
CSV <- "pdac_spot_level.csv.gz"
if (!file.exists(CSV)) {
  alt <- path.expand("~/Downloads/pdac_spot_level.csv.gz")
  if (file.exists(alt)) CSV <- alt else
    stop("Cannot find pdac_spot_level.csv.gz. setwd() to the folder holding it.")
}
cat("reading:", CSV, "\n")

## ---- stand in for their 01_prepare_data.R output ---------------------------
## Their `md` is a data.table with columns patient, site, cc.
## Their site labels are PT / HM / LNM; ours are Pancreas / Liver / Lymph node.
md <- as.data.table(read.csv(gzfile(CSV), stringsAsFactors = FALSE))
md[, site := c("Pancreas"="PT", "Liver"="HM", "Lymph node"="LNM",
               "Normal Pancreas"="NP")[site]]
setnames(md, "ecotype", "cc")
cat("spots:", nrow(md), " patients:", uniqueN(md$patient), "\n\n")

## ---- THEIR CODE, lines 86-93 ------------------------------------------------
pw <- md[site %in% c("PT","HM","LNM"),
         .N, by = .(patient, site, cc)]
pw[, prop := N / sum(N), by = .(patient, site)]
propw <- dcast(pw, patient + site ~ cc, value.var = "prop", fill = 0)

cat("propw dimensions:", nrow(propw), "patient-site rows x",
    ncol(propw) - 2, "composition classes\n")
print(table(propw$site))
cat("\n")

## ---- THEIR CODE, lines 159-161 ---------------------------------------------
comp_mat <- as.matrix(propw[, -(1:2)])
rownames(comp_mat) <- paste(propw$patient, propw$site, sep = "|")
bc <- vegdist(comp_mat, method = "bray")
## THEIR CALL, verbatim:
pa <- tryCatch(adonis2(bc ~ site, data = propw, permutations = 999,
                       strata = propw$patient),
               error = function(e) { warning("PERMANOVA failed: ",
                                             conditionMessage(e)); NULL })

cat("=========== THEIR PERMANOVA, AS CODED ===========\n")
if (is.null(pa)) {
  cat("*** The adonis2 call ERRORED. Their script wraps it in tryCatch(),\n")
  cat("*** so this would have failed silently for them too. Report this.\n")
} else {
  print(pa)
  cat("\n  F  =", sprintf("%.3f", pa$F[1]),
      "  R2 =", sprintf("%.3f", pa$R2[1]),
      "  P  =", sprintf("%.4f", pa$`Pr(>F)`[1]), "\n")
  cat("  expected: F ~ 6.89, R2 ~ 0.375, P = 0.001 (the floor for 999 perms)\n")
}

## ---- P is a random draw: repeat it ------------------------------------------
cat("\n--- same call, ten different seeds ---\n")
ps <- sapply(1:10, function(s) {
  set.seed(s)
  r <- tryCatch(adonis2(bc ~ site, data = propw, permutations = 999,
                        strata = propw$patient), error = function(e) NULL)
  if (is.null(r)) NA_real_ else r$`Pr(>F)`[1]
})
print(round(ps, 4))
cat("  range", sprintf("%.4f", min(ps, na.rm = TRUE)), "-",
    sprintf("%.4f", max(ps, na.rm = TRUE)), "\n")

## ---- for comparison: our version -------------------------------------------
## Nine matched PT-HM pairs, six shared classes (CC6/9/10 normal, CC8 absent
## from every primary), permutation restricted within patient and enumerated
## in full rather than sampled.
cat("\n=========== OUR VERSION, FOR COMPARISON ===========\n")
PATS <- c("PT_2","PT_3","PT_4","PT_6","PT_8","PT_9","PT_10","PT_11","PT_12")
E6   <- c("CC1","CC2","CC3","CC4","CC5","CC7")

sub <- md[patient %in% PATS & site %in% c("PT","HM") & cc %in% E6]
pw6 <- sub[, .N, by = .(patient, site, cc)]
pw6[, prop := N / sum(N), by = .(patient, site)]
w6  <- dcast(pw6, patient + site ~ cc, value.var = "prop", fill = 0)
w6  <- w6[order(match(patient, PATS), factor(site, levels = c("HM","PT")))]

X <- as.matrix(w6[, ..E6]); grp <- w6$site
D <- as.matrix(vegdist(X, method = "bray"))
n <- nrow(X); G <- { A <- -0.5 * D^2; C <- diag(n) - 1/n; C %*% A %*% C }
pseudoF <- function(g) {
  SST <- sum(diag(G)); tot <- 0
  for (l in unique(g)) { i <- which(g == l); tot <- tot + sum(G[i,i])/length(i) }
  SSA <- tot - sum(G)/n
  c(F = (SSA/1)/((SST - SSA)/(n - 2)), R2 = SSA/SST)
}
obs <- pseudoF(grp)
hits <- 0
for (i in 0:511) {
  f <- as.integer(intToBits(i))[1:9]; g <- grp
  for (k in which(f == 1)) { a <- 2*k-1; b <- 2*k; tmp <- g[a]; g[a] <- g[b]; g[b] <- tmp }
  if (pseudoF(g)["F"] >= obs["F"] - 1e-12) hits <- hits + 1
}
cat(sprintf("  F = %.2f   R2 = %.3f   P = %d/512 = %.4f  (exact enumeration)\n",
            obs["F"], obs["R2"], hits, hits/512))
cat("  expected: F ~ 3.88, R2 ~ 0.195, P = 2/512 = 0.0039\n")

cat("\n=========== WHAT TO REPORT BACK ===========\n")
cat("  1. The adonis2 table above, and whether it errored\n")
cat("  2. The ten-seed P values\n")
cat("  3. Whether our F = 3.88 / R2 = 0.195 reproduces under vegan\n")
cat("  A disagreement with the expected values here is the finding -- it\n")
cat("  means the reported number needs correcting, not this script.\n\n")
