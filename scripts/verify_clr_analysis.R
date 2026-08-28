## verify_clr_analysis.R
##
## Reproduces the compositional (CLR) analysis from Ma & Li's own
## 03_sensitivity_analysis.R, section B -- their code, unaltered, just
## pointed at the exported spot table instead of the 11.6 GB Seurat
## object -- then repeats it restricted to the six classes that actually
## exist in both primary and metastatic arms.
##
## Expected, cross-checked against their deposited 03_results/.../clr_paired.tsv:
##    CC1  t p 0.000400  FDR 0.000933
##    CC2  t p 0.234     FDR 0.327
##    CC3  t p 0.292     FDR 0.340
##    CC4  t p 0.00301   FDR 0.00526
##    CC5  t p 5.98e-05  FDR 0.000209
##     CC7  t p 0.342     FDR 0.342
##     CC8  t p 3.34e-07  FDR 2.34e-06
## =============================================================================

library(data.table)

CSV <- "/Users/akhaliq/Desktop/rebuttle_natgen/new/pdac_spot_level.csv.gz"

md <- as.data.table(read.csv(gzfile(CSV), stringsAsFactors = FALSE))
md[, site := c("Pancreas"="PT", "Liver"="HM", "Lymph node"="LNM",
               "Normal Pancreas"="NP")[site]]
setnames(md, "ecotype", "cc")

tumor_cc <- c("CC1","CC2","CC3","CC4","CC5","CC7","CC8")   # their line 79

## ---- THEIR CODE, 03_sensitivity_analysis.R lines 43-60 ----------------------
pw <- md[site %in% c("PT","HM"), .N, by = .(patient, site, cc)]
pw[, prop := N / sum(N), by = .(patient, site)]
propw <- dcast(pw, patient + site ~ cc, value.var = "prop", fill = 0)
mat <- as.matrix(propw[, -(1:2)])
minpos <- min(mat[mat > 0]) / 2
mat_z <- mat; mat_z[mat_z == 0] <- minpos
clr <- t(scale(t(log(mat_z)), center = TRUE, scale = FALSE))   # row-wise clr
clrw <- cbind(propw[, .(patient, site)], as.data.table(clr))
resB <- rbindlist(lapply(tumor_cc, function(cc_i) {
  a <- clrw[site == "PT", .(patient, pt = get(cc_i))]
  b <- clrw[site == "HM", .(patient, hm = get(cc_i))]
  j <- merge(a, b, by = "patient")
  tt <- tryCatch(t.test(j$pt, j$hm, paired = TRUE), error = function(e) list(p.value = NA_real_))
  data.table(analysis = "B_clr_paired_t", cc = cc_i, n_pairs = nrow(j),
             mean_clr_diff = mean(j$pt - j$hm), t_p = tt$p.value)
}))
resB[, fdr := p.adjust(t_p, "BH")]

cat("=========== THEIR CLR, AS CODED (geometry = all ten classes) ===========\n")
print(resB)
cat("\n  How many classes enter the CLR geometry:", ncol(mat), "\n")
cat("  Of those, classes with zero spots in every primary:\n   ")
zero_in_PT <- names(which(colSums(mat[propw$site == "PT", , drop = FALSE]) == 0))
cat(if (length(zero_in_PT)) paste(zero_in_PT, collapse = ", ") else "none", "\n")
cat("  Each such zero was replaced by", signif(minpos, 3), "and logged.\n")

## ---- the same thing, restricted to classes present in both arms ------------
E6 <- c("CC1","CC2","CC3","CC4","CC5","CC7")
sub6  <- md[site %in% c("PT","HM") & cc %in% E6]
pw6   <- sub6[, .N, by = .(patient, site, cc)]
pw6[, prop := N / sum(N), by = .(patient, site)]
p6    <- dcast(pw6, patient + site ~ cc, value.var = "prop", fill = 0)
m6    <- as.matrix(p6[, -(1:2)])
m6[m6 == 0] <- min(m6[m6 > 0]) / 2
m6    <- m6 / rowSums(m6)
clr6  <- t(scale(t(log(m6)), center = TRUE, scale = FALSE))
c6w   <- cbind(p6[, .(patient, site)], as.data.table(clr6))
res6 <- rbindlist(lapply(E6, function(cc_i) {
  a <- c6w[site == "PT", .(patient, pt = get(cc_i))]
  b <- c6w[site == "HM", .(patient, hm = get(cc_i))]
  j <- merge(a, b, by = "patient")
  tt <- t.test(j$pt, j$hm, paired = TRUE)
  wt <- suppressWarnings(wilcox.test(j$pt, j$hm, paired = TRUE, exact = TRUE))
  data.table(cc = cc_i, n_pairs = nrow(j), mean_clr_diff = mean(j$pt - j$hm),
             t_p = tt$p.value, w_p = wt$p.value)
}))
res6[, `:=`(t_fdr = p.adjust(t_p, "BH"), w_fdr = p.adjust(w_p, "BH"))]

cat("\n=========== SAME CLR, SIX CLASSES PRESENT IN BOTH ARMS ===========\n")
print(res6)
cat("\n  expected: CC2 t FDR 0.038 and signed-rank FDR 0.023;",
    "\n            CC5 t FDR 0.003 and signed-rank FDR 0.012\n")

cat("\n=========== WHAT TO COMPARE ===========\n")
cat("  1. The first table against your deposited clr_paired.tsv — it should match.\n")
cat("  2. CC2 in the first table (FDR 0.327) against CC2 in the second.\n")
cat("  3. CC8's mean_clr_diff against every other class in the first table.\n\n")
