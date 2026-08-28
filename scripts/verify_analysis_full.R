## verify_analysis_full.R
##
## Full reproduction check for the patient-level composition analysis:
## per-class relative abundance on both denominators, the community-level
## PERMANOVA and every robustness variant of it, the ISCHIA exclusion
## criteria applied independently to this cohort, the direction-of-shift
## contrast, and the purity-vs-enrichment comparison. Everything below
## prints an OK/MISMATCH line against the value quoted in the reply text,
## so a mismatch here means the reply needs to change, not the script.

library(data.table)
library(vegan)

CSV <- "/Users/akhaliq/Desktop/rebuttle_natgen/new/pdac_spot_level.csv.gz"

## ---- pass/fail helper --------------------------------------------------------
.RESULTS <- list()
ok <- function(label, got, expected, tol = 0.01, rel = FALSE) {
  d <- if (rel) abs(got - expected) / abs(expected) else abs(got - expected)
  pass <- is.finite(got) && d <= tol
  .RESULTS[[length(.RESULTS) + 1]] <<- pass
  cat(sprintf("  %-52s %10s  (expected %s)  %s\n",
              label,
              if (is.numeric(got)) sprintf("%.4f", got) else as.character(got),
              if (is.numeric(expected)) sprintf("%.4f", expected) else as.character(expected),
              if (pass) "OK" else "*** MISMATCH ***"))
  invisible(pass)
}

md <- as.data.table(read.csv(gzfile(CSV), stringsAsFactors = FALSE))
md[, site9 := c("Pancreas"="PT", "Liver"="HM", "Lymph node"="LNM",
                "Normal Pancreas"="NP")[site]]
setnames(md, "ecotype", "cc")
cat("spots:", nrow(md), " patients:", uniqueN(md$patient), "\n")

PATS <- c("PT_2","PT_3","PT_4","PT_6","PT_8","PT_9","PT_10","PT_11","PT_12")
E6   <- c("CC1","CC2","CC3","CC4","CC5","CC7")
E7   <- c(E6, "CC8")
ALL10 <- c("CC1","CC2","CC3","CC4","CC5","CC6","CC7","CC8","CC9","CC10")

## ---- shared builder: patient x site composition, chosen classes -------------
build_prop <- function(pats = PATS, classes = E6, sites = c("PT","HM"), renorm = TRUE) {
  s <- md[patient %in% pats & site9 %in% sites]
  pw <- s[, .N, by = .(patient, site9, cc)]
  pw[, prop := N / sum(N), by = .(patient, site9)]
  w <- dcast(pw, patient + site9 ~ cc, value.var = "prop", fill = 0)
  keep <- intersect(classes, names(w))
  X <- as.matrix(w[, ..keep])
  if (renorm) X <- X / rowSums(X)
  # compute the sort index ONCE, on the pre-sort w, then apply identically to
  # X and to the returned site/patient vectors so they stay aligned. (A prior
  # version computed this index on an already-sorted copy of w, which
  # silently returned X in its original order while site/patient were
  # sorted -- a real misalignment bug, not a data issue. Confirmed by
  # diffing against the original scripts' build() output on step A.)
  idx <- order(match(w$patient, pats), factor(w$site9, levels = c("HM","PT")))
  list(X = X[idx, , drop = FALSE], site = w$site9[idx], patient = w$patient[idx])
}

## ---- exact restricted permutation on N matched pairs (2^N arrangements) -----
## N.B. must scale with however many pairs are actually present -- hardcoding
## 9 pairs/512 perms breaks silently (out-of-bounds -> NA) whenever this is
## called on a reduced set, e.g. leave-one-patient-out (8 pairs).
exact_permanova <- function(X, method = "bray") {
  D <- as.matrix(vegdist(X, method = method))
  n <- nrow(D)
  npairs <- n / 2
  stopifnot(npairs == round(npairs))
  nperm <- 2^npairs
  A <- -0.5 * D^2; C <- diag(n) - 1/n; G <- C %*% A %*% C
  pseudoF <- function(g) {
    SST <- sum(diag(G)); tot <- 0
    for (l in unique(g)) { i <- which(g == l); tot <- tot + sum(G[i,i]) / length(i) }
    SSA <- tot - sum(G) / n
    c(F = (SSA/1) / ((SST - SSA)/(n - 2)), R2 = SSA/SST)
  }
  grp0 <- rep(c("HM","PT"), npairs)  # matches build_prop's HM,PT ordering per patient
  obs <- pseudoF(grp0); hits <- 0
  for (i in 0:(nperm - 1)) {
    f <- as.integer(intToBits(i))[1:npairs]; g <- grp0
    for (k in which(f == 1)) { a <- 2*k-1; b <- 2*k; tmp <- g[a]; g[a] <- g[b]; g[b] <- tmp }
    if (pseudoF(g)["F"] >= obs["F"] - 1e-12) hits <- hits + 1
  }
  c(obs, P = hits/nperm, hits = hits, nperm = nperm)
}

clr_transform <- function(X) {
  Z <- X; Z[Z == 0] <- min(Z[Z > 0]) / 2; Z <- Z / rowSums(Z)
  t(apply(log(Z), 1, function(r) r - mean(r)))
}

cat("\n==============================================================================\n")
cat("SECTION 1 -- Per-class relative abundance (exact paired Wilcoxon, BH across classes)\n")
cat("==============================================================================\n")

## ---- six-habitat denominator -------------------------------------------------
b6 <- build_prop(classes = E6)
w6 <- data.table(patient = b6$patient, site = b6$site, b6$X)
wide <- dcast(w6, patient ~ site, value.var = E6)
res1 <- rbindlist(lapply(E6, function(cc) {
  hm <- wide[[paste0(cc, "_HM")]] * 100; pt <- wide[[paste0(cc, "_PT")]] * 100
  wt <- suppressWarnings(wilcox.test(hm, pt, paired = TRUE, exact = TRUE))
  data.table(cc = cc, delta = mean(hm - pt), pos = sum(hm > pt), P = wt$p.value)
}))
res1[, BH := p.adjust(P, "BH")]
print(res1)

ok("CC1 mean shift (pp)", res1[cc=="CC1", delta], -16.4, tol = 0.1)
ok("CC1 BH-FDR",          res1[cc=="CC1", BH],    0.0117, tol = 0.001)
ok("CC5 mean shift (pp)", res1[cc=="CC5", delta], -17.7, tol = 0.1)
ok("CC5 BH-FDR",          res1[cc=="CC5", BH],    0.0117, tol = 0.001)
ok("CC2 mean shift (pp)", res1[cc=="CC2", delta],  19.9, tol = 0.15)
ok("CC2 count of 9",      res1[cc=="CC2", pos],    8, tol = 0)
ok("CC2 P",               res1[cc=="CC2", P],      0.0078, tol = 0.001)
ok("CC2 BH-FDR",          res1[cc=="CC2", BH],     0.0156, tol = 0.001)

## ---- all-ten denominator, seven tumour classes tested ------------------------
b10 <- build_prop(classes = E7, renorm = FALSE)  # their denominator: proportion of ALL spots
w10 <- data.table(patient = b10$patient, site = b10$site, b10$X)
wide10 <- dcast(w10, patient ~ site, value.var = E7)
res10 <- rbindlist(lapply(E7, function(cc) {
  hm <- wide10[[paste0(cc, "_HM")]] * 100; pt <- wide10[[paste0(cc, "_PT")]] * 100
  wt <- suppressWarnings(wilcox.test(hm, pt, paired = TRUE, exact = TRUE))
  data.table(cc = cc, P = wt$p.value)
}))
res10[, BH := p.adjust(P, "BH")]
print(res10)

ok("CC1 P (all-ten denom)", res10[cc=="CC1", P],  0.0039, tol = 0.001)
ok("CC1 BH (all-ten denom)", res10[cc=="CC1", BH], 0.0091, tol = 0.001)
ok("CC5 P (all-ten denom)", res10[cc=="CC5", P],  0.0039, tol = 0.001)
ok("CC5 BH (all-ten denom)", res10[cc=="CC5", BH], 0.0091, tol = 0.001)
ok("CC8 P (all-ten denom)", res10[cc=="CC8", P],  0.0039, tol = 0.001)
ok("CC8 BH (all-ten denom)", res10[cc=="CC8", BH], 0.0091, tol = 0.001)
ok("CC4 P (all-ten denom)", res10[cc=="CC4", P],  0.0742, tol = 0.005)
ok("CC4 BH (all-ten denom)", res10[cc=="CC4", BH], 0.1299, tol = 0.005)
ok("CC7 P (all-ten denom)", res10[cc=="CC7", P],  0.7344, tol = 0.01)
ok("CC2 BH (all-ten denom, theirs)", res10[cc=="CC2", BH], 0.9102, tol = 0.01)

cat("\n==============================================================================\n")
cat("SECTION 2 -- PERMANOVA (community-level), six-habitat denominator + robustness\n")
cat("==============================================================================\n")

main <- exact_permanova(b6$X, "bray")
ok("Bray-Curtis F", main["F"], 3.88, tol = 0.02)
ok("Bray-Curtis R2", main["R2"], 0.195, tol = 0.005)
ok("Bray-Curtis P", main["P"], 0.0039, tol = 0.001)

## leave-one-class-out
cat("\n  leave-one-class-out:\n")
loco <- sapply(E6, function(drop) {
  bb <- build_prop(classes = setdiff(E6, drop))
  exact_permanova(bb$X, "bray")
})
print(round(t(loco), 4))
ok("LOCO: number at P<=0.0117", sum(loco["P",] <= 0.0117 + 1e-4), 6, tol = 0)
ok("LOCO: F range min", min(loco["F",]), 2.38, tol = 0.1)
ok("LOCO: F range max", max(loco["F",]), 4.69, tol = 0.1)

## leave-one-patient-out
cat("\n  leave-one-patient-out:\n")
loo <- sapply(PATS, function(drop) {
  bb <- build_prop(pats = setdiff(PATS, drop), classes = E6)
  exact_permanova(bb$X, "bray")
})
print(round(t(loo), 4))
ok("LOO: all nine at P=0.0078", sum(abs(loo["P",] - 0.0078) < 5e-4), 9, tol = 0)
ok("LOO: F range min", min(loo["F",]), 2.80, tol = 0.1)
ok("LOO: F range max", max(loo["F",]), 4.09, tol = 0.1)

## partial PERMANOVA conditioning on tumour-epithelial fraction (uses tumour_frac col)
tf <- md[patient %in% PATS & site9 %in% c("PT","HM") & cc %in% E6,
         .(tf = mean(tumour_frac)), by = .(patient, site9)]
tf <- tf[order(match(patient, PATS), factor(site9, levels = c("HM","PT")))]
Xb <- b6$X; grp <- b6$site
Db <- as.matrix(vegdist(Xb, "bray"))
n <- nrow(Db)
## adonis2's reported p-value for this test doesn't match full enumeration,
## even when it's given the exact restricted permutation matrix by hand --
## checked this against a McArdle-Anderson decomposition computed directly
## over all 512 arrangements. The F it returns is right, the p-value isn't.
## Computing the sequential SS decomposition directly instead of trusting
## adonis2's own significance test.
tfc <- matrix(tf$tf - mean(tf$tf), ncol = 1)
Htf <- tfc %*% solve(t(tfc) %*% tfc) %*% t(tfc)
A <- -0.5 * Db^2; C <- diag(n) - 1/n; G <- C %*% A %*% C
SS_tf <- sum(diag(Htf %*% G))
SS_total <- sum(diag(G))
seq_partial_F <- function(g) {
  gd <- as.numeric(factor(g)) - mean(as.numeric(factor(g)))
  gd <- matrix(gd, ncol = 1)
  Xfull <- cbind(tfc, gd)
  Hfull <- Xfull %*% solve(t(Xfull) %*% Xfull) %*% t(Xfull)
  SS_full <- sum(diag(Hfull %*% G))
  SS_grp <- SS_full - SS_tf
  SS_resid <- SS_total - SS_full
  c(F = (SS_grp / 1) / (SS_resid / (n - 3)), R2 = SS_grp / SS_total)
}
grp0 <- grp
obs <- seq_partial_F(grp0)
hits <- 0L
for (i in 0:511) {
  bits <- as.integer(intToBits(i))[1:9]; g <- grp0
  for (k in which(bits == 1)) { a <- 2*k-1; b <- 2*k; tmp <- g[a]; g[a] <- g[b]; g[b] <- tmp }
  if (seq_partial_F(g)["F"] >= obs["F"] - 1e-8) hits <- hits + 1L
}
Fpartial <- obs["F"]; R2partial <- obs["R2"]; Ppartial <- hits / 512
cat(sprintf("partial PERMANOVA (exact enumeration): F=%.4f  R2=%.4f  P=%d/512=%.4f\n",
            Fpartial, R2partial, hits, Ppartial))
ok("partial PERMANOVA F (site | tf)", Fpartial, 6.83, tol = 0.5)
ok("partial PERMANOVA P (site | tf)", Ppartial, 0.0039, tol = 0.01)

## Aitchison
clr <- clr_transform(b6$X)
# reuse exact_permanova machinery but on Euclidean distance of CLR coords directly
aitch_permanova <- function(clr_mat) {
  D <- as.matrix(dist(clr_mat))
  n <- nrow(D); npairs <- n/2; nperm <- 2^npairs
  A <- -0.5*D^2; C <- diag(n)-1/n; G <- C %*% A %*% C
  pseudoF <- function(g) {
    SST <- sum(diag(G)); tot <- 0
    for (l in unique(g)) { i <- which(g==l); tot <- tot + sum(G[i,i])/length(i) }
    SSA <- tot - sum(G)/n
    c(F=(SSA/1)/((SST-SSA)/(n-2)), R2=SSA/SST)
  }
  grp0 <- rep(c("HM","PT"), npairs)
  obs <- pseudoF(grp0); hits <- 0
  for (i in 0:(nperm-1)) {
    f <- as.integer(intToBits(i))[1:npairs]; g <- grp0
    for (k in which(f==1)) { a<-2*k-1; b<-2*k; tmp<-g[a]; g[a]<-g[b]; g[b]<-tmp }
    if (pseudoF(g)["F"] >= obs["F"] - 1e-12) hits <- hits+1
  }
  c(obs, P = hits/nperm)
}
obsA <- aitch_permanova(clr)
ok("Aitchison F", obsA["F"], 3.86, tol = 0.05)
ok("Aitchison R2", obsA["R2"], 0.194, tol = 0.01)
ok("Aitchison P", obsA["P"], 0.0039, tol = 0.001)

cat("\n==============================================================================\n")
cat("SECTION 3 -- Their code, run verbatim\n")
cat("==============================================================================\n")
cat("  (Run RUN_THEIR_PERMANOVA.R and RUN_THEIR_CLR.R separately -- their full\n")
cat("   console output is long and already prints its own OK-style comparisons.\n")
cat("   Expected: PERMANOVA F=6.894 R2=0.375 P=0.001 (10 seeds: 0.001 x9, 0.002 x1)\n")
cat("             CLR ten-class: CC2 FDR 0.327; CC8 FDR 2.3e-6; CC8 mean diff 6.62\n")
cat("             CLR six-class: CC2 FDR 0.038(t)/0.023(w); CC5 FDR 0.003(t)/0.012(w)\n")
cat("             (corrected 27 Aug 2026 from RUN_THEIR_CLR_output.txt; the letter\n")
cat("              now quotes these. The old 0.050/0.008/0.023 line was wrong.)\n\n")

cat("==============================================================================\n")
cat("SECTION 4 -- ISCHIA-criterion denominator (independent derivation)\n")
cat("==============================================================================\n")
crit_s <- md[patient %in% PATS & site9 %in% c("PT","HM")]
crit <- crit_s[, {
  byp <- sort(table(patient), decreasing = TRUE)
  .(spots = .N, n_PT = uniqueN(patient[site9=="PT"]), n_HM = uniqueN(patient[site9=="HM"]),
    top_patient_pct = round(100*byp[1]/.N, 1),
    tumour = round(mean(tumour_frac),3), hepatocyte = round(mean(hep_frac),3))
}, by = cc][order(cc)]
print(crit)

cat("\n  ISCHIA's stated exclusion rule is a documented judgment call (anatomy +\n")
cat("  sample-specificity), not a single computable threshold -- so this checks\n")
cat("  the cited justification numbers themselves, not a re-derivation:\n\n")
ok("CC6 hepatocyte fraction (anatomy: normal liver)", crit[cc=="CC6", hepatocyte], 0.805, tol=0.01)
ok("CC6 n_PT (should be 0 -- absent from every primary)", crit[cc=="CC6", n_PT], 0, tol=0)
ok("CC8 hepatocyte fraction (anatomy: tumour-liver interface)", crit[cc=="CC8", hepatocyte], 0.426, tol=0.01)
ok("CC8 n_PT (should be 0)", crit[cc=="CC8", n_PT], 0, tol=0)
ok("CC10 total spots (anatomy: near-absent)", crit[cc=="CC10", spots], 2, tol=0)
ok("CC9 top-patient share (sample-specificity)", crit[cc=="CC9", top_patient_pct], 74.5, tol=1)
ok("CC9 n_PT (6 of 9 -- present, but sample-specific)", crit[cc=="CC9", n_PT], 6, tol=0)

EXCLUDE <- c("CC6","CC8","CC10","CC9")
retained <- setdiff(ALL10, EXCLUDE)
cat("\n  excluded on ISCHIA's stated criteria:", paste(EXCLUDE, collapse=", "), "\n")
cat("  retained:", paste(sort(retained), collapse=", "), "\n")
ok("retained set == E6 (six shared classes)", as.integer(setequal(retained, E6)), 1, tol=0)

cat("\n==============================================================================\n")
cat("SECTION 4b -- CC8 ladder, A through D\n")
cat("==============================================================================\n")
bA <- build_prop(pats = unique(md$patient), classes = ALL10, sites = c("PT","HM"), renorm = FALSE)
DA <- vegdist(bA$X, "bray")
set.seed(1)
pA <- adonis2(DA ~ site, data = data.frame(site = bA$site), permutations = 999, strata = bA$patient)
ok("A: all aggregates, 10 classes, F", pA$F[1], 8.46, tol = 0.1)
ok("A: R2", pA$R2[1], 0.308, tol = 0.01)

rB <- exact_permanova(build_prop(classes = ALL10, renorm = FALSE)$X, "bray")
ok("B: 9 pairs, 10 classes, F", rB["F"], 6.95, tol = 0.05)
ok("B: R2", rB["R2"], 0.303, tol = 0.005)

rC <- exact_permanova(build_prop(classes = E7)$X, "bray")
ok("C: 9 pairs, 7 classes incl CC8, F", rC["F"], 6.00, tol = 0.05)
ok("C: R2", rC["R2"], 0.273, tol = 0.005)

rD <- main  # already computed above, six shared classes
ok("D: 9 pairs, 6 shared classes, F", rD["F"], 3.88, tol = 0.02)
ok("D: R2", rD["R2"], 0.195, tol = 0.005)

cat("\n==============================================================================\n")
cat("SECTION 5 -- Direction of the shift (CC1+CC5 myCAF vs CC2+CC3 tumour)\n")
cat("==============================================================================\n")
myCAF <- c("CC1","CC5"); tumC <- c("CC2","CC3")
w5 <- data.table(patient = b6$patient, site = b6$site, b6$X)  # b6$X already sums to 1 over E6
wide5 <- dcast(w5, patient ~ site, value.var = E6)
wide5[, myCAF_HM := CC1_HM + CC5_HM]; wide5[, myCAF_PT := CC1_PT + CC5_PT]
wide5[, tum_HM := CC2_HM + CC3_HM];   wide5[, tum_PT := CC2_PT + CC3_PT]
wide5[, contrast_HM := tum_HM - myCAF_HM]; wide5[, contrast_PT := tum_PT - myCAF_PT]

wm <- suppressWarnings(wilcox.test(wide5$myCAF_HM, wide5$myCAF_PT, paired=TRUE, exact=TRUE))
wt <- suppressWarnings(wilcox.test(wide5$tum_HM, wide5$tum_PT, paired=TRUE, exact=TRUE))
wc <- suppressWarnings(wilcox.test(wide5$contrast_HM, wide5$contrast_PT, paired=TRUE, exact=TRUE))

ok("myCAF mean HM (%)", mean(wide5$myCAF_HM)*100, 7.7, tol = 1)
ok("myCAF mean PT (%)", mean(wide5$myCAF_PT)*100, 41.7, tol = 1)
ok("myCAF: 0 of 9 rise", sum(wide5$myCAF_HM > wide5$myCAF_PT), 0, tol=0)
ok("myCAF P", wm$p.value, 0.0039, tol = 0.001)

ok("tumour-dominated mean HM (%)", mean(wide5$tum_HM)*100, 57.0, tol = 1)
ok("tumour-dominated mean PT (%)", mean(wide5$tum_PT)*100, 26.9, tol = 1)
ok("tumour-dominated: 8 of 9 rise", sum(wide5$tum_HM > wide5$tum_PT), 8, tol=0)
ok("tumour-dominated P", wt$p.value, 0.0078, tol = 0.001)

ok("contrast: 9 of 9 move", sum(wide5$contrast_HM > wide5$contrast_PT), 9, tol=0)
ok("contrast P", wc$p.value, 0.0039, tol = 0.001)

## excluding PT_9
wide5_no9 <- wide5[patient != "PT_9"]
wc9 <- suppressWarnings(wilcox.test(wide5_no9$contrast_HM, wide5_no9$contrast_PT, paired=TRUE, exact=TRUE))
ok("excl. PT_9: contrast mean shift (pp)", mean(wide5_no9$contrast_HM - wide5_no9$contrast_PT)*100, 61.9, tol=1)
ok("excl. PT_9: 8 of 8 move", sum(wide5_no9$contrast_HM > wide5_no9$contrast_PT), 8, tol=0)
ok("excl. PT_9: P", wc9$p.value, 0.0078, tol = 0.001)

## correlation check: myCAF loss vs tumour gain should be near zero (closure warning)
loss <- wide5$myCAF_PT - wide5$myCAF_HM
gain <- wide5$tum_HM - wide5$tum_PT
ok("corr(myCAF loss, tumour gain), Pearson", cor(loss, gain), 0.08, tol = 0.15)
ok("corr(myCAF loss, tumour gain), Spearman", cor(loss, gain, method="spearman"), 0.00, tol = 0.15)

cat("\n==============================================================================\n")
cat("SECTION 6 -- Purity versus enrichment (both on six-habitat denominator)\n")
cat("==============================================================================\n")
tf6 <- md[patient %in% PATS & site9 %in% c("PT","HM") & cc %in% E6,
          .(tf = mean(tumour_frac)), by = .(patient, site9)]
tf6w <- dcast(tf6, patient ~ site9, value.var = "tf")
wtf <- suppressWarnings(wilcox.test(tf6w$HM, tf6w$PT, paired=TRUE, exact=TRUE))
ok("tumour-epithelial fraction, HM mean", mean(tf6w$HM), 0.448, tol = 0.02)
ok("tumour-epithelial fraction, PT mean", mean(tf6w$PT), 0.383, tol = 0.02)
ok("tumour fraction: 4 of 9 rise", sum(tf6w$HM > tf6w$PT), 4, tol=0)
ok("tumour fraction P", wtf$p.value, 0.65, tol = 0.05)

wshare <- data.table(patient = b6$patient, site = b6$site,
                      share = b6$X[,"CC2"] + b6$X[,"CC3"])
wshare_w <- dcast(wshare, patient ~ site, value.var = "share")
wsh <- suppressWarnings(wilcox.test(wshare_w$HM, wshare_w$PT, paired=TRUE, exact=TRUE))
ok("tumour-dominated habitat share, HM mean", mean(wshare_w$HM), 0.570, tol = 0.02)
ok("tumour-dominated habitat share, PT mean", mean(wshare_w$PT), 0.269, tol = 0.02)
ok("habitat share: 8 of 9 rise", sum(wshare_w$HM > wshare_w$PT), 8, tol=0)
ok("habitat share P", wsh$p.value, 0.0078, tol = 0.001)

## all-spot purity control (should FALL, illustrating the denominator trap)
tf_all <- md[patient %in% PATS & site9 %in% c("PT","HM"), .(tf = mean(tumour_frac)), by=.(patient, site9)]
tf_all_w <- dcast(tf_all, patient ~ site9, value.var = "tf")
ok("all-spot tumour fraction, HM mean (control, should be LOWER)", mean(tf_all_w$HM), 0.305, tol = 0.02)
ok("all-spot tumour fraction, PT mean (control)", mean(tf_all_w$PT), 0.379, tol = 0.02)

cat("\n==============================================================================\n")
cat("SUMMARY\n")
cat("==============================================================================\n")
n_fail <- sum(!unlist(.RESULTS))
n_tot <- length(.RESULTS)
cat(sprintf("  %d / %d checks passed\n", n_tot - n_fail, n_tot))
if (n_fail > 0) {
  cat(sprintf("  *** %d CHECK(S) FAILED -- a disagreement here is the finding,\n", n_fail))
  cat("      not something to quietly patch. Fix the number in the reply, not\n")
  cat("      the check.\n")
} else {
  cat("  All checks in this script passed. Still worth running separately:\n")
  cat("    - run_their_permanova.R and run_their_clr.R (their deposited code, verbatim)\n")
  cat("    - the CC6 interface control (not automated here)\n")
  cat("    - make_all_figures.R, to regenerate the main figures\n")
}
cat("\n  REMINDER: 0.0039 and 0.0078 are permutation FLOORS (2/512, 4/512), not\n")
cat("  measures of strength. Do not describe them as 'overwhelming'.\n")
