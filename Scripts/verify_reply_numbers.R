## verify_reply_numbers.R
##
## Reproduces every statistic quoted in the reply text against the raw
## spot-level data, block by block, and prints the computed value beside
## the quoted one -- OK if they match, MISMATCH if they don't.
##
## Needs vegan for the independent PERMANOVA cross-check in block M;
## everything else is base R. Runs fine without vegan, just skips that
## one block.
##
## Input: pdac_spot_level.csv.gz -- one row per QC-passed spot, columns
## spot, patient, section, site, ecotype, tumour_frac, hep_frac. Site
## levels: "Pancreas", "Liver", "Lymph node", "Normal Pancreas".

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
CSV <- find_spot_table()
cat("reading:", CSV, "\n")
spot <- read.csv(gzfile(CSV), stringsAsFactors = FALSE)

PATS  <- c("PT_2","PT_3","PT_4","PT_6","PT_8","PT_9","PT_10","PT_11","PT_12")
ALL10 <- c("CC1","CC2","CC3","CC4","CC5","CC6","CC7","CC8","CC9","CC10")
E7    <- c("CC1","CC2","CC3","CC4","CC5","CC7","CC8")   # tumour-associated
E6    <- c("CC1","CC2","CC3","CC4","CC5","CC7")         # six shared classes
DUP   <- "IU_PDA_HM2_2"                                  # duplicate hepatic section

fails <- 0L
ok <- function(label, got, want, tol = 1e-3) {
  hit <- is.finite(got) && abs(got - want) <= tol
  if (!hit) fails <<- fails + 1L
  cat(sprintf("  %-52s %12.4f  (letter %.4f)  %s\n",
              label, got, want, if (hit) "OK" else "*** MISMATCH ***"))
  invisible(hit)
}
hdr <- function(x) cat("\n", strrep("=", 78), "\n", x, "\n", strrep("=", 78), "\n", sep = "")

## ---------------------------------------------------------------- helpers ---

## Composition (percent) of one patient-site aggregate over `keep`,
## with the denominator taken over `src`.
comp <- function(d, p, s, keep, src = keep) {
  x <- d[d$patient == p & d$site == s & d$ecotype %in% src, ]
  v <- table(factor(x$ecotype, levels = src))
  as.numeric(v[keep]) / sum(v) * 100
}

## Paired difference matrix. dir = "HM-PT" or "PT-HM".
paired <- function(d, keep, src = keep, pats = PATS, dir = "HM-PT") {
  m <- t(sapply(pats, function(p) {
    a <- comp(d, p, "Liver", keep, src); b <- comp(d, p, "Pancreas", keep, src)
    if (dir == "HM-PT") a - b else b - a
  }))
  colnames(m) <- keep; rownames(m) <- pats; m
}

## Exact two-sided paired Wilcoxon signed-rank; zeros dropped.
wx <- function(x) {
  x <- x[x != 0]; n <- length(x)
  if (n == 0) return(1)
  r <- rank(abs(x)); W <- sum(r[x > 0]); tot <- sum(r); cnt <- 0
  for (i in 0:(2^n - 1)) {
    s <- as.integer(intToBits(i))[1:n]
    if (abs(sum(r[s == 1]) - tot/2) >= abs(W - tot/2) - 1e-9) cnt <- cnt + 1
  }
  cnt / 2^n
}

## Two-sided t tail, df = v.
tsf2 <- function(t, v) 2 * pt(-abs(t), v)

## OLS intercept test: y ~ covariate, report intercept (effect at covariate = 0).
site_effect <- function(y, tv) {
  f <- summary(lm(y ~ tv))$coefficients
  c(effect = f[1,1], t = f[1,3], p = f[1,4])
}

## Multiplicative zero replacement (Martin-Fernandez et al. 2003), as in
## RUN_analysis.R: zeros -> m/n, non-zeros rescaled to keep the sum at 1.
mult_repl <- function(M, nvec, m) {
  for (i in seq_len(nrow(M))) {
    r <- M[i, ]; delta <- m / nvec[i]; z <- r == 0
    if (any(z)) { r[z] <- delta; r[!z] <- r[!z] * (1 - sum(z) * delta) }
    M[i, ] <- r / sum(r)
  }
  M
}
clr <- function(M) log(M) - rowMeans(log(M))

## All 2^9 = 512 within-patient swaps of the site label.
## Rows of the data are ordered Liver_1, Pancreas_1, Liver_2, Pancreas_2, ...
perm_matrix <- function(npair = 9) {
  n <- 2 * npair
  t(sapply(0:(2^npair - 1), function(i) {
    f <- as.integer(intToBits(i))[1:npair]; idx <- 1:n
    for (k in which(f == 1)) { a <- 2*k - 1; b <- 2*k; idx[c(a,b)] <- c(b,a) }
    idx
  }))
}

## Gower-centred distance matrix.
gower <- function(D) { n <- nrow(D); A <- -0.5 * D^2; C <- diag(n) - 1/n; C %*% A %*% C }

bray_mat <- function(X) {
  n <- nrow(X); D <- matrix(0, n, n)
  for (i in 1:n) for (j in 1:n)
    D[i,j] <- sum(abs(X[i,] - X[j,])) / (sum(X[i,]) + sum(X[j,]))
  D
}
euclid_mat <- function(X) as.matrix(dist(X))

## Exact restricted PERMANOVA. If `cov` is supplied, the site term is tested
## sequentially after it (Type I SS), and only the site column is permuted.
permanova_exact <- function(D, site, cov = NULL, P = perm_matrix(length(site)/2)) {
  G <- gower(D); n <- nrow(D); SST <- sum(diag(G))
  trH <- function(M) { Q <- qr.Q(qr(M)); sum(diag(crossprod(Q, G %*% Q))) }
  one <- matrix(1, n, 1)
  stat <- function(sv) {
    S <- matrix(sv, n, 1)
    if (is.null(cov)) {
      ssA <- trH(cbind(one, S)) - trH(one); ssR <- SST - trH(cbind(one, S)); df <- n - 2
    } else {
      Cv <- matrix(cov, n, 1)
      ssA <- trH(cbind(one, Cv, S)) - trH(cbind(one, Cv))
      ssR <- SST - trH(cbind(one, Cv, S)); df <- n - 3
    }
    c(F = (ssA/1)/(ssR/df), R2 = ssA/SST)
  }
  obs <- stat(site)
  Fs  <- apply(P, 1, function(idx) stat(site[idx])["F"])
  c(obs, P = mean(Fs >= obs["F"] - 1e-12))
}

## Build the 18-row design used by every PERMANOVA below.
build_design <- function(keep, src = keep, pats = PATS, d = spot) {
  X <- do.call(rbind, lapply(pats, function(p)
    rbind(comp(d, p, "Liver", keep, src), comp(d, p, "Pancreas", keep, src))))
  X <- X / rowSums(X)
  ## FIX (audit 27 Aug 2026): the covariate must use the SAME six-habitat
  ## denominator as X. Averaging tumour_frac over ALL spots pulls in normal
  ## liver (CC6/CC8), which is 96.6% of PT_8's hepatic section, and yields
  ## F=7.61 / R2=0.240 / tf-alone R2=0.286 / Aitchison 9.05 -- the four
  ## MISMATCHes in the 26 Aug console log. Restricting to `keep` reproduces
  ## the reply text exactly: 6.83 / 0.157 / 0.498 / 11.20.
  tf <- unlist(lapply(pats, function(p) c(
    mean(d$tumour_frac[d$patient == p & d$site == "Liver"    & d$ecotype %in% keep]),
    mean(d$tumour_frac[d$patient == p & d$site == "Pancreas" & d$ecotype %in% keep]))))
  list(X = X, site = rep(c(1, 0), length(pats)), tf = tf)
}

cat("\nREPRODUCE_LETTER.R —", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n")
cat("spots:", nrow(spot), " sections:", length(unique(spot$section)), "\n")

## =========================================================== A. SIX-CLASS ===
hdr("A. Six-class paired Wilcoxon (CC8 removed) — their test, our denominator")

D6 <- paired(spot, E6)
res <- data.frame(
  cc   = E6,
  pos  = apply(D6, 2, function(v) sum(v > 0)),
  P    = apply(D6, 2, wx),
  mean = apply(D6, 2, mean),
  med  = apply(D6, 2, median))
res$BH <- p.adjust(res$P, "BH")
print(res, row.names = FALSE, digits = 4)
cat("\n")
ok("CC2 hepatic-higher count (of 9)",        res$pos[res$cc=="CC2"],  8,      0)
ok("CC2 exact paired Wilcoxon P",            res$P[res$cc=="CC2"],    0.0078, 5e-4)
ok("CC2 BH-FDR across six",                  res$BH[res$cc=="CC2"],   0.0156, 1e-3)
ok("CC2 mean shift (pp)",                    res$mean[res$cc=="CC2"], 19.9,   0.1)
ok("CC2 median shift (pp)",                  res$med[res$cc=="CC2"],  12.6,   0.1)
ok("CC1 BH-FDR",                             res$BH[res$cc=="CC1"],   0.0117, 1e-3)
ok("CC5 BH-FDR",                             res$BH[res$cc=="CC5"],   0.0117, 1e-3)
ok("CC1 mean shift (pp)",                    res$mean[res$cc=="CC1"], -16.4,  0.1)
ok("CC5 mean shift (pp)",                    res$mean[res$cc=="CC5"], -17.7,  0.1)
ok("CC3 count (letter: five of nine)",       res$pos[res$cc=="CC3"],  5,      0)
ok("CC3 P (letter: not significant)",        res$P[res$cc=="CC3"],    0.5703, 1e-3)
ok("CC3 mean shift (pp)",                    res$mean[res$cc=="CC3"], 10.2,   0.1)
ok("CC7 count (letter: five of nine)",       res$pos[res$cc=="CC7"],  5,      0)
ok("CC7 mean shift (pp)",                    res$mean[res$cc=="CC7"], 12.9,   0.1)

## Cross-check the exact Wilcoxon against base R.
p_base <- suppressWarnings(wilcox.test(
  sapply(PATS, function(p) comp(spot, p, "Liver", E6)[2]),
  sapply(PATS, function(p) comp(spot, p, "Pancreas", E6)[2]),
  paired = TRUE, exact = TRUE)$p.value)
ok("CC2 P from base wilcox.test (agreement)", p_base, res$P[res$cc=="CC2"], 5e-4)

## ================================================ B. THEIR SENSITIVITIES ====
hdr("B. Their three pre-specified sensitivity analyses, on the six classes")

loo <- sapply(1:9, function(i) wx(D6[-i, "CC2"]))
names(loo) <- PATS; print(round(loo, 4))
ok("CC2 leave-one-out minimum", min(loo), 0.0078, 5e-4)
ok("CC2 leave-one-out maximum", max(loo), 0.0156, 5e-4)
ok("CC2 leave-one-out: all significant", sum(loo < 0.05), 9, 0)

v <- D6[PATS != "PT_9", "CC2"]
ok("no-neoadjuvant (drop PT_9): count of 8", sum(v > 0),  8,      0)
ok("no-neoadjuvant: P",                      wx(v),       0.0078, 5e-4)
ok("no-neoadjuvant: mean (pp)",              mean(v),     22.4,   0.1)

D6s <- paired(spot[spot$section != DUP, ], E6)
ok("single hepatic lesion: count of 9",      sum(D6s[,"CC2"] > 0), 8,      0)
ok("single hepatic lesion: P",               wx(D6s[,"CC2"]),      0.0078, 5e-4)
ok("single hepatic lesion: mean (pp)",       mean(D6s[,"CC2"]),    19.6,   0.1)

## ============================================== C. SIZE-THRESHOLD SERIES ====
hdr("C. Minimum hepatic-microenvironment size")

sz <- sapply(PATS, function(p)
  sum(spot$patient == p & spot$site == "Liver" & spot$ecotype %in% E6))
print(sz)
ok("smallest hepatic TME (spots)", min(sz), 90,   0)
ok("largest hepatic TME (spots)",  max(sz), 3209, 0)
for (thr in c(100, 200)) {
  k <- which(sz >= thr); v <- D6[k, "CC2"]
  cat(sprintf("  >= %d spots: n=%d  %d/%d  P=%.4f  mean=%+.1f\n",
              thr, length(k), sum(v > 0), length(k), wx(v), mean(v)))
}
k <- which(sz >= 100); ok(">=100 spots: P", wx(D6[k,"CC2"]), 0.0156, 5e-4)
ok(">=100 spots: mean (pp)", mean(D6[k,"CC2"]), 21.2, 0.15)
k <- which(sz >= 200); ok(">=200 spots: P", wx(D6[k,"CC2"]), 0.0312, 1e-3)
ok(">=200 spots: mean (pp)", mean(D6[k,"CC2"]), 16.6, 0.15)

## ================================================== D. ABUNDANCE (FIG 1G) ===
hdr("D. Patient-level relative abundance — ISCHIA's Fig. 1G construction")

L  <- colMeans(t(sapply(PATS, function(p) comp(spot, p, "Liver",    E6))))
Pm <- colMeans(t(sapply(PATS, function(p) comp(spot, p, "Pancreas", E6))))
share <- setNames(L / (L + Pm), E6)
print(round(share, 3))
.expected_share <- c(CC7 = 0.74, CC2 = 0.69, CC3 = 0.66, CC4 = 0.39, CC5 = 0.16, CC1 = 0.14)
for (.cc_name in names(.expected_share)) {
  ok(sprintf("%s hepatic share", .cc_name), share[[.cc_name]], .expected_share[[.cc_name]], 0.006)
}

## ======================================================= E. CC8 STRUCTURE ===
hdr("E. CC8 as a structural zero")

c8h <- mean(sapply(PATS, function(p) comp(spot, p, "Liver",    E7, E7)[which(E7=="CC8")]))
c8p <- max( sapply(PATS, function(p) comp(spot, p, "Pancreas", E7, E7)[which(E7=="CC8")]))
ok("CC8 mean share of hepatic TME (%)", c8h, 31.8, 0.1)
ok("CC8 maximum share in ANY primary (%)", c8p, 0, 1e-6)

## ============================================= F. BRAY-CURTIS PERMANOVA =====
hdr("F. PERMANOVA on Bray-Curtis distances, exact restricted permutation")

des <- build_design(E6)
Db  <- bray_mat(des$X)
r1  <- permanova_exact(Db, des$site)
cat(sprintf("  site alone:  F=%.2f  R2=%.3f  P=%.4f\n", r1["F"], r1["R2"], r1["P"]))
ok("Bray-Curtis F",  r1[["F"]],  3.88,   0.02)
ok("Bray-Curtis R2", r1[["R2"]], 0.195,  0.002)
ok("Bray-Curtis P",  r1[["P"]],  0.0039, 5e-4)

cat("\n  leave-one-class-out:\n")
loco <- sapply(E6, function(c_) {
  dd <- build_design(setdiff(E6, c_))
  permanova_exact(bray_mat(dd$X), dd$site)
})
print(round(t(loco), 4))
ok("LOCO: number at P = 0.0039", sum(abs(loco["P",] - 0.0039) < 5e-4), 5, 0)
ok("LOCO: P when CC5 dropped",   loco["P","CC5"], 0.0117, 1e-3)

cat("\n  leave-one-patient-out:\n")
loop <- sapply(PATS, function(p_) {
  dd <- build_design(E6, pats = setdiff(PATS, p_))
  permanova_exact(bray_mat(dd$X), dd$site, P = perm_matrix(8))
})
print(round(t(loop), 4))
ok("LOO: all nine at P = 0.0078", sum(abs(loop["P",] - 0.0078) < 5e-4), 9, 0)
ok("LOO: minimum F", min(loop["F",]), 2.80, 0.03)
ok("LOO: maximum F", max(loop["F",]), 4.09, 0.03)

## ============================================== G. PARTIAL PERMANOVA ========
hdr("G. Partial PERMANOVA — site conditioned on tumour-epithelial fraction")

r2p <- permanova_exact(Db, des$site, cov = des$tf)
cat(sprintf("  site | tumour fraction:  F=%.2f  partial R2=%.3f  P=%.4f\n",
            r2p["F"], r2p["R2"], r2p["P"]))
ok("partial PERMANOVA F",  r2p[["F"]],  6.83,   0.05)
ok("partial PERMANOVA R2", r2p[["R2"]], 0.157,  0.003)
ok("partial PERMANOVA P",  r2p[["P"]],  0.0039, 5e-4)

Gb <- gower(Db); trH <- function(M) { Q <- qr.Q(qr(M)); sum(diag(crossprod(Q, Gb %*% Q))) }
one <- matrix(1, 18, 1)
r2tf <- (trH(cbind(one, matrix(des$tf, 18, 1))) - trH(one)) / sum(diag(Gb))
ok("tumour fraction alone, R2", r2tf, 0.498, 0.005)

## ============================================ H. AITCHISON PERMANOVA ========
hdr("H. PERMANOVA on Aitchison distance (CLR coordinates)")
cat("  NOTE: the reply text quotes the SIMPLE-replacement values. The multiplicative\n")
cat("  values are computed here too. If they differ materially, report both.\n\n")

Pmat <- des$X
nvec <- unlist(lapply(PATS, function(p) c(
  sum(spot$patient == p & spot$site == "Liver"    & spot$ecotype %in% E6),
  sum(spot$patient == p & spot$site == "Pancreas" & spot$ecotype %in% E6))))

# simple replacement: half the smallest positive proportion
Ms <- Pmat; Ms[Ms == 0] <- min(Pmat[Pmat > 0]) / 2; Ms <- Ms / rowSums(Ms)
Da_s <- euclid_mat(clr(Ms))
a1 <- permanova_exact(Da_s, des$site)
a2 <- permanova_exact(Da_s, des$site, cov = des$tf)
cat(sprintf("  simple repl   — site alone: F=%.2f R2=%.3f P=%.4f | site|tf: F=%.2f R2=%.3f P=%.4f\n",
            a1["F"], a1["R2"], a1["P"], a2["F"], a2["R2"], a2["P"]))
ok("Aitchison (simple) site F",       a1[["F"]],  3.86,   0.05)
ok("Aitchison (simple) site P",       a1[["P"]],  0.0039, 5e-4)
ok("Aitchison (simple) site|tf F",    a2[["F"]],  11.20,  0.20)
ok("Aitchison (simple) site|tf P",    a2[["P"]],  0.0039, 5e-4)

# multiplicative replacement, as used in RUN_analysis.R blocks 6-9
Mm <- mult_repl(Pmat, nvec, 0.65)
Da_m <- euclid_mat(clr(Mm))
b1 <- permanova_exact(Da_m, des$site)
b2 <- permanova_exact(Da_m, des$site, cov = des$tf)
cat(sprintf("  mult repl m=.65— site alone: F=%.2f R2=%.3f P=%.4f | site|tf: F=%.2f R2=%.3f P=%.4f\n",
            b1["F"], b1["R2"], b1["P"], b2["F"], b2["R2"], b2["P"]))
cat("  --> P must be 0.0039 in both. If the F values differ from the simple-\n")
cat("      replacement ones by more than ~0.3, report them and the reply text will\n")
cat("      be changed to quote the multiplicative version.\n")
ok("Aitchison (mult) site P",    b1[["P"]], 0.0039, 5e-4)
ok("Aitchison (mult) site|tf P", b2[["P"]], 0.0039, 5e-4)

## ================================================ I. ADJUSTED PER-CLASS =====
hdr("I. Adjusted per-class table (seven ecotypes, PT minus HM)")
cat("  Covariate = paired difference in mean tumour_frac over ALL spots, PT-HM.\n\n")

D7 <- paired(spot, E7, dir = "PT-HM")
dtf <- sapply(PATS, function(p)
  mean(spot$tumour_frac[spot$patient == p & spot$site == "Pancreas"]) -
  mean(spot$tumour_frac[spot$patient == p & spot$site == "Liver"]))
tab <- t(sapply(E7, function(e) {
  f <- site_effect(D7[, e], dtf); c(unadj = mean(D7[, e]), f)
}))
tab <- as.data.frame(tab); tab$BH <- p.adjust(tab$p, "BH")
print(round(tab[, c("unadj","effect","p","BH")], 4))
cat("\n")
ok("dtf minimum", min(dtf), -0.467, 0.002)
ok("dtf maximum", max(dtf),  0.506, 0.002)
for (e in list(c("CC5",20.18,0.0001,0.0004), c("CC2",-9.91,0.0066,0.0227),
               c("CC3",-13.16,0.0101,0.0227), c("CC8",-25.64,0.0130,0.0227),
               c("CC1",21.53,0.0183,0.0256),  c("CC4",13.06,0.1150,0.1350),
               c("CC7",-6.07,0.5710,0.5710))) {
  cc <- e[1]
  ok(paste(cc, "adjusted effect (pp)"), tab[cc,"effect"], as.numeric(e[2]), 0.02)
  ok(paste(cc, "adjusted P"),           tab[cc,"p"],      as.numeric(e[3]), 5e-4)
  ok(paste(cc, "BH-FDR"),               tab[cc,"BH"],     as.numeric(e[4]), 1e-3)
}

## ==================================================== J. SAMPLING TABLE =====
hdr("J. What each hepatic specimen contains")

samp <- t(sapply(PATS, function(p) {
  h <- spot[spot$patient == p & spot$site == "Liver", ]
  v <- table(factor(h$ecotype, levels = ALL10)) / nrow(h) * 100
  c(CC6_CC8 = sum(v[c("CC6","CC8")]), TME = sum(v[E6]), tumour = mean(h$tumour_frac))
}))
print(round(samp[order(-samp[,"TME"]), ], 3))
cat("\n")
ok("P8  CC6+CC8 (%)",  samp["PT_8","CC6_CC8"],  96.6, 0.1)
ok("P8  TME (%)",      samp["PT_8","TME"],       3.4, 0.1)
ok("P3  CC6+CC8 (%)",  samp["PT_3","CC6_CC8"],  92.1, 0.1)
ok("P10 CC6+CC8 (%)",  samp["PT_10","CC6_CC8"], 76.9, 0.1)
ok("P12 TME (%)",      samp["PT_12","TME"],     98.0, 0.1)
ok("P8  mean tumour fraction", samp["PT_8","tumour"],  0.096, 0.002)
ok("P12 mean tumour fraction", samp["PT_12","tumour"], 0.665, 0.002)

hp8 <- spot[spot$patient == "PT_8" & spot$site == "Liver" & spot$ecotype %in% E6, ]
ok("P8: CC2 as % of its captured TME", mean(hp8$ecotype == "CC2") * 100, 86, 1.5)

## ============================================ K. REPRODUCING LI & MA ========
hdr("K. Reproducing Li & Ma's own FDRs, and the denominator/family mismatch")

D10 <- paired(spot, E7, src = ALL10)          # proportions over all ten
p10 <- apply(D10, 2, wx)
q_family7  <- p.adjust(p10, "BH")             # family = seven  -> their numbers
cat("  denominator = all ten, family = seven (as their script does):\n")
print(round(rbind(P = p10, BH = q_family7), 4))
ok("CC4 FDR (letter: 0.130)", q_family7[["CC4"]], 0.130, 2e-3)
ok("CC7 FDR (letter: 0.910)", q_family7[["CC7"]], 0.910, 2e-3)

D10all <- paired(spot, ALL10)
q10 <- p.adjust(apply(D10all, 2, wx), "BH")
ok("CC4 FDR if corrected across all ten", q10[["CC4"]], 0.148, 2e-3)
ok("CC7 FDR if corrected across all ten", q10[["CC7"]], 1.000, 2e-3)

D7d <- paired(spot, E7)
q7 <- p.adjust(apply(D7d, 2, wx), "BH")
ok("CC4 FDR on the seven-ecotype denominator", q7[["CC4"]], 0.171, 2e-3)
ok("CC7 FDR on the seven-ecotype denominator", q7[["CC7"]], 0.820, 2e-3)

## ============================================= L. LYMPH NODE FLOOR =========
hdr("L. The four primary-lymph node pairs")

lnp <- sort(intersect(unique(spot$patient[spot$site == "Lymph node"]),
                      unique(spot$patient[spot$site == "Pancreas"])))
cat("  matched primary-LNM patients:", paste(lnp, collapse = ", "), "\n")
ok("number of primary-LNM pairs", length(lnp), 4, 0)
ok("smallest attainable two-sided exact P at n=4", 2/2^4, 0.125, 1e-9)
Dl <- t(sapply(lnp, function(p) {
  x <- spot[spot$ecotype %in% E6, ]
  a <- x[x$patient == p & x$site == "Lymph node", ]; b <- x[x$patient == p & x$site == "Pancreas", ]
  ta <- table(factor(a$ecotype, levels = E6))/nrow(a)*100
  tb <- table(factor(b$ecotype, levels = E6))/nrow(b)*100
  as.numeric(ta - tb)
}))
colnames(Dl) <- E6
ok("CC1 concordant in 4 of 4", sum(Dl[,"CC1"] < 0), 4, 0)
ok("CC5 concordant in 4 of 4", sum(Dl[,"CC5"] < 0), 4, 0)
ok("CC1 P (cannot beat 0.125)", wx(Dl[,"CC1"]), 0.125, 1e-6)

## ========================================== M. CROSS-CHECK WITH vegan =======
hdr("M. Independent cross-check of the PERMANOVAs against vegan::adonis2")

if (requireNamespace("vegan", quietly = TRUE)) {
  df <- data.frame(site = factor(des$site), tf = des$tf)
  bc <- as.dist(Db)
  P512 <- perm_matrix(9)
  v1 <- vegan::adonis2(bc ~ site, data = df, permutations = P512)
  cat("  adonis2(bc ~ site):\n"); print(v1)
  ok("vegan F agrees with ours",  v1$F[1],  r1[["F"]],  0.02)
  ok("vegan R2 agrees with ours", v1$R2[1], r1[["R2"]], 0.002)
  ## NOTE on the P value: vegan reports (count + 1) / (nperm + 1), so with all
  ## 512 restricted permutations it returns 3/513 = 0.00585 where the exact
  ## enumeration gives 2/512 = 0.0039. Both are the floor of the same test and
  ## neither is wrong; the reply text quotes the exact-enumeration convention.
  cat(sprintf("  vegan P = %.5f  (its (count+1)/(nperm+1) convention)\n",
              v1$`Pr(>F)`[1]))
  cat(sprintf("  ours    = %.5f  (exact enumeration, count/512)\n", r1[["P"]]))
  ok("vegan P is at its own floor", v1$`Pr(>F)`[1], 3/513, 1e-4)

  v2 <- vegan::adonis2(bc ~ tf + site, data = df, permutations = P512, by = "terms")
  cat("\n  adonis2(bc ~ tf + site, by = 'terms'):\n"); print(v2)
  ok("vegan partial F agrees with ours", v2$F[2], r2p[["F"]], 0.05)
  cat("  (vegan permutes all terms; the site row is the one to compare.)\n")
} else {
  cat("  vegan not installed - block skipped. Please install and rerun:\n")
  cat("      install.packages('vegan')\n")
  cat("  This is the only independent check on the PERMANOVA implementation.\n")
}

## ======================================================== SUMMARY ==========
hdr("SUMMARY")
if (fails == 0) {
  cat("  All checks passed. Every number quoted in the reply text is reproduced.\n")
} else {
  cat(sprintf("  *** %d CHECK(S) FAILED ***\n", fails))
  cat("  Send this console output back before the reply text is sent.\n")
}
cat("\n  Still outstanding, not covered here (they live in RUN_analysis.R\n")
cat("  blocks 6-9 and verify_ateeq_checks.R block J):\n")
cat("    - seven-part CLR sweep: CC2 0.5645 -> 0.1138, CC3 0.4137 -> 0.1172\n")
cat("    - six-part CLR: CC2 0.0076, CC3 0.0072 at m = 0.65 conservative\n")
cat("  Those are the only figures in the reply text not checked by this script.\n\n")
