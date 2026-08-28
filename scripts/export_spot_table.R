## export_spot_table.R
##
## Pulls a flat spot-level table out of the Seurat object for downstream
## analysis: one row per spot, with patient, section, anatomical site,
## ecotype assignment, and the RCTD-derived tumour/hepatocyte/myCAF
## fractions the reply's statistics run on.
##
## ecotype is read from cc_ischia_10 -- the k=10 clustering, matching the
## published elbow-method choice. CompositionCluster_CC looks like the
## obvious column to use but isn't: it's a relabeling of cc_ischia_18,
## an 18-cluster run, and would silently give the wrong ecotype for
## every spot if used here instead.
##
## site comes straight from Origin -- level counts match Fig. 1a exactly
## (Liver 28,520; Lymph node 17,698; Normal Pancreas 9,820; Pancreas
## 35,458), so no recoding needed.

library(Seurat)

pdac <- readRDS("/Users/akhaliq/Desktop/rebuttle_natgen/Pdac_allres_most_updated_V2.rds")

md <- pdac@meta.data

# --- pull RCTD proportions for the two features verify_reply_numbers.R needs ---
rctd_mat <- GetAssayData(pdac, assay = "rctd_fullfinal", slot = "data")
stopifnot(all(c("Tumor Epithelial cells", "Hepatocytes") %in% rownames(rctd_mat)))

# colnames(rctd_mat) carries an extra `names` attribute (visible in str() as
# "Named chr") that identical() treats as a mismatch even when the barcode
# strings themselves are in identical order. Strip it before comparing.
rctd_colnames <- unname(colnames(rctd_mat))
md_rownames   <- unname(rownames(md))

if (!identical(rctd_colnames, md_rownames)) {
  # genuine order/content mismatch -- do NOT proceed silently, diagnose instead
  cat("Length rctd_mat cols:", length(rctd_colnames), " Length md rows:", length(md_rownames), "\n")
  cat("Set difference (in rctd, not in md):", length(setdiff(rctd_colnames, md_rownames)), "\n")
  cat("Set difference (in md, not in rctd):", length(setdiff(md_rownames, rctd_colnames)), "\n")
  cat("First few rctd colnames:\n"); print(head(rctd_colnames))
  cat("First few md rownames:\n"); print(head(md_rownames))
  stop("Spot order/content mismatch between rctd_fullfinal and meta.data -- inspect output above before proceeding.")
} else {
  cat("Spot order confirmed identical after stripping names attribute: OK\n")
}

tumour_frac <- as.numeric(rctd_mat["Tumor Epithelial cells", ])
hep_frac    <- as.numeric(rctd_mat["Hepatocytes", ])
myCAF_frac  <- as.numeric(rctd_mat["myCAF", ])
## Added 2026-08-27 per Ashiq's audit: two numbers in the letter (CC1/CC5 at
## 0.63/0.48 myCAF, CC7 at 0.20 myCAF in the CC6->CC8->CC7 continuum) and
## make_all_figures.py both read this column. treatment and iCAF_frac are
## NOT needed -- no script reads them and no letter number depends on them
## (PT_9 is excluded by patient ID in code, not via the treatment column).

spot_export <- data.frame(
  spot        = rownames(md),
  patient     = md$patient,
  section     = md$orig.ident,
  site        = as.character(md$Origin),
  ecotype     = as.character(md$cc_ischia_10),
  tumour_frac = tumour_frac,
  hep_frac    = hep_frac,
  myCAF_frac  = myCAF_frac,
  stringsAsFactors = FALSE
)

## --- sanity checks before writing ---
cat("Rows:", nrow(spot_export), "\n")
cat("Site table:\n"); print(table(spot_export$site))
cat("Ecotype table:\n"); print(table(spot_export$ecotype))
cat("Unique patients:", length(unique(spot_export$patient)), "\n")
print(sort(unique(spot_export$patient)))

# these nine must all be present for verify_reply_numbers.R to run
PATS <- c("PT_2","PT_3","PT_4","PT_6","PT_8","PT_9","PT_10","PT_11","PT_12")
missing_pats <- setdiff(PATS, unique(spot_export$patient))
if (length(missing_pats) > 0) {
  stop("Missing expected matched patients: ", paste(missing_pats, collapse = ", "))
} else {
  cat("All 9 matched primary-hepatic patients present: OK\n")
}

# check for NAs in the two numeric columns -- would break comp()/wx() silently
cat("NAs in tumour_frac:", sum(is.na(spot_export$tumour_frac)), "\n")
cat("NAs in hep_frac:", sum(is.na(spot_export$hep_frac)), "\n")
cat("NAs in myCAF_frac:", sum(is.na(spot_export$myCAF_frac)), "\n")

## --- write ---
write.csv(spot_export, gzfile("pdac_spot_level.csv.gz"), row.names = FALSE)
cat("\nWritten: pdac_spot_level.csv.gz\n")
