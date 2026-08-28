## figure_two_denominators.R
##
## Same test, same nine patients, three denominators: proportion of every
## spot in the specimen (all ten classes present), a seven-class
## subcomposition that keeps CC8 but drops the classes with no counterpart
## in a primary (CC6, CC9, CC10), and the six-class subcomposition shared
## by both arms. Everything computed from the spot table directly, not
## hardcoded.

library(data.table)
library(ggplot2)

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
md <- as.data.table(read.csv(gzfile(CSV), stringsAsFactors = FALSE))
setnames(md, "ecotype", "cc")
md[, site := c("Pancreas" = "PT", "Liver" = "HM")[site]]

patients <- c("PT_2", "PT_3", "PT_4", "PT_6", "PT_8", "PT_9", "PT_10", "PT_11", "PT_12")
test_classes <- c("CC1", "CC2", "CC3", "CC4", "CC5", "CC7", "CC8")

## denom_classes = NULL means normalize against every spot in the specimen
## (all ten classes); otherwise normalize against just the listed classes
## (a renormalized subcomposition).
paired_test <- function(denom_classes, label) {
  s <- md[patient %in% patients & site %in% c("PT", "HM")]
  if (!is.null(denom_classes)) s <- s[cc %in% denom_classes]
  counts <- s[, .N, by = .(patient, site, cc)]
  counts[, prop := N / sum(N), by = .(patient, site)]
  wide <- dcast(counts, patient + site ~ cc, value.var = "prop", fill = 0)
  wide <- wide[order(match(patient, patients))]
  res <- rbindlist(lapply(intersect(test_classes, names(wide)), function(cl) {
    pt <- wide[site == "PT"][[cl]] * 100
    hm <- wide[site == "HM"][[cl]] * 100
    w <- suppressWarnings(wilcox.test(pt, hm, paired = TRUE, exact = TRUE))
    data.table(cc = cl, mean_diff = mean(pt - hm), P = w$p.value)
  }))
  res[, BH := p.adjust(P, "BH")]
  res[, denom := label]
  res
}

all_spots <- paired_test(NULL, "All ten classes")
seven     <- paired_test(c("CC1","CC2","CC3","CC4","CC5","CC7","CC8"), "Six shared + CC8")
six       <- paired_test(c("CC1","CC2","CC3","CC4","CC5","CC7"), "Six habitats present in both arms")

combined <- rbind(all_spots, seven, six)
combined[, denom := factor(denom, levels = c("All ten classes", "Six shared + CC8",
                                              "Six habitats present in both arms"))]
combined[, cc := factor(cc, levels = rev(c("CC1","CC5","CC4","CC2","CC3","CC7","CC8")))]
combined[, sig := BH < 0.05]

ggplot(combined, aes(x = mean_diff, y = cc, colour = sig)) +
  geom_vline(xintercept = 0, colour = "grey50") +
  geom_point(size = 3) +
  geom_text(aes(label = sprintf("FDR %.4f", BH)), hjust = -0.1, size = 2.9, colour = "black") +
  scale_colour_manual(values = c("TRUE" = "#1F3864", "FALSE" = "grey60"), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.08, 0.55))) +
  facet_wrap(~ denom, nrow = 1) +
  labs(x = "mean paired difference in ecotype proportion, PT - HM (points)",
       y = NULL,
       title = "The same test, the same patients, three denominators") +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"))

ggsave("figure_two_denominators.png", width = 13, height = 4, dpi = 300)

cat("\n--- all-ten (raw proportion of every spot) ---\n"); print(all_spots[order(cc)])
cat("\n--- six shared + CC8 ---\n"); print(seven[order(cc)])
cat("\n--- six habitats ---\n"); print(six[order(cc)])
