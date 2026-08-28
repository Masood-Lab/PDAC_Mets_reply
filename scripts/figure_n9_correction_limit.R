## figure_n9_correction_limit.R
##
## At nine pairs, the exact signed-rank test only returns a handful of
## possible p-values (0.0039, 0.0078, 0.0117, ...), and correcting across
## six or seven dependent classes can move a result across the nominal
## threshold depending on how many other classes also reach significance.
## This plots each class's exact p-value against the Benjamini-Hochberg
## step-up threshold for the six-class family, computed directly from data.

library(data.table)
library(ggplot2)

CSV <- "/Users/akhaliq/Desktop/rebuttle_natgen/new/pdac_spot_level.csv.gz"
md <- as.data.table(read.csv(gzfile(CSV), stringsAsFactors = FALSE))
setnames(md, "ecotype", "cc")
md[, site := c("Pancreas" = "PT", "Liver" = "HM")[site]]

patients <- c("PT_2", "PT_3", "PT_4", "PT_6", "PT_8", "PT_9", "PT_10", "PT_11", "PT_12")
classes  <- c("CC1", "CC2", "CC3", "CC4", "CC5", "CC7")

s <- md[patient %in% patients & site %in% c("PT","HM") & cc %in% classes]
counts <- s[, .N, by = .(patient, site, cc)]
counts[, prop := N / sum(N), by = .(patient, site)]
wide <- dcast(counts, patient + site ~ cc, value.var = "prop", fill = 0)
wide <- wide[order(match(patient, patients))]

res <- rbindlist(lapply(classes, function(cl) {
  pt <- wide[site == "PT"][[cl]] * 100
  hm <- wide[site == "HM"][[cl]] * 100
  w <- suppressWarnings(wilcox.test(pt, hm, paired = TRUE, exact = TRUE))
  data.table(cc = cl, P = w$p.value)
}))
res <- res[order(P)]
res[, rank := .I]
n <- nrow(res)
res[, bh_threshold := (rank / n) * 0.05]
res[, BH := p.adjust(P, "BH")]

ggplot(res, aes(x = rank, y = P)) +
  geom_step(aes(y = bh_threshold), colour = "grey50", linetype = "dashed", direction = "hv") +
  geom_point(size = 3, colour = "#C0504D") +
  geom_text(aes(label = sprintf("%s\nBH %.4f", cc, BH)), vjust = -0.7, size = 2.8) +
  scale_y_log10(expand = expansion(mult = c(0.05, 0.35))) +
  labs(title = "At nine pairs, the correction decides the verdict as much as the data",
       subtitle = "exact signed-rank P against the Benjamini-Hochberg step-up threshold, six-class family",
       x = "rank of the class among the six tested",
       y = "P value (log scale)") +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank())

ggsave("figure_n9_correction_limit.png", width = 8, height = 5.5, dpi = 300)

print(res)
