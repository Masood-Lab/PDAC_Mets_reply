## figure_purity_vs_enrichment.R
##
## Two ways of describing the same nine matched pairs: how much tumour is
## present (mean tumour-epithelial fraction across the six shared
## habitats), versus how the tissue is organized (share of the composition
## occupied by the two tumour-dominated habitats, CC2+CC3). The first is
## flat; the second roughly doubles. Computed directly from the spot table.

library(data.table)
library(ggplot2)
library(patchwork)

CSV <- "/Users/akhaliq/Desktop/rebuttle_natgen/new/pdac_spot_level.csv.gz"
md <- as.data.table(read.csv(gzfile(CSV), stringsAsFactors = FALSE))
setnames(md, "ecotype", "cc")
md[, site := c("Pancreas" = "PT", "Liver" = "HM")[site]]

patients <- c("PT_2", "PT_3", "PT_4", "PT_6", "PT_8", "PT_9", "PT_10", "PT_11", "PT_12")
classes  <- c("CC1", "CC2", "CC3", "CC4", "CC5", "CC7")

sub6 <- md[patient %in% patients & site %in% c("PT","HM") & cc %in% classes]

## panel 1: mean tumour-epithelial fraction within the six shared habitats
purity <- sub6[, .(tf = mean(tumour_frac)), by = .(patient, site)]
purity_wide <- dcast(purity, patient ~ site, value.var = "tf")
w1 <- suppressWarnings(wilcox.test(purity_wide$PT, purity_wide$HM, paired = TRUE, exact = TRUE))

## panel 2: share of the six-habitat composition occupied by CC2+CC3
comp <- function(p, s) {
  s2 <- sub6[patient == p & site == s]
  tab <- table(factor(s2$cc, levels = classes))
  as.numeric(tab) / sum(tab)
}
share <- rbindlist(lapply(patients, function(p) {
  pt <- comp(p, "PT"); hm <- comp(p, "HM")
  data.table(patient = p, PT = pt[2] + pt[3], HM = hm[2] + hm[3])
}))
w2 <- suppressWarnings(wilcox.test(share$PT, share$HM, paired = TRUE, exact = TRUE))

long <- function(dt, valcol_pt = "PT", valcol_hm = "HM") {
  rbind(data.table(patient = dt$patient, site = "primary", value = dt[[valcol_pt]]),
        data.table(patient = dt$patient, site = "hepatic\nmetastasis", value = dt[[valcol_hm]]))
}

d1 <- long(purity_wide); d1$site <- factor(d1$site, levels = c("primary", "hepatic\nmetastasis"))
d2 <- long(share); d2$site <- factor(d2$site, levels = c("primary", "hepatic\nmetastasis"))

p1 <- ggplot(d1, aes(site, value, group = patient)) +
  geom_line(colour = "grey70") + geom_point(colour = "#1F3864", size = 2.4) +
  labs(title = "How much tumour is present",
       subtitle = sprintf("mean tumour-epithelial fraction\n%d of 9 higher in metastasis, P = %.4f",
                           sum(purity_wide$HM > purity_wide$PT), w1$p.value),
       x = NULL, y = "tumour-epithelial fraction") +
  theme_minimal(base_size = 10) +
  theme(plot.subtitle = element_text(size = 9))

p2 <- ggplot(d2, aes(site, value, group = patient)) +
  geom_line(colour = "grey70") + geom_point(colour = "#C0504D", size = 2.4) +
  labs(title = "How the tissue is organized",
       subtitle = sprintf("tumour-dominated habitat share (CC2+CC3)\n%d of 9 higher, P = %.4f",
                           sum(share$HM > share$PT), w2$p.value),
       x = NULL, y = "share of composition") +
  theme_minimal(base_size = 10) +
  theme(plot.subtitle = element_text(size = 9))

combined <- (p1 | p2) + plot_annotation(
  title = "The seed is much the same; the soil is not"
)
ggsave("figure_purity_vs_enrichment.png", combined, width = 9, height = 4.3, dpi = 300)

cat("purity PT/HM means:", round(mean(purity_wide$PT),4), round(mean(purity_wide$HM),4),
    " P =", w1$p.value, "\n")
cat("share PT/HM means:", round(mean(share$PT),4), round(mean(share$HM),4),
    " P =", w2$p.value, "\n")
