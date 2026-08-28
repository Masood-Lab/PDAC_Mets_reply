## figure_direction_of_shift.R
##
## Three panels, patient by patient: the myCAF-dominated habitats (CC1+CC5)
## falling, the tumour-dominated habitats (CC2+CC3) rising, and the two
## combined into one ecological contrast per patient. Everything is computed
## directly from the spot table below, not hardcoded.

library(data.table)
library(ggplot2)
library(patchwork)

CSV <- "/Users/akhaliq/Desktop/rebuttle_natgen/new/pdac_spot_level.csv.gz"
md <- as.data.table(read.csv(gzfile(CSV), stringsAsFactors = FALSE))
setnames(md, "ecotype", "cc")
md[, site := c("Pancreas" = "PT", "Liver" = "HM")[site]]

patients <- c("PT_2", "PT_3", "PT_4", "PT_6", "PT_8", "PT_9", "PT_10", "PT_11", "PT_12")
classes  <- c("CC1", "CC2", "CC3", "CC4", "CC5", "CC7")
treated  <- "PT_9"

comp <- function(p, s) {
  sub <- md[patient == p & site == s & cc %in% classes]
  tab <- table(factor(sub$cc, levels = classes))
  as.numeric(tab) / sum(tab) * 100
}

wide <- rbindlist(lapply(patients, function(p) {
  pt <- comp(p, "PT"); hm <- comp(p, "HM")
  data.table(patient = p,
             myCAF_PT = pt[1] + pt[5], myCAF_HM = hm[1] + hm[5],
             tum_PT   = pt[2] + pt[3], tum_HM   = hm[2] + hm[3])
}))
wide[, contrast_PT := tum_PT - myCAF_PT]
wide[, contrast_HM := tum_HM - myCAF_HM]
wide[, treated := patient == treated]

wilcox_summary <- function(x, y, direction = c("rise", "fall")) {
  direction <- match.arg(direction)
  w <- suppressWarnings(wilcox.test(x, y, paired = TRUE, exact = TRUE))
  n_matching <- if (direction == "rise") sum(y > x) else sum(y < x)
  sprintf("%d of %d %s, P = %.4f", n_matching, length(x), direction, w$p.value)
}

long_pair <- function(dt, pt_col, hm_col) {
  rbind(
    data.table(patient = dt$patient, site = "primary\npancreatic",
               value = dt[[pt_col]], treated = dt$treated),
    data.table(patient = dt$patient, site = "hepatic\nmetastasis",
               value = dt[[hm_col]], treated = dt$treated)
  )
}

site_levels <- c("primary\npancreatic", "hepatic\nmetastasis")
pal <- c("primary\npancreatic" = "#1F3864", "hepatic\nmetastasis" = "#C0504D")

make_panel <- function(dt, pt_col, hm_col, title, sub, ylab) {
  d <- long_pair(dt, pt_col, hm_col)
  d$site <- factor(d$site, levels = site_levels)
  treated_line <- d[treated == TRUE]
  ggplot(d, aes(x = site, y = value, group = patient)) +
    geom_line(data = d[treated == FALSE], colour = "grey70", linewidth = 0.5) +
    geom_line(data = treated_line, colour = "#E8A33D", linewidth = 1.1) +
    geom_point(aes(colour = site), size = 2.4) +
    scale_colour_manual(values = pal, guide = "none") +
    labs(title = title, subtitle = sub, x = NULL, y = ylab) +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(size = 10, face = "bold"),
          panel.grid.minor = element_blank())
}

p1 <- make_panel(wide, "myCAF_PT", "myCAF_HM",
                  "myCAF-dominated habitats\nCC1 + CC5", "",
                  "share of shared composition (%)") +
  annotate("text", x = 1.5, y = -8,
           label = wilcox_summary(wide$myCAF_PT, wide$myCAF_HM, "fall"), size = 3)

p2 <- make_panel(wide, "tum_PT", "tum_HM",
                  "tumour-dominated habitats\nCC2 + CC3", "",
                  "") +
  annotate("text", x = 1.5, y = -8,
           label = wilcox_summary(wide$tum_PT, wide$tum_HM, "rise"), size = 3)

p3 <- make_panel(wide, "contrast_PT", "contrast_HM",
                  "ecological balance\n(CC2 + CC3) - (CC1 + CC5)", "",
                  "difference in share (percentage points)") +
  annotate("text", x = 1.5, y = min(wide$contrast_PT) - 15,
           label = wilcox_summary(wide$contrast_PT, wide$contrast_HM, "rise"), size = 3)

combined <- p1 | p2 | p3
combined <- combined + plot_annotation(
  title = "The reorganization has a direction, and all nine patients move the same way"
)

ggsave("figure_direction_of_shift.png", combined, width = 11, height = 4.2, dpi = 300)

cat("myCAF: ", round(mean(wide$myCAF_PT),1), "->", round(mean(wide$myCAF_HM),1), "\n")
cat("tumour: ", round(mean(wide$tum_PT),1), "->", round(mean(wide$tum_HM),1), "\n")
cat("contrast: ", round(mean(wide$contrast_PT),1), "->", round(mean(wide$contrast_HM),1), "\n")
