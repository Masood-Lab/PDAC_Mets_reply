## figure_permanova_robustness.R
##
## One figure summarizing every version of the PERMANOVA test we ran on
## the primary-vs-metastasis comparison: the full denominator ladder
## (unpaired ten-class down to the paired six-class subcomposition),
## leave-one-class-out, leave-one-patient-out, conditioning on
## tumour-epithelial fraction, and Aitchison distance instead of
## Bray-Curtis. Numbers come straight from the verification output, not
## recomputed here.

library(ggplot2)

results <- data.frame(
  group = c(
    "Main result",
    rep("Denominator ladder", 4),
    rep("Leave-one-class-out", 6),
    rep("Leave-one-patient-out", 2),
    rep("Alternate conditioning", 2)
  ),
  label = c(
    "Six shared classes (n = 9 pairs)",
    "A: unpaired, all 10 classes",
    "B: 9 pairs, all 10 classes",
    "C: 9 pairs, 7 classes (CC8 retained)",
    "D: 9 pairs, 6 shared classes (reported)",
    "drop CC1", "drop CC2", "drop CC3", "drop CC4", "drop CC5", "drop CC7",
    "minimum F (PT_12)", "maximum F (PT_8)",
    "conditioned on tumour fraction",
    "Aitchison distance"
  ),
  F = c(3.88,
        8.46, 6.95, 6.00, 3.88,
        3.89, 3.69, 4.69, 4.50, 2.38, 3.50,
        2.80, 4.09,
        6.83, 3.86),
  P = c(0.0039,
        0.002, 0.0039, 0.0039, 0.0039,
        0.0039, 0.0039, 0.0039, 0.0039, 0.0117, 0.0039,
        0.0078, 0.0078,
        0.0039, 0.0039)
)

results$group <- factor(results$group, levels = rev(unique(results$group)))
results$label <- factor(results$label, levels = rev(results$label))
results$at_floor <- ifelse(results$P <= 0.0039 + 1e-9, "at permutation floor",
                     ifelse(results$P <= 0.0117, "significant, above floor", "n.s. at this n"))

p <- ggplot(results, aes(x = label, y = F, fill = at_floor)) +
  geom_col(width = 0.65, colour = "black", linewidth = 0.3) +
  geom_text(aes(label = sprintf("F = %.2f, P = %.4f", F, P)),
            hjust = -0.05, size = 3.1) +
  coord_flip(clip = "off") +
  facet_grid(rows = vars(group), scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_manual(values = c("at permutation floor" = "#1E7B34",
                                "significant, above floor" = "#4C72B0",
                                "n.s. at this n" = "#C0504D")) +
  scale_y_continuous(limits = c(0, 11), expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "PERMANOVA pseudo-F", fill = NULL,
       title = "PERMANOVA site effect under every perturbation applied") +
  theme_minimal(base_size = 11) +
  theme(
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    plot.margin = margin(5, 60, 5, 5)
  )

ggsave("figure_permanova_robustness.png", p, width = 9, height = 7, dpi = 300)
