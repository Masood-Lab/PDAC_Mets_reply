# Spatial transcriptomic analysis of primary and metastatic pancreatic cancers highlights tumor microenvironmental heterogeneity

Code supporting Khaliq et al., *Spatial
transcriptomic analysis of primary and metastatic pancreatic cancers
highlights tumor microenvironmental heterogeneity*, Nature Genetics 56,
2455-2465 (2024).

## Reproducing the analysis

Run `export_spot_table.R` once against the Seurat object to produce
`pdac_spot_level.csv.gz`. Every other script reads from that file and can
then be run independently, in any order.

```r
install.packages(c("data.table", "vegan", "ggplot2", "patchwork", "Seurat"))
```

## Verification scripts

- **verify_reply_numbers.R** -- reproduces every statistic in the reply
  text block by block: the paired Wilcoxon on both denominators, PERMANOVA
  under exact restricted permutation, the CC8 denominator ladder, the
  adjusted per-class regression, and an independent cross-check against
  `vegan::adonis2`.
- **verify_analysis_full.R** -- broader sweep covering per-class abundance,
  community-level PERMANOVA and every robustness variant (leave-one-
  class-out, leave-one-patient-out, conditioning on tumour fraction,
  Aitchison distance), the ISCHIA exclusion criteria applied independently
  to this cohort, the direction-of-shift contrast, and purity vs.
  enrichment.
- **verify_clr_analysis.R** -- the compositional (CLR) analysis, run first
  exactly as coded in the deposited reanalysis, then restricted to the six
  classes shared between primary and metastatic tissue.
- **permdisp_dispersion_test.R** -- the dispersion check reported in the
  reply text: confirms the PERMANOVA result reflects a genuine
  compositional shift rather than unequal within-group spread.

## Figure 1 (the figure enclosed with the reply)

**figure_direction_of_shift.R** produces the only titled, enclosed figure
in the reply text: the myCAF-dominated and tumour-dominated habitat
shares, and their contrast, patient by patient.

## Supplementary figures

Not individually titled or enclosed in the reply text, but each
visualizes an argument made in it and is included so the full analysis is
inspectable, not just the one submitted figure.

- **figure_two_denominators.R** -- the same paired test under three
  denominators (all ten classes, six shared classes plus CC8, six shared
  classes only), showing which per-class results are denominator-dependent
  and which are not.
- **figure_pcoa_two_geometries.R** -- ordination under Bray-Curtis and
  Aitchison distance, with the PERMANOVA result for each.
- **figure_purity_vs_enrichment.R** -- tumour-epithelial content (flat)
  against tumour-dominated habitat share (roughly doubles) -- the
  "purity vs. enrichment" distinction discussed in the reply.
- **figure_n9_correction_limit.R** -- the achievable exact p-values at
  nine matched pairs against the Benjamini-Hochberg step-up threshold.
- **figure_permanova_robustness.R** -- every version of the PERMANOVA test
  reported in the reply (the denominator ladder, leave-one-out variants,
  conditioning on tumour fraction, Aitchison distance) in one figure.

## A note on two numbers that took more than one attempt to get right

The partial PERMANOVA (site conditioned on tumour-epithelial fraction)
and the PERMDISP dispersion test both have the same issue: `vegan`'s own
reported p-value doesn't match full enumeration under a restricted
permutation scheme, even when given the exact permutation matrix by hand.
The F-statistics `vegan` returns are correct in both cases; the p-values
are not. Both scripts bypass `vegan`'s built-in significance test for
that specific line and compute the exact p-value directly by refitting
the model across all 512 valid arrangements and counting.
