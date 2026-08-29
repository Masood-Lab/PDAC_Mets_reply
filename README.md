# Spatial transcriptomic analysis of primary and metastatic pancreatic cancers highlights tumor microenvironmental heterogeneity

Code supporting the analysis published in Nature Genetics 56, 2455-2465
(2024).

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

## Why we report exact enumeration rather than vegan's permutation p-value

One result in this repository -- the community-level PERMANOVA (`bc ~
site`) -- is reported from full enumeration of the restricted permutation
set rather than from `vegan::adonis2`'s own significance test. The
F-statistic is identical either way; only the p-value convention differs,
and the difference is fully accounted for.

**The permutation floor at nine matched pairs.** Permuting site labels
within patient gives 2^9 = 512 distinct arrangements. Exactly two of them
reach the observed F: the identity arrangement and its complete
complement. Full enumeration returns 2/512 = 0.0039, the smallest value
this design admits.

**Why `adonis2` reports a different number.** `adonis2` uses the
convention `(count + 1) / (nperm + 1)`. Supplied the full 512-row
permutation matrix (identity included), it returns 3/513 = 0.00585 --
because the matrix already contains the identity arrangement and
`adonis2` adds its own "+1" for the observed statistic on top of that,
counting the identity row twice. This is demonstrated directly in
`verify_reply_numbers.R`, Section M:

    adonis2(bc ~ site, permutations = P512)   ->  P = 0.00585  (3/513)
    exact enumeration (this repository)       ->  P = 0.00391  (2/512)

Both are the floor of the same test under two different conventions;
neither is wrong. Section M's own comment states this directly, and the
check `ok("vegan P is at its own floor", ..., 3/513, ...)` asserts it in
code, so this is reproducible by running that script, not just claimed
here.

**A note on scope.** We looked into whether a similarly precise account
could be given for the partial PERMANOVA (site conditioned on
tumour-epithelial fraction), which this repository also reports via exact
enumeration rather than `adonis2`'s own test. An earlier version of this
section made a specific claim about that case; on review, the claim was
not fully supported by what is actually checked in the deposited scripts,
so it has been removed rather than repeated. The partial-model p-value
reported in the reply (0.0039) comes from direct enumeration in
`verify_analysis_full.R`, independent of `adonis2` entirely -- that
computation stands on its own and does not depend on the withdrawn
claim.

**PERMDISP.** `permutest.betadisper` also does not return the same
p-value as full enumeration under this restricted permutation scheme,
even when handed the exact permutation matrix directly -- the F-statistic
agrees, the p-value does not. Unlike the partial-PERMANOVA case above,
this is demonstrated directly in `permdisp_dispersion_test.R` (not just
asserted): the script includes a direct `vegan::permutest.betadisper`
call on the same 512-arrangement set as a side-by-side comparison. The
p-values used in the reply come from refitting `betadisper()` under each
of the 512 valid arrangements directly and counting, not from
`permutest()`.
