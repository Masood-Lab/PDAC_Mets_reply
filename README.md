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

## Why we report exact enumeration rather than vegan's permutation p-value

Two results in the reply -- the partial PERMANOVA (site conditioned on
tumour-epithelial fraction) and the PERMDISP dispersion test -- are reported
from full enumeration of the restricted permutation set rather than from
`vegan`'s own significance test. The F-statistics are identical either way.
The reason is a difference in bookkeeping convention, not an error in
`vegan` -- but the exact mechanism is more specific than "supply 511 instead
of 512" or vice versa. It goes the *opposite* way for the two PERMANOVA
models used in this reply, and both directions are worth stating precisely
so anyone re-deriving this lands on the same numbers.

**The permutation floor at nine matched pairs.** Permuting site labels
within patient gives 2^9 = 512 distinct arrangements. Exactly two of them
reach the observed F for the site-only model: the identity arrangement and
its complete complement. Full enumeration returns 2/512 = 0.0039, the
smallest value this design admits.

**Simple model (`bc ~ site`).** Supplying `adonis2` the full 512-row matrix,
identity included, gives P = 0.005848 (3/513) -- `adonis2` adds its own "+1"
for the observed statistic on top of a matrix that already contains it,
double-counting the identity row. Supplying only the 511 non-identity
arrangements removes that double count and returns 2/512 = 0.003906 exactly,
matching hand enumeration. Confirmed directly:

    adonis2(D ~ site, permutations = P512)  ->  P = 0.005848
    adonis2(D ~ site, permutations = P511)  ->  P = 0.003906  (correct)

**Partial model (`bc ~ tf + site, by = "terms"`) -- the opposite fix.**
Tested the same way, both directions:

    adonis2(D ~ tf + site, permutations = P512, by="terms")  ->  P = 0.003899  (correct)
    adonis2(D ~ tf + site, permutations = P511, by="terms")  ->  P = 0.001953  (wrong)

For this model the full 512-row matrix, identity included, gives the exact
answer; excluding identity undercounts by one hit and returns 1/512 instead
of 2/512. This is the reverse of the simple-model case. `adonis2`'s
row-recycling for a sequential (`by="terms"`) two-term model evidently
handles the observed-vs-permuted accounting differently once a covariate
term precedes the term of interest, and empirically that difference exactly
cancels the "+1" convention that requires the 511-only fix for the single-
term model. We have not traced this to a specific line in `vegan`'s source;
what's stated above is what direct testing shows, not a claim about why.

Both directions are independently reproducible from this repository
(`verify_analysis_full.R`, `verify_reply_numbers.R` Section M) and were
confirmed on a second machine before being written up here. **If you plan to
fix this by feeding `adonis2` a reduced permutation matrix, check which of
the two directions above applies to your specific formula** -- it is not
the same fix for a single-term and a sequential two-term model.

**PERMDISP.** `permutest.betadisper` permutes distances-to-centroid with
group centroids held fixed from the observed fit; `permdisp_dispersion_test.R`
instead refits `betadisper()` under each permuted labelling, so centroids
are recomputed each time. This is the more conservative of the two and is
what we report. Both approaches agree on the substantive conclusion here:
all four tests (Bray-Curtis and Aitchison, centroid and spatial median) are
non-significant.
