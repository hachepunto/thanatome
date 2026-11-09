# thanatome

Donor death signatures for post-mortem control tissue.

`thanatome` annotates each sample of a post-mortem tissue cohort (GTEx or similar) with its position along two transcriptional axes of donor death derived from GTEx v8 in ten tissues, **ventilator** (perfused versus post-mortem) and prolonged **agony**, plus the **ischemia** component of the ventilator axis. For every sample it returns a signature score per axis, a position scaled between the two reference poles of the same tissue, a state label and technical contamination flags, and optionally (joint regression) the fraction of the sample's deviation from the tissue centroid that lies along the axes. Samples can then be plotted, selected with explicit criteria, and the scores used as covariates in differential-expression designs or batch-correction tools. `compare_controls()` repeats a case-versus-control differential expression with each candidate control subset against random subsets of the same size; `annotate_de()` and `annotate_mra()` annotate existing tables with the axes their genes fall in.

## Installation

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("thanatome")
## development version
BiocManager::install("hachepunto/thanatome")
```

## Quick start

```r
library(thanatome)
data(thanatome_example)                       # 74 GTEx ovary controls + 40 TCGA-OV tumours (recount3)
ctl <- thanatome_example[, thanatome_example$group == "control"]
sc <- thanatome_score(ctl, reference = "ovary")
plot_thanatome(sc, hardy = ctl$hardy)
sel <- select_samples(sc, agony = c(-Inf, 0.5))   # explicit criterion, no default
cv <- death_covariates(sc)                    # for limma::removeBatchEffect or a design matrix
cc <- compare_controls(thanatome_example, group = "group", control = "control",
                       criteria = list(no_agony = list(preset = "no_agony")), reference = "ovary")
audit_report(sc, selection = sel, comparison = cc)
```

## What the analyses behind the package say

- The two axes reproduce across tissues and predict, for each sample, whether the donor died on a ventilator (AUC 0.87-0.94 with signatures derived without that tissue) and whether the death was slow (AUC 0.69-0.90).
- Seen from outside GTEx, the ventilator axis is **perfused versus post-mortem**: living colon biopsies fall beyond the ventilated pole and accident autopsies without intensive care fall at the non-ventilated pole. Ventilated donors are the post-mortem samples closest to living and to surgical tissue; "no ICU" is not the clean control.
- Ischemia is not an independent axis: its signature is mostly the ventilator signature with the opposite sign, and it tracks ischemic time only among perfused donors (partial correlation 0.20-0.50 within GTEx Hardy 0).
- The axes explain a small fraction of each sample's variance (under 1 percent of the whole transcriptome in median with the signature genes) and do not explain the difference between GTEx and tumour-adjacent normal tissue. The package annotates and quantifies; it does not promise to correct death.
- There is no default "acceptable" criterion: in GTEx no subset of post-mortem controls is closer to surgical tissue than the full set except, slightly, the ventilated donors, and excluding samples costs power. `select_samples()` makes the criteria explicit and reports their cost.

The packaged references were built from recount3 (Monorail, GENCODE v26) and apply to data quantified the same way; for other pipelines pass your own cohort as reference. `inst/scripts/` documents how every data object was built, with the source of each dataset; the analysis repository is `hachepunto/thanatome-manuscript`.

## Provenance

The R code, tests, vignette and data scripts of this package were written with the assistance of Claude Code (Anthropic), under the direction and review of the authors, following the Bioconductor policy on AI-assisted contributions. Each of those files carries a note at its top. The analyses that define the signatures and the numbers quoted here were run by the authors in the analysis repository.
