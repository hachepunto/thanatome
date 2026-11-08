# thanatome

Donor death signatures for post-mortem control tissue.

`thanatome` annotates each sample of a post-mortem tissue cohort (GTEx or similar) with its position along three transcriptional axes of donor death derived from GTEx v8 in ten tissues: **ventilator** (perfused versus post-mortem), post-mortem **ischemia** and prolonged **agony**. For every sample it returns a signature score per axis, a position scaled between the two reference poles of the same tissue, the fraction of the sample's deviation from the tissue centroid that lies along the death axes, a state label and technical contamination flags. Samples can then be plotted, selected or excluded with explicit criteria, and the scores can be used as covariates in differential-expression designs or batch-correction tools. A second group of functions annotates an existing differential-expression table with the axes its genes fall in.

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
data(thanatome_example)                       # 67 GTEx v8 ovary samples
sc <- thanatome_score(thanatome_example, reference = "ovary")
plot_thanatome(sc, hardy = thanatome_example$hardy)
sel <- select_samples(sc, agony = c(-Inf, 0.5))   # explicit criterion, no default
cv <- death_covariates(sc)                    # for limma::removeBatchEffect or a design matrix
```

## What the analyses behind the package say

- The three axes reproduce across tissues and predict, for each sample, whether the donor died on a ventilator (AUC 0.87-0.94 with signatures derived without that tissue) and whether the death was slow (AUC 0.69-0.90).
- Seen from outside GTEx, the ventilator axis is **perfused versus post-mortem**: living colon biopsies fall beyond the ventilated pole and accident autopsies without intensive care fall at the non-ventilated pole. Ventilated donors are the post-mortem samples closest to living and to surgical tissue; "no ICU" is not the clean control.
- The axes explain a small fraction of each sample's variance (under 1 percent of the whole transcriptome in median with the signature genes) and do not explain the difference between GTEx and tumour-adjacent normal tissue. The package annotates and quantifies; it does not promise to correct death.
- There is no default "acceptable" criterion: no subset of post-mortem controls is closer to surgical tissue than the full set except, slightly, the ventilated donors, and excluding samples costs power. `select_samples()` makes the criteria explicit and reports their cost.

The analysis repository is `hachepunto/gtex-controls` (private until publication); `data-raw/` documents how every data object was built from it.

## Provenance

The code, documentation and tests of this package were written with the assistance of Claude Code (Anthropic), under the direction and review of the authors, following the Bioconductor policy on AI-assisted contributions. Every file produced with assistance carries a note at its top. The analyses that define the signatures and the numbers quoted here were run by the authors in the analysis repository.
