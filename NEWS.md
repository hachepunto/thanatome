# thanatome 0.99.0

- First version submitted to Bioconductor.
- `thanatome_score()`: per-sample signature score, scaled position, state label and technical flags along the ventilator and agony axes, plus the ischemia component of the ventilator axis, with a GTEx tissue (recount3 quantification) or a cohort as reference; `method = "regression"` adds the joint coefficients and the fraction of deviation explained.
- `thanatome_score()` leaves out the genes a cohort does not express (`drop_unexpressed`, cohorts of at least 5 samples) and can quantile-normalise external cohorts to the reference tissue (`normalize = "quantile"`); `thanatome_references` carries the statistics of the quantile-normalised GTEx samples (`qn`).
- `thanatome_score()` and `contamination_flags()` take `tissue`, which switches off the flag that is physiological in that tissue (surfactant in lung, pancreatic enzymes in pancreas), also with a cohort reference. Cohort references check their poles (`min_pole`, default 8 samples per pole, and an internal AUC >= 0.6); unusable poles give `NA` positions with a warning and are recorded in `metadata()$poles`.
- `compare_controls()` takes `design` (covariates of the user's DE, the group term added last), `BPPARAM` (BiocParallel; random subsets are drawn before fitting, so results do not depend on the backend) and a `DGEList`.
- `select_samples()`: `preset = "no_agony"` excludes only samples labelled `postmortem_agony` and warns when the agony axis is unavailable; the `agony` interval applies only to non-perfused samples; intervals are validated.
- `annotate_de()` maps HGNC symbols with org.Hs.eg.db and reports the matches; `annotate_mra()` defines the sign of `direction_<axis>` per axis and takes `padj_threshold`.
- `thanatome_score()` accepts a limma `EList`, ignores the normalisation factors of a `DGEList` (as the references do), rejects duplicated sample names, and `drop_unexpressed` now leaves out only genes absent from the cohort (zero counts in most samples).
- `add_thanatome_score()` returns the `SummarizedExperiment` with the score in its `colData` and the score metadata in `metadata(x)$thanatome`.
- `death_gene_sets()` exports the signatures as signed gene sets (named list or GSEABase `GeneSetCollection`) for gene-set tools.
- `select_samples()`, `plot_thanatome()` and `death_covariates()` for sample selection, display and adjustment; `contamination_flags()` and `hardy_layers()` as helpers.
- `annotate_de()` and `as_de_table()` for annotating an existing differential-expression table.
- Data: `thanatome_signatures` (consensus signatures), `thanatome_references` (per-tissue references and poles, from recount3) and `thanatome_example` (74 GTEx ovary controls and 40 TCGA-OV tumours from recount3); accessors `death_signatures()` and `death_reference()`; build scripts in `inst/scripts/`.
- `compare_controls()` repeats the DE with explicit control subsets against random subsets of the same size; `annotate_mra()` annotates regulator tables; `audit_report()` writes the summary. The DE inside `compare_controls()` follows the edgeR 4 quasi-likelihood pipeline (`glmQLFit()` without `estimateDisp()`).
