# thanatome 0.99.0

- First version submitted to Bioconductor.
- `thanatome_score()`: per-sample signature score, scaled position, state label and technical flags along the ventilator and agony axes, plus the ischemia component of the ventilator axis, with a GTEx tissue (recount3 quantification) or a cohort as reference; `method = "regression"` adds the joint coefficients and the fraction of deviation explained.
- `add_thanatome_score()` returns the `SummarizedExperiment` with the score in its `colData` and the score metadata in `metadata(x)$thanatome`.
- `death_gene_sets()` exports the signatures as signed gene sets (list, GSEABase `GeneSetCollection` or `BiocSet`) for gene-set tools.
- `select_samples()`, `plot_thanatome()` and `death_covariates()` for sample selection, display and adjustment; `contamination_flags()` and `hardy_layers()` as helpers.
- `annotate_de()` and `as_de_table()` for annotating an existing differential-expression table.
- Data: `thanatome_signatures` (consensus signatures), `thanatome_references` (per-tissue references and poles, from recount3) and `thanatome_example` (74 GTEx ovary controls and 40 TCGA-OV tumours from recount3); accessors `death_signatures()` and `death_reference()`; build scripts in `inst/scripts/`.
- `compare_controls()` repeats the DE with explicit control subsets against random subsets of the same size; `annotate_mra()` annotates regulator tables; `audit_report()` writes the summary. The DE inside `compare_controls()` follows the edgeR 4 quasi-likelihood pipeline (`glmQLFit()` without `estimateDisp()`).
