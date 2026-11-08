# thanatome 0.99.0

- First version submitted to Bioconductor.
- `thanatome_score()`: per-sample signature score, scaled position, state label and technical flags along the ventilator, ischemia and agony axes, with a GTEx tissue or a cohort as reference; `method = "regression"` adds the joint coefficients and the fraction of deviation explained.
- `select_samples()`, `plot_thanatome()` and `death_covariates()` for sample selection, display and adjustment.
- `annotate_de()` and `as_de_table()` for annotating an existing differential-expression table.
- Data: `thanatome_firmas` (consensus signatures), `thanatome_referencias` (per-tissue references and poles) and `thanatome_example` (74 GTEx ovary controls and 40 TCGA-OV tumours from recount3).
- `compare_controls()` repeats the DE with explicit control subsets against random subsets of the same size; `annotate_mra()` annotates regulator tables; `audit_report()` writes the summary.
