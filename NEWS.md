# thanatome 0.99.0

- First version submitted to Bioconductor.
- `thanatome_score()`: per-sample signature score, scaled position, state label and technical flags along the ventilator, ischemia and agony axes, with a GTEx tissue or a cohort as reference; `method = "regression"` adds the joint coefficients and the fraction of deviation explained.
- `select_samples()`, `plot_thanatome()` and `death_covariates()` for sample selection, display and adjustment.
- `annotate_de()` and `as_de_table()` for annotating an existing differential-expression table.
- Data: `thanatome_firmas` (consensus signatures), `thanatome_referencias` (per-tissue references and poles) and `thanatome_example` (67 GTEx ovary samples).
- Planned for a later version: `annotate_mra()` for regulon tables, `compare_controls()` to repeat a DE with explicit control subsets, and `audit_report()`.
