# data-raw: how the package data are built

Nothing in this folder ships in the package tarball (it is listed in `.Rbuildignore`); it is here so that every object in `data/` can be rebuilt from its source and so that a reviewer can see where each number came from.

- `build_firmas.R` reads the consensus tables of the three death axes produced by the analysis repository `hachepunto/gtex-controls` (`results/firmas/{ventilador,isquemia,agonia}_consenso.tsv`, commit recorded in `sessioninfo/build_firmas.txt`) and writes `data/thanatome_firmas.rda`: the signature genes per axis (concordant sign in at least 5 of 10 tissues, |median log2 fold change| > 0.5), the consensus log fold change, the Stouffer z and the number of concordant tissues, plus the technical flag gene sets.
- `build_referencias.R` reads the processed GTEx v8 counts of the ten tissues from the same repository (`data/processed/<tissue>/counts.rds` and `metadata.tsv`) and writes `data/thanatome_referencias.rda`: per tissue, the mean and standard deviation of log2-CPM of every signature gene across the GTEx samples (to centre and scale new samples) and the poles of each axis (median score of the ventilated and non-ventilated donors, and of rapid and slow deaths) that define the scaled position.
- `build_example.R` writes `data/thanatome_example.rda`: a `SummarizedExperiment` with raw counts of the signature genes for 67 GTEx v8 ovary samples stratified by Hardy scale, used in the examples, the tests and the vignette. GTEx v8 gene read counts are open access.
- `record_session.R` writes `sessioninfo/<script>.txt` at the end of each script.

Run them from the package root with the environment variable `GTEX_CONTROLS` pointing at a checkout of `hachepunto/gtex-controls` that already contains the results (`GTEX_CONTROLS=~/gtex-proyect Rscript data-raw/build_firmas.R`).

Assisted by Claude Code (Anthropic); see the Provenance section of README.md.
