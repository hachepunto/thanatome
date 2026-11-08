## Asistido por Claude Code (Anthropic); ver README.md, sección Provenance.

#' Consensus death signatures from GTEx v8
#'
#' Genes with a concordant sign of differential expression in at least 5 of
#' 10 GTEx v8 tissues and an absolute median log2 fold change above 0.5
#' (ventilator, agony) or 0.25 (ischemia, whose contrast is per standard
#' deviation of ischemic time) for each of three contrasts within the GTEx donors: `ventilator` (Hardy scale
#' 0, death on a ventilator, versus Hardy 1-4), `ischemia` (post-mortem
#' ischemic time as a continuous covariate, adjusted for ventilator) and
#' `agony` (slow deaths, Hardy 3-4, versus fast deaths, Hardy 1-2, among
#' non-ventilated donors). Built by `data-raw/build_firmas.R` from the
#' analysis repository `hachepunto/gtex-controls`.
#'
#' @format A list with three elements:
#' \describe{
#'   \item{signatures}{data.frame with one row per gene and axis: `ensembl`
#'     (without version), `symbol`, `axis`, `lfc_consensus` (median log2 fold
#'     change across tissues), `z_stouffer`, `n_tissues` (concordant tissues)
#'     and `direction` (`up` or `down`).}
#'   \item{flags}{data.frame of technical flag genes: pancreatic enzymes
#'     (contamination of neighbouring tissue) and surfactant proteins.}
#'   \item{thresholds}{the selection thresholds `min_tissues` and `min_lfc` (per axis).}
#' }
#' @source GTEx v8 (dbGaP phs000424.v8), ten tissues: ovary, breast, thyroid,
#'   lung, sigmoid colon, transverse colon, liver, prostate, uterus and skin
#'   (lower leg).
#' @examples
#' data(thanatome_firmas)
#' table(thanatome_firmas$signatures$axis, thanatome_firmas$signatures$direction)
"thanatome_firmas"

#' Per-tissue GTEx references for centring new samples
#'
#' For each of the ten GTEx v8 tissues used to derive the signatures: the mean
#' and standard deviation of log2-CPM of every signature and flag gene across
#' the GTEx samples of that tissue, and the poles of each axis (median
#' signature score of the ventilated and non-ventilated donors for
#' `ventilator`; of fast and slow deaths for `agony`) that define the scaled
#' position returned by [thanatome_score()]. Built by
#' `data-raw/build_referencias.R`.
#'
#' @format A named list, one element per tissue (`ovary`, `breast`,
#'   `thyroid`, `lung`, `colon_sigmoid`, `colon_transverse`, `liver`,
#'   `prostate`, `uterus`, `skin_leg`), each a list with `tissue`,
#'   `gtex_name`, `n`, `hardy` (table of donors by Hardy scale), `mean`,
#'   `sd` (named numeric vectors over Ensembl ids), `poles` (data.frame) and
#'   `ischemic_range` (minutes).
#' @source GTEx v8 gene read counts.
#' @examples
#' data(thanatome_referencias)
#' names(thanatome_referencias)
#' thanatome_referencias$ovary$poles
"thanatome_referencias"

#' Example data: 67 GTEx v8 ovary samples
#'
#' A `SummarizedExperiment` with raw read counts of the signature and flag
#' genes for 67 GTEx v8 ovary samples stratified by Hardy scale, with the
#' full library size (`lib_size`, needed because only a subset of genes is
#' kept), donor, Hardy scale, ischemic time, RIN, sex and age in `colData`. GTEx v8 gene
#' read counts are open access. Built by `data-raw/build_example.R`.
#'
#' @format A `SummarizedExperiment` with 1522 genes and 67 samples.
#' @source GTEx v8 gene read counts.
#' @examples
#' data(thanatome_example)
#' table(thanatome_example$hardy)
"thanatome_example"
