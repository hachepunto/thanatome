## Asistido por Claude Code (Anthropic); ver README.md, sección Provenance.

#' Consensus death signatures from GTEx v8
#'
#' Genes with a concordant sign of differential expression in at least 5 of
#' 10 GTEx v8 tissues for each of three contrasts within the GTEx donors:
#' `ventilator` (Hardy scale 0, death on a ventilator, versus Hardy 1-4;
#' |median log2 fold change| > 0.5), `agony` (slow deaths, Hardy 3-4, versus
#' fast deaths, Hardy 1-2, among non-ventilated donors; > 0.5) and
#' `ischemia` (post-mortem ischemic time as a continuous covariate, adjusted
#' for ventilator; > 0.25 because the contrast is per standard deviation of
#' ischemic time). Ischemia is not a third axis: 61 of its 86 genes belong
#' to the ventilator signature with the opposite sign (`in_ventilator`) and
#' only 10-23 were expressed per tissue in the analysis (15-29 in the
#' packaged recount3 references); it is the component of the
#' ventilator axis sensitive to ischemic time. Built by
#' `inst/scripts/build_signatures.R` (see `system.file("scripts",
#' package = "thanatome")`).
#'
#' @format A list with four elements:
#' \describe{
#'   \item{signatures}{data.frame with one row per gene and signature:
#'     `ensembl` (without version), `symbol`, `axis`, `lfc_consensus`
#'     (median log2 fold change across tissues), `z_stouffer`, `n_tissues`
#'     (concordant tissues), `direction` (`up` or `down`) and
#'     `in_ventilator`.}
#'   \item{flags}{data.frame of technical flag genes: pancreatic enzymes
#'     (contamination by neighbouring tissue) and surfactant proteins.}
#'   \item{thresholds}{`min_tissues` and the per-axis `min_lfc`.}
#'   \item{source}{one line describing the origin.}
#' }
#' @source GTEx v8 gene read counts from the GTEx portal
#'   (\url{https://gtexportal.org}, dbGaP phs000424.v8), ten tissues: ovary,
#'   breast, thyroid, lung, sigmoid colon, transverse colon, liver, prostate,
#'   uterus and skin (lower leg); DESeq2 per tissue and Stouffer consensus,
#'   analysis repository \url{https://github.com/hachepunto/thanatome-manuscript}.
#' @references GTEx Consortium (2020) Science 369:1318-1330
#'   \doi{10.1126/science.aaz1776}; Ferreira PG et al. (2018) Nature
#'   Communications 9:490 \doi{10.1038/s41467-017-02772-x}; Nieuwenhuis TO
#'   et al. (2020) Nature Communications 11:1933
#'   \doi{10.1038/s41467-020-15821-9}.
#' @usage data(thanatome_signatures)
#' @examples
#' data(thanatome_signatures)
#' table(thanatome_signatures$signatures$axis, thanatome_signatures$signatures$direction)
"thanatome_signatures"

#' Per-tissue GTEx references for centring new samples
#'
#' For each of the ten GTEx v8 tissues used to derive the signatures: the
#' mean and standard deviation of log2-CPM (prior count 1, full library
#' size) of every signature and flag gene across the GTEx samples of that
#' tissue, which of those genes are expressed there (mean log2-CPM above 0;
#' only these are used for scoring), and the poles of each axis (median
#' signature score of the ventilated and non-ventilated donors for
#' `ventilator`; of fast and slow deaths for `agony`) that define the scaled
#' position returned by [thanatome_score()]. Counts come from recount3
#' (Monorail pipeline, GENCODE v26); the references apply to data
#' quantified the same way. Built by `inst/scripts/build_references.R`.
#'
#' @format A named list, one element per tissue (`ovary`, `breast`,
#'   `thyroid`, `lung`, `colon_sigmoid`, `colon_transverse`, `liver`,
#'   `prostate`, `uterus`, `skin_leg`), each a list with `tissue`,
#'   `gtex_name`, `pipeline`, `n_samples`, `hardy` (table of samples by
#'   Hardy scale), `mean`, `sd`, `expressed` (named vectors over Ensembl
#'   ids), `poles` (data.frame) and `ischemic_time_range` (minutes).
#' @source recount3 (\url{https://rna.recount.bio}; Wilks C et al. (2021)
#'   Genome Biology 22:323 \doi{10.1186/s13059-021-02533-6}), GTEx projects
#'   OVARY, BREAST, THYROID, LUNG, COLON, LIVER, PROSTATE, UTERUS and SKIN,
#'   RNASEQ freeze samples, one run per sample.
#' @usage data(thanatome_references)
#' @examples
#' data(thanatome_references)
#' names(thanatome_references)
#' thanatome_references$ovary$poles
"thanatome_references"

#' Example data: GTEx ovary controls and TCGA-OV tumours
#'
#' A `SummarizedExperiment` with raw read counts of the signature and flag
#' genes for 74 GTEx v8 ovary samples (stratified by Hardy scale) and 40
#' TCGA-OV primary tumours, both from recount3 (same Monorail pipeline),
#' with the full library size (`lib_size`, needed because only a subset of
#' genes is kept), `group` (`control` or `tumor`), donor, Hardy scale,
#' ischemic time and RIN in `colData`. The GTEx samples entered the
#' consensus, the references and the poles, so any accuracy measured on
#' them is in-sample. Built by `inst/scripts/build_example.R`.
#'
#' @format A `SummarizedExperiment` with 1663 genes and 114 samples.
#' @source recount3 projects OVARY (GTEx v8, RNASEQ freeze) and OV (TCGA
#'   primary tumours), gene read counts, GENCODE v26.
#' @usage data(thanatome_example)
#' @examples
#' data(thanatome_example)
#' table(thanatome_example$group, thanatome_example$hardy, useNA = "ifany")
"thanatome_example"
