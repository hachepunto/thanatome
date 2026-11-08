## Asistido por Claude Code (Anthropic); ver README.md, sección Provenance.

#' thanatome: donor death signatures for post-mortem control tissue
#'
#' Annotates each sample of a post-mortem tissue cohort with its position
#' along three transcriptional axes of donor death derived from GTEx v8 in
#' ten tissues (ventilator / perfusion, post-mortem ischemia and prolonged
#' agony), lets the user plot and select samples with explicit criteria,
#' exports the per-sample scores as covariates, and annotates existing
#' differential-expression tables with the axes their genes fall in.
#'
#' The main entry point is [thanatome_score()]. See the vignette
#' `vignette("thanatome")` for a worked example.
#'
#' @keywords internal
#' @importFrom methods is
#' @importFrom stats median sd setNames lm.fit quantile
#' @importFrom utils head
#' @importFrom ggplot2 .data
"_PACKAGE"
