## Assisted by Claude Code (Anthropic); see README.md, section Provenance.

#' Per-sample death scores as covariates for designs and batch correction
#'
#' Returns a numeric matrix (one row per sample, one column per axis) ready
#' to be passed as `covariates` to `limma::removeBatchEffect()`, added to a
#' `model.matrix()` for limma, edgeR or DESeq2, or given to `sva` as known
#' surrogate variables. `ComBat` takes only categorical batches; the `state`
#' label can be used there, but a continuous adjustment is the honest use of
#' these scores.
#'
#' Two cautions are printed on each call unless `quiet = TRUE`. First, the
#' axes vary only among the post-mortem samples: in a tumour-versus-control
#' design every tumour sits at the same (undefined) value, so adjusting across
#' groups partially confounds the death axes with the group effect; adjust
#' within controls, or include the covariates in a design where the group term
#' is kept.
#' Second, the adjustment centres the controls at their own mean, not at the
#' living pole: it removes the variation between donors, not the distance
#' between post-mortem and living tissue.
#'
#' @param scores The `DataFrame` returned by [thanatome_score()].
#' @param axes Columns to export (default the two axes). The ischemia
#'   component can be added explicitly, but it is nearly collinear with the
#'   ventilator score (about -0.8 in GTEx) and only informative among
#'   perfused samples.
#' @param center Logical; centre each column (default `TRUE`).
#' @param quiet Logical; `TRUE` suppresses the cautions.
#' @return A numeric matrix with the samples in rows.
#' @examples
#' data(thanatome_example)
#' ctl <- thanatome_example[, thanatome_example$group == "control"]
#' sc <- thanatome_score(ctl, reference = "ovary")
#' cv <- death_covariates(sc)
#' head(cv)
#' ## e.g. limma::removeBatchEffect(logcpm, covariates = cv)
#' @export
death_covariates <- function(scores, axes = c("ventilator", "agony"),
                            center = TRUE, quiet = FALSE) {
    missing_axes <- setdiff(axes, colnames(scores))
    if (length(missing_axes)) {
        warning("not in 'scores': ", paste(missing_axes, collapse = ", "))
    }
    axes <- intersect(axes, colnames(scores))
    if (!length(axes)) stop("no requested axis found in 'scores'")
    m <- as.matrix(as.data.frame(scores[, axes, drop = FALSE]))
    rownames(m) <- rownames(scores)
    if (center) m <- scale(m, center = TRUE, scale = FALSE)
    if (!quiet) {
        message(
            "death_covariates: adjust within post-mortem controls or keep ",
            "the group term in the design; ",
            "the adjustment centres controls at their mean, not at the ",
            "living pole"
        )
    }
    m[, , drop = FALSE]
}
