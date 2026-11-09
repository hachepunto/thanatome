## Assisted by Claude Code (Anthropic); see README.md, section Provenance.

#' Add the death-axis annotation to the colData of a SummarizedExperiment
#'
#' Convenience wrapper around [thanatome_score()] that returns `x` itself,
#' with the per-sample annotation appended as `colData` columns, so that the
#' scores travel with the counts into `DESeqDataSet`, `DGEList` or any other
#' downstream object. The metadata of the score (reference, method, poles
#' and the internal check against Hardy) is stored in
#' `metadata(x)$thanatome`.
#'
#' @param x A `SummarizedExperiment` with raw counts (see [thanatome_score()]).
#' @param reference,... Passed to [thanatome_score()].
#' @param prefix Character string prepended to the new column names, to
#'   avoid clashes with existing `colData` columns. Default `""`.
#' @return `x` with the columns of [thanatome_score()] added to its
#'   `colData` (existing columns with the same name are replaced, with a
#'   message) and `metadata(x)$thanatome`.
#' @examples
#' data(thanatome_example)
#' ctl <- thanatome_example[, thanatome_example$group == "control"]
#' ctl <- add_thanatome_score(ctl, reference = "ovary")
#' table(ctl$state)
#' @seealso [thanatome_score()]
#' @export
add_thanatome_score <- function(x, reference, ..., prefix = "") {
    if (!is(x, "SummarizedExperiment")) {
        stop("'x' must be a SummarizedExperiment; ",
            "use thanatome_score() for matrices and DGEList objects")
    }
    sc <- thanatome_score(x, reference = reference, ...)
    colnames(sc) <- paste0(prefix, colnames(sc))
    cd <- SummarizedExperiment::colData(x)
    replaced <- intersect(colnames(sc), colnames(cd))
    if (length(replaced)) {
        message("replacing colData columns: ", paste(replaced, collapse = ", "))
        cd <- cd[, setdiff(colnames(cd), replaced), drop = FALSE]
    }
    SummarizedExperiment::colData(x) <- cbind(cd, sc[colnames(x), ])
    S4Vectors::metadata(x)$thanatome <- S4Vectors::metadata(sc)
    x
}
