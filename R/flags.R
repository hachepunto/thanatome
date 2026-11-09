## Assisted by Claude Code (Anthropic); see README.md, section Provenance.

#' Technical contamination flags per sample
#'
#' Flags samples whose mean log2-CPM over the pancreatic enzyme genes
#' (PNLIP, CELA3A, PRSS1, CPA1, CELA2A, CTRB1; contamination by neighbouring
#' tissue, Nieuwenhuis et al. 2020) or over the surfactant genes (SFTPC,
#' SFTPB, SFTPA1, SFTPA2) exceeds an absolute threshold. The threshold is a
#' heuristic on absolute expression, not a z-score against the reference,
#' because the reference tissues do not express these genes and a z-score
#' against a near-zero standard deviation is meaningless. A flag is not
#' raised for the tissue that expresses its genes physiologically: the
#' surfactant flag for lung and the pancreatic flag for pancreas, when
#' `tissue` says so.
#'
#' @param x A `SummarizedExperiment` or matrix of raw counts (or log2-CPM with
#'   `log = TRUE`), as in [thanatome_score()].
#' @param threshold Mean log2-CPM above which the flag is raised. Default 4
#'   (about 16 CPM).
#' @param tissue Optional tissue of the samples (case-insensitive; names
#'   starting with `"lung"` disable the surfactant flag and names starting
#'   with `"pancrea"` the pancreatic flag). [thanatome_score()] passes the
#'   tissue of the packaged reference unless told otherwise.
#' @param assay,lib_size,log Passed to the normalisation; see
#'   [thanatome_score()].
#' @return A `data.frame` with logical columns `pancreas` and `surfactant`,
#'   the mean log2-CPM of each set (`pancreas_logcpm`, `surfactant_logcpm`)
#'   and the number of genes found (`pancreas_n`, `surfactant_n`). When a set
#'   has no gene in `x` the flag is `NA` and a message says so.
#' @examples
#' data(thanatome_example)
#' fl <- contamination_flags(thanatome_example)
#' table(fl$pancreas, useNA = "ifany")
#' @export
contamination_flags <- function(x, threshold = 4, tissue = NULL,
                                assay = "counts", lib_size = NULL,
                                log = FALSE) {
    m <- .as_logcpm(x, assay = assay, log = log, lib_size = lib_size)
    .flags_from_logcpm(m, threshold = threshold, tissue = tissue)
}
## TRUE if 'tissue' names a tissue whose name starts with 'prefix'
.is_tissue <- function(tissue, prefix) {
    length(tissue) == 1 && !is.na(tissue) &&
        startsWith(tolower(tissue), prefix)
}
.flags_from_logcpm <- function(m, threshold = 4, tissue = NULL) {
    flags <- .signatures()$flags
    one <- function(f) {
        g <- intersect(flags$ensembl[flags$flag == f], rownames(m))
        if (!length(g)) {
            message(
                "contamination_flags: no ", f,
                " gene found in the data; flag set to NA"
            )
            return(list(v = rep(NA_real_, ncol(m)), n = 0L))
        }
        list(v = colMeans(m[g, , drop = FALSE]), n = length(g))
    }
    pz <- one("pancreas")
    sz <- one("surfactant")
    off <- rep(FALSE, ncol(m))
    data.frame(
        pancreas = if (.is_tissue(tissue, "pancrea")) off else pz$v > threshold,
        surfactant = if (.is_tissue(tissue, "lung")) off else sz$v > threshold,
        pancreas_logcpm = pz$v, surfactant_logcpm = sz$v,
        pancreas_n = pz$n, surfactant_n = sz$n,
        row.names = colnames(m), stringsAsFactors = FALSE
    )
}
