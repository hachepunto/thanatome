## Asistido por Claude Code (Anthropic); ver README.md, sección Provenance.

#' Repeat a differential-expression analysis with explicit control subsets
#'
#' Scores the controls, applies each selection criterion, repeats the
#' case-versus-control differential expression (edgeR quasi-likelihood) with
#' the controls each criterion keeps, and compares every result with the one
#' obtained with all controls and with `n_random` random subsets of the same
#' size. The random subsets are the yardstick: a criterion that loses or
#' gains no more genes than a random subset of the same size is acting
#' through power, not through the death axes. The function is the packaged
#' form of the sensitivity analysis behind the package and replaces the
#' tempting loop "filter, look at the DE, filter again", which selects on the
#' outcome. Decide the candidate criteria first, run them all once, and
#' choose with their cost in view.
#'
#' @param x A `SummarizedExperiment` with raw counts of cases and controls.
#' @param group Name of the `colData` column with the group, or a vector.
#' @param control Value of `group` that marks the controls.
#' @param criteria Named list of criteria; each element is a list of
#'   arguments for [select_samples()] (for example
#'   `list(preset = "perfused")` or `list(agony = c(-Inf, 0.5))`). The
#'   comparison with all controls is always added as `all`.
#' @param scores Optional `DataFrame` from [thanatome_score()] for the
#'   controls; computed from `reference` when `NULL`.
#' @param reference Reference for [thanatome_score()] when `scores` is
#'   `NULL` (a tissue name or the controls themselves).
#' @param n_random Number of random control subsets of the same size per
#'   criterion. Default 10.
#' @param padj,lfc Thresholds that define a differentially expressed gene.
#' @param min_controls Minimum controls a criterion must keep. Default 8.
#' @param mra Optional function taking the DE table of one run (a
#'   `DataFrame` with `gene`, `log2FoldChange`, `stat`, `padj`) and returning
#'   a named numeric vector of regulator activities (for example a wrapper
#'   around viper or NaRnEA). When given, the activities of each criterion
#'   are compared with those of `all` by Spearman correlation and overlap of
#'   the top 25 regulators.
#' @param assay,lib_size,log Passed to the normalisation; see
#'   [thanatome_score()].
#' @return A list with `summary` (one row per criterion: controls kept, DEGs,
#'   lost and gained with respect to `all`, the same for the random subsets,
#'   Jaccard, standard deviation of the change in log fold change, fraction
#'   of DEGs that are signature genes, sign flips, and the MRA agreement when
#'   `mra` is given), `de` (the DE table of each criterion), `kept` (the
#'   control names kept by each criterion), `scores` and `mra`.
#' @examples
#' data(thanatome_example)
#' cc <- compare_controls(thanatome_example, group = "group", control = "control",
#'                        criteria = list(no_agony = list(preset = "no_agony")),
#'                        reference = "ovary", n_random = 2)
#' cc$summary
#' @export
compare_controls <- function(x, group, control = "control", criteria = list(), scores = NULL, reference = NULL,
                             n_random = 10, padj = 0.05, lfc = 1, min_controls = 8, mra = NULL,
                             assay = "counts", lib_size = NULL, log = FALSE) {
    g <- if (length(group) == 1 && group %in% colnames(SummarizedExperiment::colData(x))) x[[group]] else group
    is_ctrl <- g == control
    if (sum(is_ctrl) < min_controls || sum(!is_ctrl) < 2) stop("too few controls or cases")
    if (is.null(scores)) {
        if (is.null(reference)) stop("give 'scores' or a 'reference' to score the controls")
        scores <- thanatome_score(x[, is_ctrl], reference = reference, assay = assay, lib_size = lib_size, log = log)
    }
    ctrl_names <- colnames(x)[is_ctrl]
    if (!all(rownames(scores) %in% ctrl_names)) stop("'scores' must be the scores of the control samples of 'x'")
    m <- .as_logcpm(x, assay = assay, log = TRUE, lib_size = lib_size)    # conteos crudos (log = TRUE evita normalizar aquí)
    libs <- if (!is.null(lib_size)) lib_size else if ("lib_size" %in% colnames(SummarizedExperiment::colData(x))) x$lib_size else colSums(m)
    names(libs) <- colnames(m)
    run_de <- function(keep) .de_edger(m, libs, cases = colnames(x)[!is_ctrl], controls = keep)
    kept <- c(list(all = ctrl_names), lapply(criteria, function(k) do.call(select_samples, c(list(scores = scores), k))$keep))
    too_few <- vapply(kept, length, integer(1)) < min_controls
    if (any(too_few)) { warning("criteria dropped for keeping fewer than ", min_controls, " controls: ", paste(names(kept)[too_few], collapse = ", ")); kept <- kept[!too_few] }
    de <- lapply(kept, run_de)
    deg <- function(d) d$gene[!is.na(d$padj) & d$padj <= padj & abs(d$log2FoldChange) >= lfc]
    base <- de$all; deg_base <- deg(base)
    sig_genes <- unique(.firmas()$signatures$ensembl)
    res <- lapply(names(kept), function(k) .criterion_metrics(k, de, kept, ctrl_names, run_de, deg, deg_base, sig_genes, n_random))
    names(res) <- names(kept)
    summary <- S4Vectors::DataFrame(do.call(rbind, lapply(res, function(r) r$row)))
    rownames(summary) <- NULL
    mra_out <- NULL
    if (!is.null(mra)) {
        act <- lapply(de, mra); act_base <- act$all
        agree <- function(a) { common <- intersect(names(a), names(act_base)); top <- function(v) names(sort(v, decreasing = TRUE))[seq_len(min(25, length(v)))]
            c(mra_spearman = stats::cor(a[common], act_base[common], method = "spearman"), mra_top25_overlap = length(intersect(top(a), top(act_base))) / 25) }
        ag <- t(vapply(act, agree, numeric(2)))
        rnd_ag <- t(vapply(names(kept), function(k) { r <- res[[k]]$rnd; if (!length(r)) return(c(NA_real_, NA_real_)); apply(vapply(r, function(d) agree(mra(d)), numeric(2)), 1, stats::median) }, numeric(2)))
        colnames(rnd_ag) <- paste0(colnames(ag), "_random")
        summary <- cbind(summary, S4Vectors::DataFrame(ag, rnd_ag))
        mra_out <- act
    }
    list(summary = summary, de = de, kept = kept, scores = scores, mra = mra_out, thresholds = c(padj = padj, lfc = lfc))
}
## Métricas de un criterio frente a todos los controles y frente a subconjuntos al azar del mismo n
.criterion_metrics <- function(k, de, kept, ctrl_names, run_de, deg, deg_base, sig_genes, n_random) {
    base <- de$all; d <- de[[k]]; dk <- deg(d); n <- length(kept[[k]])
    rnd <- if (k == "all" || n >= length(ctrl_names)) list() else lapply(seq_len(n_random), function(i) run_de(sample(ctrl_names, n)))
    rnd_deg <- lapply(rnd, deg)
    delta <- d$log2FoldChange[match(base$gene, d$gene)] - base$log2FoldChange
    rnd_delta <- vapply(rnd, function(r) sd(r$log2FoldChange[match(base$gene, r$gene)] - base$log2FoldChange, na.rm = TRUE), numeric(1))
    both <- intersect(dk, deg_base)
    flips <- sum(sign(d$log2FoldChange[match(both, d$gene)]) != sign(base$log2FoldChange[match(both, base$gene)]))
    row <- data.frame(criterion = k, n_controls = n, n_deg = length(dk),
        lost = length(setdiff(deg_base, dk)), gained = length(setdiff(dk, deg_base)),
        lost_random = .med(vapply(rnd_deg, function(r) length(setdiff(deg_base, r)), numeric(1))),
        gained_random = .med(vapply(rnd_deg, function(r) length(setdiff(r, deg_base)), numeric(1))),
        jaccard = .jaccard(dk, deg_base), jaccard_random = .med(vapply(rnd_deg, .jaccard, numeric(1), y = deg_base)),
        sd_delta_lfc = sd(delta, na.rm = TRUE), sd_delta_lfc_random = .med(rnd_delta),
        frac_deg_in_axes = if (length(dk)) mean(dk %in% sig_genes) else NA_real_, sign_flips = flips, stringsAsFactors = FALSE)
    list(row = row, rnd = rnd)
}
## edgeR cuasi-verosimilitud, casos frente a controles, sin covariables
.de_edger <- function(counts, libs, cases, controls) {
    cols <- c(cases, controls); y <- edgeR::DGEList(counts[, cols, drop = FALSE], lib.size = libs[cols])
    keep <- edgeR::filterByExpr(y, group = rep(c("case", "control"), c(length(cases), length(controls))))
    y <- edgeR::normLibSizes(y[keep, , keep.lib.sizes = TRUE])
    design <- stats::model.matrix(~ rep(c(1, 0), c(length(cases), length(controls))))
    y <- edgeR::estimateDisp(y, design)
    fit <- edgeR::glmQLFit(y, design); tt <- edgeR::glmQLFTest(fit, coef = 2)$table
    S4Vectors::DataFrame(gene = rownames(tt), log2FoldChange = tt$logFC, stat = sign(tt$logFC) * sqrt(tt$F),
                         pvalue = tt$PValue, padj = stats::p.adjust(tt$PValue, "BH"))
}
.jaccard <- function(x, y) { u <- length(union(x, y)); if (u == 0) NA_real_ else length(intersect(x, y)) / u }
.med <- function(v) if (length(v)) stats::median(v) else NA_real_
