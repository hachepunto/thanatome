## Assisted by Claude Code (Anthropic); see README.md, section Provenance.

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
#' @param flags Logical; `TRUE` (default) removes the controls with a
#'   technical flag from every subset, including `all` and the random ones,
#'   so that all comparisons share the same pool.
#' @param assay,lib_size Passed to the normalisation; see
#'   [thanatome_score()]. `x` must hold raw counts (edgeR needs them). With
#'   fewer than 5000 genes in `x`, the library sizes are used as given and
#'   no TMM factors are computed, because TMM assumes that most genes do not
#'   change.
#' @return A list with `summary` (one row per criterion: controls kept, DEGs,
#'   lost and gained with respect to `all`, the same for the random subsets,
#'   Jaccard, standard deviation of the change in log fold change, fraction
#'   of DEGs that are signature genes, sign flips, and the MRA agreement when
#'   `mra` is given), `de` (the DE table of each criterion), `kept` (the
#'   control names kept by each criterion), `scores` and `mra`.
#' @examples
#' data(thanatome_example)
#' cc <- compare_controls(thanatome_example, group = "group",
#'                        control = "control",
#'                        criteria = list(no_agony = list(preset = "no_agony")),
#'                        reference = "ovary", n_random = 2)
#' cc$summary
#' @export
compare_controls <- function(x, group, control = "control", criteria = list(),
                            scores = NULL, reference = NULL, n_random = 10,
                            padj = 0.05, lfc = 1, min_controls = 8,
                            mra = NULL, flags = TRUE, assay = "counts",
                            lib_size = NULL) {
    .check_criteria(criteria)
    is_ctrl <- .group_vector(x, group) == control
    if (sum(is_ctrl) < min_controls || sum(!is_ctrl) < 2) {
        stop("too few controls or cases")
    }
    scores <- .control_scores(x, is_ctrl, scores, reference, assay, lib_size)
    ## Common set: controls without technical flags (if flags = TRUE);
    ## 'all' and the random subsets are drawn from it
    ctrl_names <- select_samples(scores, flags = flags)$keep
    m <- .raw_counts(x, assay)
    libs <- .lib_sizes(x, m, lib_size)
    cases <- colnames(x)[!is_ctrl]
    run_de <- function(keep) .de_edger(m, libs, cases = cases, controls = keep)
    kept <- .kept_controls(criteria, scores, flags, ctrl_names, min_controls)
    de <- lapply(kept, run_de)
    deg <- function(d) {
        d$gene[!is.na(d$padj) & d$padj <= padj & abs(d$log2FoldChange) >= lfc]
    }
    sig_genes <- unique(.signatures()$signatures$ensembl)
    res <- lapply(names(kept), .criterion_metrics, de = de, kept = kept,
                ctrl_names = ctrl_names, run_de = run_de, deg = deg,
                deg_base = deg(de$all), sig_genes = sig_genes,
                n_random = n_random)
    names(res) <- names(kept)
    summary <- S4Vectors::DataFrame(do.call(rbind, lapply(res, `[[`, "row")))
    rownames(summary) <- NULL
    mra_out <- NULL
    if (!is.null(mra)) {
        ag <- .mra_agreement(mra, de, res)
        summary <- cbind(summary, ag$table)
        mra_out <- ag$activities
    }
    list(summary = summary, de = de, kept = kept, scores = scores,
        mra = mra_out, thresholds = c(padj = padj, lfc = lfc))
}

## ---- internal helpers -------------------------------------------------------

.check_criteria <- function(criteria) {
    nm <- names(criteria)
    bad_names <- length(criteria) &&
        (is.null(nm) || anyDuplicated(nm) || any(nm == ""))
    if (!is.list(criteria) || bad_names) {
        stop("'criteria' must be a named list with unique names")
    }
    invisible(TRUE)
}
## Group vector from a colData column name or given directly
.group_vector <- function(x, group) {
    in_cd <- length(group) == 1 &&
        group %in% colnames(SummarizedExperiment::colData(x))
    g <- if (in_cd) x[[group]] else group
    if (anyNA(g)) stop("'group' has missing values")
    if (is.null(colnames(x))) stop("'x' needs sample names")
    g
}
## Scores of the controls: given, or computed against 'reference'
.control_scores <- function(x, is_ctrl, scores, reference, assay, lib_size) {
    if (is.null(scores)) {
        if (is.null(reference)) {
            stop("give 'scores' or a 'reference' to score the controls")
        }
        scores <- thanatome_score(x[, is_ctrl], reference = reference,
                                assay = assay, lib_size = lib_size[is_ctrl])
    }
    if (!all(rownames(scores) %in% colnames(x)[is_ctrl])) {
        stop("'scores' must be the scores of the control samples of 'x'")
    }
    scores
}
## Raw count matrix with unversioned ids; refuses normalised values
.raw_counts <- function(x, assay) {
    ## log = TRUE only skips the CPM step: the values stay raw counts
    m <- .as_logcpm(x, assay = assay, log = TRUE)
    if (any(m < 0, na.rm = TRUE) || any(m != round(m), na.rm = TRUE)) {
        stop("compare_controls needs raw counts in 'x'")
    }
    m
}
.lib_sizes <- function(x, m, lib_size) {
    cd <- SummarizedExperiment::colData(x)
    libs <- if (!is.null(lib_size)) {
        lib_size
    } else if ("lib_size" %in% colnames(cd)) {
        x$lib_size
    } else {
        colSums(m)
    }
    names(libs) <- colnames(m)
    libs
}
## Controls kept by 'all' and by each criterion; drops criteria left with
## fewer than 'min_controls'
.kept_controls <- function(criteria, scores, flags, ctrl_names, min_controls) {
    one <- function(k) {
        args <- utils::modifyList(list(scores = scores, flags = flags),
                                as.list(k))
        intersect(do.call(select_samples, args)$keep, ctrl_names)
    }
    kept <- c(list(all = ctrl_names), lapply(criteria, one))
    too_few <- lengths(kept) < min_controls
    if (any(too_few)) {
        warning("criteria dropped for keeping fewer than ", min_controls,
                " controls: ", paste(names(kept)[too_few], collapse = ", "))
        kept <- kept[!too_few]
    }
    kept
}
## Metrics of one criterion against all controls and against random subsets
## of the same n
.criterion_metrics <- function(k, de, kept, ctrl_names, run_de, deg, deg_base,
                                sig_genes, n_random) {
    base <- de$all
    d <- de[[k]]
    dk <- deg(d)
    n <- length(kept[[k]])
    rnd <- if (k == "all" || n >= length(ctrl_names)) {
        list()
    } else {
        lapply(seq_len(n_random), function(i) run_de(sample(ctrl_names, n)))
    }
    rnd_deg <- lapply(rnd, deg)
    rnd_n <- function(f) .med(vapply(rnd_deg, f, numeric(1)))
    row <- data.frame(
        criterion = k, n_controls = n, n_deg = length(dk),
        lost = length(setdiff(deg_base, dk)),
        gained = length(setdiff(dk, deg_base)),
        lost_random = rnd_n(function(r) length(setdiff(deg_base, r))),
        gained_random = rnd_n(function(r) length(setdiff(r, deg_base))),
        jaccard = .jaccard(dk, deg_base),
        jaccard_random = rnd_n(function(r) .jaccard(r, deg_base)),
        sd_delta_lfc = .sd_delta_lfc(d, base),
        sd_delta_lfc_random = .med(vapply(rnd, .sd_delta_lfc, numeric(1),
                                        base = base)),
        frac_deg_in_axes = if (length(dk)) mean(dk %in% sig_genes)
                            else NA_real_,
        sign_flips = .sign_flips(d, base, intersect(dk, deg_base)),
        stringsAsFactors = FALSE)
    list(row = row, rnd = rnd)
}
## SD of the change in lfc of 'd' with respect to 'base', over base genes
.sd_delta_lfc <- function(d, base) {
    delta <- d$log2FoldChange[match(base$gene, d$gene)] - base$log2FoldChange
    sd(delta, na.rm = TRUE)
}
.sign_flips <- function(d, base, genes) {
    s_d <- sign(d$log2FoldChange[match(genes, d$gene)])
    s_b <- sign(base$log2FoldChange[match(genes, base$gene)])
    sum(s_d != s_b)
}
## Agreement of regulator activities with those of 'all': Spearman and
## overlap of the top 25; for the random subsets, the median over subsets
.mra_agreement <- function(mra, de, res) {
    act <- lapply(de, mra)
    agree <- function(a) .mra_agree(a, act$all)
    ag <- t(vapply(act, agree, numeric(2)))
    rnd_one <- function(k) {
        r <- res[[k]]$rnd
        if (!length(r)) return(c(NA_real_, NA_real_))
        apply(vapply(r, function(d) agree(mra(d)), numeric(2)), 1,
                stats::median)
    }
    rnd_ag <- t(vapply(names(res), rnd_one, numeric(2)))
    colnames(rnd_ag) <- paste0(colnames(ag), "_random")
    list(table = S4Vectors::DataFrame(ag, rnd_ag), activities = act)
}
.mra_agree <- function(a, b) {
    common <- intersect(names(a), names(b))
    k <- min(25, length(common))
    top <- function(v) names(sort(v[common], decreasing = TRUE))[seq_len(k)]
    c(mra_spearman = stats::cor(a[common], b[common], method = "spearman"),
        mra_top25_overlap = length(intersect(top(a), top(b))) / k)
}
## edgeR quasi-likelihood, cases versus controls, no covariates
.de_edger <- function(counts, libs, cases, controls) {
    cols <- c(cases, controls)
    grp <- rep(c("case", "control"), c(length(cases), length(controls)))
    y <- edgeR::DGEList(counts[, cols, drop = FALSE], lib.size = libs[cols])
    keep <- edgeR::filterByExpr(y, group = grp)
    y <- y[keep, , keep.lib.sizes = TRUE]
    ## TMM only with a reasonably complete transcriptome
    if (nrow(counts) >= 5000) y <- edgeR::normLibSizes(y)
    design <- stats::model.matrix(~ as.numeric(grp == "case"))
    ## edgeR >= 4: glmQLFit estimates the dispersions itself (no estimateDisp)
    fit <- edgeR::glmQLFit(y, design)
    tt <- edgeR::glmQLFTest(fit, coef = 2)$table
    S4Vectors::DataFrame(gene = rownames(tt), log2FoldChange = tt$logFC,
                        stat = sign(tt$logFC) * sqrt(tt$F),
                        pvalue = tt$PValue,
                        padj = stats::p.adjust(tt$PValue, "BH"))
}
.jaccard <- function(x, y) {
    u <- length(union(x, y))
    if (u == 0) NA_real_ else length(intersect(x, y)) / u
}
.med <- function(v) if (length(v)) stats::median(v) else NA_real_
