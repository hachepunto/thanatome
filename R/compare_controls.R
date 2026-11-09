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
#' @param x A `SummarizedExperiment` or an edgeR `DGEList` with raw counts
#'   of cases and controls. For a `DGEList` the counts are `x$counts`, the
#'   library sizes `x$samples$lib.size` (`norm.factors` are ignored: TMM
#'   is recomputed in every fit), `group` may name a column of `x$samples`
#'   and `design` refers to the columns of `x$samples`.
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
#'   criterion. Default 10. The subsets are drawn with [sample()], so the
#'   results depend on the random number generator: call [set.seed()]
#'   beforehand to reproduce them.
#' @param padj,lfc Thresholds that define a differentially expressed gene.
#' @param min_controls Minimum controls a criterion must keep. Default 8.
#' @param mra Optional function taking the DE table of one run (a
#'   `DataFrame` with `gene`, `log2FoldChange`, `stat`, `padj`) and returning
#'   a named numeric vector of regulator activities (for example a wrapper
#'   around viper or NaRnEA). When given, the activities of each criterion
#'   are compared with those of `all` by Spearman correlation and overlap of
#'   the top 25 regulators. It is always called serially.
#' @param flags Logical; `TRUE` (default) removes the controls with a
#'   technical flag from every subset, including `all` and the random ones,
#'   so that all comparisons share the same pool.
#' @param assay,lib_size Passed to the normalisation; see
#'   [thanatome_score()]. `x` must hold raw counts (edgeR needs them). With
#'   fewer than 5000 genes in `x`, the library sizes are used as given and
#'   no TMM factors are computed, because TMM assumes that most genes do not
#'   change.
#' @param design Optional one-sided formula over the sample annotation
#'   (`colData(x)`, or `x$samples` for a `DGEList`) with the covariates to
#'   include in every DE fit, for example `~ batch + sex`. The case/control
#'   indicator is added as the last term and is the tested coefficient;
#'   `NULL` (default) fits `~ group`. The covariates should be the same as
#'   in your real DE analysis, so that the comparison reproduces it. The
#'   function stops if a variable is missing, has missing values, or the
#'   design matrix of a fit is not of full column rank (for example a
#'   covariate confounded with the group, or constant in the samples of a
#'   fit).
#' @param BPPARAM A [BiocParallel::BiocParallelParam-class] object. The
#'   edgeR fits (one per criterion plus the random-subset fits) are run with
#'   [BiocParallel::bplapply()]. The random subsets are always drawn first,
#'   in the main process, so the results do not depend on the backend.
#'   Default [BiocParallel::SerialParam()].
#' @return A list with `summary` (one row per criterion: controls kept, DEGs,
#'   lost and gained with respect to `all`, the same for the random subsets,
#'   Jaccard, standard deviation of the change in log fold change, fraction
#'   of DEGs that are signature genes, sign flips, and the MRA agreement when
#'   `mra` is given), `de` (the DE table of each criterion), `kept` (the
#'   control names kept by each criterion), `scores` and `mra`.
#' @examples
#' data(thanatome_example)
#' set.seed(1)
#' cc <- compare_controls(thanatome_example, group = "group",
#'                        control = "control",
#'                        criteria = list(no_agony = list(preset = "no_agony")),
#'                        reference = "ovary", n_random = 2)
#' cc$summary
#'
#' ## Covariates of the real DE in every fit, and the fits spread over
#' ## workers (MulticoreParam(4), for example, outside Windows)
#' x <- thanatome_example
#' x$batch <- factor(rep(c("a", "b"), length.out = ncol(x)))
#' set.seed(1)
#' cc2 <- compare_controls(x, group = "group", control = "control",
#'     criteria = list(no_agony = list(preset = "no_agony")),
#'     reference = "ovary", n_random = 2, design = ~ batch,
#'     BPPARAM = BiocParallel::SerialParam())
#' cc2$summary[, c("criterion", "n_controls", "n_deg", "lost")]
#' @export
compare_controls <- function(x, group, control = "control", criteria = list(),
                            scores = NULL, reference = NULL, n_random = 10,
                            padj = 0.05, lfc = 1, min_controls = 8,
                            mra = NULL, flags = TRUE, assay = "counts",
                            lib_size = NULL, design = NULL,
                            BPPARAM = BiocParallel::SerialParam()) {
    .check_criteria(criteria)
    cd <- .sample_data(x)
    .check_design(design, cd)
    is_ctrl <- .group_vector(x, group, cd) == control
    if (sum(is_ctrl) < min_controls || sum(!is_ctrl) < 2) {
        stop("too few controls or cases")
    }
    scores <- .control_scores(x, is_ctrl, scores, reference, assay, lib_size)
    ## Common set: controls without technical flags (if flags = TRUE);
    ## 'all' and the random subsets are drawn from it
    ctrl_names <- select_samples(scores, flags = flags)$keep
    kept <- .kept_controls(criteria, scores, flags, ctrl_names, min_controls)
    ## Draw every random subset first (main process, fixed RNG order)
    rnd_sets <- .draw_subsets(kept, ctrl_names, n_random)
    fits <- .fit_all(x, cd, assay, lib_size, colnames(x)[!is_ctrl], design,
                    c(kept, unlist(rnd_sets, recursive = FALSE)), BPPARAM)
    de <- fits[seq_along(kept)]
    names(de) <- names(kept)
    rnd_de <- split(fits[-seq_along(kept)],
                    factor(rep(names(kept), lengths(rnd_sets)), names(kept)))
    .assemble(de, rnd_de, kept, ctrl_names, mra, padj, lfc, scores)
}

## ---- internal helpers -------------------------------------------------------

## Metrics, summary table and MRA agreement from the fitted DE tables
.assemble <- function(de, rnd_de, kept, ctrl_names, mra, padj, lfc, scores) {
    deg <- function(d) {
        d$gene[!is.na(d$padj) & d$padj <= padj & abs(d$log2FoldChange) >= lfc]
    }
    sig_genes <- unique(.signatures()$signatures$ensembl)
    one <- function(k) {
        rnd <- unname(rnd_de[[k]])
        .criterion_metrics(k, de[[k]], de$all, length(kept[[k]]), rnd, deg,
                            sig_genes)
    }
    res <- lapply(names(kept), one)
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
## Sample annotation: colData of a SummarizedExperiment, samples of a DGEList
.sample_data <- function(x) {
    if (is(x, "DGEList")) x$samples else SummarizedExperiment::colData(x)
}
## Random subsets of the same size as each criterion (none for 'all' or for
## criteria that keep the whole pool); order: criterion, then replicate
.draw_subsets <- function(kept, ctrl_names, n_random) {
    one <- function(k) {
        n <- length(kept[[k]])
        if (k == "all" || n >= length(ctrl_names)) return(list())
        lapply(seq_len(n_random), function(i) sample(ctrl_names, n))
    }
    lapply(setNames(names(kept), names(kept)), one)
}
## DE fit of every control set (list of control names), in parallel
.fit_all <- function(x, cd, assay, lib_size, cases, design, sets, BPPARAM) {
    m <- .raw_counts(x, assay)
    libs <- .lib_sizes(x, m, lib_size, cd)
    jobs <- lapply(sets, function(s) {
        list(controls = s, design = .fit_design(design, cd, cases, s))
    })
    BiocParallel::bplapply(jobs, .fit_job, counts = m, libs = libs,
                        cases = cases, BPPARAM = BPPARAM)
}
.fit_job <- function(job, counts, libs, cases) {
    .de_edger(counts, libs, cases, job$controls, job$design)
}
.check_design <- function(design, cd) {
    if (is.null(design)) return(invisible(TRUE))
    if (!inherits(design, "formula") || length(design) != 2L) {
        stop("'design' must be a one-sided formula, e.g. ~ batch + sex")
    }
    miss <- setdiff(all.vars(design), colnames(cd))
    if (length(miss)) {
        stop("design variables not found in the sample annotation: ",
            paste(miss, collapse = ", "))
    }
    invisible(TRUE)
}
## Model matrix of one fit: covariates first, case indicator as last column
.fit_design <- function(design, cd, cases, controls) {
    cols <- c(cases, controls)
    vars <- if (is.null(design)) character() else all.vars(design)
    dat <- as.data.frame(cd[match(cols, rownames(cd)), vars, drop = FALSE],
                        optional = TRUE)
    if (anyNA(dat)) stop("missing values in design variables")
    dat$.thanatome_case <- as.numeric(seq_along(cols) <= length(cases))
    tt <- if (is.null(design)) NULL else stats::terms(design)
    f <- stats::reformulate(c(attr(tt, "term.labels"), ".thanatome_case"),
                            intercept = is.null(tt) || attr(tt, "intercept"))
    mm <- stats::model.matrix(f, dat)
    .check_rank(mm, length(cases), length(controls))
}
## Full column rank and residual degrees of freedom, or stop
.check_rank <- function(mm, n_cases, n_controls) {
    if (qr(mm)$rank < ncol(mm) || nrow(mm) <= ncol(mm)) {
        stop("design matrix is not of full column rank with residual ",
            "degrees of freedom (qr rank < number of columns) in a fit ",
            "with ", n_cases, " cases and ", n_controls, " controls")
    }
    mm
}


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
.group_vector <- function(x, group, cd) {
    in_cd <- length(group) == 1 && group %in% colnames(cd)
    g <- if (in_cd) cd[[group]] else group
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
.lib_sizes <- function(x, m, lib_size, cd) {
    libs <- if (!is.null(lib_size)) {
        lib_size
    } else if (is(x, "DGEList")) {
        x$samples$lib.size
    } else if ("lib_size" %in% colnames(cd)) {
        cd$lib_size
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
## of the same n ('rnd': DE tables of the random subsets)
.criterion_metrics <- function(k, d, base, n, rnd, deg, sig_genes) {
    deg_base <- deg(base)
    dk <- deg(d)
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
## edgeR quasi-likelihood, cases versus controls; the tested coefficient is
## the last column of 'design' (rows: cases, then controls)
.de_edger <- function(counts, libs, cases, controls, design) {
    cols <- c(cases, controls)
    grp <- rep(c("case", "control"), c(length(cases), length(controls)))
    y <- edgeR::DGEList(counts[, cols, drop = FALSE], lib.size = libs[cols])
    keep <- edgeR::filterByExpr(y, group = grp)
    y <- y[keep, , keep.lib.sizes = TRUE]
    ## TMM only with a reasonably complete transcriptome
    if (nrow(counts) >= 5000) y <- edgeR::normLibSizes(y)
    ## edgeR >= 4: glmQLFit estimates the dispersions itself (no estimateDisp)
    fit <- edgeR::glmQLFit(y, design)
    tt <- edgeR::glmQLFTest(fit, coef = ncol(design))$table
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
