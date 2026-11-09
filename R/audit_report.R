## Assisted by Claude Code (Anthropic); see README.md, section Provenance.

#' Plain-text audit report of the control samples
#'
#' Summarises, in text and tables that can be pasted into a supplement, the
#' per-sample annotation ([thanatome_score()]), the selection applied
#' ([select_samples()]) and the sensitivity of the differential expression to
#' the control subset ([compare_controls()]).
#'
#' @param scores The `DataFrame` from [thanatome_score()].
#' @param selection Optional result of [select_samples()].
#' @param comparison Optional result of [compare_controls()].
#' @param file Optional path; the report is written there as well.
#' @return An object of class `thanatome_report`: the character vector with
#'   the report lines, printed line by line when shown.
#' @examples
#' data(thanatome_example)
#' ctl <- thanatome_example[, thanatome_example$group == "control"]
#' sc <- thanatome_score(ctl, reference = "ovary")
#' audit_report(sc, selection = select_samples(sc, preset = "no_agony"))
#' @export
audit_report <- function(scores, selection = NULL, comparison = NULL,
                        file = NULL) {
    lines <- .report_scores(scores)
    if (!is.null(selection)) {
        lines <- c(lines, .report_selection(selection))
    }
    if (!is.null(comparison)) {
        lines <- c(lines, .report_comparison(comparison))
    }
    lines <- c(lines, "", paste0(
        "In GTEx v8, no subset of post-mortem controls was closer to ",
        "surgical or living tissue than the full set except, slightly, ",
        "the ventilated donors; the package annotates and quantifies, ",
        "it does not correct death."
    ))
    if (!is.null(file)) writeLines(lines, file)
    structure(lines, class = "thanatome_report")
}

# Header, states, positions, flags and internal check of the report
.report_scores <- function(scores) {
    meta <- S4Vectors::metadata(scores)
    fmt <- function(x, d = 2) formatC(x, digits = d, format = "f")
    states <- c(
        "perfused", "postmortem", "postmortem_rapid",
        "postmortem_agony"
    )
    st <- table(
        factor(scores$state, levels = intersect(states, unique(scores$state))),
        useNA = "ifany"
    )
    q <- function(v, p) fmt(quantile(v, p, na.rm = TRUE))
    pv <- scores$position_ventilator
    pa <- scores$position_agony
    lines <- c(
        "thanatome audit report", "======================",
        sprintf(
            "Reference: %s; method: %s; %d samples.",
            meta$reference, meta$method, nrow(scores)
        ),
        sprintf(
            "States (threshold %s): %s.", fmt(meta$state_threshold, 2),
            paste(names(st), st, sep = " = ", collapse = ", ")
        ),
        sprintf(
            paste0(
                "Median position on the ventilator axis %s (IQR %s to ",
                "%s); on the agony axis %s (IQR %s to %s)."
            ),
            fmt(median(pv, na.rm = TRUE)), q(pv, 0.25), q(pv, 0.75),
            fmt(median(pa, na.rm = TRUE)), q(pa, 0.25), q(pa, 0.75)
        ),
        sprintf(
            "Technical flags: pancreas %d, surfactant %d.",
            sum(scores$flag_pancreas, na.rm = TRUE),
            sum(scores$flag_surfactant, na.rm = TRUE)
        )
    )
    if (!is.null(meta$auc)) {
        lines <- c(lines, sprintf(
            paste0(
                "Internal check against the Hardy label: AUC ventilator ",
                "%s, agony %s; state agreement %s."
            ),
            fmt(meta$auc[["ventilator"]]), fmt(meta$auc[["agony"]]),
            fmt(meta$state_agreement)
        ))
    }
    lines
}

# Selection section of the report
.report_selection <- function(selection) {
    s <- selection$summary
    c(
        "", "Selection", "---------",
        sprintf(
            "%d of %d samples kept. Excluded by criterion: %s.",
            s[["kept"]], s[["n"]],
            paste(names(s)[-c(1, 2)], s[-c(1, 2)],
                sep = " = ",
                collapse = ", "
            )
        )
    )
}

# Sensitivity section of the report
.report_comparison <- function(comparison) {
    d <- as.data.frame(comparison$summary)
    cols <- c(
        "criterion", "n_controls", "n_deg", "lost", "lost_random",
        "gained", "gained_random", "jaccard", "jaccard_random",
        "sd_delta_lfc", "sd_delta_lfc_random", "frac_deg_in_axes",
        "sign_flips", "mra_spearman", "mra_spearman_random"
    )
    c(
        "", "Sensitivity of the differential expression to the control subset",
        "----------------------------------------------------------------",
        sprintf(
            paste0(
                "DEG defined by padj <= %s and |log2FC| >= %s; random ",
                "subsets of the same size are the yardstick."
            ),
            comparison$thresholds[["padj"]],
            comparison$thresholds[["lfc"]]
        ),
        utils::capture.output(print(d[, intersect(cols, colnames(d))],
            digits = 3, row.names = FALSE
        )),
        paste0(
            "A criterion whose losses and gains do not exceed those of the ",
            "random subsets acts through sample size, not through the death ",
            "axes."
        )
    )
}

#' @rdname audit_report
#' @param x A `thanatome_report`.
#' @param ... Ignored.
#' @export
print.thanatome_report <- function(x, ...) {
    writeLines(unclass(x))
    invisible(x)
}
