## Assisted by Claude Code (Anthropic); see README.md, section Provenance.

#' Select samples by explicit criteria on the death axes
#'
#' Returns the samples whose scaled positions, ischemia component, state
#' label and technical flags satisfy the criteria given, together with a table
#' of what was excluded and why. No criterion on the death axes is applied by
#' default (technical flags are excluded unless `flags = FALSE`): the analyses
#' behind this package found that no subset of post-mortem controls is closer to
#' surgical tissue than the full set except, slightly, the ventilated donors,
#' and that excluding samples costs power. Named presets are offered for the two
#' criteria the evidence does support as explicit choices; each is documented
#' with its cost.
#'
#' @param scores The `DataFrame` returned by [thanatome_score()].
#' @param ventilator,agony Numeric length-2 intervals (closed) on
#'   `position_ventilator` and `position_agony`; `NULL` means no criterion.
#'   The agony axis is defined only among non-perfused donors, so the
#'   `agony` interval is applied to non-perfused samples only; samples
#'   labelled `perfused` pass it.
#' @param ischemia Numeric length-2 interval on the raw ischemia component
#'   (meaningful among perfused samples only).
#' @param state Character vector of allowed state labels.
#' @param flags Logical; `TRUE` excludes samples with any technical flag.
#' @param preset `"none"` (default), `"perfused"` (keep only samples labelled
#'   `perfused`: closest to living tissue, but in the GTEx references it
#'   discards the non-ventilated donors, 22-62 percent of the samples
#'   depending on the tissue) or `"no_agony"` (exclude samples labelled
#'   `postmortem_agony`: the only axis that displaces controls coherently
#'   with respect to adjacent normal tissue, at a cost of the Hardy 3-4
#'   donors, 9-17 percent of the samples, with a label that, built from the
#'   joint-regression positions in tissue not seen by the signatures, was
#'   right for 61-84 percent of samples). `"no_agony"` keeps every other
#'   sample (including those without a label) and warns, excluding nothing
#'   on agony, when no sample carries an agony label (agony axis missing).
#'   A preset only acts when `state` is `NULL`; the other criteria are
#'   combined with it.
#' @return A list with `keep` (sample names retained), `excluded` (a
#'   `DataFrame` with the sample and the first criterion it failed) and
#'   `summary` (counts per criterion).
#' @examples
#' data(thanatome_example)
#' ctl <- thanatome_example[, thanatome_example$group == "control"]
#' sc <- thanatome_score(ctl, reference = "ovary")
#' sel <- select_samples(sc, agony = c(-Inf, 0.5))
#' sel$summary
#' head(sel$excluded)
#' @export
select_samples <- function(scores, ventilator = NULL, agony = NULL,
                            ischemia = NULL, state = NULL, flags = TRUE,
                            preset = c("none", "perfused", "no_agony")) {
    preset <- match.arg(preset)
    if (preset == "perfused" && is.null(state)) state <- "perfused"
    no_agony <- preset == "no_agony" && is.null(state)
    if (no_agony) .check_agony_labels(scores)
    crit <- list(
        ventilator = .in_interval(scores$position_ventilator, ventilator),
        agony = .agony_criterion(scores, agony),
        ischemia = .in_interval(scores$ischemia, ischemia),
        state = .state_criterion(scores, state, no_agony),
        flags = .flag_criterion(scores, flags)
    )
    if (flags && (anyNA(scores$flag_pancreas) ||
        anyNA(scores$flag_surfactant))) {
        message(
            "select_samples: some technical flags are NA (genes absent) ",
            "and are treated as not raised"
        )
    }
    if (is.null(rownames(scores))) stop("'scores' has no sample names")
    pass <- Reduce(`&`, crit)
    first_fail <- apply(!do.call(cbind, crit), 1, function(r) {
        if (any(r)) names(crit)[which(r)[1]] else NA_character_
    })
    excluded <- S4Vectors::DataFrame(
        sample = rownames(scores)[!pass],
        reason = first_fail[!pass]
    )
    summary <- vapply(crit, function(k) sum(!k), integer(1))
    list(
        keep = rownames(scores)[pass], excluded = excluded,
        summary = c(n = nrow(scores), kept = sum(pass), summary)
    )
}
## Agony interval applies to non-perfused samples only
.agony_criterion <- function(scores, agony) {
    ok <- .in_interval(scores$position_agony, agony)
    ok | scores$state %in% "perfused"
}
.state_criterion <- function(scores, state, no_agony) {
    if (no_agony) {
        return(!(scores$state %in% "postmortem_agony"))
    }
    if (is.null(state)) {
        return(rep(TRUE, nrow(scores)))
    }
    !is.na(scores$state) & scores$state %in% state
}
.flag_criterion <- function(scores, flags) {
    if (!flags) {
        return(rep(TRUE, nrow(scores)))
    }
    !(scores$flag_pancreas %in% TRUE | scores$flag_surfactant %in% TRUE)
}
## Warn when no sample has an agony label (axis unavailable)
.check_agony_labels <- function(scores) {
    labelled <- scores$state %in% c("postmortem_rapid", "postmortem_agony")
    if (!any(labelled) || all(is.na(scores$position_agony))) {
        warning(
            "select_samples: the agony axis is unavailable (no sample has ",
            "an agony label); preset 'no_agony' excludes nothing on agony",
            call. = FALSE
        )
    }
}
.in_interval <- function(v, interval) {
    if (is.null(interval)) {
        return(rep(TRUE, length(v)))
    }
    if (!is.numeric(interval) || length(interval) != 2 ||
        anyNA(interval) || interval[1] > interval[2]) {
        stop(
            "an interval must be a numeric vector c(lower, upper) with ",
            "lower <= upper",
            call. = FALSE
        )
    }
    !is.na(v) & v >= interval[1] & v <= interval[2]
}
