## Assisted by Claude Code (Anthropic); see README.md, section Provenance.

#' Select samples by explicit criteria on the death axes
#'
#' Returns the samples whose scaled positions, ischemia component, state
#' label and technical flags satisfy the criteria given, together with a table
#' of what was excluded and why. There is no default criterion: the analyses
#' behind this package found that no subset of post-mortem controls is closer to
#' surgical tissue than the full set except, slightly, the ventilated donors,
#' and that excluding samples costs power. Named presets are offered for the two
#' criteria the evidence does support as explicit choices; each is documented
#' with its cost.
#'
#' @param scores The `DataFrame` returned by [thanatome_score()].
#' @param ventilator,agony Numeric length-2 intervals (closed) on
#'   `position_ventilator` and `position_agony`; `NULL` means no criterion.
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
#'   right for 61-84 percent of samples). A preset only sets `state`
#'   when `state` is `NULL`; the other criteria are combined with it.
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
    if (preset == "no_agony" && is.null(state)) {
        state <- c("perfused", "postmortem_rapid")
    }
    crit <- list(
        ventilator = .in_interval(scores$position_ventilator, ventilator),
        agony = .in_interval(scores$position_agony, agony),
        ischemia = .in_interval(scores$ischemia, ischemia),
        state = if (is.null(state)) {
            rep(TRUE, nrow(scores))
        } else {
            !is.na(scores$state) & scores$state %in% state
        },
        flags = if (!flags) {
            rep(TRUE, nrow(scores))
        } else {
            !(scores$flag_pancreas %in% TRUE |
                scores$flag_surfactant %in% TRUE)
        }
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
.in_interval <- function(v, interval) {
    if (is.null(interval)) {
        return(rep(TRUE, length(v)))
    }
    stopifnot(length(interval) == 2)
    !is.na(v) & v >= interval[1] & v <= interval[2]
}
