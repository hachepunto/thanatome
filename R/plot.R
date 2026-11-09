## Assisted by Claude Code (Anthropic); see README.md, section Provenance.

#' Plot the per-sample death scores
#'
#' `type = "scatter"` draws the samples on two axes (by default the scaled
#' positions on the ventilator and agony axes) with the reference poles at
#' 0 and 1, coloured by state label; `type = "ranked"` draws every sample
#' as a point on each axis, ordered by position, with the state threshold,
#' so that individual samples can be found. Both return a `ggplot` object
#' that can be modified further.
#'
#' @param scores The `DataFrame` returned by [thanatome_score()].
#' @param type `"scatter"` or `"ranked"`.
#' @param x,y Column names used by the scatter plot.
#' @param colour Column used for colour (default `state`); `NULL` for none.
#' @param hardy Optional vector of Hardy scale values (same order as
#'   `scores`) drawn as point shape, to see the agreement with the label.
#' @param label Logical; `TRUE` writes the sample names next to the points
#'   in the scatter plot (useful for small cohorts); ignored by `"ranked"`.
#' @return A `ggplot` object.
#' @examples
#' data(thanatome_example)
#' ctl <- thanatome_example[, thanatome_example$group == "control"]
#' sc <- thanatome_score(ctl, reference = "ovary")
#' plot_thanatome(sc, hardy = ctl$hardy)
#' plot_thanatome(sc, type = "ranked")
#' @export
plot_thanatome <- function(scores, type = c("scatter", "ranked"),
                            x = "position_ventilator", y = "position_agony",
                            colour = "state", hardy = NULL, label = FALSE) {
    type <- match.arg(type)
    d <- as.data.frame(
        scores[, setdiff(colnames(scores), "n_genes"), drop = FALSE]
    )
    d$sample <- rownames(d)
    if (!is.null(hardy)) d$hardy <- factor(hardy)
    thr <- S4Vectors::metadata(scores)$state_threshold
    if (type == "scatter") {
        return(.plot_scatter(d, x, y, colour, hardy, label, thr))
    }
    .plot_ranked(d, colour, hardy, thr)
}

# Scatter of two axes with reference lines at 0, 1 and the state threshold
.plot_scatter <- function(d, x, y, colour, hardy, label, thr) {
    p <- ggplot2::ggplot(d, ggplot2::aes(x = .data[[x]], y = .data[[y]])) +
        ggplot2::geom_vline(
            xintercept = c(0, 1), linetype = "dotted",
            colour = "grey50"
        ) +
        ggplot2::geom_hline(
            yintercept = c(0, 1), linetype = "dotted",
            colour = "grey50"
        )
    if (!is.null(thr) && grepl("^position_", x)) {
        p <- p + ggplot2::geom_vline(xintercept = thr, linetype = "dashed")
    }
    if (!is.null(thr) && grepl("^position_", y)) {
        p <- p + ggplot2::geom_hline(yintercept = thr, linetype = "dashed")
    }
    p <- p + .points(colour, hardy) +
        ggplot2::labs(x = .axis_label(x), y = .axis_label(y))
    if (label) {
        p <- p + ggplot2::geom_text(ggplot2::aes(label = .data$sample),
            size = 2.5, vjust = -0.7
        )
    }
    p + ggplot2::theme_bw()
}

# One facet per axis with the samples ranked by value
.plot_ranked <- function(d, colour, hardy, thr) {
    axes <- intersect(
        c("position_ventilator", "position_agony", "ischemia"),
        colnames(d)
    )
    long <- do.call(rbind, lapply(axes, function(a) {
        data.frame(d[, c("sample", colour, if (!is.null(hardy)) "hardy")],
            axis = .axis_label(a), value = d[[a]], stringsAsFactors = FALSE
        )
    }))
    long <- long[!is.na(long$value), ]
    long$rank <- stats::ave(long$value, long$axis, FUN = rank)
    p <- ggplot2::ggplot(long, ggplot2::aes(x = .data$rank, y = .data$value)) +
        ggplot2::facet_wrap(~axis, scales = "free") +
        .points(colour, hardy) +
        ggplot2::labs(x = "sample (ranked)", y = "score / position") +
        ggplot2::theme_bw()
    if (!is.null(thr)) {
        thr_axes <- .axis_label(c("position_ventilator", "position_agony"))
        p <- p + ggplot2::geom_hline(
            data = data.frame(axis = thr_axes, thr = thr),
            ggplot2::aes(yintercept = .data$thr), linetype = "dashed"
        )
    }
    p
}
.points <- function(colour, hardy) {
    if (is.null(colour) && is.null(hardy)) {
        return(ggplot2::geom_point(alpha = 0.8))
    }
    if (is.null(hardy)) {
        return(ggplot2::geom_point(ggplot2::aes(colour = .data[[colour]]),
            alpha = 0.8
        ))
    }
    if (is.null(colour)) {
        return(ggplot2::geom_point(ggplot2::aes(shape = .data$hardy),
            alpha = 0.8
        ))
    }
    ggplot2::geom_point(
        ggplot2::aes(colour = .data[[colour]], shape = .data$hardy),
        alpha = 0.8
    )
}
.axis_label <- function(k) {
    lab <- c(
        position_ventilator = paste0(
            "ventilator axis (0 = perfused pole, ",
            "1 = post-mortem pole)"
        ),
        position_agony = "agony axis (0 = fast deaths, 1 = slow deaths)",
        ischemia = paste0(
            "ischemia component (informative among perfused ",
            "samples)"
        ),
        ventilator = "ventilator score", agony = "agony score"
    )
    ifelse(k %in% names(lab), lab[k], k)
}
