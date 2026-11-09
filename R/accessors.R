## Assisted by Claude Code (Anthropic); see README.md, section Provenance.

#' Access the consensus death signatures
#'
#' @param axes Which signatures to return: `"ventilator"`, `"agony"` and/or
#'   `"ischemia"` (the ischemia component of the ventilator axis; see
#'   [thanatome_signatures]).
#' @return A `data.frame` with one row per gene and signature (columns
#'   `ensembl`, `symbol`, `axis`, `lfc_consensus`, `z_stouffer`,
#'   `n_tissues`, `direction`, `in_ventilator`).
#' @examples
#' head(death_signatures("agony"))
#' @export
death_signatures <- function(axes = c("ventilator", "agony", "ischemia")) {
    axes <- match.arg(axes, several.ok = TRUE)
    s <- .signatures()$signatures
    s[s$axis %in% axes, , drop = FALSE]
}

#' Access the GTEx reference of one tissue
#'
#' @param tissue One of `names(thanatome_references)`.
#' @return The reference list of that tissue (`mean`, `sd`, `expressed`,
#'   `poles`, `hardy`, ..., and `qn` with the statistics after quantile
#'   normalisation); see [thanatome_references].
#' @examples
#' death_reference("ovary")$poles
#' @export
death_reference <- function(tissue) {
    refs <- .references()
    tissue <- match.arg(tissue, names(refs))
    refs[[tissue]]
}

#' Map the GTEx Hardy scale to the death axes
#'
#' @param hardy Vector of Hardy scale values (0-4; 0 = ventilator case,
#'   1 = violent and fast death, 2 = fast death of natural causes,
#'   3 = intermediate, 4 = slow death). Integers, numerics, characters or a
#'   factor with those labels are accepted.
#' @return A `data.frame` with `hardy` as integer, the logical labels
#'   `ventilator` (Hardy 0) and `agony` (Hardy 3-4 versus 1-2; `NA` for
#'   ventilator cases, which were not used for that contrast) and the `state`
#'   label (`perfused`, `postmortem_rapid`, `postmortem_agony`).
#' @examples
#' hardy_layers(c(0, 1, 2, 3, 4, NA))
#' hardy_layers(factor(c("0", "4")))
#' @export
hardy_layers <- function(hardy) {
    hardy <- .as_hardy(hardy)
    ventilator <- hardy == 0L
    agony <- ifelse(hardy >= 3L, TRUE, ifelse(hardy >= 1L, FALSE, NA))
    state <- ifelse(ventilator, "perfused",
        ifelse(agony, "postmortem_agony", "postmortem_rapid")
    )
    data.frame(
        hardy = hardy, ventilator = ventilator, agony = agony,
        state = state, stringsAsFactors = FALSE
    )
}

## Hardy as integer 0-4 (accepts a factor or character with those labels)
.as_hardy <- function(hardy) {
    h <- as.character(hardy)
    if (!all(is.na(h) | grepl("^[0-4]$", h))) {
        stop("'hardy' must contain values 0-4 (or NA)")
    }
    as.integer(h)
}
## Lazy loading of the package data (LazyData: false) without touching the
## global environment
.signatures <- function() {
    e <- new.env(parent = emptyenv())
    utils::data("thanatome_signatures", package = "thanatome", envir = e)
    e$thanatome_signatures
}
.references <- function() {
    e <- new.env(parent = emptyenv())
    utils::data("thanatome_references", package = "thanatome", envir = e)
    e$thanatome_references
}
