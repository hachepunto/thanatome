## Asistido por Claude Code (Anthropic); ver README.md, sección Provenance.

#' Access the consensus signatures
#'
#' @param axes Which axes to return; any of `"ventilator"`, `"ischemia"`,
#'   `"agony"`.
#' @return A `data.frame` with the signature genes of the requested axes
#'   (see [thanatome_firmas]).
#' @examples
#' head(thanatome_firmas(axes = "agony"))
#' @export
thanatome_firmas <- function(axes = c("ventilator", "ischemia", "agony")) {
    axes <- match.arg(axes, several.ok = TRUE)
    s <- .firmas()$signatures
    s[s$axis %in% axes, , drop = FALSE]
}

#' Access the GTEx reference of one tissue
#'
#' @param tissue One of the ten tissue names in [thanatome_referencias]
#'   (`names(thanatome_referencias)`).
#' @return The reference list of that tissue (`mean`, `sd`, `poles`, ...).
#' @examples
#' thanatome_referencia("ovary")$poles
#' @export
thanatome_referencia <- function(tissue) {
    refs <- .referencias()
    tissue <- match.arg(tissue, names(refs))
    refs[[tissue]]
}

#' Map the GTEx Hardy scale to the death axes
#'
#' @param hardy Integer vector of Hardy scale values (0-4; 0 = ventilator
#'   case, 1 = violent and fast death, 2 = fast death of natural causes,
#'   3 = intermediate, 4 = slow death).
#' @return A `data.frame` with the logical labels used to derive the axes:
#'   `ventilator` (Hardy 0), `agony` (Hardy 3-4 versus 1-2; `NA` for
#'   ventilator cases, which were not used for this contrast) and the
#'   `state` label (`perfused`, `postmortem_rapid`, `postmortem_agony`).
#' @examples
#' hardy_layers(c(0, 1, 2, 3, 4, NA))
#' @export
hardy_layers <- function(hardy) {
    hardy <- as.integer(hardy)
    ventilator <- hardy == 0L
    agony <- ifelse(hardy >= 3L, TRUE, ifelse(hardy >= 1L, FALSE, NA))
    state <- ifelse(ventilator, "perfused", ifelse(agony, "postmortem_agony", "postmortem_rapid"))
    data.frame(hardy = hardy, ventilator = ventilator, agony = agony, state = state, stringsAsFactors = FALSE)
}

## Carga perezosa de los datos del paquete (LazyData: false) sin tocar el entorno global
.firmas <- function() {
    e <- new.env(); utils::data("thanatome_firmas", package = "thanatome", envir = e); e$thanatome_firmas
}
.referencias <- function() {
    e <- new.env(); utils::data("thanatome_referencias", package = "thanatome", envir = e); e$thanatome_referencias
}
