## Assisted by Claude Code (Anthropic); see README.md, section Provenance.

#' Annotate a master-regulator table with the death axes
#'
#' For each regulator, reports which fraction of its regulon falls in each
#' death axis, whether that fraction is larger than expected (Fisher test
#' against the regulon universe) and a death-direction score: the mean, over
#' the targets that belong to the axis, of the product of the regulon mode
#' (sign of `tfmode`) and the sign of the consensus log fold change. A score
#' near +1 means the regulator's inferred activity rises along the axis
#' (its activated targets go up with the death process), near -1 that it
#' falls; regulators with a significant enrichment are flagged as ones whose
#' differential activity may reflect the state of the controls rather than
#' the cases.
#'
#' @param nes A named numeric vector of regulator activities (NES from viper,
#'   NaRnEA or similar), or a `data.frame` with columns `regulator` and
#'   `nes`.
#' @param regulons A list named by regulator. Each element is either a list
#'   with a named numeric `tfmode` (viper format; `likelihood` is ignored),
#'   a named numeric vector of modes, or a character vector of target gene
#'   ids (mode +1). Gene ids are Ensembl, with or without version; HGNC
#'   symbols are mapped to Ensembl with `org.Hs.eg.db` (Suggests) when the
#'   ids do not look like Ensembl, with a message, and the universe is
#'   mapped the same way.
#' @param axes Axes to annotate (default the two axes; `"ischemia"` adds the
#'   component).
#' @param min_targets Regulators with fewer targets in the universe are
#'   dropped. Default 10.
#' @param universe Gene universe for the Fisher test; by default the union
#'   of all regulon targets.
#' @return A `DataFrame` with one row per regulator: `nes`, `n_targets`,
#'   and per axis `frac_<axis>`, `fisher_padj_<axis>`, `direction_<axis>`,
#'   plus `flag_control_state`.
#' @examples
#' sig <- death_signatures("ventilator")
#' up <- sig$ensembl[sig$direction == "up"][1:30]
#' other <- paste0("ENSG", seq_len(30))
#' regs <- list(TF_vent = list(tfmode = setNames(rep(1, 30), up)),
#'              TF_other = list(tfmode = setNames(rep(1, 30), other)))
#' annotate_mra(c(TF_vent = 3.2, TF_other = -1.1), regs)
#' @export
annotate_mra <- function(nes, regulons, axes = c("ventilator", "agony"),
                        min_targets = 10, universe = NULL) {
    axes <- match.arg(axes, c("ventilator", "agony", "ischemia"),
                        several.ok = TRUE)
    if (is.data.frame(nes)) nes <- setNames(nes$nes, nes$regulator)
    modes <- lapply(regulons, .regulon_modes)
    all_ids <- unique(unlist(lapply(modes, names)))
    ## HGNC symbols: mapped with org.Hs.eg.db (regulons and universe)
    if (mean(grepl("^ENSG", all_ids)) < 0.5) {
        message("annotate_mra: target ids do not look like Ensembl; ",
                "mapping symbols with org.Hs.eg.db")
        modes <- lapply(modes, .modes_to_ensembl)
        if (!is.null(universe)) {
            universe <- .symbols_to_ensembl(universe)
            universe <- universe[!is.na(universe)]
        }
    }
    if (is.null(universe)) universe <- unique(unlist(lapply(modes, names)))
    universe <- sub("\\..*$", "", universe)
    modes <- lapply(modes, function(tf) tf[names(tf) %in% universe])
    modes <- modes[lengths(modes) >= min_targets]
    if (!length(modes)) {
        stop("no regulon has at least ", min_targets,
            " targets in the universe (after id mapping)")
    }
    sig <- death_signatures(axes)
    ann <- do.call(cbind, lapply(axes, function(a) {
        .axis_enrichment(modes, sig[sig$axis == a, ], universe, a)
    }))
    res <- S4Vectors::DataFrame(regulator = names(modes),
                                nes = unname(nes[names(modes)]),
                                n_targets = lengths(modes), ann)
    padj_cols <- grep("^fisher_padj_", colnames(res), value = TRUE)
    padj <- as.matrix(as.data.frame(res[, padj_cols, drop = FALSE]))
    res$flag_control_state <- apply(padj <= 0.05, 1, any)
    rownames(res) <- NULL
    res
}

## ---- internal helpers -------------------------------------------------------

## Named vector of target modes from any of the accepted regulon formats
.regulon_modes <- function(r) {
    tf <- if (is.list(r)) {
        r$tfmode
    } else if (is.numeric(r) && !is.null(names(r))) {
        r
    } else {
        setNames(rep(1, length(r)), r)
    }
    names(tf) <- sub("\\..*$", "", names(tf))
    tf
}
.modes_to_ensembl <- function(tf) {
    e <- .symbols_to_ensembl(names(tf))
    keep <- !is.na(e)
    setNames(tf[keep], e[keep])
}
## Per regulator: fraction of targets in the axis, Fisher enrichment against
## the universe and death-direction score
.axis_enrichment <- function(modes, s, universe, a) {
    n_axis <- sum(universe %in% s$ensembl)
    n_out <- length(universe) - n_axis
    rows <- lapply(modes, function(tf) {
        hit <- names(tf) %in% s$ensembl
        k <- sum(hit)
        tab <- matrix(c(k, length(tf) - k, n_axis - k,
                        n_out - (length(tf) - k)), 2)
        p <- stats::fisher.test(tab, alternative = "greater")$p.value
        lfc_sign <- sign(s$lfc_consensus[match(names(tf)[hit], s$ensembl)])
        dir <- if (k) mean(sign(tf[hit]) * lfc_sign) else NA_real_
        c(frac = k / length(tf), p = p, direction = dir)
    })
    out <- S4Vectors::DataFrame(do.call(rbind, rows))
    out$padj <- stats::p.adjust(out$p, "BH")
    out$p <- NULL
    colnames(out) <- paste0(c("frac_", "direction_", "fisher_padj_"), a)
    out
}
## HGNC symbols to Ensembl with org.Hs.eg.db (clear error if not installed)
.symbols_to_ensembl <- function(symbols) {
    if (!requireNamespace("org.Hs.eg.db", quietly = TRUE) ||
        !requireNamespace("AnnotationDbi", quietly = TRUE)) {
        stop("gene ids are not Ensembl; install org.Hs.eg.db to map ",
            "symbols, or provide Ensembl ids")
    }
    db <- org.Hs.eg.db::org.Hs.eg.db
    keys <- intersect(unique(symbols),
                        AnnotationDbi::keys(db, keytype = "SYMBOL"))
    if (!length(keys)) return(rep(NA_character_, length(symbols)))
    e <- suppressMessages(AnnotationDbi::mapIds(db, keys = keys,
                                                column = "ENSEMBL",
                                                keytype = "SYMBOL",
                                                multiVals = "first"))
    unname(e[symbols])
}
