## Asistido por Claude Code (Anthropic); ver README.md, sección Provenance.

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
#'   symbols are mapped to Ensembl through the signature table when no
#'   Ensembl id is found, with a message.
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
#' regs <- list(TF_vent = list(tfmode = setNames(rep(1, 30), up)),
#'              TF_other = list(tfmode = setNames(rep(1, 30), paste0("ENSG", seq_len(30)))))
#' annotate_mra(c(TF_vent = 3.2, TF_other = -1.1), regs)
#' @export
annotate_mra <- function(nes, regulons, axes = c("ventilator", "agony"), min_targets = 10, universe = NULL) {
    axes <- match.arg(axes, c("ventilator", "agony", "ischemia"), several.ok = TRUE)
    if (is.data.frame(nes)) nes <- setNames(nes$nes, nes$regulator)
    modes <- lapply(regulons, function(r) {
        tf <- if (is.list(r)) r$tfmode else if (is.numeric(r) && !is.null(names(r))) r else setNames(rep(1, length(r)), r)
        names(tf) <- sub("\\..*$", "", names(tf)); tf
    })
    all_ids <- unique(unlist(lapply(modes, names)))
    if (mean(grepl("^ENSG", all_ids)) < 0.5) {   # símbolos HGNC: se mapean por la tabla de firmas
        message("annotate_mra: target ids do not look like Ensembl; mapping symbols through the signature table")
        full <- .signatures()$signatures; map <- setNames(full$ensembl, full$symbol)
        modes <- lapply(modes, function(tf) { e <- unname(map[names(tf)]); keep <- !is.na(e); setNames(tf[keep], e[keep]) })
    }
    if (is.null(universe)) universe <- unique(unlist(lapply(modes, names)))
    universe <- sub("\\..*$", "", universe)
    modes <- lapply(modes, function(tf) tf[names(tf) %in% universe])
    modes <- modes[vapply(modes, length, integer(1)) >= min_targets]
    if (!length(modes)) stop("no regulon has at least ", min_targets, " targets in the universe (after id mapping)")
    sig <- death_signatures(axes)
    per_axis <- function(a) {
        s <- sig[sig$axis == a, ]; in_axis <- universe %in% s$ensembl
        rows <- lapply(modes, function(tf) {
            hit <- names(tf) %in% s$ensembl; k <- sum(hit)
            p <- stats::fisher.test(matrix(c(k, length(tf) - k, sum(in_axis) - k, sum(!in_axis) - (length(tf) - k)), 2), alternative = "greater")$p.value
            dir <- if (k) mean(sign(tf[hit]) * sign(s$lfc_consensus[match(names(tf)[hit], s$ensembl)])) else NA_real_
            c(frac = k / length(tf), p = p, direction = dir)
        })
        out <- S4Vectors::DataFrame(do.call(rbind, rows))
        out$padj <- stats::p.adjust(out$p, "BH"); out$p <- NULL
        colnames(out) <- paste0(c("frac_", "direction_", "fisher_padj_"), a)
        out
    }
    ann <- do.call(cbind, lapply(axes, per_axis))
    res <- S4Vectors::DataFrame(regulator = names(modes), nes = unname(nes[names(modes)]),
                                n_targets = vapply(modes, length, integer(1)), ann)
    padj_cols <- grep("^fisher_padj_", colnames(res), value = TRUE)
    res$flag_control_state <- apply(as.matrix(as.data.frame(res[, padj_cols, drop = FALSE])) <= 0.05, 1, any)
    rownames(res) <- NULL
    res
}
