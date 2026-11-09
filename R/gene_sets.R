## Assisted by Claude Code (Anthropic); see README.md, section Provenance.

#' Export the death signatures as gene sets
#'
#' Returns the consensus signatures as gene sets for the gene-set tools of
#' Bioconductor (for example `fgsea`, `limma::camera()`, `GSVA`,
#' `singscore` or `clusterProfiler`). By default every axis is split by
#' direction, because the signatures are signed: `ventilator_up` holds the
#' genes higher in ventilated donors (Hardy 0) than in the rest, `agony_up`
#' the genes higher in slow deaths (Hardy 3-4) than in fast deaths
#' (Hardy 1-2), and `ischemia_up` the genes that rise with post-mortem
#' ischemic time. With `split = FALSE` each axis becomes a single unsigned
#' set, suitable only for tests that ignore direction.
#'
#' @param axes Axes to export: `"ventilator"`, `"agony"` and/or
#'   `"ischemia"`. Default the two axes.
#' @param id Gene identifiers: `"ensembl"` (unversioned, default) or
#'   `"symbol"` (genes without a symbol are dropped).
#' @param split Logical; `TRUE` (default) gives one set per axis and
#'   direction (`<axis>_up`, `<axis>_down`), `FALSE` one set per axis.
#' @param format `"list"` (default, a named list of character vectors),
#'   `"GeneSetCollection"` (package GSEABase) or `"BiocSet"` (package
#'   BiocSet).
#' @return The gene sets in the requested format.
#' @examples
#' gs <- death_gene_sets()
#' lengths(gs)
#' death_gene_sets("ischemia", id = "symbol")$ischemia_up[1:5]
#' if (requireNamespace("GSEABase", quietly = TRUE)) {
#'     death_gene_sets(format = "GeneSetCollection")
#' }
#' @seealso [death_signatures()], [annotate_de()]
#' @export
death_gene_sets <- function(axes = c("ventilator", "agony"),
                            id = c("ensembl", "symbol"), split = TRUE,
                            format = c("list", "GeneSetCollection",
                                        "BiocSet")) {
    axes <- match.arg(axes, c("ventilator", "agony", "ischemia"),
                        several.ok = TRUE)
    id <- match.arg(id)
    format <- match.arg(format)
    s <- death_signatures(axes)
    s$gene <- s[[id]]
    s <- s[!is.na(s$gene) & s$gene != "", , drop = FALSE]
    key <- if (split) paste(s$axis, s$direction, sep = "_") else s$axis
    lv <- if (split) {
        paste(rep(axes, each = 2), c("up", "down"), sep = "_")
    } else {
        axes
    }
    sets <- lapply(split(s$gene, factor(key, levels = lv)), unique)
    sets <- sets[lengths(sets) > 0]
    switch(format,
        list = sets,
        GeneSetCollection = .as_gene_set_collection(sets, id),
        BiocSet = .as_biocset(sets)
    )
}

## ---- internal helpers -------------------------------------------------------

.as_gene_set_collection <- function(sets, id) {
    if (!requireNamespace("GSEABase", quietly = TRUE)) {
        stop("format = \"GeneSetCollection\" needs the GSEABase package")
    }
    id_type <- if (id == "ensembl") {
        GSEABase::ENSEMBLIdentifier()
    } else {
        GSEABase::SymbolIdentifier()
    }
    gs <- lapply(names(sets), function(n) {
        GSEABase::GeneSet(sets[[n]], setName = n, geneIdType = id_type,
                        shortDescription = paste("thanatome", n))
    })
    GSEABase::GeneSetCollection(gs)
}
.as_biocset <- function(sets) {
    if (!requireNamespace("BiocSet", quietly = TRUE)) {
        stop("format = \"BiocSet\" needs the BiocSet package")
    }
    do.call(BiocSet::BiocSet, sets)
}
