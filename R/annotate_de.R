## Asistido por Claude Code (Anthropic); ver README.md, sección Provenance.

#' Coerce a differential-expression result to a standard table
#'
#' Accepts a `DESeqResults` or any `data.frame`-like object from DESeq2,
#' `limma::topTable()` or `edgeR::topTags()` and returns a `DataFrame` with
#' the columns `gene` (Ensembl id without version), `log2FoldChange` and
#' `padj`, keeping the other columns.
#'
#' @param de The result object or table.
#' @return A `DataFrame`.
#' @examples
#' de <- data.frame(logFC = c(1.2, -0.8), adj.P.Val = c(0.01, 0.2),
#'                  row.names = c("ENSG00000166189.7", "ENSG00000172409.5"))
#' as_de_table(de)
#' @export
as_de_table <- function(de) {
    if (is(de, "TopTags")) de <- de$table
    d <- S4Vectors::DataFrame(as.data.frame(de))
    gene <- if ("gene" %in% colnames(d)) d$gene else rownames(d)
    if (is.null(gene)) stop("gene ids must be in the row names or in a 'gene' column")
    d$gene <- sub("\\..*$", "", gene)
    lfc <- intersect(c("log2FoldChange", "logFC"), colnames(d))
    adj <- intersect(c("padj", "adj.P.Val", "FDR"), colnames(d))
    if (!length(lfc) || !length(adj)) stop("could not find the log fold change or adjusted p-value column")
    d$log2FoldChange <- d[[lfc[1]]]; d$padj <- d[[adj[1]]]
    d
}

#' Annotate a differential-expression table with the death axes
#'
#' For each gene, adds the axis and direction in which it falls in the
#' consensus signatures, the consensus log fold change, and the change in
#' the tumour-versus-control log fold change expected when the control
#' subset moves from non-ventilated to ventilated donors (ventilator axis)
#' or from fast to slow deaths (agony axis). Within GTEx that change equals
#' the contrast between the two subsets of controls, which the consensus of
#' the other tissues predicts with a slope close to -1; the prediction is
#' therefore `-lfc_consensus`. The flag `sign_depends_on_control` marks the
#' genes whose predicted change exceeds their own log fold change.
#'
#' @param de A DE result accepted by [as_de_table()].
#' @param axes Axes to annotate.
#' @param min_tissues Minimum number of concordant tissues for a gene to be
#'   annotated (default: the signature threshold, 5).
#' @return The `DataFrame` from [as_de_table()] with the extra columns
#'   `axis_<axis>`, `lfc_consensus_<axis>`, `delta_lfc_<axis>` and
#'   `sign_depends_on_control`.
#' @examples
#' de <- data.frame(logFC = c(1.2, -0.8, 0.3), adj.P.Val = c(0.01, 0.2, 0.5),
#'                  row.names = c("ENSG00000166189", "ENSG00000172409", "ENSG00000000003"))
#' annotate_de(de)
#' @export
annotate_de <- function(de, axes = c("ventilator", "ischemia", "agony"), min_tissues = NULL) {
    axes <- match.arg(axes, several.ok = TRUE)
    d <- as_de_table(de)
    sig <- thanatome_firmas(axes)
    if (!is.null(min_tissues)) sig <- sig[sig$n_tissues >= min_tissues, ]
    cols <- lapply(axes, function(a) {
        s <- sig[sig$axis == a, ]; i <- match(d$gene, s$ensembl)
        S4Vectors::DataFrame(axis = s$direction[i], lfc_consensus = s$lfc_consensus[i], delta_lfc = -s$lfc_consensus[i])
    })
    names(cols) <- axes
    ann <- do.call(cbind, lapply(axes, function(a) { x <- cols[[a]]; colnames(x) <- paste0(colnames(x), "_", a); x }))
    d <- cbind(d, ann)
    dm <- do.call(cbind, lapply(cols, function(x) x$delta_lfc)); dm[is.na(dm)] <- 0
    d$sign_depends_on_control <- apply(abs(dm) > abs(d$log2FoldChange), 1, any)
    d
}
