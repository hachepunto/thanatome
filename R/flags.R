## Asistido por Claude Code (Anthropic); ver README.md, sección Provenance.

#' Technical contamination flags per sample
#'
#' Flags samples whose mean z-score over the pancreatic enzyme genes
#' (PNLIP, CELA3A, PRSS1, CPA1, CELA2A, CTRB1; contamination of neighbouring
#' tissue, Nieuwenhuis et al. 2020) or over the surfactant genes (SFTPC,
#' SFTPB, SFTPA1, SFTPA2) exceeds `threshold`. The surfactant flag is not
#' raised for lung, where those genes are expressed physiologically.
#'
#' @param z Matrix of per-gene z-scores (genes in rows, Ensembl ids), as
#'   computed inside [thanatome_score()].
#' @param flags The flag gene table (default `thanatome_firmas$flags`).
#' @param threshold Mean z above which the flag is raised. Default 2.
#' @param lung Logical; `TRUE` disables the surfactant flag.
#' @return A `data.frame` with logical columns `pancreas` and `surfactant`
#'   and the mean z of each set (`pancreas_z`, `surfactant_z`).
#' @examples
#' data(thanatome_firmas)
#' fl <- thanatome_firmas$flags
#' z <- matrix(rnorm(5 * nrow(fl)), nrow(fl), 5, dimnames = list(fl$ensembl, NULL))
#' contamination_flags(z, flags = fl)
#' @export
contamination_flags <- function(z, flags = NULL, threshold = 2, lung = FALSE) {
    if (is.null(flags)) flags <- .firmas()$flags
    one <- function(f) {
        g <- intersect(flags$ensembl[flags$flag == f], rownames(z))
        if (!length(g)) return(rep(NA_real_, ncol(z)))
        colMeans(z[g, , drop = FALSE])
    }
    pz <- one("pancreas"); sz <- one("surfactant")
    data.frame(pancreas = !is.na(pz) & pz > threshold, surfactant = !lung & !is.na(sz) & sz > threshold,
               pancreas_z = pz, surfactant_z = sz, row.names = colnames(z), stringsAsFactors = FALSE)
}
