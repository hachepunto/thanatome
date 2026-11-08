## Asistido por Claude Code (Anthropic); ver README.md, sección Provenance.

#' Annotate each sample with its position along the death axes
#'
#' For every sample, computes a signature score per death axis (mean of the
#' per-gene z-scores, with the sign of the consensus), the position scaled
#' between the two poles of the reference tissue (0 = ventilated / fast
#' deaths, 1 = non-ventilated / slow deaths; values outside 0-1 are allowed
#' and meaningful: living tissue falls below 0 on the ventilator axis), a
#' state label, technical contamination flags and, with
#' `method = "regression"`, the coefficients of the joint regression of the
#' sample's z profile on the consensus vectors and the fraction of the
#' sample's deviation from the tissue centroid explained by the axes.
#'
#' A reference of the same tissue is mandatory: a sample has no position on
#' its own, since z-scores are computed per gene against a centroid. Use the
#' name of one of the ten GTEx tissues for GTEx-like data, or your own cohort
#' (a `SummarizedExperiment` or count matrix with at least 20 samples) when
#' the platform or protocol differs; in that case the poles are estimated
#' from the cohort if its `colData` carries a `hardy` (or `DTHHRDY`) column,
#' and are `NA` otherwise.
#'
#' Ischemia is always reported as the marginal signature score, never as a
#' joint regression coefficient: in GTEx ventilator cases and ischemic time
#' are confounded and the joint coefficient loses its meaning. The ischemia
#' score tracks ischemic time only among perfused (ventilated) donors.
#'
#' @param x A `SummarizedExperiment` with raw counts in `assay(x, assay)`,
#'   or a numeric matrix of raw counts (genes in rows, Ensembl ids with or
#'   without version). Set `log = TRUE` if the values are already
#'   log2-normalised.
#' @param reference Either the name of a GTEx tissue (see
#'   `names(thanatome_referencias)`), or a `SummarizedExperiment` / count
#'   matrix of the same tissue to be used as cohort reference.
#' @param axes Axes to score: any of `"ventilator"`, `"ischemia"`,
#'   `"agony"`.
#' @param method `"signature"` (default) or `"regression"`. The regression
#'   adds the columns `beta_*`, `r2_signature` and `r2_total` (the latter
#'   only with a cohort reference, which covers all genes); positions and
#'   labels always come from the signature score.
#' @param assay Name of the assay with counts when `x` is a
#'   `SummarizedExperiment`.
#' @param log Logical; `TRUE` if `x` already holds log2-normalised values.
#' @param lib_size Optional library sizes (one per sample) for the CPM
#'   normalisation, needed when `x` holds only a subset of the genes (as the
#'   example data do). If `NULL`, a `lib_size` column of `colData(x)` is used
#'   when present, and the column sums otherwise.
#' @param state_threshold Position above which a sample is labelled
#'   non-perfused (ventilator axis) or in agony (agony axis). Default 0.5.
#' @return A `DataFrame` (one row per sample) with the signature scores
#'   (`ventilator`, `ischemia`, `agony`), `position_ventilator`,
#'   `position_agony`, `state`, `flag_pancreas`, `flag_surfactant`, the
#'   regression columns when requested, and `n_genes` used per axis. Its
#'   `metadata()` records the reference, the method and, when the input
#'   carries a `hardy` column, the AUC of each score against the Hardy label
#'   as an internal check.
#' @examples
#' data(thanatome_example)
#' sc <- thanatome_score(thanatome_example, reference = "ovary")
#' head(sc)
#' S4Vectors::metadata(sc)$auc
#' @seealso [select_samples()], [plot_thanatome()], [death_covariates()]
#' @export
thanatome_score <- function(x, reference, axes = c("ventilator", "ischemia", "agony"),
                            method = c("signature", "regression"), assay = "counts",
                            log = FALSE, lib_size = NULL, state_threshold = 0.5) {
    if (missing(reference) || is.null(reference)) {
        stop("'reference' is required: a tissue name from names(thanatome_referencias) or a cohort of the same tissue")
    }
    axes <- match.arg(axes, several.ok = TRUE); method <- match.arg(method)
    cd <- .coldata(x)
    m <- .as_logcpm(x, assay = assay, log = log, lib_size = lib_size)
    ref <- .resolve_reference(reference, assay = assay, log = log, lib_size = lib_size)
    firmas <- .firmas()
    ## z por gen frente a la referencia, en los genes comunes
    g <- intersect(rownames(m), names(ref$mean))
    if (length(g) < 50) stop("fewer than 50 reference genes found in 'x'; check that row names are Ensembl ids")
    z <- (m[g, , drop = FALSE] - ref$mean[g]) / ref$sd[g]
    z[!is.finite(z)] <- 0
    sig <- firmas$signatures
    scores <- vapply(axes, function(a) .signature_score(z, sig[sig$axis == a, ]), numeric(ncol(z)))
    n_genes <- vapply(axes, function(a) sum(sig$ensembl[sig$axis == a] %in% g), integer(1))
    out <- S4Vectors::DataFrame(matrix(scores, ncol = length(axes), dimnames = list(colnames(m), axes)))
    ## Posición escalada por los polos de la referencia
    pos <- .positions(out, ref$poles)
    out$position_ventilator <- pos$ventilator; out$position_agony <- pos$agony
    out$state <- .state(out$position_ventilator, out$position_agony, state_threshold)
    fl <- contamination_flags(z, firmas$flags, lung = identical(ref$tissue, "lung"))
    out$flag_pancreas <- fl$pancreas; out$flag_surfactant <- fl$surfactant
    if (method == "regression") {
        rg <- .regression(z, sig[sig$axis %in% axes, ], total_ss = if (ref$cohort) .total_ss(m, ref) else NULL)
        out <- cbind(out, rg)
    }
    out$n_genes <- matrix(n_genes, nrow = nrow(out), ncol = length(axes), byrow = TRUE, dimnames = list(NULL, axes))
    meta <- list(reference = ref$label, method = method, poles = ref$poles, thresholds = firmas$thresholds,
                 state_threshold = state_threshold)
    hardy <- .hardy_column(cd)
    if (!is.null(hardy)) {
        hl <- hardy_layers(hardy)
        meta$auc <- c(ventilator = if ("ventilator" %in% axes) .auc(hl$ventilator, out$ventilator) else NA,
                      agony = if ("agony" %in% axes) .auc(hl$agony, out$agony) else NA)
        meta$state_agreement <- mean(out$state == hl$state, na.rm = TRUE)
    }
    S4Vectors::metadata(out) <- meta
    out
}

## ---- auxiliares internos ----------------------------------------------------

## Matriz log2-CPM (prior 1) con Ensembl sin versión
.as_logcpm <- function(x, assay = "counts", log = FALSE, lib_size = NULL) {
    m <- if (is(x, "SummarizedExperiment")) SummarizedExperiment::assay(x, assay) else as.matrix(x)
    if (is.null(rownames(m))) stop("'x' needs gene ids (Ensembl) as row names")
    rownames(m) <- sub("\\..*$", "", rownames(m))
    m <- m[!duplicated(rownames(m)), , drop = FALSE]
    if (is.null(lib_size)) { cd <- .coldata(x); if (!is.null(cd) && "lib_size" %in% colnames(cd)) lib_size <- cd$lib_size }
    if (!log) m <- edgeR::cpm(m, log = TRUE, prior.count = 1, lib.size = lib_size)
    m
}
.coldata <- function(x) if (is(x, "SummarizedExperiment")) SummarizedExperiment::colData(x) else NULL
.hardy_column <- function(cd) {
    if (is.null(cd)) return(NULL)
    k <- intersect(c("hardy", "DTHHRDY"), colnames(cd))
    if (length(k)) cd[[k[1]]] else NULL
}
## Referencia: tejido del paquete o cohorte del usuario (media y SD por gen; polos si hay Hardy)
.resolve_reference <- function(reference, assay, log, lib_size = NULL) {
    if (is.character(reference)) {
        r <- thanatome_referencia(reference)
        return(list(mean = r$mean, sd = r$sd, poles = r$poles, tissue = r$tissue, label = paste("GTEx", r$tissue), cohort = FALSE))
    }
    m <- .as_logcpm(reference, assay = assay, log = log, lib_size = if (identical(nrow(reference), length(lib_size))) NULL else lib_size)
    if (ncol(m) < 20) warning("cohort reference with fewer than 20 samples: centroid and poles will be noisy")
    mu <- rowMeans(m); s <- apply(m, 1, sd); s[!is.finite(s) | s == 0] <- NA
    poles <- NULL
    hardy <- .hardy_column(.coldata(reference))
    if (!is.null(hardy)) {
        firmas <- .firmas(); sig <- firmas$signatures
        z <- (m - mu) / s; z[!is.finite(z)] <- 0
        sv <- .signature_score(z, sig[sig$axis == "ventilator", ]); sa <- .signature_score(z, sig[sig$axis == "agony", ])
        h <- as.integer(hardy)
        poles <- data.frame(axis = c("ventilator", "agony"),
            pole0 = c(median(sv[h == 0], na.rm = TRUE), median(sa[h %in% c(1L, 2L)], na.rm = TRUE)),
            pole1 = c(median(sv[h >= 1], na.rm = TRUE), median(sa[h %in% c(3L, 4L)], na.rm = TRUE)),
            n_pole0 = c(sum(h == 0, na.rm = TRUE), sum(h %in% c(1L, 2L))), n_pole1 = c(sum(h >= 1, na.rm = TRUE), sum(h %in% c(3L, 4L))),
            stringsAsFactors = FALSE)
    } else {
        message("cohort reference without a 'hardy' column: positions and state labels will be NA")
    }
    list(mean = mu, sd = s, poles = poles, tissue = NA_character_, label = sprintf("cohort of %d samples", ncol(m)), cohort = TRUE)
}
## Score de firma: media de z con el signo del consenso
.signature_score <- function(z, s) {
    g <- intersect(s$ensembl, rownames(z))
    if (length(g) < 5) return(rep(NA_real_, ncol(z)))
    colMeans(z[g, , drop = FALSE] * sign(s$lfc_consensus[match(g, s$ensembl)]))
}
.positions <- function(out, poles) {
    one <- function(axis) {
        if (is.null(poles) || !axis %in% colnames(out)) return(rep(NA_real_, nrow(out)))
        p <- poles[poles$axis == axis, ]
        (out[[axis]] - p$pole0) / (p$pole1 - p$pole0)
    }
    list(ventilator = one("ventilator"), agony = one("agony"))
}
.state <- function(pv, pa, thr) {
    ifelse(is.na(pv), NA_character_, ifelse(pv < thr, "perfused", ifelse(is.na(pa) | pa < thr, "postmortem_rapid", "postmortem_agony")))
}
## Regresión conjunta del perfil z sobre los vectores de consenso (sin intercepto)
.regression <- function(z, s, total_ss = NULL) {
    axes <- unique(s$axis); g <- intersect(unique(s$ensembl), rownames(z))
    B <- vapply(axes, function(a) { b <- s$lfc_consensus[s$axis == a][match(g, s$ensembl[s$axis == a])]; b[is.na(b)] <- 0; b }, numeric(length(g)))
    B <- matrix(B, ncol = length(axes), dimnames = list(g, axes))
    zg <- z[g, , drop = FALSE]
    fit <- lm.fit(B, zg)
    beta <- t(matrix(fit$coefficients, nrow = length(axes), dimnames = list(paste0("beta_", axes), colnames(z))))
    ss_res <- colSums(fit$residuals^2); ss_sig <- colSums(zg^2)
    out <- S4Vectors::DataFrame(beta, r2_signature = 1 - ss_res / ss_sig,
                                r2_total = if (is.null(total_ss)) NA_real_ else (ss_sig - ss_res) / total_ss)
    out$beta_ischemia <- NULL   # isquemia solo como score marginal (ver detalles)
    out
}
## Suma de cuadrados de z sobre todos los genes de la cohorte (para r2_total)
.total_ss <- function(m, ref) {
    g <- intersect(rownames(m), names(ref$mean))
    z <- (m[g, , drop = FALSE] - ref$mean[g]) / ref$sd[g]; z[!is.finite(z)] <- 0
    colSums(z^2)
}
## AUC por Mann-Whitney (sin dependencias)
.auc <- function(y, s) {
    ok <- !is.na(y) & !is.na(s); y <- as.logical(y[ok]); s <- s[ok]
    n1 <- sum(y); n0 <- sum(!y)
    if (n1 < 2 || n0 < 2) return(NA_real_)
    r <- rank(s)
    (sum(r[y]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
