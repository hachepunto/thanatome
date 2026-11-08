## Asistido por Claude Code (Anthropic); ver README.md, sección Provenance.

#' Annotate each sample with its position along the death axes
#'
#' For every sample, computes a signature score per death axis (mean of the
#' per-gene z-scores against the reference, with the sign of the consensus),
#' the position scaled between the two poles of the reference tissue (0 =
#' ventilated donors / fast deaths, 1 = non-ventilated donors / slow deaths;
#' values outside 0-1 are allowed and meaningful: living tissue falls below 0
#' on the ventilator axis), a state label and technical contamination flags.
#' With `method = "regression"` it adds the coefficients of the joint
#' regression of the sample's z profile on the consensus vectors and the
#' fraction of the sample's deviation from the tissue centroid explained by
#' the axes.
#'
#' **Reference.** A reference of the same tissue is mandatory: a sample has
#' no position on its own, since z-scores are computed per gene against a
#' centroid. The packaged references (`names(thanatome_references)`) were
#' built from GTEx v8 counts quantified by the recount3 Monorail pipeline
#' (GENCODE v26) and apply to data quantified the same way; for another
#' pipeline, platform or protocol pass your own cohort (a
#' `SummarizedExperiment` or count matrix, at least 20 samples) as
#' `reference`: the centroid and, if its `colData` carries a `hardy` (or
#' `DTHHRDY`) column, the poles are estimated from it. Without Hardy the
#' positions and labels are `NA`. Scores are rank-stable across pipelines;
#' positions are not.
#'
#' **Ischemia.** Ischemia is not a third axis. Its signature (86 genes) is
#' mostly a subset of the ventilator signature with the opposite sign, and
#' only 10-23 of its genes are expressed per tissue. It is reported as the
#' column `ischemia`, the component of the ventilator axis sensitive to
#' post-mortem ischemic time, and it tracks ischemic time only among perfused
#' (ventilated) donors (partial Spearman correlation with ischemic time given
#' the ventilator score 0.20-0.50 within GTEx Hardy 0). It never enters the
#' joint regression as a coefficient.
#'
#' **Reliability.** With signatures derived without the scored tissue, the
#' ventilator score recovers the Hardy 0 label with AUC 0.87-0.94 and the
#' agony score the Hardy 3-4 label with AUC 0.69-0.90; the state label is
#' right for 61-84 percent of samples in unseen donors, against 45-80 percent
#' for the majority class, and its weak side is rapid versus agony.
#'
#' @param x A `SummarizedExperiment` with raw counts in `assay(x, assay)`, or
#'   a numeric matrix of raw counts (genes in rows, Ensembl ids with or
#'   without version, samples with names). Set `log = TRUE` if the values are
#'   already log2-CPM (prior count 1, full library size); VST or rlog values
#'   are not comparable with the reference.
#' @param reference Either the name of a GTEx tissue (see
#'   `names(thanatome_references)`), or a `SummarizedExperiment` / count
#'   matrix of the same tissue to be used as cohort reference.
#' @param axes Axes to score: `"ventilator"`, `"agony"` or both (default).
#'   The ischemia component is always reported.
#' @param method `"signature"` (default) or `"regression"`. The regression
#'   adds `beta_ventilator`, `beta_agony`, `r2_signature` (over the signature
#'   genes) and `r2_total` (over all genes of the reference; `NA` unless the
#'   reference is a cohort with at least 10 000 genes). Positions and labels
#'   always come from the signature score.
#' @param assay Name of the assay with counts when `x` is a
#'   `SummarizedExperiment`.
#' @param log Logical; `TRUE` if `x` already holds log2-CPM values.
#' @param lib_size Optional library sizes of the samples of `x` (one per
#'   column) for the CPM normalisation, needed when `x` holds only a subset
#'   of the genes (as the example data do). If `NULL`, a `lib_size` column of
#'   `colData(x)` is used when present, and the column sums otherwise; a
#'   warning is given when `x` has fewer than 10 000 genes and no library
#'   sizes. A cohort reference always uses its own `lib_size` column or its
#'   column sums.
#' @param state_threshold Position above which a sample is labelled
#'   non-perfused (ventilator axis) or in agony (agony axis). Default 0.5.
#' @param flag_threshold Mean log2-CPM above which a contamination flag is
#'   raised; see [contamination_flags()].
#' @return A `DataFrame` (one row per sample) with the signature scores
#'   (`ventilator`, `agony`, `ischemia`), `position_ventilator`,
#'   `position_agony`, `state` (`perfused`, `postmortem_rapid`,
#'   `postmortem_agony`, or `postmortem` when the agony axis is not
#'   available), `flag_pancreas`, `flag_surfactant`, the regression columns
#'   when requested, and `n_genes` used per score. The raw `ventilator` score
#'   grows towards the perfused pole; the position inverts it so that 1 is
#'   post-mortem. `metadata()` records the reference, the method, the poles
#'   and, when `x` carries a `hardy` column, the AUC of each score against
#'   the Hardy label and the agreement of the state label, as an internal
#'   check (in-sample if the samples are GTEx).
#' @examples
#' data(thanatome_example)
#' ctl <- thanatome_example[, thanatome_example$group == "control"]
#' sc <- thanatome_score(ctl, reference = "ovary")
#' head(sc)
#' S4Vectors::metadata(sc)$auc
#' @seealso [select_samples()], [plot_thanatome()], [death_covariates()],
#'   [compare_controls()]
#' @export
thanatome_score <- function(x, reference, axes = c("ventilator", "agony"), method = c("signature", "regression"),
                            assay = "counts", log = FALSE, lib_size = NULL, state_threshold = 0.5, flag_threshold = 4) {
    if (missing(reference) || is.null(reference)) {
        stop("'reference' is required: a tissue name from names(thanatome_references) or a cohort of the same tissue")
    }
    axes <- match.arg(axes, several.ok = TRUE); method <- match.arg(method)
    cd <- .coldata(x)
    m <- .as_logcpm(x, assay = assay, log = log, lib_size = lib_size)
    if (is.null(colnames(m))) stop("'x' needs sample names (column names)")
    ref <- .resolve_reference(reference, assay = assay, log = log)
    sig <- .signatures()$signatures
    ## z por gen frente a la referencia, solo en los genes expresados en ella
    g <- intersect(rownames(m), names(ref$mean)[ref$expressed])
    if (length(g) < 50) stop("fewer than 50 reference genes found in 'x'; check that row names are Ensembl ids")
    if (length(g) < 0.9 * sum(ref$expressed)) warning("only ", length(g), " of ", sum(ref$expressed), " reference genes are in 'x'; positions may be biased")
    z <- (m[g, , drop = FALSE] - ref$mean[g]) / ref$sd[g]
    z[!is.finite(z)] <- 0
    cols <- c(axes, "ischemia")
    scores <- vapply(cols, function(a) .signature_score(z, sig[sig$axis == a, ]), numeric(ncol(z)))
    n_genes <- vapply(cols, function(a) sum(sig$ensembl[sig$axis == a] %in% g), integer(1))
    out <- S4Vectors::DataFrame(matrix(scores, ncol = length(cols), dimnames = list(colnames(m), cols)))
    pos <- .positions(out, ref$poles)
    out$position_ventilator <- pos$ventilator; out$position_agony <- pos$agony
    out$state <- .state(out$position_ventilator, out$position_agony, state_threshold)
    fl <- .flags_from_logcpm(m, threshold = flag_threshold, lung = identical(ref$tissue, "lung"))
    out$flag_pancreas <- fl$pancreas; out$flag_surfactant <- fl$surfactant
    if (method == "regression") {
        rg <- .regression(z, sig[sig$axis %in% c(axes, "ischemia"), ], total_ss = if (ref$cohort && length(ref$mean) >= 10000) .total_ss(m, ref) else NULL)
        out <- cbind(out, rg)
    }
    out$n_genes <- matrix(n_genes, nrow = nrow(out), ncol = length(cols), byrow = TRUE, dimnames = list(NULL, cols))
    meta <- list(reference = ref$label, pipeline = ref$pipeline, method = method, poles = ref$poles,
                 thresholds = .signatures()$thresholds, state_threshold = state_threshold)
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
    if (log) return(m)
    if (is.null(lib_size)) { cd <- .coldata(x); if (!is.null(cd) && "lib_size" %in% colnames(cd)) lib_size <- cd$lib_size }
    if (!is.null(lib_size) && length(lib_size) != ncol(m)) stop("'lib_size' must have one value per sample of 'x'")
    if (is.null(lib_size) && nrow(m) < 10000) warning("'x' has only ", nrow(m), " genes and no library sizes: CPM will be computed over that subset and will not be comparable with the reference")
    edgeR::cpm(m, log = TRUE, prior.count = 1, lib.size = lib_size)
}
.coldata <- function(x) if (is(x, "SummarizedExperiment")) SummarizedExperiment::colData(x) else NULL
.hardy_column <- function(cd) {
    if (is.null(cd)) return(NULL)
    k <- intersect(c("hardy", "DTHHRDY"), colnames(cd))
    if (length(k)) cd[[k[1]]] else NULL
}
## Referencia: tejido del paquete o cohorte del usuario (media y SD por gen; polos si hay Hardy)
.resolve_reference <- function(reference, assay, log) {
    if (is.character(reference)) {
        r <- death_reference(reference)
        return(list(mean = r$mean, sd = r$sd, expressed = r$expressed, poles = r$poles, tissue = r$tissue,
                    pipeline = r$pipeline, label = paste("GTEx", r$tissue), cohort = FALSE))
    }
    m <- .as_logcpm(reference, assay = assay, log = log)
    if (ncol(m) < 20) warning("cohort reference with fewer than 20 samples: centroid and poles will be noisy")
    mu <- rowMeans(m); s <- apply(m, 1, sd); s[!is.finite(s) | s == 0] <- NA
    expressed <- mu > 0 & !is.na(s)
    poles <- NULL
    hardy <- .hardy_column(.coldata(reference))
    if (!is.null(hardy)) {
        sig <- .signatures()$signatures
        z <- (m[expressed, , drop = FALSE] - mu[expressed]) / s[expressed]; z[!is.finite(z)] <- 0
        poles <- .poles_from_scores(.signature_score(z, sig[sig$axis == "ventilator", ]), .signature_score(z, sig[sig$axis == "agony", ]), .as_hardy(hardy))
    } else {
        message("cohort reference without a 'hardy' column: positions and state labels will be NA")
    }
    list(mean = mu, sd = s, expressed = expressed, poles = poles, tissue = NA_character_, pipeline = "cohort",
         label = sprintf("cohort of %d samples", ncol(m)), cohort = TRUE)
}
.poles_from_scores <- function(sv, sa, h) {
    data.frame(axis = c("ventilator", "agony"),
        pole0 = c(median(sv[h == 0], na.rm = TRUE), median(sa[h %in% c(1L, 2L)], na.rm = TRUE)),
        pole1 = c(median(sv[h >= 1], na.rm = TRUE), median(sa[h %in% c(3L, 4L)], na.rm = TRUE)),
        n_pole0 = c(sum(h == 0, na.rm = TRUE), sum(h %in% c(1L, 2L))), n_pole1 = c(sum(h >= 1, na.rm = TRUE), sum(h %in% c(3L, 4L))),
        stringsAsFactors = FALSE)
}
## Score de firma: media de z con el signo del consenso (NA con menos de 5 genes)
.signature_score <- function(z, s) {
    g <- intersect(s$ensembl, rownames(z))
    if (length(g) < 5) return(rep(NA_real_, ncol(z)))
    colMeans(z[g, , drop = FALSE] * sign(s$lfc_consensus[match(g, s$ensembl)]))
}
.positions <- function(out, poles) {
    one <- function(axis) {
        if (is.null(poles) || !axis %in% colnames(out)) return(rep(NA_real_, nrow(out)))
        p <- poles[poles$axis == axis, ]
        if (!is.finite(p$pole0) || !is.finite(p$pole1) || p$pole1 == p$pole0) return(rep(NA_real_, nrow(out)))
        (out[[axis]] - p$pole0) / (p$pole1 - p$pole0)
    }
    list(ventilator = one("ventilator"), agony = one("agony"))
}
.state <- function(pv, pa, thr) {
    ifelse(is.na(pv), NA_character_, ifelse(pv < thr, "perfused", ifelse(is.na(pa), "postmortem", ifelse(pa < thr, "postmortem_rapid", "postmortem_agony"))))
}
## Regresión conjunta del perfil z sobre los vectores de consenso (sin intercepto; isquemia entra como columna pero no se reporta)
.regression <- function(z, s, total_ss = NULL) {
    axes <- unique(s$axis); g <- intersect(unique(s$ensembl), rownames(z))
    B <- vapply(axes, function(a) { b <- s$lfc_consensus[s$axis == a][match(g, s$ensembl[s$axis == a])]; b[is.na(b)] <- 0; b }, numeric(length(g)))
    B <- matrix(B, ncol = length(axes), dimnames = list(g, axes))
    zg <- z[g, , drop = FALSE]
    fit <- lm.fit(B, zg)
    coefs <- matrix(fit$coefficients, nrow = length(axes), dimnames = list(axes, colnames(z)))
    res <- as.matrix(fit$residuals)
    ss_res <- colSums(res^2); ss_sig <- colSums(zg^2)
    out <- S4Vectors::DataFrame(t(coefs[setdiff(axes, "ischemia"), , drop = FALSE]))
    colnames(out) <- paste0("beta_", setdiff(axes, "ischemia"))
    out$r2_signature <- 1 - ss_res / ss_sig
    out$r2_total <- if (is.null(total_ss)) NA_real_ else (ss_sig - ss_res) / total_ss
    out
}
## Suma de cuadrados de z sobre todos los genes de la cohorte (para r2_total)
.total_ss <- function(m, ref) {
    g <- intersect(rownames(m), names(ref$mean)[ref$expressed])
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
