## Assisted by Claude Code (Anthropic); see README.md, section Provenance.

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
#' **External cohorts.** Against a packaged reference, a cohort quantified
#' or sequenced differently carries a batch effect. Two safeguards apply.
#' Genes absent from `x` (zero counts in most of its samples) are left out
#' by default, since an absent gene is not a down-regulated
#' one (`drop_unexpressed`). When `x` holds the whole transcriptome, the
#' samples are first quantile-normalised to the distribution of the
#' reference tissue and scored against the statistics and poles of the
#' equally normalised GTEx samples (`normalize`, by default). Without it,
#' autopsies of external cohorts could land on the perfused side of the
#' ventilator axis; within GTEx it barely changes the scores. In external
#' GTEx-independent cohorts (living colon biopsies, accident autopsies) the
#' side of the poles on which a quantile-normalised cohort falls was stable
#' across gene filters and batch handling, whereas its distance beyond the
#' poles was not: read positions outside 0-1 as a direction, not a
#' magnitude.
#'
#' **Ischemia.** Ischemia is not a third axis. Its signature (86 genes) is
#' mostly a subset of the ventilator signature with the opposite sign, and
#' only 10-23 of its genes were expressed per tissue in the analysis (15-29
#' in the packaged references). It is reported as the
#' column `ischemia`, the component of the ventilator axis sensitive to
#' post-mortem ischemic time, and it tracks ischemic time only among perfused
#' (ventilated) donors (partial Spearman correlation with ischemic time given
#' the ventilator score 0.20-0.50 within GTEx Hardy 0). It never enters the
#' joint regression as a coefficient.
#'
#' **Reliability.** With signatures derived without the scored tissue, the
#' ventilator score recovers the Hardy 0 label with AUC 0.87-0.94 and the
#' agony score the Hardy 3-4 label with AUC 0.69-0.90. The state label built
#' from the joint-regression positions, in tissue not seen by the
#' signatures, was right for 61-84 percent of samples against 45-80 percent
#' for the majority class; its weak side is rapid versus agony.
#'
#' **Extreme genes.** The per-gene z-scores are not truncated. Clipping
#' them at +/-5 or +/-3, or replacing the mean by the median, left the
#' leave-one-tissue-out AUC unchanged within GTEx (maximum |z| 7.4) and did
#' not improve the discrimination of external cohorts (|z| up to 22);
#' clipping at 3 only compressed the positions beyond the poles, and the
#' median lost part of the ischemia signal.
#'
#' @param x A `SummarizedExperiment` with raw counts in `assay(x, assay)`, an
#'   edgeR `DGEList` (its library sizes are used; normalisation factors are
#'   ignored, as in the references), a limma `EList` (its log2-CPM `E` is
#'   used) or a numeric matrix of raw counts (genes in rows, Ensembl ids with
#'   or without version, unique sample names). Set `log = TRUE` if the values
#'   are already log2-CPM (`edgeR::cpm(log = TRUE, prior.count = 1)` on the
#'   full library size; edgeR scales the prior count by the library size, so
#'   lowly expressed genes depend slightly on the depth of the cohort); VST
#'   or rlog values are not comparable with the reference.
#' @param reference Either the name of a GTEx tissue (see
#'   `names(thanatome_references)`), or a `SummarizedExperiment` / count
#'   matrix of the same tissue to be used as cohort reference.
#' @param axes Axes to score: `"ventilator"`, `"agony"` or both (default).
#'   The ischemia component is always reported.
#' @param method `"signature"` (default) or `"regression"`. The regression
#'   adds `beta_ventilator`, `beta_agony`, `r2_signature` (over the signature
#'   genes) and `r2_total` (over all genes; `NA` unless both the reference,
#'   as a cohort, and `x` carry at least 10 000 genes). Positions and labels
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
#' @param normalize `"auto"` (default), `"quantile"` or `"none"`.
#'   `"quantile"` quantile-normalises the samples to the reference
#'   distribution before scoring (see *External cohorts*) and requires at
#'   least 10 000 genes in `x`; `"auto"` uses it when `x` has that many
#'   genes and `"none"` otherwise (subsets of genes, such as the example
#'   data). The choice made is stored in `metadata(result)$normalize`.
#' @param drop_unexpressed Logical; `TRUE` (default) leaves out the genes
#'   absent from `x`: zero counts in more than half of its samples (for log
#'   input, log2-CPM below log2(0.1) in more than half), when `x` has at
#'   least 5 samples (with fewer, an absent gene cannot be told apart from a
#'   lowly expressed one and all genes are kept). The criterion is far from
#'   the reference's own expression cut, so a sample's score does not depend
#'   on the other samples through it. The genes left out are listed in
#'   `metadata()$unexpressed_signature_genes`.
#' @param tissue Optional tissue of the samples, used only to switch off
#'   the contamination flag that is physiological there (surfactant in
#'   lung, pancreatic enzymes in pancreas); by default the tissue of the
#'   packaged reference, and none for a cohort reference.
#' @param min_pole For a cohort reference, minimum number of samples in each
#'   pole. Poles with fewer samples, or whose score separates the two Hardy
#'   groups with an internal AUC below 0.6, are not used: positions and
#'   state labels of that axis are `NA`, with a warning, and
#'   `metadata()$poles` records the counts, the AUC and `ok`. Default 8.
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
thanatome_score <- function(x, reference, axes = c("ventilator", "agony"),
                            method = c("signature", "regression"),
                            assay = "counts", log = FALSE, lib_size = NULL,
                            normalize = c("auto", "quantile", "none"),
                            drop_unexpressed = TRUE, tissue = NULL,
                            min_pole = 8, state_threshold = 0.5,
                            flag_threshold = 4) {
    if (missing(reference) || is.null(reference)) {
        stop("'reference' is required: a tissue name from ",
            "names(thanatome_references) or a cohort of the same tissue")
    }
    axes <- match.arg(axes, several.ok = TRUE)
    method <- match.arg(method)
    normalize <- match.arg(normalize)
    m <- .as_logcpm(x, assay = assay, log = log, lib_size = lib_size)
    if (is.null(colnames(m))) stop("'x' needs sample names (column names)")
    normalize <- .resolve_normalize(normalize, m)
    ref <- .resolve_reference(reference, assay = assay, log = log,
                                normalize = normalize, min_pole = min_pole)
    sig <- .signatures()$signatures
    mz <- .prepare_for_z(m, ref, normalize, drop_unexpressed)
    z <- .z_against_reference(mz, ref, sig)
    cols <- c(axes, "ischemia")
    out <- .axis_scores(z, sig, cols)
    pos <- .positions(out, ref$poles)
    out$position_ventilator <- pos$ventilator
    out$position_agony <- pos$agony
    out$state <- .state(out$position_ventilator, out$position_agony,
                        state_threshold)
    fl <- .flags_from_logcpm(m, threshold = flag_threshold,
                            tissue = tissue %||% ref$tissue)
    out$flag_pancreas <- fl$pancreas
    out$flag_surfactant <- fl$surfactant
    if (method == "regression") {
        full <- ref$cohort && length(ref$mean) >= 10000 && nrow(mz) >= 10000
        rg <- .regression(z, sig[sig$axis %in% cols, ],
                        total_ss = if (full) .total_ss(mz, ref) else NULL)
        out <- cbind(out, rg)
    }
    out$n_genes <- .n_genes(z, sig, cols, nrow(out))
    meta <- .score_metadata(out, ref, method, axes, state_threshold,
                            .coldata(x))
    meta$normalize <- normalize
    meta$unexpressed_signature_genes <- attr(mz, "dropped_signature")
    S4Vectors::metadata(out) <- meta
    out
}

## ---- internal helpers -------------------------------------------------------

## Signature genes used per axis, as a matrix column of the result
.n_genes <- function(z, sig, cols, n) {
    k <- vapply(cols, function(a) {
        sum(sig$ensembl[sig$axis == a] %in% rownames(z))
    }, integer(1))
    matrix(k, nrow = n, ncol = length(cols), byrow = TRUE,
        dimnames = list(NULL, cols))
}
## "auto": quantile normalisation only for whole transcriptomes
.resolve_normalize <- function(normalize, m) {
    if (normalize != "auto") return(normalize)
    if (nrow(m) >= 10000) "quantile" else "none"
}
## Matrix used for the z-scores: optionally quantile-normalised to the
## reference distribution, without the genes 'x' does not express
.prepare_for_z <- function(m, ref, normalize, drop_unexpressed) {
    ## with fewer than 5 samples an absent gene and a lowly expressed one
    ## cannot be told apart: all genes are kept
    drop <- drop_unexpressed && ncol(m) >= 5
    low <- if (drop) attr(m, "absent") else character()
    if (normalize == "quantile") {
        if (nrow(m) < 10000) {
            stop("normalize = \"quantile\" needs the whole transcriptome ",
                "(at least 10 000 genes) in 'x'")
        }
        g <- intersect(rownames(m), ref$qn_genes)
        m <- .qn_apply(m[g, , drop = FALSE], ref$qn_target)
    }
    m <- m[!rownames(m) %in% low, , drop = FALSE]
    sig_ref <- intersect(.signatures()$signatures$ensembl,
                        names(ref$mean)[ref$expressed])
    attr(m, "dropped_signature") <- intersect(low, sig_ref)
    m
}
## Per-gene z against the reference, only for genes expressed in it
.z_against_reference <- function(m, ref, sig) {
    expressed <- names(ref$mean)[ref$expressed]
    g <- intersect(rownames(m), expressed)
    if (length(g) < 50) {
        stop("fewer than 50 reference genes found in 'x'; ",
            "check that row names are Ensembl ids")
    }
    sig_ref <- intersect(sig$ensembl, expressed)
    n_found <- length(intersect(g, sig_ref))
    n_low <- length(attr(m, "dropped_signature"))
    if (n_found < 0.9 * length(sig_ref)) {
        warning("only ", n_found, " of ", length(sig_ref),
                " signature genes of the reference are used (", n_low,
                " not expressed in 'x', ",
                length(sig_ref) - n_found - n_low,
                " absent); positions may be biased")
    }
    z <- (m[g, , drop = FALSE] - ref$mean[g]) / ref$sd[g]
    z[!is.finite(z)] <- 0
    z
}
## Target distribution: 2000 quantiles of the mean sorted profile
.qn_quantiles <- function(l, n_q = 2000) {
    o <- sort(rowMeans(apply(l, 2, sort)))
    stats::approx(seq(0, 1, length.out = length(o)), o,
                seq(0, 1, length.out = n_q))$y
}
## Quantile normalisation of each column to the target (ties averaged)
.qn_apply <- function(m, target_q) {
    n <- nrow(m)
    tg <- stats::approx(seq(0, 1, length.out = length(target_q)), target_q,
                        seq(0, 1, length.out = n))$y
    out <- apply(m, 2, function(v) {
        stats::approx(seq_len(n), tg, xout = rank(v, ties.method = "average"))$y
    })
    dimnames(out) <- dimnames(m)
    out
}
## Signature score per axis, one row per sample
.axis_scores <- function(z, sig, cols) {
    one <- function(a) .signature_score(z, sig[sig$axis == a, ])
    scores <- vapply(cols, one, numeric(ncol(z)))
    S4Vectors::DataFrame(matrix(scores, ncol = length(cols),
                                dimnames = list(colnames(z), cols)))
}
## Reference, method, poles and, if 'x' carries Hardy, the internal check
.score_metadata <- function(out, ref, method, axes, state_threshold, cd) {
    meta <- list(reference = ref$label, pipeline = ref$pipeline,
                method = method, poles = ref$poles,
                thresholds = .signatures()$thresholds,
                state_threshold = state_threshold)
    hardy <- .hardy_column(cd)
    if (is.null(hardy)) return(meta)
    hl <- hardy_layers(hardy)
    no_agony <- !"agony" %in% axes || all(is.na(out$position_agony))
    if (no_agony) hl$state[hl$state != "perfused"] <- "postmortem"
    auc_axis <- function(a) if (a %in% axes) .auc(hl[[a]], out[[a]]) else NA
    meta$auc <- c(ventilator = auc_axis("ventilator"),
                agony = auc_axis("agony"))
    ok <- !is.na(out$state) & !is.na(hl$state)
    meta$state_agreement <- if (any(ok)) {
        mean(out$state[ok] == hl$state[ok])
    } else {
        NA_real_
    }
    meta
}
## log2-CPM matrix (edgeR prior count 1) with unversioned Ensembl ids;
## attribute "absent": genes with zero counts in more than half of the
## samples (log input: log2-CPM below log2(0.1) in more than half)
.as_logcpm <- function(x, assay = "counts", log = FALSE, lib_size = NULL) {
    if (is(x, "EList")) {
        message("using the log2-CPM of the EList (voom uses a prior count ",
                "of 0.5, not 1)")
        x <- x$E
        log <- TRUE
    }
    if (is(x, "DGEList")) {
        ## library sizes only: the references were built without TMM
        if (is.null(lib_size)) lib_size <- x$samples$lib.size
        x <- x$counts
    }
    m <- .gene_matrix(x, assay)
    absent <- if (log) {
        rowMeans(m < log2(0.1)) > 0.5
    } else {
        rowMeans(m == 0) > 0.5
    }
    if (!log) m <- .cpm_from_counts(m, x, lib_size)
    attr(m, "absent") <- rownames(m)[absent]
    m
}
## Count or log matrix with unversioned, unique Ensembl ids and unique
## sample names
.gene_matrix <- function(x, assay) {
    m <- if (is(x, "SummarizedExperiment")) {
        SummarizedExperiment::assay(x, assay)
    } else {
        as.matrix(x)
    }
    if (is.null(rownames(m))) stop("'x' needs gene ids (Ensembl) as row names")
    if (anyDuplicated(colnames(m))) stop("'x' has duplicated sample names")
    rownames(m) <- sub("\\..*$", "", rownames(m))
    m[!duplicated(rownames(m)), , drop = FALSE]
}
.cpm_from_counts <- function(m, x, lib_size) {
    if (is.null(lib_size)) {
        cd <- .coldata(x)
        if (!is.null(cd) && "lib_size" %in% colnames(cd)) {
            lib_size <- cd$lib_size
        }
    }
    if (!is.null(lib_size) && length(lib_size) != ncol(m)) {
        stop("'lib_size' must have one value per sample of 'x'")
    }
    if (is.null(lib_size) && nrow(m) < 10000) {
        warning("'x' has only ", nrow(m), " genes and no library sizes: ",
                "CPM will be computed over that subset and will not be ",
                "comparable with the reference")
    }
    edgeR::cpm(m, log = TRUE, prior.count = 1, lib.size = lib_size)
}
.coldata <- function(x) {
    if (!is(x, "SummarizedExperiment")) return(NULL)
    SummarizedExperiment::colData(x)
}
.hardy_column <- function(cd) {
    if (is.null(cd)) return(NULL)
    k <- intersect(c("hardy", "DTHHRDY"), colnames(cd))
    if (length(k)) cd[[k[1]]] else NULL
}
## Reference: package tissue or user cohort (per-gene mean and SD; poles if
## Hardy is present)
.resolve_reference <- function(reference, assay, log, normalize = "none",
                                min_pole = 8) {
    if (is.character(reference)) return(.tissue_reference(reference, normalize))
    m <- .as_logcpm(reference, assay = assay, log = log)
    if (ncol(m) < 20) {
        warning("cohort reference with fewer than 20 samples: ",
                "centroid and poles will be noisy")
    }
    qn_genes <- qn_target <- NULL
    if (normalize == "quantile") {
        qn_genes <- rownames(m)[rowMeans(m) > 0]
        qn_target <- .qn_quantiles(m[qn_genes, , drop = FALSE])
        m <- .qn_apply(m[qn_genes, , drop = FALSE], qn_target)
    }
    mu <- rowMeans(m)
    s <- apply(m, 1, sd)
    s[!is.finite(s) | s == 0] <- NA
    expressed <- mu > 0 & !is.na(s)
    hardy <- .hardy_column(.coldata(reference))
    poles <- if (is.null(hardy)) {
        message("cohort reference without a 'hardy' column: ",
                "positions and state labels will be NA")
        NULL
    } else {
        .cohort_poles(m[expressed, , drop = FALSE], mu[expressed],
                        s[expressed], hardy, min_pole)
    }
    list(mean = mu, sd = s, expressed = expressed, poles = poles,
        tissue = NA_character_, pipeline = "cohort",
        label = sprintf("cohort of %d samples", ncol(m)), cohort = TRUE,
        qn_genes = qn_genes, qn_target = qn_target)
}
## Packaged GTEx reference; with normalize = "quantile" the statistics of
## the quantile-normalised GTEx samples
.tissue_reference <- function(tissue, normalize) {
    r <- death_reference(tissue)
    st <- r
    if (normalize == "quantile") {
        if (is.null(r$qn)) stop("this reference has no quantile statistics")
        st <- r$qn
    }
    list(mean = st$mean, sd = st$sd, expressed = st$expressed,
        poles = st$poles, tissue = r$tissue, pipeline = r$pipeline,
        label = paste("GTEx", r$tissue), cohort = FALSE,
        qn_genes = r$qn$genes, qn_target = r$qn$target)
}
## Poles of a cohort reference from its own Hardy labels
.cohort_poles <- function(m, mu, s, hardy, min_pole = 8) {
    sig <- .signatures()$signatures
    z <- (m - mu) / s
    z[!is.finite(z)] <- 0
    sv <- .signature_score(z, sig[sig$axis == "ventilator", ])
    sa <- .signature_score(z, sig[sig$axis == "agony", ])
    h <- .as_hardy(hardy)
    .check_poles(.poles_from_scores(sv, sa, h), sv, sa, h, min_pole)
}
## Poles of a cohort are kept only with at least 'min_pole' samples each and
## an internal AUC >= 0.6 in the expected direction (Hardy 0 above the rest
## on the ventilator axis, slow above fast deaths on the agony axis);
## otherwise they are set to NA, so positions and labels are NA
.check_poles <- function(poles, sv, sa, h, min_pole) {
    vent <- ifelse(is.na(h), NA, h == 0L)
    slow <- ifelse(h %in% c(1L, 2L, 3L, 4L), h >= 3L, NA)
    poles$auc <- c(.auc(vent, sv), .auc(slow, sa))
    big <- poles$n_pole0 >= min_pole & poles$n_pole1 >= min_pole
    poles$ok <- big & !is.na(poles$auc) & poles$auc >= 0.6
    if (any(!poles$ok)) {
        bad <- poles[!poles$ok, ]
        warning("cohort reference: poles not usable for ",
                paste0(bad$axis, " (n = ", bad$n_pole0, "/", bad$n_pole1,
                        ", AUC = ", round(bad$auc, 2), ")", collapse = "; "),
                "; positions and state labels of that axis are NA ",
                "(needs >= ", min_pole, " samples per pole and AUC >= 0.6)")
        poles$pole0[!poles$ok] <- NA_real_
        poles$pole1[!poles$ok] <- NA_real_
    }
    poles
}
.poles_from_scores <- function(sv, sa, h) {
    fast <- h %in% c(1L, 2L)
    slow <- h %in% c(3L, 4L)
    data.frame(axis = c("ventilator", "agony"),
                pole0 = c(median(sv[h == 0], na.rm = TRUE),
                        median(sa[fast], na.rm = TRUE)),
                pole1 = c(median(sv[h >= 1], na.rm = TRUE),
                        median(sa[slow], na.rm = TRUE)),
                n_pole0 = c(sum(h == 0, na.rm = TRUE), sum(fast)),
                n_pole1 = c(sum(h >= 1, na.rm = TRUE), sum(slow)),
                stringsAsFactors = FALSE)
}
## Signature score: mean z with the consensus sign (NA with fewer than 5 genes)
.signature_score <- function(z, s) {
    g <- intersect(s$ensembl, rownames(z))
    if (length(g) < 5) return(rep(NA_real_, ncol(z)))
    colMeans(z[g, , drop = FALSE] * sign(s$lfc_consensus[match(g, s$ensembl)]))
}
.positions <- function(out, poles) {
    one <- function(axis) {
        na <- rep(NA_real_, nrow(out))
        if (is.null(poles) || !axis %in% colnames(out)) return(na)
        p <- poles[poles$axis == axis, ]
        ok <- is.finite(p$pole0) && is.finite(p$pole1) && p$pole1 != p$pole0
        if (!ok) return(na)
        (out[[axis]] - p$pole0) / (p$pole1 - p$pole0)
    }
    list(ventilator = one("ventilator"), agony = one("agony"))
}
.state <- function(pv, pa, thr) {
    post <- ifelse(is.na(pa), "postmortem",
                    ifelse(pa < thr, "postmortem_rapid", "postmortem_agony"))
    ifelse(is.na(pv), NA_character_, ifelse(pv < thr, "perfused", post))
}
## Joint regression of the z profile on the consensus vectors (no intercept;
## ischemia enters as a column but is not reported)
.regression <- function(z, s, total_ss = NULL) {
    axes <- unique(s$axis)
    g <- intersect(unique(s$ensembl), rownames(z))
    B <- vapply(axes, function(a) {
        sa <- s[s$axis == a, ]
        b <- sa$lfc_consensus[match(g, sa$ensembl)]
        b[is.na(b)] <- 0
        b
    }, numeric(length(g)))
    B <- matrix(B, ncol = length(axes), dimnames = list(g, axes))
    zg <- z[g, , drop = FALSE]
    fit <- lm.fit(B, zg)
    coefs <- matrix(fit$coefficients, nrow = length(axes),
                    dimnames = list(axes, colnames(z)))
    ss_res <- colSums(as.matrix(fit$residuals)^2)
    ss_sig <- colSums(zg^2)
    reported <- setdiff(axes, "ischemia")
    out <- S4Vectors::DataFrame(t(coefs[reported, , drop = FALSE]))
    colnames(out) <- paste0("beta_", reported)
    out$r2_signature <- 1 - ss_res / ss_sig
    out$r2_total <- if (is.null(total_ss)) {
        NA_real_
    } else {
        (ss_sig - ss_res) / total_ss
    }
    out
}
## Sum of squares of z over all cohort genes (for r2_total)
.total_ss <- function(m, ref) {
    g <- intersect(rownames(m), names(ref$mean)[ref$expressed])
    z <- (m[g, , drop = FALSE] - ref$mean[g]) / ref$sd[g]
    z[!is.finite(z)] <- 0
    colSums(z^2)
}
## Mann-Whitney AUC (no dependencies)
.auc <- function(y, s) {
    ok <- !is.na(y) & !is.na(s)
    y <- as.logical(y[ok])
    s <- s[ok]
    n1 <- sum(y)
    n0 <- sum(!y)
    if (n1 < 2 || n0 < 2) return(NA_real_)
    r <- rank(s)
    (sum(r[y]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
