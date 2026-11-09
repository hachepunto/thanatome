## Assisted by Claude Code (Anthropic); see README.md, section Provenance.
data(thanatome_example)
ctl <- thanatome_example[, thanatome_example$group == "control"]

test_that("thanatome_score needs a reference and returns one row per sample", {
    expect_error(thanatome_score(ctl), "reference")
    sc <- thanatome_score(ctl, reference = "ovary")
    expect_s4_class(sc, "DataFrame")
    expect_equal(nrow(sc), ncol(ctl))
    expect_true(all(c(
        "ventilator", "ischemia", "agony",
        "position_ventilator", "position_agony", "state"
    ) %in%
        colnames(sc)))
    expect_true(all(sc$state %in%
        c("perfused", "postmortem_rapid", "postmortem_agony")))
})

test_that("the ventilator score separates Hardy 0 from the rest", {
    sc <- thanatome_score(ctl, reference = "ovary")
    auc <- S4Vectors::metadata(sc)$auc
    expect_gt(auc[["ventilator"]], 0.8)
    expect_gt(auc[["agony"]], 0.6)
    ## ventilated samples sit near pole 0 and the rest near pole 1
    h <- ctl$hardy
    expect_lt(median(sc$position_ventilator[h == 0]), 0.5)
    expect_gt(median(sc$position_ventilator[h >= 1]), 0.5)
})

test_that("regression adds beta and r2 columns but keeps ischemia marginal", {
    sc <- thanatome_score(ctl, reference = "ovary", method = "regression")
    expect_true(all(c(
        "beta_ventilator", "beta_agony", "r2_signature",
        "r2_total"
    ) %in% colnames(sc)))
    expect_false("beta_ischemia" %in% colnames(sc))
    expect_true("ischemia" %in% colnames(sc))
    expect_true(all(sc$r2_signature >= 0 & sc$r2_signature <= 1))
    expect_true(all(is.na(sc$r2_total)))
    expect_gt(cor(sc$beta_ventilator, sc$ventilator), 0.8)
})

test_that("a cohort reference works, with poles when hardy is present", {
    sc <- thanatome_score(ctl, reference = ctl, method = "regression")
    expect_false(any(is.na(sc$position_ventilator)))
    # the example cohort does not carry the whole transcriptome
    expect_true(all(is.na(sc$r2_total)))
    expect_true(all(sc$r2_signature > 0 & sc$r2_signature < 1))
    m <- SummarizedExperiment::assay(ctl)
    expect_message(
        sc2 <- suppressWarnings(
            thanatome_score(m, reference = m, lib_size = ctl$lib_size)
        ),
        "hardy"
    )
    expect_true(all(is.na(sc2$state)))
})

test_that("matrix input and already-log input give the same scores", {
    m <- SummarizedExperiment::assay(ctl)
    a <- thanatome_score(m, reference = "ovary", lib_size = ctl$lib_size)
    lcpm <- edgeR::cpm(m, lib.size = ctl$lib_size, log = TRUE, prior.count = 1)
    b <- thanatome_score(lcpm, reference = "ovary", log = TRUE)
    expect_equal(a$ventilator, b$ventilator)
})

test_that("add_thanatome_score appends the score to colData", {
    sc <- thanatome_score(ctl, reference = "ovary")
    ann <- add_thanatome_score(ctl, reference = "ovary")
    expect_s4_class(ann, "SummarizedExperiment")
    expect_identical(dim(ann), dim(ctl))
    expect_equal(ann$ventilator, sc$ventilator)
    expect_identical(ann$state, sc$state)
    expect_identical(S4Vectors::metadata(ann)$thanatome$auc,
                    S4Vectors::metadata(sc)$auc)
    pre <- add_thanatome_score(ctl, reference = "ovary", prefix = "tht_")
    expect_true("tht_state" %in% colnames(SummarizedExperiment::colData(pre)))
    expect_message(add_thanatome_score(ann, reference = "ovary"), "replacing")
    expect_error(add_thanatome_score(assay(ctl), reference = "ovary"),
                "SummarizedExperiment")
})

test_that("genes not expressed in x are left out and reported", {
    x <- ctl
    g <- death_signatures("ventilator")$ensembl[1:3]
    g <- intersect(g, rownames(x))
    assay(x)[g, ] <- 0L
    sc <- suppressWarnings(thanatome_score(x, reference = "ovary"))
    expect_true(all(g %in%
        S4Vectors::metadata(sc)$unexpressed_signature_genes))
    keep <- suppressWarnings(thanatome_score(x, reference = "ovary",
                                            drop_unexpressed = FALSE))
    expect_false(isTRUE(all.equal(sc$ventilator, keep$ventilator)))
})

## Whole-transcriptome stand-in: signature genes plus background genes
whole <- function(x, n_bg = Inf, seed = 1) {
    set.seed(seed)
    bg <- setdiff(death_reference("ovary")$qn$genes, rownames(x))
    bg <- bg[seq_len(min(n_bg, length(bg)))]
    rbind(assay(x), matrix(rpois(length(bg) * ncol(x), 50), ncol = ncol(x),
        dimnames = list(bg, colnames(x))))
}

test_that("quantile normalisation needs the whole transcriptome", {
    expect_error(thanatome_score(ctl, reference = "ovary",
                                normalize = "quantile"), "10 000")
    m <- whole(ctl)
    sc <- thanatome_score(m, reference = "ovary", normalize = "quantile",
                        lib_size = ctl$lib_size)
    expect_equal(nrow(sc), ncol(ctl))
    expect_identical(S4Vectors::metadata(sc)$normalize, "quantile")
    expect_true(all(is.finite(sc$position_ventilator)))
})

test_that(".qn_apply maps every sample onto the target distribution", {
    set.seed(2)
    m <- matrix(rnorm(3000), ncol = 3,
                dimnames = list(paste0("g", 1:1000), c("a", "b", "c")))
    tq <- thanatome:::.qn_quantiles(matrix(rexp(5000), ncol = 5))
    q <- thanatome:::.qn_apply(m, tq)
    tg <- stats::approx(seq(0, 1, length.out = length(tq)), tq,
                        seq(0, 1, length.out = 1000))$y
    expect_equal(unname(sort(q[, "a"])), tg)
    expect_identical(order(q[, "b"]), order(m[, "b"]))
    expect_identical(dimnames(q), dimnames(m))
})

test_that("normalize = 'auto' quantile-normalises whole transcriptomes only", {
    sub <- thanatome_score(ctl, reference = "ovary")
    expect_identical(S4Vectors::metadata(sub)$normalize, "none")
    m <- whole(ctl)
    expect_message(auto <- thanatome_score(m, reference = "ovary",
                                        lib_size = ctl$lib_size), "auto")
    qn <- thanatome_score(m, reference = "ovary", lib_size = ctl$lib_size,
                        normalize = "quantile")
    expect_identical(S4Vectors::metadata(auto)$normalize, "quantile")
    expect_equal(auto$ventilator, qn$ventilator)
    none <- thanatome_score(m, reference = "ovary", lib_size = ctl$lib_size,
                        normalize = "none")
    expect_identical(S4Vectors::metadata(none)$normalize, "none")
    expect_false(isTRUE(all.equal(none$ventilator, auto$ventilator)))
})

test_that("normalisation counts the genes shared with the reference", {
    m <- whole(ctl)
    ## filtered matrix: enough genes, too little of the reference
    f <- m[c(rownames(ctl), setdiff(rownames(m), rownames(ctl))[1:12000]), ]
    expect_message(sc <- thanatome_score(f, reference = "ovary",
                                        lib_size = ctl$lib_size), "filtered")
    expect_identical(S4Vectors::metadata(sc)$normalize, "none")
    expect_warning(thanatome_score(f, reference = "ovary",
                                lib_size = ctl$lib_size,
                                normalize = "quantile"), "filtered")
    ## many rows, few Ensembl ids
    junk <- m
    rownames(junk)[-seq_len(nrow(ctl))] <- paste0("X", seq_len(nrow(m) -
                                                            nrow(ctl)))
    sc <- thanatome_score(junk, reference = "ovary", lib_size = ctl$lib_size)
    expect_identical(S4Vectors::metadata(sc)$normalize, "none")
    ## cohort reference with the signature genes only
    expect_error(thanatome_score(m, reference = ctl, normalize = "quantile",
                                lib_size = ctl$lib_size), "10 000")
    sc <- suppressWarnings(thanatome_score(m, reference = ctl,
                                        lib_size = ctl$lib_size))
    expect_identical(S4Vectors::metadata(sc)$normalize, "none")
})

test_that("a cohort reference works with quantile normalisation", {
    se <- SummarizedExperiment::SummarizedExperiment(
        list(counts = whole(ctl)), colData = SummarizedExperiment::colData(ctl))
    se$lib_size <- NULL
    sc <- thanatome_score(se, reference = se, normalize = "quantile")
    expect_true(all(S4Vectors::metadata(sc)$poles$ok))
    expect_gt(S4Vectors::metadata(sc)$auc[["ventilator"]], 0.8)
    rg <- thanatome_score(se, reference = se, method = "regression")
    expect_true(all(is.finite(rg$r2_total)))
    expect_true(all(rg$r2_total < 0.5))
    rn <- thanatome_score(se, reference = se, method = "regression",
                        normalize = "none")
    expect_true(all(is.finite(rn$r2_total)))
})

test_that("input containers and sample names are checked", {
    dup <- ctl
    colnames(dup)[2] <- colnames(dup)[1]
    expect_error(thanatome_score(dup, reference = "ovary"), "duplicated")
    y <- edgeR::DGEList(assay(ctl), lib.size = ctl$lib_size)
    y$samples$norm.factors <- seq(0.8, 1.2, length.out = ncol(ctl))
    sc <- thanatome_score(ctl, reference = "ovary")
    sy <- thanatome_score(y, reference = "ovary")
    expect_equal(unname(sy$ventilator), unname(sc$ventilator))
    el <- new("EList", list(E = edgeR::cpm(assay(ctl), log = TRUE,
        prior.count = 1, lib.size = ctl$lib_size)))
    expect_message(se <- thanatome_score(el, reference = "ovary"), "EList")
    expect_equal(unname(se$ventilator), unname(sc$ventilator))
})

test_that("state agreement is NA without Hardy labels", {
    x <- ctl
    x$hardy <- NULL
    sc <- thanatome_score(x, reference = "ovary")
    expect_null(S4Vectors::metadata(sc)$state_agreement)
})

test_that("audit_report writes to a file", {
    sc <- thanatome_score(ctl, reference = "ovary")
    f <- tempfile(fileext = ".txt")
    audit_report(sc, file = f)
    expect_true(file.exists(f))
    expect_gt(length(readLines(f)), 3)
})
