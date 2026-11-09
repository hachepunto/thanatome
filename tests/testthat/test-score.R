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

test_that("quantile normalisation needs the whole transcriptome", {
    expect_error(thanatome_score(ctl, reference = "ovary",
                                normalize = "quantile"), "transcriptome")
    ## whole-transcriptome stand-in: signature genes plus background genes
    set.seed(1)
    bg <- setdiff(death_reference("ovary")$qn$genes, rownames(ctl))[1:11000]
    m <- rbind(assay(ctl), matrix(rpois(length(bg) * ncol(ctl), 50),
                                ncol = ncol(ctl),
                                dimnames = list(bg, colnames(ctl))))
    sc <- thanatome_score(m, reference = "ovary", normalize = "quantile",
                        lib_size = ctl$lib_size)
    expect_equal(nrow(sc), ncol(ctl))
    expect_identical(S4Vectors::metadata(sc)$normalize, "quantile")
    expect_true(all(is.finite(sc$position_ventilator)))
})
