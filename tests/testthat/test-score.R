## Asistido por Claude Code (Anthropic); ver README.md, sección Provenance.
data(thanatome_example)

test_that("thanatome_score needs a reference and returns one row per sample", {
    expect_error(thanatome_score(thanatome_example), "reference")
    sc <- thanatome_score(thanatome_example, reference = "ovary")
    expect_s4_class(sc, "DataFrame")
    expect_equal(nrow(sc), ncol(thanatome_example))
    expect_true(all(c("ventilator", "ischemia", "agony", "position_ventilator", "position_agony", "state") %in% colnames(sc)))
    expect_true(all(sc$state %in% c("perfused", "postmortem_rapid", "postmortem_agony")))
})

test_that("the ventilator score separates Hardy 0 from the rest in the example", {
    sc <- thanatome_score(thanatome_example, reference = "ovary")
    auc <- S4Vectors::metadata(sc)$auc
    expect_gt(auc[["ventilator"]], 0.8)
    expect_gt(auc[["agony"]], 0.6)
    ## los ventilados quedan cerca del polo 0 y los demás cerca del 1
    h <- thanatome_example$hardy
    expect_lt(median(sc$position_ventilator[h == 0]), 0.5)
    expect_gt(median(sc$position_ventilator[h >= 1]), 0.5)
})

test_that("regression adds beta and r2 columns but keeps ischemia marginal", {
    sc <- thanatome_score(thanatome_example, reference = "ovary", method = "regression")
    expect_true(all(c("beta_ventilator", "beta_agony", "r2_signature", "r2_total") %in% colnames(sc)))
    expect_false("beta_ischemia" %in% colnames(sc))
    expect_true(all(sc$r2_signature >= 0 & sc$r2_signature <= 1))
    expect_true(all(is.na(sc$r2_total)))
    expect_gt(cor(sc$beta_ventilator, sc$ventilator), 0.8)
})

test_that("a cohort reference works, with poles when hardy is present", {
    sc <- thanatome_score(thanatome_example, reference = thanatome_example, method = "regression")
    expect_false(any(is.na(sc$position_ventilator)))
    expect_true(all(sc$r2_total > 0 & sc$r2_total < 1))
    m <- SummarizedExperiment::assay(thanatome_example)
    expect_message(sc2 <- thanatome_score(m, reference = m), "hardy")
    expect_true(all(is.na(sc2$state)))
})

test_that("matrix input and already-log input give the same scores", {
    m <- SummarizedExperiment::assay(thanatome_example)
    a <- thanatome_score(m, reference = "ovary")
    b <- thanatome_score(edgeR::cpm(m, log = TRUE, prior.count = 1), reference = "ovary", log = TRUE)
    expect_equal(a$ventilator, b$ventilator)
})
