## Assisted by Claude Code (Anthropic); see README.md, section Provenance.
data(thanatome_example)
ctl <- thanatome_example[, thanatome_example$group == "control"]
sc <- thanatome_score(ctl, reference = "ovary")

test_that("select_samples keeps everything without criteria and filters with them", {
    all <- select_samples(sc, flags = FALSE)
    expect_equal(length(all$keep), nrow(sc))
    sel <- select_samples(sc, ventilator = c(-Inf, 0.5), flags = FALSE)
    expect_true(all(sc[sel$keep, "position_ventilator"] <= 0.5))
    expect_equal(length(sel$keep) + nrow(sel$excluded), nrow(sc))
    expect_true(all(sel$excluded$reason == "ventilator"))
    pre <- select_samples(sc, preset = "perfused")
    expect_true(all(sc[pre$keep, "state"] == "perfused"))
    expect_equal(unname(pre$summary["kept"]), length(pre$keep))
})

test_that("plots are ggplot objects", {
    expect_s3_class(plot_thanatome(sc), "ggplot")
    expect_s3_class(plot_thanatome(sc, hardy = ctl$hardy, label = TRUE), "ggplot")
    expect_s3_class(plot_thanatome(sc, type = "ranked"), "ggplot")
})

test_that("death_covariates returns a centred matrix", {
    cv <- death_covariates(sc, quiet = TRUE)
    expect_equal(dim(cv), c(nrow(sc), 2))
    expect_warning(death_covariates(sc, axes = c("ventilator", "nope"), quiet = TRUE), "nope")
    expect_equal(unname(colMeans(cv)), c(0, 0), tolerance = 1e-12)
    expect_message(death_covariates(sc), "within")
})

test_that("annotate_de annotates signature genes and flags sign dependence", {
    sig <- death_signatures("ventilator")
    de <- data.frame(logFC = c(0.1, 3), adj.P.Val = c(0.01, 0.01), row.names = c(sig$ensembl[1], "ENSG00000000003"))
    a <- annotate_de(de)
    expect_equal(a$axis_ventilator[1], sig$direction[1])
    expect_true(is.na(a$axis_ventilator[2]))
    expect_true(a$sign_depends_on_control[1]); expect_false(a$sign_depends_on_control[2])
    expect_equal(hardy_layers(c(0, 4))$state, c("perfused", "postmortem_agony"))
})
