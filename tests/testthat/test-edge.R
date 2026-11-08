## Asistido por Claude Code (Anthropic); ver README.md, sección Provenance.
data(thanatome_example)
ctl <- thanatome_example[, thanatome_example$group == "control"]

test_that("hardy as factor or character gives the same result as integer", {
    expect_equal(hardy_layers(factor(c("0", "4")))$state, c("perfused", "postmortem_agony"))
    expect_equal(hardy_layers(c("0", "1", NA))$hardy, c(0L, 1L, NA))
    expect_error(hardy_layers(c(0, 7)), "0-4")
    f <- ctl; f$hardy <- factor(f$hardy)
    expect_equal(S4Vectors::metadata(thanatome_score(f, reference = "ovary"))$auc, S4Vectors::metadata(thanatome_score(ctl, reference = "ovary"))$auc)
})

test_that("a single sample works with both methods", {
    one <- ctl[, 1]
    sc <- thanatome_score(one, reference = "ovary", method = "regression")
    expect_equal(nrow(sc), 1)
    expect_true(is.finite(sc$beta_ventilator))
})

test_that("inputs without names or with wrong lib_size fail clearly", {
    m <- SummarizedExperiment::assay(ctl); colnames(m) <- NULL
    expect_error(thanatome_score(m, reference = "ovary", lib_size = ctl$lib_size), "sample names")
    expect_error(thanatome_score(ctl, reference = "ovary", lib_size = 1:3), "one value per sample")
    expect_warning(thanatome_score(SummarizedExperiment::assay(ctl), reference = "ovary"), "library sizes")
})

test_that("lib_size of x is not applied to a cohort reference", {
    sub <- ctl[, 1:30]
    a <- thanatome_score(sub, reference = ctl, lib_size = sub$lib_size)
    b <- thanatome_score(sub, reference = ctl)
    expect_equal(a$ventilator, b$ventilator)
})

test_that("state is 'postmortem' when the agony axis is not requested", {
    sc <- thanatome_score(ctl, reference = "ovary", axes = "ventilator")
    expect_true(all(sc$state %in% c("perfused", "postmortem")))
    expect_true(all(is.na(sc$position_agony)))
})

test_that("contamination flags work on counts and report gene numbers", {
    fl <- contamination_flags(ctl)
    expect_equal(nrow(fl), ncol(ctl))
    expect_true(all(fl$pancreas_n > 0))
    expect_type(fl$pancreas, "logical")
})

test_that("as_de_table handles limma and edgeR style tables", {
    tt <- data.frame(logFC = 1, adj.P.Val = 0.1, row.names = "ENSG00000000003.1")
    expect_equal(as_de_table(tt)$gene, "ENSG00000000003")
    tags <- structure(list(table = data.frame(logFC = -1, FDR = 0.2, row.names = "ENSG1")), class = "TopTags")
    expect_equal(as_de_table(tags)$padj, 0.2)
    expect_error(as_de_table(data.frame(x = 1, row.names = "g")), "log fold change")
})

test_that("annotate_mra maps symbols and fails clearly without regulons", {
    sig <- death_signatures("ventilator")
    regs <- list(TF = setNames(rep(1, 30), sig$symbol[seq_len(30)]))
    expect_message(a <- annotate_mra(c(TF = 1), regs), "symbols")
    expect_equal(a$frac_ventilator, 1)
    expect_error(annotate_mra(c(TF = 1), list(TF = c("a", "b"))), "at least")
})

test_that("compare_controls rejects missing groups and respects flags pooling", {
    bad <- thanatome_example; bad$group[1] <- NA
    expect_error(compare_controls(bad, group = "group", control = "control", reference = "ovary"), "missing")
})
