## Assisted by Claude Code (Anthropic); see README.md, section Provenance.

test_that("annotate_de maps HGNC symbols", {
    skip_if_not_installed("org.Hs.eg.db")
    sig <- death_signatures("ventilator")
    sym <- death_gene_sets("ventilator", id = "symbol")$ventilator_up[1:5]
    de <- data.frame(
        logFC = rep(1, 5), adj.P.Val = rep(0.01, 5), row.names = sym
    )
    expect_message(out <- annotate_de(de), "mapping symbols")
    expect_equal(out$gene, sym)
    expect_true(any(!is.na(out$axis_ventilator)))
})

test_that("annotate_de warns when no gene matches", {
    de <- data.frame(
        logFC = c(1, 2), adj.P.Val = c(0.1, 0.2),
        row.names = c("ENSG99999990", "ENSG99999991")
    )
    expect_warning(out <- annotate_de(de), "no gene matched")
    expect_true(all(is.na(out$axis_ventilator)))
})

.toy_regs <- function() {
    sig <- death_signatures("ventilator")
    up <- sig$ensembl[sig$direction == "up"][1:30]
    other <- paste0("ENSG", seq_len(30))
    list(
        TF_vent = list(tfmode = setNames(rep(1, 30), up)),
        TF_other = list(tfmode = setNames(rep(1, 30), other))
    )
}

test_that("annotate_mra accepts universe, data.frame nes, padj_threshold", {
    regs <- .toy_regs()
    nes <- data.frame(regulator = c("TF_vent", "TF_other"), nes = c(3, -1))
    universe <- unique(unlist(lapply(regs, function(r) names(r$tfmode))))
    out <- annotate_mra(nes, regs, axes = "ventilator", universe = universe)
    expect_equal(out$regulator, c("TF_vent", "TF_other"))
    expect_equal(out$nes, c(3, -1))
    expect_true(out$flag_control_state[1])
    strict <- annotate_mra(
        nes, regs, axes = "ventilator", universe = universe,
        padj_threshold = 1e-300
    )
    expect_false(any(strict$flag_control_state))
})
