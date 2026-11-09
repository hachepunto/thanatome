## Assisted by Claude Code (Anthropic); see README.md, section Provenance.

test_that("death_gene_sets splits each axis by direction", {
    gs <- death_gene_sets()
    expect_type(gs, "list")
    expect_named(gs, c("ventilator_up", "ventilator_down",
                        "agony_up", "agony_down"))
    sg <- death_signatures(c("ventilator", "agony"))
    up <- sg$ensembl[sg$axis == "ventilator" & sg$direction == "up"]
    expect_setequal(gs$ventilator_up, up)
    expect_length(death_gene_sets(split = FALSE), 2)
    sym <- death_gene_sets("ischemia", id = "symbol")
    expect_false(any(grepl("^ENSG", unlist(sym))))
})

test_that("death_gene_sets exports to GSEABase", {
    skip_if_not_installed("GSEABase")
    gsc <- death_gene_sets(format = "GeneSetCollection")
    expect_s4_class(gsc, "GeneSetCollection")
    expect_identical(names(gsc), names(death_gene_sets()))
})
