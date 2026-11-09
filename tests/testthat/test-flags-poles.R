## Assisted by Claude Code (Anthropic); see README.md, section Provenance.
data(thanatome_example)
ctl <- thanatome_example[, thanatome_example$group == "control"]

## Controls with the genes of one flag raised to about 1000 CPM
with_flag <- function(x, f) {
    fl <- thanatome:::.signatures()$flags
    g <- intersect(fl$ensembl[fl$flag == f], rownames(x))
    assay(x)[g, ] <- matrix(as.integer(round(x$lib_size * 1e-3)),
                            nrow = length(g), ncol = ncol(x), byrow = TRUE)
    x
}

test_that("the surfactant flag is off for lung, also with a cohort", {
    x <- with_flag(ctl, "surfactant")
    expect_true(all(contamination_flags(x)$surfactant))
    expect_false(any(contamination_flags(x, tissue = "Lung")$surfactant))
    sc <- suppressWarnings(thanatome_score(x, reference = "lung"))
    expect_false(any(sc$flag_surfactant))
    coh <- suppressWarnings(thanatome_score(x, reference = x, tissue = "lung"))
    expect_false(any(coh$flag_surfactant))
    coh2 <- suppressWarnings(thanatome_score(x, reference = x))
    expect_true(all(coh2$flag_surfactant))
})

test_that("the pancreatic flag is off for pancreas", {
    x <- with_flag(ctl, "pancreas")
    expect_true(all(contamination_flags(x)$pancreas))
    expect_false(any(contamination_flags(x, tissue = "pancreas")$pancreas))
    coh <- suppressWarnings(thanatome_score(x, reference = x,
                                            tissue = "Pancreas"))
    expect_false(any(coh$flag_pancreas))
})

test_that("cohort poles need enough samples and separation", {
    ok <- suppressWarnings(thanatome_score(ctl, reference = ctl))
    expect_true(all(S4Vectors::metadata(ok)$poles$ok))
    ## agony pole with fewer than min_pole samples
    few <- ctl
    slow <- which(few$hardy >= 3)
    few$hardy[slow[-(1:3)]] <- 2L
    expect_warning(sc <- thanatome_score(ctl, reference = few),
                    "poles not usable for agony")
    expect_true(all(is.na(sc$position_agony)))
    expect_false(all(is.na(sc$position_ventilator)))
    p <- S4Vectors::metadata(sc)$poles
    expect_false(p$ok[p$axis == "agony"])
    ## no separation: Hardy labels unrelated to the samples
    shuf <- ctl
    set.seed(3)
    shuf$hardy <- sample(shuf$hardy)
    expect_warning(sc2 <- thanatome_score(ctl, reference = shuf),
                    "poles not usable")
    p2 <- S4Vectors::metadata(sc2)$poles
    expect_true(all(p2$auc < 0.6 & !p2$ok))
    expect_true(all(is.na(sc2$position_ventilator)))
    expect_true(all(is.na(sc2$state)))
})
