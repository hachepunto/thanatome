## Assisted by Claude Code (Anthropic); see README.md, section Provenance.
data(thanatome_example)

set.seed(1)
test_that("compare_controls runs criteria against random subsets + MRA", {
    toy_mra <- function(de) {
        s <- setNames(de$stat, de$gene)
        c(
            up = mean(s[s > 0], na.rm = TRUE),
            down = mean(s[s < 0], na.rm = TRUE),
            all = mean(s, na.rm = TRUE)
        )
    }
    cc <- compare_controls(thanatome_example,
        group = "group", control = "control",
        criteria = list(
            perfused = list(preset = "perfused"),
            no_agony = list(preset = "no_agony")
        ),
        reference = "ovary", n_random = 2, mra = toy_mra
    )
    s <- cc$summary
    expect_equal(s$criterion, c("all", "perfused", "no_agony"))
    # controls without technical flags
    expect_equal(s$n_controls[1], length(select_samples(cc$scores)$keep))
    expect_lte(s$n_controls[1], sum(thanatome_example$group == "control"))
    expect_true(all(s$n_controls[-1] < s$n_controls[1]))
    expect_equal(s$lost[1], 0)
    expect_equal(s$gained[1], 0)
    expect_true(all(!is.na(s$lost_random[-1])))
    expect_true(all(c("mra_spearman", "mra_spearman_random") %in%
        colnames(s)))
    expect_equal(s$mra_spearman[1], 1)
    expect_equal(names(cc$de), s$criterion)
    expect_error(
        compare_controls(thanatome_example,
            group = "group",
            control = "control"
        ),
        "reference"
    )
})

test_that("annotate_mra flags regulons made of signature genes", {
    sig <- death_signatures("ventilator")
    up <- sig$ensembl[sig$direction == "up"][seq_len(30)]
    down <- sig$ensembl[sig$direction == "down"][seq_len(30)]
    regs <- list(
        TF_up = list(tfmode = setNames(rep(1, 30), up)),
        TF_down = list(tfmode = setNames(rep(1, 30), down)),
        TF_other = setNames(rep(1, 40), paste0("ENSG0000000", seq_len(40)))
    )
    a <- annotate_mra(c(TF_up = 2, TF_down = -2, TF_other = 0.1), regs,
        axes = "ventilator"
    )
    expect_equal(a$frac_ventilator[a$regulator == "TF_up"], 1)
    expect_equal(a$direction_ventilator[a$regulator == "TF_up"], 1)
    expect_equal(a$direction_ventilator[a$regulator == "TF_down"], -1)
    expect_true(a$flag_control_state[a$regulator == "TF_up"])
    expect_false(a$flag_control_state[a$regulator == "TF_other"])
})

test_that("audit_report prints and returns lines", {
    ctl <- thanatome_example[, thanatome_example$group == "control"]
    sc <- thanatome_score(ctl, reference = "ovary")
    r <- audit_report(sc, selection = select_samples(sc, preset = "no_agony"))
    out <- capture.output(print(r))
    expect_true(any(grepl("Selection", out)))
    expect_s3_class(r, "thanatome_report")
})
