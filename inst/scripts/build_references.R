## Assisted by Claude Code (Anthropic); see README.md, section Provenance.
## Builds data/thanatome_references.rda from recount3 (GTEx v8, gene,
## GENCODE v26):
## per tissue, mean and SD of log2-CPM of the signature and flag genes, which
## genes are expressed (mean log2-CPM > 0), and the poles of each axis;
## element `qn`: the same statistics after quantile normalisation of the
## samples to the tissue distribution (for thanatome_score(normalize =
## "quantile")), with the gene universe and the 2000-quantile target.
source("inst/scripts/record_session.R")
## the quantile helpers of the package, so that build and use share them
pkgload::load_all(".", export_all = TRUE, quiet = TRUE)
suppressPackageStartupMessages({
    library(SummarizedExperiment)
    library(edgeR)
})
load("data/thanatome_signatures.rda")
cache <- Sys.getenv("THANATOME_RECOUNT3_DIR", unset = "")
tissues <- data.frame(
    tissue = c(
        "ovary", "breast", "thyroid", "lung", "colon_sigmoid",
        "colon_transverse", "liver", "prostate", "uterus", "skin_leg"
    ),
    project = c(
        "OVARY", "BREAST", "THYROID", "LUNG", "COLON", "COLON", "LIVER",
        "PROSTATE", "UTERUS", "SKIN"
    ),
    smtsd = c(
        "Ovary", "Breast - Mammary Tissue", "Thyroid", "Lung",
        "Colon - Sigmoid", "Colon - Transverse", "Liver", "Prostate",
        "Uterus", "Skin - Sun Exposed (Lower leg)"
    ),
    stringsAsFactors = FALSE
)
get_project <- function(project) {
    f <- file.path(cache, paste0(project, "_gtex.rds"))
    if (nzchar(cache) && file.exists(f)) {
        return(readRDS(f))
    }
    hp <- recount3::available_projects()
    pi <- hp[hp$project == project & hp$file_source == "gtex" &
        hp$project_type == "data_sources", ]
    rse <- recount3::create_rse(pi, type = "gene", annotation = "gencode_v26")
    assay(rse, "counts") <- recount3::compute_read_counts(rse)
    assay(rse, "raw_counts") <- NULL
    rse
}
sig <- thanatome_signatures$signatures
genes_ref <- unique(c(sig$ensembl, thanatome_signatures$flags$ensembl))
score_axis <- function(z, axis) {
    s <- sig[sig$axis == axis, ]
    g <- intersect(s$ensembl, rownames(z))
    sgn <- sign(s$lfc_consensus[match(g, s$ensembl)])
    colMeans(z[g, , drop = FALSE] * sgn)
}
projects <- lapply(unique(tissues$project), get_project)
names(projects) <- unique(tissues$project)
one <- function(i) {
    rse <- projects[[tissues$project[i]]]
    cd <- colData(rse)
    keep <- cd$gtex.smtsd == tissues$smtsd[i] &
        cd$gtex.smafrze == "RNASEQ" & !duplicated(cd$gtex.sampid)
    x <- assay(rse, "counts")[, keep, drop = FALSE]
    lib <- colSums(x)
    rownames(x) <- sub("\\..*$", "", rownames(x))
    x <- x[!duplicated(rownames(x)), ]
    l_all <- edgeR::cpm(x, lib.size = lib, log = TRUE, prior.count = 1)
    h <- cd$gtex.dthhrdy[keep]
    st <- ref_stats(l_all[intersect(genes_ref, rownames(l_all)), ], h)
    u <- rownames(l_all)[rowMeans(l_all) > 0]
    target <- .qn_quantiles(l_all[u, , drop = FALSE])
    lq <- .qn_apply(l_all[u, , drop = FALSE], target)
    qn <- c(list(genes = u, target = target),
            ref_stats(lq[intersect(genes_ref, u), ], h))
    list(
        tissue = tissues$tissue[i], gtex_name = tissues$smtsd[i],
        pipeline = "recount3 Monorail, GENCODE v26",
        n_samples = sum(keep),
        hardy = table(factor(h, levels = 0:4), useNA = "ifany"),
        mean = st$mean, sd = st$sd, expressed = st$expressed,
        poles = st$poles,
        ischemic_time_range = range(cd$gtex.smtsisch[keep], na.rm = TRUE),
        qn = qn
    )
}
## Mean, SD, expressed genes and poles of a log2-CPM matrix of the
## reference genes
ref_stats <- function(l, h) {
    mu <- rowMeans(l)
    s <- apply(l, 1, sd)
    expressed <- mu > 0
    z <- (l - mu) / s
    z[!is.finite(z)] <- 0
    z <- z[expressed, , drop = FALSE]
    sv <- score_axis(z, "ventilator")
    sa <- score_axis(z, "agony")
    poles <- data.frame(
        axis = c("ventilator", "agony"),
        pole0 = c(
            median(sv[h == 0], na.rm = TRUE),
            median(sa[h %in% c(1L, 2L)], na.rm = TRUE)
        ),
        pole1 = c(
            median(sv[h >= 1], na.rm = TRUE),
            median(sa[h %in% c(3L, 4L)], na.rm = TRUE)
        ),
        n_pole0 = c(sum(h == 0, na.rm = TRUE), sum(h %in% c(1L, 2L))),
        n_pole1 = c(sum(h >= 1, na.rm = TRUE), sum(h %in% c(3L, 4L))),
        stringsAsFactors = FALSE
    )
    list(mean = mu, sd = s, expressed = expressed, poles = poles)
}
thanatome_references <- lapply(seq_len(nrow(tissues)), one)
names(thanatome_references) <- tissues$tissue
print(t(vapply(thanatome_references, function(r) {
    c(
        n = r$n_samples, genes = length(r$mean),
        expressed = sum(r$expressed),
        vent0 = r$poles$pole0[1], vent1 = r$poles$pole1[1]
    )
}, numeric(5))))
save(thanatome_references,
    file = "data/thanatome_references.rda", compress = "xz"
)
record_session("build_references")
