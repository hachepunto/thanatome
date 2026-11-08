## Asistido por Claude Code (Anthropic); ver README.md, sección Provenance.
## Construye data/thanatome_referencias.rda: por tejido, media y SD del log2-CPM
## de los genes de firma y de bandera en GTEx v8, y los polos de cada eje.
source("data-raw/record_session.R")
suppressPackageStartupMessages(library(edgeR))
repo <- Sys.getenv("GTEX_CONTROLS", unset = "~/gtex-proyect")
load("data/thanatome_firmas.rda")
tissues <- c("ovary", "breast", "thyroid", "lung", "colon_sigmoid", "colon_transverse", "liver", "prostate", "uterus", "skin_leg")
genes_ref <- unique(c(thanatome_firmas$signatures$ensembl, thanatome_firmas$flags$ensembl))
sig <- thanatome_firmas$signatures
score_axis <- function(z, axis) {
    s <- sig[sig$axis == axis, ]; g <- intersect(s$ensembl, rownames(z))
    colMeans(z[g, , drop = FALSE] * sign(s$lfc_consensus[match(g, s$ensembl)]))
}
una <- function(t) {
    x <- readRDS(file.path(repo, "data", "processed", t, "counts.rds"))
    md <- read.delim(file.path(repo, "data", "processed", t, "metadata.tsv"), stringsAsFactors = FALSE)
    stopifnot(identical(colnames(x), md$SAMPID))
    l <- edgeR::cpm(x, log = TRUE, prior.count = 1)
    rownames(l) <- sub("\\..*$", "", rownames(l)); l <- l[!duplicated(rownames(l)), ]
    expressed <- rownames(l)[rowMeans(l) > 0]
    g <- intersect(genes_ref, expressed); l <- l[g, , drop = FALSE]
    mu <- rowMeans(l); sd <- apply(l, 1, sd)
    z <- (l - mu) / sd; z[is.na(z)] <- 0
    sc <- sapply(c("ventilator", "ischemia", "agony"), function(a) score_axis(z, a))
    h <- md$DTHHRDY
    poles <- data.frame(axis = c("ventilator", "agony"),
        pole0 = c(median(sc[h == 0, "ventilator"]), median(sc[h %in% 1:2, "agony"])),
        pole1 = c(median(sc[h >= 1, "ventilator"]), median(sc[h %in% 3:4, "agony"])),
        n_pole0 = c(sum(h == 0), sum(h %in% 1:2)), n_pole1 = c(sum(h >= 1), sum(h %in% 3:4)), stringsAsFactors = FALSE)
    list(tissue = t, gtex_name = md$SMTSD[1], n = ncol(l), hardy = table(factor(h, levels = 0:4)),
         mean = mu, sd = sd, poles = poles, ischemia_range = range(md$SMTSISCH, na.rm = TRUE))
}
thanatome_referencias <- lapply(tissues, una); names(thanatome_referencias) <- tissues
print(t(sapply(thanatome_referencias, function(r) c(n = r$n, genes = length(r$mean), vent0 = r$poles$pole0[1], vent1 = r$poles$pole1[1]))))
save(thanatome_referencias, file = "data/thanatome_referencias.rda", compress = "xz")
record_session("build_referencias", repo)
