## Asistido por Claude Code (Anthropic); ver README.md, sección Provenance.
## Construye data/thanatome_example.rda: SummarizedExperiment con conteos crudos
## de los genes de firma y bandera para 67 controles de ovario de GTEx v8 y
## 40 tumores primarios de ovario de TCGA-OV, ambos de recount3 (mismo
## pipeline Monorail, GENCODE v26), con el tamaño de biblioteca completo en
## colData (necesario para que los CPM sean comparables con la referencia).
## Los conteos de GTEx y TCGA en recount3 son de acceso abierto.
source("data-raw/record_session.R")
suppressPackageStartupMessages(library(SummarizedExperiment))
repo <- Sys.getenv("GTEX_CONTROLS", unset = "~/gtex-proyect")
load("data/thanatome_firmas.rda")
set.seed(2026)
genes_keep <- unique(c(thanatome_firmas$signatures$ensembl, thanatome_firmas$flags$ensembl))
preparar <- function(rse, idx) {
    m <- assay(rse, "counts")[, idx, drop = FALSE]; lib <- colSums(m)
    rownames(m) <- sub("\\..*$", "", rownames(m)); m <- m[!duplicated(rownames(m)), ]
    list(counts = m[rownames(m) %in% genes_keep, ], lib_size = lib)
}
## Controles GTEx (freeze RNASEQ), proporcionales por Hardy con al menos 6 por clase
rg <- readRDS(file.path(repo, "data", "raw", "recount3", "ov_gtex.rds")); cg <- colData(rg)
ok <- cg$gtex.smafrze == "RNASEQ" & !duplicated(cg$gtex.sampid)   # algunas muestras tienen dos corridas en recount3
h <- factor(cg$gtex.dthhrdy, levels = 0:4)
cuota <- pmax(round(67 * prop.table(table(h[ok]))), 6)
ig <- sort(unlist(lapply(levels(h), function(k) { i <- which(ok & h == k); sample(i, min(cuota[k], length(i))) })))
g <- preparar(rg, ig)
## Tumores primarios TCGA-OV, una alícuota por paciente
rt <- readRDS(file.path(repo, "data", "raw", "recount3", "ov_tcga.rds")); ct <- colData(rt)
tum <- which(ct$tcga.cgc_sample_sample_type == "Primary Tumor" & !duplicated(substr(ct$tcga.tcga_barcode, 1, 12)))
it <- sort(sample(tum, 40)); t <- preparar(rt, it)
common <- intersect(rownames(g$counts), rownames(t$counts))
counts <- cbind(g$counts[common, ], t$counts[common, ])
cd <- DataFrame(group = rep(c("control", "tumor"), c(length(ig), length(it))),
                lib_size = c(g$lib_size, t$lib_size),
                donor = c(cg$gtex.subjid[ig], substr(ct$tcga.tcga_barcode[it], 1, 12)),
                hardy = c(cg$gtex.dthhrdy[ig], rep(NA_integer_, length(it))),
                ischemic_time = c(cg$gtex.smtsisch[ig], rep(NA_real_, length(it))),
                rin = c(cg$gtex.smrin[ig], rep(NA_real_, length(it))),
                row.names = c(cg$gtex.sampid[ig], ct$tcga.tcga_barcode[it]))
colnames(counts) <- rownames(cd)
genes <- read.delim(file.path(repo, "data", "processed", "ovary", "genes.tsv"), stringsAsFactors = FALSE)
genes$ensembl <- sub("\\..*$", "", genes$ensembl)
se <- SummarizedExperiment(assays = list(counts = counts), colData = cd,
    rowData = DataFrame(symbol = genes$Symbol[match(common, genes$ensembl)], row.names = common))
metadata(se) <- list(tissue = "ovary", source = "recount3 (GTEx v8 Ovary, RNASEQ freeze; TCGA-OV primary tumours), gene read counts")
thanatome_example <- se
print(se); print(table(se$group, se$hardy, useNA = "ifany"))
save(thanatome_example, file = "data/thanatome_example.rda", compress = "xz")
record_session("build_example", repo)
