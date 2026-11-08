## Asistido por Claude Code (Anthropic); ver README.md, sección Provenance.
## Construye data/thanatome_example.rda: SummarizedExperiment con conteos crudos de
## los genes de firma y bandera para 74 controles de ovario de GTEx v8 (proyecto
## OVARY de recount3, freeze RNASEQ, una corrida por muestra) y 40 tumores primarios
## de TCGA-OV (proyecto OV de recount3), con el tamaño de biblioteca completo en
## colData (necesario para que los CPM sean comparables con la referencia).
source("inst/scripts/record_session.R")
suppressPackageStartupMessages(library(SummarizedExperiment))
load("data/thanatome_signatures.rda")
cache <- Sys.getenv("THANATOME_RECOUNT3_DIR", unset = "")
get_rse <- function(project, source) {
    f <- file.path(cache, paste0(project, "_", source, ".rds"))
    if (nzchar(cache) && file.exists(f)) return(readRDS(f))
    hp <- recount3::available_projects()
    pi <- subset(hp, project == !!project & file_source == source & project_type == "data_sources")
    rse <- recount3::create_rse(pi, type = "gene", annotation = "gencode_v26")
    assay(rse, "counts") <- recount3::compute_read_counts(rse); assay(rse, "raw_counts") <- NULL
    rse
}
set.seed(2026)
genes_keep <- unique(c(thanatome_signatures$signatures$ensembl, thanatome_signatures$flags$ensembl))
prepare <- function(rse, idx) {
    m <- assay(rse, "counts")[, idx, drop = FALSE]; lib <- colSums(m)
    rownames(m) <- sub("\\..*$", "", rownames(m)); m <- m[!duplicated(rownames(m)), ]
    list(counts = m[rownames(m) %in% genes_keep, ], lib_size = lib)
}
rg <- get_rse("OVARY", "gtex"); cg <- colData(rg)
ok <- cg$gtex.smtsd == "Ovary" & cg$gtex.smafrze == "RNASEQ" & !duplicated(cg$gtex.sampid)
h <- factor(cg$gtex.dthhrdy, levels = 0:4)
quota <- pmax(round(67 * prop.table(table(h[ok]))), 6)   # proporcional por Hardy, al menos 6 por clase
ig <- sort(unlist(lapply(levels(h), function(k) { i <- which(ok & h == k); sample(i, min(quota[k], length(i))) })))
g <- prepare(rg, ig)
rt <- get_rse("OV", "tcga"); ct <- colData(rt)
tum <- which(ct$tcga.cgc_sample_sample_type == "Primary Tumor" & !duplicated(substr(ct$tcga.tcga_barcode, 1, 12)))
it <- sort(sample(tum, 40)); tt <- prepare(rt, it)
common <- intersect(rownames(g$counts), rownames(tt$counts))
counts <- cbind(g$counts[common, ], tt$counts[common, ])
cd <- DataFrame(group = rep(c("control", "tumor"), c(length(ig), length(it))),
                lib_size = c(g$lib_size, tt$lib_size),
                donor = c(cg$gtex.subjid[ig], substr(ct$tcga.tcga_barcode[it], 1, 12)),
                hardy = c(cg$gtex.dthhrdy[ig], rep(NA_integer_, length(it))),
                ischemic_time = c(cg$gtex.smtsisch[ig], rep(NA_real_, length(it))),
                rin = c(cg$gtex.smrin[ig], rep(NA_real_, length(it))),
                row.names = c(cg$gtex.sampid[ig], ct$tcga.tcga_barcode[it]))
colnames(counts) <- rownames(cd)
symbols <- thanatome_signatures$signatures$symbol[match(common, thanatome_signatures$signatures$ensembl)]
symbols[is.na(symbols)] <- thanatome_signatures$flags$symbol[match(common[is.na(symbols)], thanatome_signatures$flags$ensembl)]
se <- SummarizedExperiment(assays = list(counts = counts), colData = cd, rowData = DataFrame(symbol = symbols, row.names = common))
metadata(se) <- list(tissue = "ovary", source = "recount3 (GTEx v8 Ovary, RNASEQ freeze; TCGA-OV primary tumours), gene read counts, GENCODE v26")
thanatome_example <- se
print(se); print(table(se$group, se$hardy, useNA = "ifany"))
save(thanatome_example, file = "data/thanatome_example.rda", compress = "xz")
record_session("build_example")
