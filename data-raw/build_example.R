## Asistido por Claude Code (Anthropic); ver README.md, sección Provenance.
## Construye data/thanatome_example.rda: SummarizedExperiment con conteos crudos
## de los genes de firma y bandera para 60 muestras de ovario de GTEx v8
## (estratificadas por Hardy). Los conteos por gen de GTEx v8 son de acceso abierto.
source("data-raw/record_session.R")
suppressPackageStartupMessages(library(SummarizedExperiment))
repo <- Sys.getenv("GTEX_CONTROLS", unset = "~/gtex-proyect")
load("data/thanatome_firmas.rda")
set.seed(2026)
x <- readRDS(file.path(repo, "data", "processed", "ovary", "counts.rds"))
md <- read.delim(file.path(repo, "data", "processed", "ovary", "metadata.tsv"), stringsAsFactors = FALSE)
genes <- read.delim(file.path(repo, "data", "processed", "ovary", "genes.tsv"), stringsAsFactors = FALSE)
rownames(x) <- sub("\\..*$", "", rownames(x)); x <- x[!duplicated(rownames(x)), ]
genes$ensembl <- sub("\\..*$", "", genes$ensembl)
keep <- rownames(x) %in% unique(c(thanatome_firmas$signatures$ensembl, thanatome_firmas$flags$ensembl))
## 60 muestras: proporcional por Hardy, con al menos 6 de cada clase presente
h <- factor(md$DTHHRDY, levels = 0:4)
cuota <- pmax(round(60 * prop.table(table(h))), 6); cuota[table(h) == 0] <- 0
idx <- unlist(lapply(levels(h), function(k) { i <- which(h == k); sample(i, min(cuota[k], length(i))) }))
idx <- sort(idx)
se <- SummarizedExperiment(assays = list(counts = as.matrix(x[keep, idx])),
    colData = DataFrame(lib_size = colSums(x[, idx]), donor = md$SUBJID[idx], hardy = md$DTHHRDY[idx], ischemic_time = md$SMTSISCH[idx],
                        rin = md$SMRIN[idx], sex = md$SEX[idx], age = md$AGE[idx], row.names = md$SAMPID[idx]),
    rowData = DataFrame(symbol = genes$Symbol[match(rownames(x)[keep], genes$ensembl)], row.names = rownames(x)[keep]))
metadata(se) <- list(tissue = "ovary", source = "GTEx v8 gene read counts (open access), 67 of 177 Ovary samples")
thanatome_example <- se
print(se); print(table(se$hardy))
save(thanatome_example, file = "data/thanatome_example.rda", compress = "xz")
record_session("build_example", repo)
