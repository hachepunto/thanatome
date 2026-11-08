## Asistido por Claude Code (Anthropic); ver README.md, sección Provenance.
## Construye data/thanatome_signatures.rda a partir de las tablas de consenso del
## repositorio de análisis (results/firmas/*_consenso.tsv). Ver README.md de esta carpeta.
source("inst/scripts/record_session.R")
repo <- Sys.getenv("GTEX_CONTROLS", unset = "~/gtex-proyect")
min_tissues <- 5
## |lfc mediano|: 0.5 para ventilador y agonía (umbral de la validación); 0.25 para la componente de
## isquemia, cuyo contraste es por desviación estándar de tiempo isquémico (con 0.5 quedan 21 genes
## casi no expresados). Isquemia no es un tercer eje: 61 de sus 86 genes están en la firma de ventilador
## con signo opuesto; se reporta como componente del eje ventilador sensible al tiempo isquémico.
min_lfc <- c(ventilator = 0.5, ischemia = 0.25, agony = 0.5)
axes <- c(ventilator = "ventilador", ischemia = "isquemia", agony = "agonia")
read_axis <- function(axis) {
    d <- read.delim(file.path(repo, "results", "firmas", paste0(axes[[axis]], "_consenso.tsv")), stringsAsFactors = FALSE)
    d <- d[d$n_deg_concordantes >= min_tissues & abs(d$lfc_mediano) > min_lfc[[axis]] & !grepl("_PAR_Y$", d$ensembl), ]
    data.frame(ensembl = sub("\\..*$", "", d$ensembl), symbol = d$Symbol, axis = axis,
               lfc_consensus = d$lfc_mediano, z_stouffer = d$z_stouffer, n_tissues = d$n_deg_concordantes,
               direction = ifelse(d$lfc_mediano > 0, "up", "down"), stringsAsFactors = FALSE)
}
signatures <- do.call(rbind, lapply(names(axes), read_axis))
rownames(signatures) <- NULL
vent <- signatures$ensembl[signatures$axis == "ventilator"]
signatures$in_ventilator <- signatures$axis == "ischemia" & signatures$ensembl %in% vent
## Conjuntos de banderas técnicas (contaminación de tejido vecino; Nieuwenhuis et al. 2020)
genes <- read.delim(file.path(repo, "data", "processed", "ovary", "genes.tsv"), stringsAsFactors = FALSE)
genes$ensembl <- sub("\\..*$", "", genes$ensembl)
flag_sets <- list(pancreas = c("PNLIP", "CELA3A", "PRSS1", "CPA1", "CELA2A", "CTRB1"),
                  surfactant = c("SFTPC", "SFTPB", "SFTPA1", "SFTPA2"))
flags <- do.call(rbind, lapply(names(flag_sets), function(f) {
    s <- flag_sets[[f]]; e <- genes$ensembl[match(s, genes$Symbol)]
    data.frame(ensembl = e, symbol = s, flag = f, stringsAsFactors = FALSE)[!is.na(e), ]
}))
rownames(flags) <- NULL
thanatome_signatures <- list(signatures = signatures, flags = flags,
                             thresholds = list(min_tissues = min_tissues, min_lfc = min_lfc),
                             source = "GTEx v8 (portal gene read counts), ten tissues; DESeq2 per tissue and Stouffer consensus; hachepunto/thanatome-manuscript results/firmas")
print(table(signatures$axis, signatures$direction)); print(table(signatures$in_ventilator[signatures$axis == "ischemia"]))
save(thanatome_signatures, file = "data/thanatome_signatures.rda", compress = "xz")
record_session("build_signatures", repo)
