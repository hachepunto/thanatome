## Asistido por Claude Code (Anthropic); ver README.md, sección Provenance.
## Construye data/thanatome_firmas.rda a partir de los consensos de
## hachepunto/gtex-controls (results/firmas/*_consenso.tsv).
source("data-raw/record_session.R")
repo <- Sys.getenv("GTEX_CONTROLS", unset = "~/gtex-proyect")
min_tissues <- 5
## Umbral de |lfc mediano|: 0.5 para ventilador y agonía (validación); 0.25 para isquemia, cuyo
## contraste es por desviación estándar de tiempo isquémico y con 0.5 deja 21 genes casi no expresados
min_lfc <- c(ventilator = 0.5, ischemia = 0.25, agony = 0.5)
axes <- c(ventilator = "ventilador", ischemia = "isquemia", agony = "agonia")
leer <- function(axis) {
    d <- read.delim(file.path(repo, "results", "firmas", paste0(axes[[axis]], "_consenso.tsv")), stringsAsFactors = FALSE)
    d <- d[d$n_deg_concordantes >= min_tissues & abs(d$lfc_mediano) > min_lfc[[axis]] & !grepl("_PAR_Y$", d$ensembl), ]
    data.frame(ensembl = sub("\\..*$", "", d$ensembl), symbol = d$Symbol, axis = axis,
               lfc_consensus = d$lfc_mediano, z_stouffer = d$z_stouffer, n_tissues = d$n_deg_concordantes,
               direction = ifelse(d$lfc_mediano > 0, "up", "down"), stringsAsFactors = FALSE)
}
firmas <- do.call(rbind, lapply(names(axes), leer))
rownames(firmas) <- NULL
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
thanatome_firmas <- list(signatures = firmas, flags = flags,
                         thresholds = list(min_tissues = min_tissues, min_lfc = min_lfc),
                         source = "GTEx v8, ten tissues; hachepunto/gtex-controls results/firmas")
print(table(firmas$axis, firmas$direction))
save(thanatome_firmas, file = "data/thanatome_firmas.rda", compress = "xz")
record_session("build_firmas", repo)
