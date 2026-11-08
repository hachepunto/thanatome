## Asistido por Claude Code (Anthropic); ver README.md, sección Provenance.
## Escribe data-raw/sessioninfo/<name>.txt con el commit del repositorio de
## análisis y el sessionInfo() de la corrida que produjo los datos.
record_session <- function(name, repo = Sys.getenv("GTEX_CONTROLS")) {
    dir <- file.path("data-raw", "sessioninfo")
    dir.create(dir, showWarnings = FALSE, recursive = TRUE)
    commit <- tryCatch(system2("git", c("-C", repo, "rev-parse", "HEAD"), stdout = TRUE), error = function(e) NA_character_)
    con <- file(file.path(dir, paste0(name, ".txt")), "wt")
    on.exit(close(con))
    writeLines(c(paste("script:", name), paste("date:", format(Sys.time(), "%Y-%m-%d %H:%M %Z")),
                 paste("gtex-controls commit:", commit), ""), con)
    writeLines(capture.output(print(sessionInfo())), con)
}
