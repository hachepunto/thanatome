## Asistido por Claude Code (Anthropic); ver README.md, sección Provenance.
## Escribe inst/scripts/sessioninfo/<name>.txt con la fecha, el commit del
## repositorio de análisis (si aplica) y el sessionInfo() de la corrida.
record_session <- function(name, repo = "") {
    dir <- file.path("inst", "scripts", "sessioninfo")
    dir.create(dir, showWarnings = FALSE, recursive = TRUE)
    commit <- if (nzchar(repo)) tryCatch(system2("git", c("-C", repo, "rev-parse", "HEAD"), stdout = TRUE), error = function(e) NA_character_) else NA_character_
    con <- file(file.path(dir, paste0(name, ".txt")), "wt")
    on.exit(close(con))
    writeLines(c(paste("script:", name), paste("date:", format(Sys.time(), "%Y-%m-%d %H:%M %Z")),
                 paste("analysis repository commit:", commit), ""), con)
    writeLines(capture.output(print(sessionInfo())), con)
}
