args <- commandArgs(trailingOnly=TRUE)
seed <- as.integer(args[[1]])
output_dir <- args[[2]]
dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)

script_arg <- grep("^--file=", commandArgs(trailingOnly=FALSE), value=TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[[1]])))
Sys.setenv(SIMULATION_DEFINITIONS_ONLY="1")
source(file.path(script_dir, "simulation_fig_comparison.R"))

write_matrix <- function(X, path) {
  connection <- file(path, "wb")
  on.exit(close(connection))
  writeBin(as.double(t(X)), connection, size=8, endian="little")
}

P1 <- get_P1(n1_train, seed)[, -1]
sources <- lapply(seq_along(deltas), function(i) {
  get_random_shifts(deltas[[i]], nk_vec[[i]], seed)[, -1]
})
reference <- rbind(P1, sources[[1]])
write_matrix(reference, file.path(output_dir, "reference.bin"))

manifest <- data.frame(
  candidate=seq_len(length(deltas) - 1),
  rows=nk_vec[-1],
  columns=ncol(reference)
)
for (i in seq_len(nrow(manifest))) {
  write_matrix(
    sources[[i + 1]],
    file.path(output_dir, sprintf("candidate_%02d.bin", i))
  )
}
write.csv(manifest, file.path(output_dir, "manifest.csv"), row.names=FALSE)
