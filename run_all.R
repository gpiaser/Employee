## =====================================================================
## Employee Ownership and Entrenchment Strategies
## Runs the three scripts and writes their outputs to results/.
## Usage, from the root of the repository:  Rscript run_all.R
## Only base R is required. Script S1 takes a few minutes.
## =====================================================================

root <- getwd()
if (!dir.exists(file.path(root, "R"))) {
  stop("Run this file from the root of the repository.")
}
dir.create(file.path(root, "results"), showWarnings = FALSE)

scripts <- c("A4_mixed_equilibrium_example.R",   # Script S2
             "A5_overgranting_range.R",          # Script S3
             "S6_continuous_effort.R")           # Script S1

old <- setwd(file.path(root, "results"))
on.exit(setwd(old))
for (s in scripts) {
  cat("\n==================== ", s, " ====================\n", sep = "")
  t0 <- Sys.time()
  source(file.path(root, "R", s), echo = FALSE, local = new.env())
  cat(sprintf("Done in %.1f seconds\n",
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}
cat("\nOutputs written to", file.path(root, "results"), "\n")
print(sessionInfo())
