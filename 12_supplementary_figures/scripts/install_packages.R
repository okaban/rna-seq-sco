#!/usr/bin/env Rscript
# Install required packages

packages_to_check <- c("clusterProfiler", "enrichplot")

for (pkg in packages_to_check) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing", pkg, "...\n")
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  } else {
    cat(pkg, "is already installed\n")
  }
}

cat("Package check completed\n")
