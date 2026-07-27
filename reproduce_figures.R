#!/usr/bin/env Rscript
#
# Reproduce the single-cell figures from Erdem et al. 2026 (Diabetologia).
#
# Prerequisites
#   1. Place 060425_scRNA_v3.3.rds in data/ (see data/README.md)
#   2. Install the packages listed in DESCRIPTION
#
# Usage
#   Rscript reproduce_figures.R
#
# Writes fig_1j.svg, fig_1k.svg and esm_fig_2a.svg to figures/.

library(targets)

DATA_FILE <- "data/060425_scRNA_v3.3.rds"

REQUIRED <- c(
    "targets", "Seurat", "UCell", "BiocParallel", "dplyr", "tibble",
    "ggplot2", "lme4", "lmerTest", "glmmTMB", "broom.mixed"
)

# ---------------------------------------------------------------------------
# Pre-flight checks
# ---------------------------------------------------------------------------

missing <- REQUIRED[!vapply(REQUIRED, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
    stop(
        "Required packages not installed: ", paste(missing, collapse = ", "), "\n",
        "See DESCRIPTION for the full dependency list.",
        call. = FALSE
    )
}

if (!file.exists(DATA_FILE)) {
    stop(
        "Input data not found at ", DATA_FILE, "\n",
        "Download it first - see data/README.md.",
        call. = FALSE
    )
}

message("Input data: ", DATA_FILE,
        " (", round(file.size(DATA_FILE) / 1e9, 1), " GB)")

# UCell scoring of the full 448,935-cell object is the memory-limiting step. It
# forks one worker per core and each holds a copy of the object, so peak memory
# scales with core count - a reference run peaked at 252 GB on 4 cores. Warn
# against the requirement for the cores actually available.
cores <- {
    n <- Sys.getenv("SLURM_CPUS_PER_TASK", unset = "")
    if (n == "") n <- Sys.getenv("SLURM_CPUS_ON_NODE", unset = "")
    min(if (n != "") as.integer(n) else parallel::detectCores(), 8)
}
needed_gb <- 60 * cores

mem_gb <- tryCatch(
    as.numeric(system("free -g | awk '/Mem:/ {print $2}'", intern = TRUE)),
    warning = function(w) NA_real_,
    error = function(e) NA_real_
)
if (!is.na(mem_gb) && mem_gb < needed_gb) {
    warning(
        "Detected ", mem_gb, " GB RAM but UCell will fork ", cores,
        " workers, needing roughly ", needed_gb, " GB. Reduce the cores made ",
        "available to this job, or UCell scoring may be killed for running ",
        "out of memory.",
        call. = FALSE
    )
}

# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------

message("\nRunning pipeline (roughly an hour on first run)...\n")
started <- Sys.time()

tar_make()

elapsed <- round(as.numeric(difftime(Sys.time(), started, units = "mins")), 1)

# ---------------------------------------------------------------------------
# Report the key estimates so they can be checked against the paper
# ---------------------------------------------------------------------------

hurdle <- tar_read(hurdle_model)
slopes <- tar_read(random_slopes)

message("\nCompleted in ", elapsed, " min. Figures written to figures/\n")

message("Reported in the paper vs this run:")
message(sprintf(
    "  Part 1 odds ratio     paper 17.59 [5.53, 56]      this run %.2f [%.2f, %.2f]",
    hurdle$odds_ratio$estimate,
    hurdle$odds_ratio$ci_lower,
    hurdle$odds_ratio$ci_upper
))
message(sprintf(
    "  Part 2 beta*          paper 0.45 [0.28, 0.62]     this run %.2f [%.2f, %.2f]",
    hurdle$beta_std$estimate,
    hurdle$beta_std$ci_lower,
    hurdle$beta_std$ci_upper
))
message(sprintf(
    "  Cells / donors        paper 20,605 / 86           this run %s / %d",
    format(hurdle$n_cells_total, big.mark = ","),
    hurdle$n_donors
))
message(sprintf(
    "  Double-positive       paper 108 cells / 33 donors this run %d / %d",
    hurdle$n_cells_apc_pos,
    hurdle$n_donors_apc_pos
))
message(sprintf(
    "  Random slopes LRT     paper p < 0.001             this run p = %.3g (%d donors)",
    slopes$lrt$p_value,
    nrow(slopes$donor_slopes)
))
