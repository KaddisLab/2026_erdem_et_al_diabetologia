# PanKBase Filtering Functions
#
# Functions to filter the PanKBase Seurat object for the analysis cohort:
# - filter_pankbase_main():   control/vehicle treatments only, with a curated
#                             diabetes status (excludes experimentally
#                             perturbed cells)
# - filter_pankbase_ductal(): ductal cell types only
#
# Applied in sequence these yield the analysis cohort reported in the paper:
# 20,605 ductal cells from 86 donors.

#' Filter PanKBase for Main Inflammation/APC Analysis
#'
#' Filters to control/vehicle treatments only with curated diabetes status.
#' Excludes:
#' - Cytokine-treated cells (experimental perturbation)
#' - Drug-treated cells
#' - ER stress-induced cells
#' - Cells without curated diabetes status (NODM, AABP, or T1DM)
#'
#' Uses CONTROL_TREATMENTS constant from constants.R to define valid treatments.
#'
#' @param seurat_obj PanKBase Seurat object with mapped metadata
#' @return Filtered Seurat object (control/vehicle cells with valid diabetes status)
#' @export
filter_pankbase_main <- function(seurat_obj) {
    message("Filtering PanKBase for main analysis (control/vehicle treatments only)...")

    meta <- seurat_obj@meta.data
    n_start <- ncol(seurat_obj)

    # Filter criteria: treatment-based filtering
    keep <- (
        # Only control/vehicle treatments
        meta$treatments %in% CONTROL_TREATMENTS &
            # Must have curated diabetes status (NODM, AABP, or T1DM)
            meta$diabetes_status %in% c("NODM", "AABP", "T1DM") &
            # Must have valid donor identifier
            !is.na(meta$orig.ident) & meta$orig.ident != "" &
            # Must have valid sex and age for covariates
            !is.na(meta$sample_sex) &
            !is.na(meta$sample_age)
    )

    # Apply filter
    filtered_obj <- seurat_obj[, keep]

    # Summary
    n_end <- ncol(filtered_obj)
    message("\nMain analysis filter summary:")
    message("  Before: ", n_start, " cells")
    message("  After:  ", n_end, " cells")
    message("  Removed: ", n_start - n_end, " cells")
    message("\nRemaining diabetes status distribution:")
    print(table(filtered_obj@meta.data$diabetes_status))
    message("\nRemaining treatment distribution:")
    print(table(filtered_obj@meta.data$treatments))
    message("\nUnique donors: ", length(unique(filtered_obj@meta.data$orig.ident)))

    return(filtered_obj)
}

#' Filter to Ductal Cells
#'
#' Filters to ductal cell types only (Ductal + MUC5B+ Ductal)
#'
#' @param seurat_obj Seurat object with cell_type column
#' @return Filtered Seurat object (ductal cells only)
#' @export
filter_pankbase_ductal <- function(seurat_obj) {
    message("Filtering to ductal cells...")

    ductal_types <- c("Ductal", "MUC5B+ Ductal")
    meta <- seurat_obj@meta.data
    n_start <- ncol(seurat_obj)

    keep <- meta$cell_type %in% ductal_types
    filtered_obj <- seurat_obj[, keep]

    n_end <- ncol(filtered_obj)
    message("\nDuctal cell filter summary:")
    message("  Before: ", n_start, " cells")
    message("  After:  ", n_end, " cells (ductal + MUC5B+ ductal)")
    message("\nDuctal cells by diabetes status:")
    print(table(filtered_obj@meta.data$diabetes_status))
    message("\nUnique donors with ductal cells: ", length(unique(filtered_obj@meta.data$orig.ident)))

    return(filtered_obj)
}
