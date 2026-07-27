# PanKBase Data Loading and Metadata Mapping
#
# Functions to load and standardize the PanKBase Seurat object
# (448K cells from HPAP + IIDP + Prodo sources)

#' Load PanKBase Seurat Object
#'
#' @param path Path to the PanKBase RDS file
#' @return Seurat object (unmodified)
#' @export
load_pankbase <- function(path) {
    message("Loading PanKBase Seurat object from: ", path)
    obj <- readRDS(path)
    message("Loaded ", ncol(obj), " cells x ", nrow(obj), " features")
    return(obj)
}

#' Map PanKBase Metadata to Standardized Columns
#'
#' Creates standardized metadata columns that match the existing HPAP-only
#' pipeline conventions. This allows PanKBase QMDs to use similar code.
#'
#' @param seurat_obj PanKBase Seurat object
#' @return Seurat object with mapped metadata columns
#' @export
map_pankbase_metadata <- function(seurat_obj) {
    message("Mapping PanKBase metadata to standardized columns...")

    meta <- seurat_obj@meta.data

    # --- Core ID mapping ---
    meta$orig.ident <- meta$aliases

    # --- Cell type (use existing Cell_Type) ---
    meta$cell_type <- meta$Cell_Type

    # --- Demographic/clinical metadata ---
    meta$sample_sex <- meta$sex
    meta$sample_age <- meta$`age_(years)`
    # BMI and HbA1c keep their original names but also create clean versions
    meta$hba1c <- meta$`hba1c_(percentage)`
    meta$c_peptide <- meta$`c_peptide_(ng/ml)`

    # --- Autoantibody columns (clean names) ---
    meta$gada_positive <- meta$aab_gada_positive
    meta$gada_value <- meta$`aab_gada_value_(unit/ml)`
    meta$ia2_positive <- meta$aab_ia2_positive
    meta$ia2_value <- meta$`aab_ia2_value_(unit/ml)`
    meta$znt8_positive <- meta$aab_znt8_positive
    meta$znt8_value <- meta$`aab_znt8_value_(unit/ml)`
    meta$iaa_positive <- meta$aab_iaa_positive
    meta$iaa_value <- meta$`aab_iaa_value_(unit/ml)`

    # --- Data source ---
    meta$data_source <- meta$source

    # --- Treatment (for Prodo samples) ---
    meta$treatment <- meta$treatments

    # --- Diabetes status mapping ---
    # Any autoantibody positive
    # Note: AAb columns are logical (TRUE/FALSE/NA)
    # IIDP and Prodo have NA (no AAb data available)
    meta$any_aab_positive <- (
        isTRUE(meta$gada_positive) |
            isTRUE(meta$ia2_positive) |
            isTRUE(meta$znt8_positive) |
            isTRUE(meta$iaa_positive)
    )
    # Vectorized version since isTRUE doesn't vectorize
    meta$any_aab_positive <- (
        meta$gada_positive %in% TRUE |
            meta$ia2_positive %in% TRUE |
            meta$znt8_positive %in% TRUE |
            meta$iaa_positive %in% TRUE
    )

    # High-risk NODM flag (potential T2D controls)
    # BMI >= 35 AND HbA1c >= 6.0
    meta$high_risk_nodm <- (
        !is.na(meta$bmi) & !is.na(meta$hba1c) &
            meta$bmi >= 35 & meta$hba1c >= 6.0
    )
    meta$high_risk_nodm[is.na(meta$high_risk_nodm)] <- FALSE

    # Map diabetes status
    meta$diabetes_status <- dplyr::case_when(
        meta$description_of_diabetes_status == "type 1 diabetes" ~ "T1DM",
        meta$description_of_diabetes_status == "type 2 diabetes" ~ "T2DM",
        meta$description_of_diabetes_status == "non-diabetic" & meta$any_aab_positive ~ "AABP",
        meta$description_of_diabetes_status == "non-diabetic" & !meta$any_aab_positive & !meta$high_risk_nodm ~ "NODM",
        meta$description_of_diabetes_status == "non-diabetic" & !meta$any_aab_positive & meta$high_risk_nodm ~ "NODM_excluded",
        meta$description_of_diabetes_status == "" & meta$any_aab_positive ~ "AABP",
        meta$description_of_diabetes_status == "" & !meta$any_aab_positive & !meta$high_risk_nodm ~ "NODM",
        meta$description_of_diabetes_status == "" & !meta$any_aab_positive & meta$high_risk_nodm ~ "NODM_excluded",
        TRUE ~ NA_character_  # Prodo treated samples and unknowns
    )

    # Set factor levels for consistent plot ordering: NODM, AABP, T1DM
    meta$diabetes_status <- factor(
        meta$diabetes_status,
        levels = c("NODM", "AABP", "T1DM", "T2DM", "NODM_excluded")
    )

    # Update Seurat metadata
    seurat_obj@meta.data <- meta

    # --- Summary statistics ---
    message("\n=== PanKBase Metadata Mapping Summary ===")
    message("Total cells: ", ncol(seurat_obj))
    message("\nDiabetes status distribution:")
    print(table(meta$diabetes_status, useNA = "ifany"))
    message("\nData source distribution:")
    print(table(meta$data_source))
    message("\nCell type distribution:")
    print(sort(table(meta$cell_type), decreasing = TRUE))

    return(seurat_obj)
}
