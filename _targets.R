# Pipeline to reproduce the single-cell figures in:
#
#   Erdem N, Arribas-Layton D, Zook HN, O'Meally D, Mares J, Quijano JC,
#   Donohue C, Ortiz JA, Jou K, Vasavada RC, Montero E, Kaddis JS,
#   Reijonen H, Ku HT. Human pancreatic ductal cells from non-diabetic donors
#   function as non-professional antigen-presenting cells upon inflammatory
#   cytokine exposure. Diabetologia 2026;69(8):2323-2337.
#   https://doi.org/10.1007/s00125-026-06746-x
#
# Outputs (written to figures/):
#   fig_1j.svg       Fig. 1j     APC score by cytokine response status
#   fig_1k.svg       Fig. 1k     Cytokine-APC dose-response
#   esm_fig_2a.svg   ESM Fig. 2a Per-donor random slopes
#
# Run with: Rscript reproduce_figures.R

library(targets)

tar_option_set(
    packages = c(
        "Seurat",
        "UCell",
        "BiocParallel",
        "dplyr",
        "tibble",
        "ggplot2",
        "lme4",
        "lmerTest",
        "glmmTMB",
        "broom.mixed"
    )
)

# Load the analysis functions in R/
tar_source("R")

list(

    # ---------------------------------------------------------------------
    # Input data
    # ---------------------------------------------------------------------
    # PanKBase Seurat object (~448,000 cells, 139 donors), archived at
    # https://pubrepo.coh.org/publication/17 - see data/README.md
    tar_target(
        seurat_file,
        "data/060425_scRNA_v3.3.rds",
        format = "file"
    ),

    tar_target(
        pankbase_raw,
        load_pankbase(seurat_file)
    ),

    # Standardise metadata column names and derive diabetes_status
    tar_target(
        pankbase_mapped,
        map_pankbase_metadata(pankbase_raw)
    ),

    # ---------------------------------------------------------------------
    # Signature scoring and cohort selection
    # ---------------------------------------------------------------------
    # Score before filtering so scores are independent of cohort definition
    tar_target(
        pankbase_scored,
        calculate_ucell_scores(
            pankbase_mapped,
            signatures,
            c("ifng_chemokine_effector", "apc_function_full")
        )
    ),

    # Control/vehicle treatments with a curated diabetes status
    tar_target(
        pankbase_cohort,
        filter_pankbase_main(pankbase_scored)
    ),

    # Ductal cells only -> the analysis cohort: 20,605 cells, 86 donors
    tar_target(
        ductal_cells,
        filter_pankbase_ductal(pankbase_cohort)
    ),

    # ---------------------------------------------------------------------
    # Models
    # ---------------------------------------------------------------------
    # Two-part hurdle model.
    #   Part 1 (binary):   all 20,605 cells / 86 donors -> OR 17.59
    #   Part 2 (positive): double-positive cells, APC > 0 AND cytokine > 0,
    #                      108 cells / 33 donors        -> beta* 0.45
    tar_target(
        hurdle_model,
        fit_hurdle_model(
            cell_data = ductal_cells@meta.data,
            outcome = "apc_function_full_UCell",
            predictor = "ifng_chemokine_effector_UCell"
        )
    ),

    # Random slopes vs random intercepts, tested by likelihood ratio test.
    # Uses APC-positive cells (APC > 0; 13,718 cells / 84 donors). This is a
    # BROADER subset than the hurdle positive component above, which also
    # requires cytokine > 0 - the two are not interchangeable.
    tar_target(
        random_slopes,
        compare_random_effects(
            cell_data = ductal_cells@meta.data %>%
                dplyr::filter(apc_function_full_UCell > 0),
            outcome = "apc_function_full_UCell",
            predictor = "ifng_chemokine_effector_UCell",
            filter_nonzero = FALSE
        )
    ),

    # ---------------------------------------------------------------------
    # Figures
    # ---------------------------------------------------------------------
    # Fig. 1j - APC score in cytokine-negative vs cytokine-positive cells,
    # annotated with the hurdle Part 1 odds ratio. Cytokine-negative cells are
    # subsampled to 1% for display only.
    tar_target(
        fig_1j,
        {
            p <- plot_cytokine_boxplot(
                cell_data = ductal_cells@meta.data,
                outcome = "apc_function_full_UCell",
                predictor = "ifng_chemokine_effector_UCell",
                diabetes_col = "diabetes_status",
                hurdle_or = hurdle_model$odds_ratio
            )
            save_manuscript_figure(p, "fig_1j", width = 6, height = 5)
            p
        }
    ),

    # Fig. 1k - dose-response among double-positive cells, annotated with the
    # hurdle Part 2 standardised effect size (beta*).
    tar_target(
        fig_1k,
        {
            p <- plot_cell_level_scatter_jittered(
                cell_data = ductal_cells@meta.data,
                outcome = "apc_function_full_UCell",
                predictor = "ifng_chemokine_effector_UCell",
                diabetes_col = "diabetes_status",
                hurdle_result = hurdle_model,
                alpha = 0.5,
                point_size = 3
            )
            save_manuscript_figure(p, "fig_1k", width = 7, height = 5)
            p
        }
    ),

    # ESM Fig. 2a - per-donor slope estimates, one row per donor
    tar_target(
        esm_fig_2a,
        {
            apc_pos <- ductal_cells@meta.data %>%
                dplyr::filter(apc_function_full_UCell > 0)

            # Per-donor annotation columns shown beside the forest plot
            cytokine_status <- apc_pos %>%
                dplyr::group_by(orig.ident) %>%
                dplyr::summarise(
                    cytokine_pos = ifelse(
                        any(ifng_chemokine_effector_UCell > 0), "Yes", "No"
                    ),
                    .groups = "drop"
                ) %>%
                dplyr::rename(donor = orig.ident)

            donor_metadata <- apc_pos %>%
                dplyr::distinct(orig.ident, diabetes_status, sample_sex, sample_age) %>%
                dplyr::rename(donor = orig.ident) %>%
                dplyr::mutate(sample_age = round(sample_age, 0)) %>%
                dplyr::left_join(cytokine_status, by = "donor")

            p <- plot_random_slopes_manuscript(
                random_effects_result = random_slopes,
                donor_metadata = donor_metadata
            )
            save_manuscript_figure(p, "esm_fig_2a", width = 10, height = 12)
            p
        }
    )
)
