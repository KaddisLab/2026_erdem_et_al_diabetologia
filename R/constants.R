# Gene signatures, palettes and plot theme
#
# Only the constants used by the published single-cell analysis are defined
# here (Fig. 1j, Fig. 1k and ESM Fig. 2a).

# ============================================================================
# Gene signatures
# ============================================================================

signatures <- list(
    # Cytokine response: IFN-gamma-inducible chemokines and effectors
    ifng_chemokine_effector = c(
        "CXCL9", "CXCL10", "CXCL11", "GBP4", "GBP5", "NOS2"
    ),

    # APC programme: HLA class II, co-stimulatory and antigen-processing genes
    apc_function_full = c(
        # HLA class II and its master regulator
        "HLA-DRA", "HLA-DRB1", "HLA-DRB5",
        "HLA-DQA1", "HLA-DQA2", "HLA-DQB1",
        "HLA-DPA1", "HLA-DPB1", "CIITA",
        # Co-stimulatory
        "CD40", "ICAM1",
        # Antigen processing and loading
        "CD74", "CTSS", "HLA-DMA", "HLA-DMB", "HLA-DOA", "TAP1", "TAP2"
    )
)

# ============================================================================
# Cohort definition
# ============================================================================

# Control/vehicle treatments retained for the analysis cohort. Cells given any
# other treatment (cytokine, drug, ER stress) are experimentally perturbed and
# are excluded, so the cohort reflects donor-intrinsic variation only.
CONTROL_TREATMENTS <- c(
    "no_treatment",
    "DMSO",
    "DMSO_Control_36Hr",
    "DMSO_Control_72Hr",
    "4Hr_DMSO",
    "24h_DMSO",
    "mock",
    "EtOH"
)

# ============================================================================
# Palette and theme
# ============================================================================

# Diabetes status colours, as used in Fig. 1j, 1k and ESM Fig. 2a
diabetes_palette <- c(
    "NODM" = "#3fa36b",  # green
    "AABP" = "#F4C542",  # yellow
    "T1DM" = "#9353b3"   # purple
)

# Arial is the journal's font. Machines without it draw with a metric-compatible
# substitute; see arial_font_alias() in manuscript_plots.R.
manuscript_theme <- ggplot2::theme_classic() +
    ggplot2::theme(
        text = ggplot2::element_text(family = "Arial"),
        plot.title = ggplot2::element_text(size = 14, face = "bold"),
        plot.subtitle = ggplot2::element_text(size = 12),
        plot.caption = ggplot2::element_text(size = 10, hjust = 0),
        axis.title = ggplot2::element_text(size = 12),
        axis.text = ggplot2::element_text(size = 10),
        axis.ticks = ggplot2::element_blank(),
        axis.line = ggplot2::element_line(linewidth = 0.5),
        legend.title = ggplot2::element_text(size = 11),
        legend.text = ggplot2::element_text(size = 10),
        legend.key.size = ggplot2::unit(5, "mm"),
        legend.position = "bottom",
        strip.text = ggplot2::element_text(size = 11, face = "bold"),
        strip.background = ggplot2::element_blank(),
        panel.spacing.x = ggplot2::unit(0.5, "lines")
    )
