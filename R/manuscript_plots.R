# Manuscript Figure Functions
# Publication-ready figures for the ductal cell LMM manuscript
# All figures use "Cytokine Response Score" and "APC Score" labels
# Uses manuscript_theme from constants.R for consistent styling
#
# All plot functions return a list with:
#   $plot   - the ggplot2 object (no embedded caption or legend)
#   $caption - text string for use in Quarto document

# ============================================================================
# Figure 1a: Cell-level boxplot (Cytokine+ vs Cytokine-)
# ============================================================================

#' Plot cell-level cytokine status boxplot (Fig 1a)
#'
#' Boxplot comparing APC scores between Cytokine+ and Cytokine- cells.
#' Cytokine+ status defined as UCell score > 0 for the cytokine signature.
#'
#' @param cell_data Cell-level data frame with UCell scores
#' @param outcome APC score column name
#' @param predictor Cytokine response score column name
#' @param diabetes_col Column for diabetes status coloring
#' @param hurdle_or Optional hurdle model OR result (list with estimate, ci_lower, ci_upper, p_value)
#' @return List with $plot (ggplot2 object) and $caption (text string)
#' @export
plot_cytokine_boxplot <- function(cell_data,
                                   outcome = "apc_function_full_UCell",
                                   predictor = "ifng_chemokine_effector_UCell",
                                   diabetes_col = "diabetes_status",
                                   hurdle_or = NULL,
                                   seed = 42) {

    # The Cytokine-ve group is subsampled to 1% for display (see below), which
    # draws on R's RNG. Seed it here so the panel is reproducible when this
    # function is called directly, not only when run under targets (which seeds
    # each target from its name).
    #
    # NOTE: this seed was added when the code was prepared for publication; the
    # published panel was produced under the targets seed. The subsample is
    # cosmetic - it affects which Cytokine-ve points are drawn, not the boxplot
    # statistics or any model estimate - so the published and regenerated panels
    # are statistically identical and differ only in the grey points shown.
    set.seed(seed)

    # Create cytokine status variable
    df <- cell_data %>%
        dplyr::mutate(
            cytokine_status = ifelse(
                .data[[predictor]] > 0,
                "Cytokine+ve",
                "Cytokine-ve"
            ),
            cytokine_status = factor(cytokine_status, levels = c("Cytokine-ve", "Cytokine+ve"))
        )

    # Calculate statistics
    n_pos <- sum(df$cytokine_status == "Cytokine+ve")
    n_neg <- sum(df$cytokine_status == "Cytokine-ve")
    mean_pos <- mean(df[[outcome]][df$cytokine_status == "Cytokine+ve"], na.rm = TRUE)
    mean_neg <- mean(df[[outcome]][df$cytokine_status == "Cytokine-ve"], na.rm = TRUE)
    fold_change <- mean_pos / mean_neg

    # Calculate y position for significance bar (above max data point)
    y_max <- max(df[[outcome]], na.rm = TRUE)
    bar_y <- y_max * 1.05
    text_y <- y_max * 1.12

    # Determine annotation text based on hurdle OR or Wilcoxon
    if (!is.null(hurdle_or)) {
        # Use hurdle model OR with significance stars
        or_val <- hurdle_or$estimate
        or_p <- hurdle_or$p_value
        signif_stars <- dplyr::case_when(
            or_p < 0.001 ~ "***",
            or_p < 0.01 ~ "**",
            or_p < 0.05 ~ "*",
            TRUE ~ ""
        )
        annot_label <- sprintf("OR = %.2f%s", or_val, signif_stars)
    } else {
        # Fall back to Wilcoxon
        wilcox_p <- stats::wilcox.test(df[[outcome]] ~ df$cytokine_status)$p.value
        signif_stars <- dplyr::case_when(
            wilcox_p < 0.001 ~ "***",
            wilcox_p < 0.01 ~ "**",
            wilcox_p < 0.05 ~ "*",
            TRUE ~ ""
        )
        annot_label <- signif_stars
    }

    # Build plot - grey90 for Cytokine-ve, light red for Cytokine+ve
    # Draw boxplots separately to control fill without conflicting with point fill
    p <- ggplot2::ggplot(df, ggplot2::aes(x = cytokine_status, y = .data[[outcome]])) +
        # Cytokine-ve boxplot (grey)
        ggplot2::geom_boxplot(
            data = dplyr::filter(df, cytokine_status == "Cytokine-ve"),
            outlier.shape = NA,
            width = 0.5,
            fill = "grey90"
        ) +
        # Cytokine+ve boxplot (deeper orange - distinct from red and AABP yellow)
        ggplot2::geom_boxplot(
            data = dplyr::filter(df, cytokine_status == "Cytokine+ve"),
            outlier.shape = NA,
            width = 0.5,
            fill = "#FFAB66"
        ) +
        # Points for Cytokine-ve (sample to avoid overplotting)
        ggplot2::geom_jitter(
            data = dplyr::filter(df, cytokine_status == "Cytokine-ve") %>%
                dplyr::slice_sample(prop = 0.01, replace = FALSE),  # Reproducible 1% sample
            ggplot2::aes(fill = .data[[diabetes_col]]),
            position = ggplot2::position_jitter(width = 0.18, height = 0, seed = 42),
            alpha = 0.7,
            size = 3,
            shape = 21,
            stroke = 0.5,
            colour = "white"
        ) +
        # Points for Cytokine+ve
        ggplot2::geom_jitter(
            data = dplyr::filter(df, cytokine_status == "Cytokine+ve"),
            ggplot2::aes(fill = .data[[diabetes_col]]),
            position = ggplot2::position_jitter(width = 0.18, height = 0, seed = 42),
            alpha = 0.7,
            size = 3,
            shape = 21,
            stroke = 0.5,
            colour = "white"
        ) +
        # Significance bar (no end ticks, like reference figure)
        ggplot2::annotate("segment", x = 1, xend = 2, y = bar_y, yend = bar_y, linewidth = 0.5) +
        # family is set explicitly: annotate() does not inherit it from the theme
        ggplot2::annotate("text", x = 1.5, y = text_y, label = annot_label,
                          size = 6, family = "Arial") +
        ggplot2::scale_fill_manual(
            values = diabetes_palette,
            guide = "none"
        ) +
        ggplot2::scale_x_discrete(
            labels = c(
                "Cytokine-ve" = "=0",
                "Cytokine+ve" = ">0"
            )
        ) +
        ggplot2::labs(
            x = "Single-cell cytokine score",
            y = "Single-cell APC score"
        ) +
        ggplot2::expand_limits(y = 0) +
        manuscript_theme

    # Build caption text
    if (!is.null(hurdle_or)) {
        caption <- sprintf(
            "Box plot comparing APC scores between Cytokine<sup>−ve</sup> and Cytokine<sup>+ve</sup> (UCell IFN-γ > 0) ductal cells. Cytokine<sup>+ve</sup> cells (n = %d) exhibited significantly higher APC scores than Cytokine<sup>−ve</sup> cells (n = %s; mean %.3f vs %.3f; %.1f-fold increase). Mixed-effects logistic regression (binary component) shows Cytokine<sup>+ve</sup> status significantly increases odds of APC activation (OR = %.2f [95%% CI: %.2f–%.2f], p %s). The 6-gene IFN-γ chemokine/effector signature comprises *CXCL9*, *CXCL10*, *CXCL11*, *GBP4*, *GBP5*, *NOS2*; APC score is an 18-gene antigen presenting cell function signature. Box plots show median, interquartile range (box), and whiskers (1.5× IQR). Individual Cytokine<sup>+ve</sup> cells shown as jittered points coloured by diabetes status: green, NODM; yellow, AABP; purple, T1DM. A random 1%% sample of Cytokine<sup>−ve</sup> cells is shown to avoid overplotting.",
            n_pos,
            format(n_neg, big.mark = ","),
            mean_pos,
            mean_neg,
            fold_change,
            hurdle_or$estimate,
            hurdle_or$ci_lower,
            hurdle_or$ci_upper,
            if (hurdle_or$p_value < 0.001) "< 0.001" else sprintf("= %.3f", hurdle_or$p_value)
        )
    } else {
        wilcox_p <- stats::wilcox.test(df[[outcome]] ~ df$cytokine_status)$p.value
        caption <- sprintf(
            "Box plot comparing APC scores between Cytokine<sup>−ve</sup> and Cytokine<sup>+ve</sup> (UCell IFN-γ > 0) ductal cells. Cytokine<sup>+ve</sup> cells (n = %d) exhibited significantly higher APC scores than Cytokine<sup>−ve</sup> cells (n = %s; mean %.3f vs %.3f; %.1f-fold increase; Wilcoxon rank-sum test, p %s). The 6-gene IFN-γ chemokine/effector signature comprises *CXCL9*, *CXCL10*, *CXCL11*, *GBP4*, *GBP5*, *NOS2*; APC score is an 18-gene antigen presenting cell function signature. Box plots show median, interquartile range (box), and whiskers (1.5× IQR). Individual Cytokine<sup>+ve</sup> cells shown as jittered points coloured by diabetes status: green, NODM; yellow, AABP; purple, T1DM. A random 1%% sample of Cytokine<sup>−ve</sup> cells is shown to avoid overplotting.",
            n_pos,
            format(n_neg, big.mark = ","),
            mean_pos,
            mean_neg,
            fold_change,
            if (wilcox_p < 0.001) "< 0.001" else sprintf("= %.3f", wilcox_p)
        )
    }

    list(plot = p, caption = caption)
}

# ============================================================================
# Figure 1b Alternative: Cell-level scatter with jittered cytokine=0 cells
# ============================================================================

#' Plot cell-level scatter with jittered cytokine=0 cells (Fig 1b variant)
#'
#' Scatter plot of cytokine vs APC scores for APC-positive cells, with
#' x-axis jittering applied only to cells with cytokine score = 0.
#' This spreads out the zero-valued cells for better visualization while
#' preserving the true positions of cytokine-positive cells.
#'
#' @param cell_data Cell-level data frame with UCell scores
#' @param outcome APC score column name
#' @param predictor Cytokine response score column name
#' @param diabetes_col Column for diabetes status coloring
#' @param hurdle_result Optional hurdle model result from fit_hurdle_model()
#' @param alpha Point transparency
#' @param point_size Point size
#' @param jitter_width Width of jitter band for cytokine=0 cells (default: 0.015)
#' @return List with $plot (ggplot2 object), $caption, sample sizes, and split counts
#' @export
plot_cell_level_scatter_jittered <- function(cell_data,
                                              outcome = "apc_function_full_UCell",
                                              predictor = "ifng_chemokine_effector_UCell",
                                              diabetes_col = "diabetes_status",
                                              hurdle_result = NULL,
                                              alpha = 0.3,
                                              point_size = 1,
                                              jitter_width = 0.015) {

    # Filter to double-positive cells (APC > 0 and cytokine > 0) to match hurdle Part 2 subset
    df <- cell_data %>%
        dplyr::filter(
            !is.na(.data[[outcome]]),
            !is.na(.data[[predictor]]),
            .data[[outcome]] > 0,
            .data[[predictor]] > 0
        )

    # All cells are cytokine > 0; use actual cytokine score as x
    df_plot <- df %>%
        dplyr::mutate(x_jittered = .data[[predictor]])

    n_cells <- nrow(df_plot)
    n_donors <- dplyr::n_distinct(df_plot$orig.ident)

    # Get stats from hurdle result
    if (!is.null(hurdle_result)) {
        beta_std <- hurdle_result$beta_std$estimate
        stat_label <- sprintf("\u03b2* = %.2f\np < 0.001", beta_std)
    } else {
        cor_result <- stats::cor.test(df[[outcome]], df[[predictor]], method = "spearman")
        stat_label <- sprintf("\u03c1 = %.2f\np < 0.001", cor_result$estimate)
    }

    # Build plot
    p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = x_jittered, y = .data[[outcome]])) +
        ggplot2::geom_point(
            ggplot2::aes(fill = .data[[diabetes_col]]),
            alpha = alpha,
            size = point_size,
            shape = 21,
            stroke = 0.3,
            colour = "white"
        ) +
        ggplot2::scale_fill_manual(values = diabetes_palette, guide = "none")

    # Add trend line from hurdle model
    if (!is.null(hurdle_result$part2_model)) {
        fe <- lme4::fixef(hurdle_result$part2_model)
        intercept <- fe["(Intercept)"]
        slope <- hurdle_result$beta_raw$estimate
        p <- p +
            ggplot2::geom_abline(
                intercept = intercept,
                slope = slope,
                colour = "black",
                linewidth = 1,
                linetype = "solid"
            )
    }

    # Add annotation
    x_max <- max(df_plot$x_jittered, na.rm = TRUE)
    y_range <- range(df_plot[[outcome]], na.rm = TRUE)
    annot_x <- x_max * 0.75
    annot_y <- y_range[1] + (y_range[2] - y_range[1]) * 0.12

    p <- p +
        ggplot2::annotate(
            "text",
            x = annot_x,
            y = annot_y,
            label = stat_label,
            hjust = 0.5,
            size = 4,
            # annotate() does not inherit the family from the theme
            family = "Arial"
        ) +
        ggplot2::labs(
            x = "Single-cell cytokine score",
            y = "Single-cell APC score"
        ) +
        manuscript_theme

    # Build caption text
    if (!is.null(hurdle_result)) {
        beta_std <- hurdle_result$beta_std$estimate
        beta_ci_lower <- hurdle_result$beta_std$ci_lower
        beta_ci_upper <- hurdle_result$beta_std$ci_upper
        caption <- sprintf(
            "Dose-response relationship between cytokine and APC scores in double-positive ductal cells. This analysis includes cells with both APC score > 0 and cytokine score > 0 (n = %s cells from %d donors), matching the double-positive subset used for the hurdle model Part 2 (positive component). The standardized effect size \u03b2* = %.2f [95%% CI: %.2f\u2013%.2f] indicates that a 1 SD increase in cytokine score is associated with a %.2f SD increase in APC score. Points coloured by diabetes status; black line shows LMM fixed effect.",
            format(n_cells, big.mark = ","),
            n_donors,
            beta_std,
            beta_ci_lower,
            beta_ci_upper,
            beta_std
        )
    } else {
        cor_result <- stats::cor.test(df[[outcome]], df[[predictor]], method = "spearman")
        caption <- sprintf(
            "Co-regulation of APC and cytokine pathways in double-positive ductal cells. This analysis includes cells with both APC score > 0 and cytokine score > 0 (n = %s cells from %d donors). Within these cells, cytokine score magnitude correlates with APC score magnitude (Spearman \u03c1 = %.3f, p < 0.001). Points coloured by diabetes status; black line shows linear regression.",
            format(n_cells, big.mark = ","),
            n_donors,
            cor_result$estimate
        )
    }

    list(
        plot = p,
        caption = caption,
        n_cells = n_cells,
        n_donors = n_donors
    )
}

# ============================================================================
# Supp Figure a: Random slopes forest plot (wrapper)
# ============================================================================

#' Plot random slopes for manuscript (Supp Fig a)
#'
#' Wrapper around plot_donor_slopes() with manuscript styling.
#'
#' @param random_effects_result Result from compare_random_effects()
#' @param donor_metadata Optional donor metadata tibble
#' @return List with $plot (ggplot2 object) and $caption (text string)
#' @export
plot_random_slopes_manuscript <- function(random_effects_result,
                                           donor_metadata = NULL) {

    # Use existing plot_donor_slopes but apply manuscript theme
    p <- plot_donor_slopes(
        random_effects_result,
        donor_metadata = donor_metadata,
        title = "Per-Donor Slopes: APC ~ Cytokine Response",
        xlab = "Effect of Cytokine Response on APC Score"
    )

    # Apply manuscript styling and remove legend
    p <- p +
        manuscript_theme +
        ggplot2::theme(legend.position = "none")

    # Get number of donors and cells from the random effects result
    n_donors <- nrow(random_effects_result$donor_slopes)
    n_cells <- random_effects_result$n_cells

    # Build caption text
    caption <- sprintf(
        "Per-donor random slopes from LMM with random intercepts and slopes among cytokine-positive cells. Each point represents a donor's individual slope for the effect of cytokine response on APC score. n = %d donors (%s cells). Points coloured by diabetes status: NODM (green), AABP (yellow), T1DM (purple). Error bars show 95%% CI.",
        n_donors,
        format(n_cells, big.mark = ",")
    )

    list(plot = p, caption = caption)
}

# ============================================================================
# Helper: Save figure in multiple formats
# ============================================================================

#' Resolve the font used to draw "Arial" in figure text
#'
#' The journal's font is Arial. Where it is unavailable, Liberation Sans is a
#' metric-compatible substitute: identical advance widths, so text occupies the
#' same space and the figure is laid out exactly as it would be with Arial.
#'
#' svglite records the name held inside the font file it draws with, which would
#' write "Liberation Sans" into the SVG and lose the request for Arial even for
#' readers who have it. Passing an explicit alias to svglite's `user_fonts`
#' keeps the substitute's metrics while recording "Arial" as the family, so
#' viewers use real Arial when they have it.
#'
#' Returns an empty list when Arial itself resolves, or when no metric-compatible
#' substitute is installed. In the latter case ggplot falls back to whatever the
#' system offers and glyph spacing may differ slightly from the published
#' figures; see README.
#'
#' @return A list suitable for svglite's `user_fonts` argument, possibly empty
#' @keywords internal
arial_font_alias <- function() {
    installed <- systemfonts::system_fonts()

    # Nothing to do if the machine has genuine Arial
    if (any(installed$family == "Arial")) {
        return(list())
    }

    # Metric-compatible with Arial, best first
    substitutes <- c("Liberation Sans", "Arimo")
    available <- substitutes[substitutes %in% installed$family]
    if (length(available) == 0) {
        return(list())
    }
    faces <- installed[installed$family == available[1], ]

    pick <- function(style) {
        hit <- faces$path[faces$style == style]
        if (length(hit) == 0) NULL else list(name = "Arial", file = hit[1])
    }

    alias <- list(
        plain = pick("Regular"),
        bold = pick("Bold"),
        italic = pick("Italic"),
        bolditalic = pick("Bold Italic")
    )
    alias <- alias[!vapply(alias, is.null, logical(1))]
    if (length(alias) == 0) list() else list(Arial = alias)
}

#' Save manuscript figure
#'
#' Saves the figure as SVG, drawing text in Arial (or a metric-compatible
#' stand-in) and recording "Arial" as the font family in the output.
#'
#' @param plot ggplot2 object or list with $plot element
#' @param filename Base filename (without extension)
#' @param path Output directory
#' @param width Figure width in inches
#' @param height Figure height in inches
#' @export
save_manuscript_figure <- function(plot,
                                    filename,
                                    path = "figures",
                                    width = 7,
                                    height = 5) {

    # Handle list input (extract $plot if present)
    if (is.list(plot) && "plot" %in% names(plot)) {
        plot <- plot$plot
    }

    # Create directory if needed
    if (!dir.exists(path)) {
        dir.create(path, recursive = TRUE)
    }

    outfile <- file.path(path, paste0(filename, ".svg"))

    ggplot2::ggsave(
        filename = outfile,
        plot = plot,
        width = width,
        height = height,
        bg = "white",
        device = svglite::svglite,
        user_fonts = arial_font_alias()
    )

    message("Saved: ", outfile)
}

