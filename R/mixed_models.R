# Linear Mixed Models for APC Function Analysis
#
# Functions for fitting LMMs to properly account for pseudoreplication
# (multiple cells per donor) in scRNA-seq data.
#
# Includes hurdle (two-part) models using glmmTMB for zero-inflated data:
# - Part 1 (logistic): P(APC > 0) ~ cytokine_binary + covariates + (1|donor)
# - Part 2 (linear): APC | APC > 0 ~ cytokine + covariates + (1|donor)

#' Compare Random Intercept vs Random Slope Models
#'
#' Tests whether the predictor~outcome relationship varies by donor using
#' a likelihood ratio test comparing random intercepts vs random slopes models.
#'
#' @param cell_data Cell-level data frame
#' @param outcome Outcome variable column name
#' @param predictor Predictor variable column name
#' @param donor_col Column identifying donors (default: "orig.ident")
#' @param covariates Optional character vector of covariate columns
#' @param filter_nonzero Filter to non-zero values for both outcome and predictor
#'
#' @return List with:
#'   - lrt: Likelihood ratio test results (chisq, df, p_value)
#'   - intercept_model: Random intercept model object
#'   - slope_model: Random slope model object
#'   - donor_slopes: Per-donor slope estimates with CIs
#'   - fixed_effects_comparison: Side-by-side fixed effects
#'   - n_cells, n_donors: Sample sizes
#'
#' @export
compare_random_effects <- function(cell_data,
                                   outcome,
                                   predictor,
                                   donor_col = "orig.ident",
                                   covariates = NULL,
                                   filter_nonzero = TRUE) {

    # Prepare data
    df <- cell_data %>%
        dplyr::select(dplyr::all_of(c(outcome, predictor, donor_col, covariates))) %>%
        dplyr::filter(!is.na(.data[[outcome]]), !is.na(.data[[predictor]]))

    if (filter_nonzero) {
        df <- df %>%
            dplyr::filter(.data[[outcome]] > 0, .data[[predictor]] > 0)
    }

    n_cells <- nrow(df)
    n_donors <- dplyr::n_distinct(df[[donor_col]])

    # Build formulas
    fixed_part <- predictor
    if (!is.null(covariates) && length(covariates) > 0) {
        fixed_part <- paste(c(predictor, covariates), collapse = " + ")
    }

    # Random intercept model
    formula_intercept <- stats::as.formula(
        paste(outcome, "~", fixed_part, "+ (1|", donor_col, ")")
    )

    # Random slope model (correlated intercepts and slopes)
    formula_slope <- stats::as.formula(
        paste(outcome, "~", fixed_part, "+ (", predictor, "|", donor_col, ")")
    )

    # Fit models with ML for LRT (not REML)
    intercept_model <- lme4::lmer(formula_intercept, data = df, REML = FALSE)

    # Random slopes can fail to converge - handle gracefully
    slope_model <- tryCatch({
        lme4::lmer(formula_slope, data = df, REML = FALSE,
                   control = lme4::lmerControl(
                       optimizer = "bobyqa",
                       optCtrl = list(maxfun = 100000)
                   ))
    }, error = function(e) {
        warning("Random slope model failed to converge: ", e$message)
        NULL
    })

    # If slope model failed, return early
    if (is.null(slope_model)) {
        return(list(
            lrt = list(chisq = NA, df = NA, p_value = NA),
            intercept_model = intercept_model,
            slope_model = NULL,
            donor_slopes = NULL,
            convergence_failed = TRUE,
            n_cells = n_cells,
            n_donors = n_donors,
            formula_intercept = as.character(formula_intercept),
            formula_slope = as.character(formula_slope)
        ))
    }

    # Likelihood ratio test
    lrt_result <- stats::anova(intercept_model, slope_model)
    lrt <- list(
        chisq = lrt_result$Chisq[2],
        df = lrt_result$Df[2],
        p_value = lrt_result$`Pr(>Chisq)`[2]
    )

    # Extract per-donor slopes (BLUPs)
    ranef_slope <- lme4::ranef(slope_model)[[donor_col]]
    fixef_slope <- lme4::fixef(slope_model)

    # Get the predictor coefficient name (might have been renamed)
    pred_coef_name <- predictor

    # Per-donor slopes = fixed effect + random effect
    donor_slopes <- tibble::tibble(
        donor = rownames(ranef_slope),
        random_intercept = ranef_slope[, "(Intercept)"],
        random_slope = ranef_slope[, pred_coef_name],
        fixed_slope = fixef_slope[pred_coef_name],
        total_slope = fixed_slope + random_slope
    )

    # Get approximate CIs for random slopes using conditional variance
    # (simplified - full inference would use bootMer)
    re_var <- as.data.frame(lme4::VarCorr(slope_model))
    slope_sd <- re_var %>%
        dplyr::filter(grp == donor_col, var1 == pred_coef_name, is.na(var2)) %>%
        dplyr::pull(sdcor)

    if (length(slope_sd) > 0) {
        donor_slopes <- donor_slopes %>%
            dplyr::mutate(
                slope_se = slope_sd,  # Approximate
                ci_lower = total_slope - 1.96 * slope_se,
                ci_upper = total_slope + 1.96 * slope_se
            )
    }

    # Fixed effects comparison
    fe_intercept <- broom.mixed::tidy(intercept_model, effects = "fixed") %>%
        dplyr::mutate(model = "Random Intercept")
    fe_slope <- broom.mixed::tidy(slope_model, effects = "fixed") %>%
        dplyr::mutate(model = "Random Slope")
    fixed_effects_comparison <- dplyr::bind_rows(fe_intercept, fe_slope)

    list(
        lrt = lrt,
        intercept_model = intercept_model,
        slope_model = slope_model,
        donor_slopes = donor_slopes,
        fixed_effects_comparison = fixed_effects_comparison,
        convergence_failed = FALSE,
        n_cells = n_cells,
        n_donors = n_donors,
        formula_intercept = as.character(formula_intercept),
        formula_slope = as.character(formula_slope),
        outcome = outcome,
        predictor = predictor
    )
}
#' Fit Hurdle Mixed Model for APC Induction and Dose-Response
#'
#' Fits a two-part hurdle model to properly analyze zero-inflated single-cell
#' gene signature data without pre-filtering on the predictor:
#'
#' - Part 1 (binomial): P(APC > 0) ~ cytokine_binary + covariates + (1|donor)
#'   Tests whether cytokine exposure increases odds of APC activation
#'
#' - Part 2 (truncated normal): APC | APC > 0 ~ cytokine + covariates + (1|donor)
#'   Tests dose-response among responding cells only
#'
#' @param cell_data Cell-level data frame with UCell scores
#' @param outcome APC signature column name (e.g., "apc_function_full_UCell")
#' @param predictor IFN-γ signature column name (e.g., "ifng_chemokine_effector_UCell")
#' @param donor_col Column identifying donors (default: "orig.ident")
#' @param covariates Optional character vector of covariate columns (e.g., "sample_age")
#'
#' @return List with:
#'   - part1_model: glmmTMB logistic model (induction)
#'   - part2_model: lme4 linear model (dose-response, APC > 0 cells only)
#'   - odds_ratio: List with estimate, ci_lower, ci_upper, p_value
#'   - beta_std: List with estimate, ci_lower, ci_upper, se, p_value
#'   - n_cells_total, n_cells_apc_pos, n_donors
#'   - outcome, predictor (for reference)
#'
#' @export
fit_hurdle_model <- function(cell_data,
                             outcome,
                             predictor,
                             donor_col = "orig.ident",
                             covariates = NULL) {

    # Prepare data - create binary variables
    df <- cell_data %>%
        dplyr::mutate(
            apc_binary = as.integer(.data[[outcome]] > 0),
            cytokine_binary = as.integer(.data[[predictor]] > 0)
        ) %>%
        dplyr::filter(!is.na(.data[[outcome]]), !is.na(.data[[predictor]]))

    # Also filter for covariates if specified
    if (!is.null(covariates) && length(covariates) > 0) {
        for (cov in covariates) {
            df <- df %>% dplyr::filter(!is.na(.data[[cov]]))
        }
    }

    n_cells_total <- nrow(df)
    n_donors <- dplyr::n_distinct(df[[donor_col]])

    # -------------------------------------------------------------------------
    # Part 1: Logistic model for P(APC > 0)
    # Binary predictor: does having any cytokine increase odds of APC activation?
    # -------------------------------------------------------------------------
    formula_p1_str <- paste(
        "apc_binary ~ cytokine_binary",
        if (!is.null(covariates) && length(covariates) > 0)
            paste("+", paste(covariates, collapse = " + "))
        else "",
        "+ (1|", donor_col, ")"
    )
    formula_p1 <- stats::as.formula(formula_p1_str)

    part1_model <- glmmTMB::glmmTMB(
        formula_p1,
        data = df,
        family = stats::binomial(link = "logit")
    )

    # Extract OR for cytokine effect
    part1_summary <- summary(part1_model)
    coef_cytokine <- part1_summary$coefficients$cond["cytokine_binary", ]

    # Get confidence intervals
    ci_matrix <- stats::confint(part1_model, parm = "cytokine_binary", method = "Wald")

    or <- exp(coef_cytokine["Estimate"])
    or_ci_lower <- exp(ci_matrix[1])
    or_ci_upper <- exp(ci_matrix[2])
    or_p <- coef_cytokine["Pr(>|z|)"]

    # -------------------------------------------------------------------------
    # Part 2: Linear model for APC | APC > 0 & cytokine > 0
    # Double-positive cells: both APC and cytokine programs active
    # -------------------------------------------------------------------------
    df_apc_pos <- df %>% dplyr::filter(.data[[outcome]] > 0, .data[[predictor]] > 0)
    n_cells_apc_pos <- nrow(df_apc_pos)
    n_donors_apc_pos <- dplyr::n_distinct(df_apc_pos[[donor_col]])

    formula_p2_str <- paste(
        outcome, "~", predictor,
        if (!is.null(covariates) && length(covariates) > 0)
            paste("+", paste(covariates, collapse = " + "))
        else "",
        "+ (1|", donor_col, ")"
    )
    formula_p2 <- stats::as.formula(formula_p2_str)

    part2_model <- lmerTest::lmer(formula_p2, data = df_apc_pos, REML = TRUE)

    # Extract and standardize beta
    part2_fe <- broom.mixed::tidy(part2_model, effects = "fixed", conf.int = TRUE)
    pred_row <- part2_fe %>% dplyr::filter(term == predictor)

    beta <- pred_row$estimate
    beta_se <- pred_row$std.error
    beta_ci_lower <- pred_row$conf.low
    beta_ci_upper <- pred_row$conf.high
    beta_p <- pred_row$p.value

    # Standardize: beta* = beta * SD_x / SD_y
    sd_x <- stats::sd(df_apc_pos[[predictor]], na.rm = TRUE)
    sd_y <- stats::sd(df_apc_pos[[outcome]], na.rm = TRUE)
    beta_std <- beta * sd_x / sd_y
    beta_std_se <- beta_se * sd_x / sd_y
    beta_std_ci_lower <- beta_ci_lower * sd_x / sd_y
    beta_std_ci_upper <- beta_ci_upper * sd_x / sd_y

    # -------------------------------------------------------------------------
    # Return results
    # -------------------------------------------------------------------------
    list(
        part1_model = part1_model,
        part2_model = part2_model,
        odds_ratio = list(
            estimate = as.numeric(or),
            ci_lower = as.numeric(or_ci_lower),
            ci_upper = as.numeric(or_ci_upper),
            p_value = as.numeric(or_p)
        ),
        beta_raw = list(
            estimate = beta,
            se = beta_se,
            ci_lower = beta_ci_lower,
            ci_upper = beta_ci_upper,
            p_value = beta_p
        ),
        beta_std = list(
            estimate = beta_std,
            se = beta_std_se,
            ci_lower = beta_std_ci_lower,
            ci_upper = beta_std_ci_upper,
            p_value = beta_p
        ),
        n_cells_total = n_cells_total,
        n_cells_apc_pos = n_cells_apc_pos,
        n_donors = n_donors,
        n_donors_apc_pos = n_donors_apc_pos,
        formula_p1 = formula_p1_str,
        formula_p2 = formula_p2_str,
        outcome = outcome,
        predictor = predictor
    )
}

#' Plot Per-Donor Slopes (Forest Plot)
#'
#' Creates a forest plot showing per-donor slope estimates from a random
#' slopes model, with confidence intervals and optional metadata annotation
#' bars on the right side of the plot.
#'
#' @param random_effects_result Result from compare_random_effects()
#' @param donor_metadata Optional tibble with donor metadata for annotation.
#'   Must have 'donor' column matching donor IDs. Other columns will be displayed
#'   as annotation bars (e.g., diabetes_status, sample_sex, sample_age).
#' @param title Plot title
#' @param xlab X-axis label
#' @param highlight_significant Highlight donors with CIs not crossing zero
#' @param show_metadata Whether to show metadata annotation bars (default: TRUE if provided)
#'
#' @return ggplot2 object with annotation bars drawn using annotate()
#'
#' @export
plot_donor_slopes <- function(random_effects_result,
                              donor_metadata = NULL,
                              title = "Per-Donor Slopes",
                              xlab = "Slope (Effect of Predictor)",
                              highlight_significant = TRUE,
                              show_metadata = TRUE) {

    if (random_effects_result$convergence_failed || is.null(random_effects_result$donor_slopes)) {
        return(
            ggplot2::ggplot() +
                ggplot2::annotate("text", x = 0.5, y = 0.5,
                                 label = "Random slopes model did not converge",
                                 size = 5) +
                ggplot2::theme_void()
        )
    }

    plot_data <- random_effects_result$donor_slopes %>%
        dplyr::arrange(total_slope) %>%
        dplyr::mutate(
            donor = factor(donor, levels = donor),
            significant = ci_lower > 0 | ci_upper < 0,
            y_num = as.numeric(donor)
        )

    # Fixed effect reference line
    fixed_slope <- random_effects_result$donor_slopes$fixed_slope[1]

    # LRT annotation
    lrt <- random_effects_result$lrt
    lrt_label <- sprintf("LRT: \u03c7\u00b2 = %.2f, df = %d, p = %.3f",
                         lrt$chisq, lrt$df, lrt$p_value)

    # Get x-axis range for positioning annotation bars
    x_range <- range(c(plot_data$ci_lower, plot_data$ci_upper), na.rm = TRUE)
    x_span <- diff(x_range)

    # Define color palettes
    status_colors <- c("NODM" = "#3fa36b", "AABP" = "#F4C542", "T1DM" = "#9353b3")
    sex_colors <- c("female" = "#E91E63", "male" = "#2196F3", "Female" = "#E91E63", "Male" = "#2196F3")

    # Build base forest plot
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = total_slope, y = donor)) +
        ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
        ggplot2::geom_errorbarh(
            ggplot2::aes(xmin = ci_lower, xmax = ci_upper),
            height = 0.3,
            color = ifelse(highlight_significant & plot_data$significant,
                          "#2E86AB", "gray60")
        ) +
        ggplot2::geom_point(
            size = 2,
            color = ifelse(highlight_significant & plot_data$significant,
                          "#2E86AB", "gray40")
        )

    # Add metadata annotation bars if provided
    if (!is.null(donor_metadata) && show_metadata) {
        # Join metadata
        plot_data <- plot_data %>%
            dplyr::left_join(donor_metadata, by = "donor")

        meta_cols <- setdiff(names(donor_metadata), "donor")
        n_meta <- length(meta_cols)

        if (n_meta > 0) {
            # Bar positioning: right of plot, each bar width = 0.08 of x_span
            bar_width <- x_span * 0.08
            bar_gap <- x_span * 0.02
            bar_start <- x_range[2] + x_span * 0.05  # Start 5% beyond plot

            for (i in seq_along(meta_cols)) {
                col <- meta_cols[i]
                col_vals <- as.character(plot_data[[col]])

                # Calculate bar x positions
                xmin <- bar_start + (i - 1) * (bar_width + bar_gap)
                xmax <- xmin + bar_width

                # Get fill colors based on column type
                if (col == "diabetes_status") {
                    fills <- status_colors[col_vals]
                    fills[is.na(fills)] <- "gray80"
                    labels <- col_vals
                    text_col <- "white"
                } else if (col == "sample_sex") {
                    fills <- sex_colors[col_vals]
                    fills[is.na(fills)] <- "gray80"
                    labels <- substr(col_vals, 1, 1)  # M/F
                    text_col <- "white"
                } else if (col == "sample_age") {
                    ages <- as.numeric(col_vals)
                    age_norm <- (ages - min(ages, na.rm = TRUE)) /
                                (max(ages, na.rm = TRUE) - min(ages, na.rm = TRUE))
                    age_norm[is.na(age_norm)] <- 0.5
                    fills <- scales::col_numeric(
                        palette = c("#2E8B57", "#8B4513"),
                        domain = c(0, 1)
                    )(age_norm)
                    labels <- col_vals
                    text_col <- "white"
                } else if (col == "cytokine_pos") {
                    # Cytokine status: Yes = red/salmon, No = gray
                    fills <- ifelse(col_vals == "Yes", "#FFCCCB", "gray80")
                    labels <- col_vals
                    text_col <- "black"
                } else {
                    fills <- rep("gray70", length(col_vals))
                    labels <- col_vals
                    text_col <- "black"
                }

                # Add rectangles for each donor
                p <- p + ggplot2::annotate(
                    "rect",
                    xmin = xmin, xmax = xmax,
                    ymin = plot_data$y_num - 0.4,
                    ymax = plot_data$y_num + 0.4,
                    fill = fills
                )

                # Add text labels
                # family is set explicitly: annotate() does not inherit it from
                # the theme applied later
                p <- p + ggplot2::annotate(
                    "text",
                    x = (xmin + xmax) / 2,
                    y = plot_data$y_num,
                    label = labels,
                    size = 2.5,
                    color = text_col,
                    fontface = "bold",
                    family = "Arial"
                )

                # Add column header
                col_label <- dplyr::case_when(
                    col == "diabetes_status" ~ "Status",
                    col == "sample_sex" ~ "Sex",
                    col == "sample_age" ~ "Age",
                    col == "cytokine_pos" ~ "Cyt>0",
                    TRUE ~ col
                )
                p <- p + ggplot2::annotate(
                    "text",
                    x = (xmin + xmax) / 2,
                    y = nrow(plot_data) + 1.2,
                    label = col_label,
                    size = 3,
                    fontface = "bold",
                    family = "Arial"
                )
            }

            # Extend x-axis to accommodate annotations
            x_max_extended <- bar_start + n_meta * (bar_width + bar_gap)
            p <- p + ggplot2::coord_cartesian(
                xlim = c(x_range[1] - x_span * 0.05, x_max_extended),
                ylim = c(0.5, nrow(plot_data) + 1.5),
                clip = "off"
            )
        }
    }

    # Add labels and theme
    p <- p +
        ggplot2::labs(
            title = title,
            subtitle = lrt_label,
            x = xlab,
            y = "Donor"
        ) +
        ggplot2::theme_minimal(base_size = 11) +
        ggplot2::theme(
            axis.text.y = ggplot2::element_text(size = 8),
            plot.title = ggplot2::element_text(face = "bold", size = 12),
            plot.margin = ggplot2::margin(10, 60, 10, 10)  # Extra right margin for bars
        )

    p
}
