################################################################################
## MS-specific F1-F4 logistic regression and full-KEGG F4 sensitivity
################################################################################


library(dplyr)
library(ggplot2)



features <- c(
  F1 = "recognized_ebv_target",
  F2 = "log_recognized_ebv_length",
  F3 = "log_recognized_ebv_host_burden",
  F4 = "max_pathway_ebv_perturbation"
)


feature_sets <- list()

for (k in 1:4) {
  
  x <- combn(
    names(features),
    k,
    simplify = FALSE
  )
  
  feature_sets <- c(
    feature_sets,
    x
  )
}

model_names <- sapply(
  feature_sets,
  function(x) paste(x, collapse = "+")
)

names(feature_sets) <- model_names


logistic_fits <- list()
logistic_results <- list()

for (model_name in names(feature_sets)) {
  
  selected_F <- feature_sets[[model_name]]
  
  selected_variables <- unname(
    features[selected_F]
  )
  
  predictors <- c(
    M0,
    selected_variables
  )
  
  model_formula <- as.formula(
    paste(
      "ms_autoantigen ~",
      paste(
        predictors,
        collapse = " + "
      )
    )
  )
  
  fit <- glm(
    model_formula,
    data = ms_model_pathway,
    family = binomial()
  )
  
  logistic_fits[[model_name]] <- fit
  
  coef_table <- summary(fit)$coefficients
  
  coef_table <- coef_table[
    selected_variables,
    ,
    drop = FALSE
  ]
  
  temp <- data.frame(
    model = model_name,
    feature = selected_F,
    variable = selected_variables,
    beta = coef_table[, "Estimate"],
    SE = coef_table[, "Std. Error"],
    z = coef_table[, "z value"],
    p_value = coef_table[, "Pr(>|z|)"],
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  
  temp$OR <- exp(temp$beta)
  
  temp$CI_low <- exp(
    temp$beta - 1.96 * temp$SE
  )
  
  temp$CI_high <- exp(
    temp$beta + 1.96 * temp$SE
  )
  
  logistic_results[[model_name]] <- temp
}




logistic_results_all <- dplyr::bind_rows(
  logistic_results
)



logistic_results_display <- logistic_results_all %>%
  dplyr::mutate(
    beta = round(beta, 4),
    SE = round(SE, 4),
    z = round(z, 3),
    OR = round(OR, 3),
    CI_low = round(CI_low, 3),
    CI_high = round(CI_high, 3),
    P = format.pval(
      p_value,
      digits = 3,
      eps = 1e-300
    )
  ) %>%
  dplyr::select(
    model,
    feature,
    beta,
    SE,
    z,
    OR,
    CI_low,
    CI_high,
    P
  )


print(
  tibble::as_tibble(logistic_results_display),
  n = Inf
)



# ============================================================
# FOREST-PLOT DATA
# ============================================================

plot_df <- logistic_results_all %>%
  dplyr::mutate(
    
    estimate = ifelse(
      feature == "F4",
      beta * 0.1,
      beta
    ),
    
    SE_plot = ifelse(
      feature == "F4",
      SE * 0.1,
      SE
    ),
    
    lower = estimate - 1.96 * SE_plot,
    
    upper = estimate + 1.96 * SE_plot,
    
    OR_plot = exp(estimate),
    
    significance = ifelse(
      p_value < 0.05,
      "p < 0.05",
      "n.s."
    ),
    
    predictor = paste0(
      model,
      " : ",
      feature
    ),
    
    OR_label = paste0(
      "OR=",
      sprintf("%.2f", OR_plot)
    )
  )



plot_df$predictor <- factor(
  plot_df$predictor,
  levels = rev(plot_df$predictor)
)



p <- ggplot(
  plot_df,
  aes(
    x = estimate,
    y = predictor,
    colour = significance
  )
) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.6,
    colour = "grey40"
  ) +
  
  geom_errorbarh(
    aes(
      xmin = lower,
      xmax = upper
    ),
    height = 0.18,
    linewidth = 0.7
  ) +
  
  geom_point(
    size = 3
  ) +
  
  geom_text(
    aes(
      x = upper,
      label = OR_label
    ),
    hjust = -0.15,
    colour = "black",
    size = 3.3
  ) +
  
  scale_colour_manual(
    values = c(
      "n.s." = "grey55",
      "p < 0.05" = "#C62828"
    ),
    breaks = c(
      "n.s.",
      "p < 0.05"
    )
  ) +
  
  labs(
    title = "MS Autoantigen Model: Logistic Regression Coefficients",
    subtitle = "Coefficients with 95% confidence intervals and odds ratios",
    x = "Estimate (Log-Odds)",
    y = "Predictor",
    colour = "Significance"
  ) +
  
  theme_minimal(
    base_size = 13
  ) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 16,
      hjust = 0.5
    ),
    
    plot.subtitle = element_text(
      size = 12,
      hjust = 0.5
    ),
    
    axis.title.x = element_text(
      face = "bold",
      size = 13
    ),
    
    axis.title.y = element_text(
      face = "bold",
      size = 13
    ),
    
    axis.text.y = element_text(
      size = 9
    ),
    
    axis.text.x = element_text(
      size = 10
    ),
    
    legend.title = element_text(
      face = "bold"
    ),
    
    legend.position = "right",
    
    panel.grid.minor = element_blank(),
    
    panel.grid.major.y = element_line(
      colour = "grey90",
      linewidth = 0.5
    ),
    
    panel.grid.major.x = element_line(
      colour = "grey90",
      linewidth = 0.5
    ),
    
    plot.margin = margin(
      t = 10,
      r = 80,
      b = 10,
      l = 10
    )
  ) +
  
  coord_cartesian(
    clip = "off"
  )

# ============================================================
# F4 SENSITIVITY ANALYSIS
# FULL HUMAN KEGG PATHWAY UNIVERSE
# ============================================================

library(dplyr)


ebv_target_genes <- ppis_std %>%
  dplyr::filter(
    !is.na(Gene),
    Gene != ""
  ) %>%
  dplyr::distinct(Gene) %>%
  dplyr::pull(Gene)


kegg_full <- kegg_all_path_genes %>%
  dplyr::filter(
    !is.na(pathway_id),
    !is.na(Gene),
    Gene != ""
  ) %>%
  dplyr::distinct(
    pathway_id,
    Gene
  )


# ============================================================
# LEAVE-ONE-OUT EBV PERTURBATION
# ============================================================

kegg_full_perturbation <- kegg_full %>%
  
  dplyr::group_by(pathway_id) %>%
  
  dplyr::mutate(
    
    pathway_n_genes = dplyr::n(),
    
    pathway_n_ebv_targets =
      sum(Gene %in% ebv_target_genes),
    
    gene_is_ebv_target =
      as.integer(Gene %in% ebv_target_genes),
    
    loo_n_genes =
      pathway_n_genes - 1,
    
    loo_n_ebv_targets =
      pathway_n_ebv_targets -
      gene_is_ebv_target,
    
    loo_ebv_perturbation =
      ifelse(
        loo_n_genes > 0,
        loo_n_ebv_targets / loo_n_genes,
        NA_real_
      )
  ) %>%
  
  dplyr::ungroup()


gene_full_kegg_F4 <- kegg_full_perturbation %>%
  
  dplyr::filter(
    !is.na(loo_ebv_perturbation)
  ) %>%
  
  dplyr::group_by(Gene) %>%
  
  dplyr::summarise(
    
    F4_full_KEGG =
      max(
        loo_ebv_perturbation,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  ) %>%
  
  dplyr::rename(
    gene = Gene
  )


# ============================================================
# ADD FULL-KEGG F4 TO MS MODEL
# ============================================================

ms_model_full_kegg <- ms_model_pathway %>%
  
  dplyr::select(
    -dplyr::any_of("F4_full_KEGG")
  ) %>%
  
  dplyr::left_join(
    gene_full_kegg_F4,
    by = "gene"
  ) %>%
  
  dplyr::mutate(
    
    F4_full_KEGG =
      dplyr::coalesce(
        F4_full_KEGG,
        0
      )
  )


# ============================================================
#  BACKGROUND-ADJUSTED LOGISTIC REGRESSION
# ============================================================

fit_F4_full_kegg <- glm(
  
  ms_autoantigen ~
    
    log_human_length +
    log_string_degree_imp +
    log_kegg_pathways +
    in_string_graph +
    
    F4_full_KEGG,
  
  data = ms_model_full_kegg,
  
  family = binomial()
)


# ============================================================
#  F4 RESULT PER 10-PERCENTAGE-POINT INCREASE
# ============================================================

x <- summary(fit_F4_full_kegg)$coefficients[
  "F4_full_KEGG",
]

beta <- x["Estimate"]
SE   <- x["Std. Error"]
P    <- x["Pr(>|z|)"]


F4_full_kegg_result <- data.frame(
  
  beta = beta,
  
  OR_per_10pp =
    exp(beta * 0.1),
  
  CI_low =
    exp(
      (beta - 1.96 * SE) * 0.1
    ),
  
  CI_high =
    exp(
      (beta + 1.96 * SE) * 0.1
    ),
  
  P = P
)


print(F4_full_kegg_result)


# ============================================================
# REVISED FOREST PLOT
# Faceted by feature (F1-F4)
# Coefficients shown on log-odds scale
# Odds ratios printed beside estimates
# ============================================================

library(dplyr)
library(ggplot2)
library(forcats)


plot_df_facet <- plot_df %>%
  mutate(
    feature = factor(
      feature,
      levels = c("F1", "F2", "F3", "F4")
    ),
    
    # Keep model order stable within each facet
    model_label = factor(
      model,
      levels = unique(model)
    ),
    
    model_label = forcats::fct_rev(model_label)
  )



x_min <- min(plot_df_facet$lower, na.rm = TRUE)
x_max <- max(plot_df_facet$upper, na.rm = TRUE)

x_range <- x_max - x_min

# Add room on right for OR text
plot_right <- x_max + 0.35 * x_range

# Fixed position for all OR labels
or_x <- x_max + 0.05 * x_range



p <- ggplot(
  plot_df_facet,
  aes(
    x = estimate,
    y = model_label,
    colour = significance
  )
) +
  
  # Null line
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.55,
    colour = "grey40"
  ) +
  
  # 95% CI
  geom_errorbarh(
    aes(
      xmin = lower,
      xmax = upper
    ),
    height = 0.15,
    linewidth = 0.65
  ) +
  
  # Coefficient estimate
  geom_point(
    size = 2.8
  ) +
  
  # OR labels at fixed horizontal position
  geom_text(
    aes(
      x = or_x,
      label = OR_label
    ),
    inherit.aes = TRUE,
    hjust = 0,
    colour = "black",
    size = 3.0,
    show.legend = FALSE
  ) +
  


facet_wrap(
  ~ feature,
  ncol = 2,
  scales = "free_y",
  labeller = as_labeller(
    c(
      F1 = "F1: Recognized-EBV target",
      F2 = "F2: Viral protein length",
      F3 = "F3: Viral host-interaction burden",
      F4 = "F4: Pathway EBV perturbation"
    )
  )
) +
  

scale_colour_manual(
  values = c(
    "p < 0.05" = "#C62828",
    "n.s." = "grey55"
  ),
  breaks = c(
    "p < 0.05",
    "n.s."
  )
) +
  
  # Keep same x-axis across all four facets
  coord_cartesian(
    xlim = c(x_min, plot_right),
    clip = "off"
  ) +
  
  labs(
    x = "Coefficient estimate (log-odds)",
    y = "Model specification",
    colour = "Significance"
  ) +
  
theme_minimal(
  base_size = 14
) +
  
  theme(
    strip.text = element_text(
      face = "bold",
      size = 13,
      hjust = 0.5
    ),
    
    strip.background = element_rect(
      fill = "grey95",
      colour = NA
    ),
    
    axis.title.x = element_text(
      face = "bold",
      size = 14,
      hjust = 0.5,
      margin = margin(t = 8)
    ),
    
    axis.title.y = element_text(
      face = "bold",
      size = 14,
      hjust = 0.5
    ),
    
    axis.text.y = element_text(
      size = 11.5
    ),
    
    axis.text.x = element_text(
      size = 11
    ),
    
    # Center legend below entire figure
    legend.position = "bottom",
    legend.justification = "center",
    legend.box.just = "center",
    
    legend.title = element_text(
      face = "bold",
      size = 13
    ),
    
    panel.grid.minor = element_blank(),
    
    panel.grid.major.y = element_line(
      colour = "grey92",
      linewidth = 0.4
    ),
    
    panel.grid.major.x = element_line(
      colour = "grey90",
      linewidth = 0.4
    ),
    
    panel.spacing = unit(
      1.2,
      "lines"
    ),
    
    plot.margin = margin(
      t = 10,
      r = 20,
      b = 10,
      l = 10
    )
  )


print(p)

