################################################################################
## BROAD PROTEIN-LEVEL PREDICTIVE MODELLING 
## Network-blocked cross-validation of M0, M1, M2, M0+M1, and M3.
################################################################################

library(dplyr)

# model56_df contains the outcome, the original M0/M1/M2 features,
# and the network-blocked fold assignment.
model_final <- model56_df

stopifnot(
  nrow(model_final) == 20596,
  !anyNA(model_final$autoantigen),
  !anyNA(model_final$fold)
)




formula_M0 <- autoantigen ~
  log_human_length +
  log_string_degree_imp +
  log_kegg_pathways +
  in_string_graph



formula_M1 <- autoantigen ~
  has_8aa_ebv_match +
  log_8aa_matches +
  longest_exact_ebv_match


formula_M0_M1 <- autoantigen ~
  log_human_length +
  log_string_degree_imp +
  log_kegg_pathways +
  in_string_graph +
  has_8aa_ebv_match +
  log_8aa_matches +
  longest_exact_ebv_match


formula_M2 <- autoantigen ~
  direct_ebv_target +
  log_n_ebv_proteins +
  ebv_neighbor_fraction_imp +
  distance_to_ebv_target_imp +
  in_string_graph +
  string_disconnected


# ------------------------------------------------------------
# M3 = M0 + M1 + M2
# ------------------------------------------------------------

formula_M3 <- autoantigen ~
  log_human_length +
  log_string_degree_imp +
  log_kegg_pathways +
  in_string_graph +
  has_8aa_ebv_match +
  log_8aa_matches +
  longest_exact_ebv_match +
  direct_ebv_target +
  log_n_ebv_proteins +
  ebv_neighbor_fraction_imp +
  distance_to_ebv_target_imp +
  string_disconnected



formulas_final <- list(
  M0 = formula_M0,
  M1 = formula_M1,
  M2 = formula_M2,
  M0_M1 = formula_M0_M1,
  M3 = formula_M3
)



pred_final <- model_final %>%
  dplyr::select(
    gene,
    autoantigen,
    fold
  )



for (model_name in names(formulas_final)) {
  
  pred_col <-
    paste0(
      "pred_",
      model_name
    )
  
  pred_final[[pred_col]] <- NA_real_
  
  
  for (f in sort(unique(model_final$fold))) {
    
    train <- model_final %>%
      dplyr::filter(
        fold != f
      )
    
    test <- model_final %>%
      dplyr::filter(
        fold == f
      )
    
    
    fit <- glm(
      formulas_final[[model_name]],
      data = train,
      family = binomial()
    )
    
    
    idx <-
      model_final$fold == f
    
    
    # KEEP THIS [[ ]] EXPRESSION ON ONE LINE
    pred_final[[pred_col]][idx] <- predict(
      fit,
      newdata = test,
      type = "response"
    )
  }
}




prediction_columns <-
  paste0(
    "pred_",
    names(formulas_final)
  )

stopifnot(
  all(
    sapply(
      pred_final[
        prediction_columns
      ],
      function(x) {
        !anyNA(x)
      }
    )
  )
)


# ============================================================
# POOLED AUPRC + AUROC
# ============================================================

results_final <- data.frame(
  
  model =
    names(formulas_final),
  
  AUPRC =
    NA_real_,
  
  AUROC =
    NA_real_,
  
  stringsAsFactors = FALSE
)


for (i in seq_len(nrow(results_final))) {
  
  model_name <-
    results_final$model[i]
  
  pred_col <-
    paste0(
      "pred_",
      model_name
    )
  
  
  results_final$AUPRC[i] <-
    calculate_auprc(
      pred_final$autoantigen,
      pred_final[[pred_col]]
    )
  
  
  results_final$AUROC[i] <-
    calculate_auroc(
      pred_final$autoantigen,
      pred_final[[pred_col]]
    )
}

