
################################################################################
## Network-blocked broad prediction + exact peptide-sharing / EBV-protein analysis
################################################################################


library(dplyr)

set.seed(20260929)

K <- 5L


features_M0 <- c(
  "log_human_length",
  "log_string_degree_imp",
  "log_kegg_pathways",
  "in_string_graph"
)

features_M1 <- c(
  "has_8aa_ebv_match",
  "log_8aa_matches",
  "longest_exact_ebv_match"
)

features_M2 <- c(
  "direct_ebv_target",
  "log_n_ebv_proteins",
  "ebv_neighbor_fraction_imp",
  "distance_to_ebv_target_imp",
  "string_disconnected"
)

model_features <- list(
  M0 = features_M0,
  M1 = features_M1,
  M2 = features_M2,
  M0_M1 = c(features_M0, features_M1),
  M3 = c(features_M0, features_M1, features_M2)
)


required_vars <- unique(
  c(
    "gene",
    "autoantigen",
    "community",
    unlist(model_features)
  )
)

missing_vars <- setdiff(
  required_vars,
  names(model_blocked)
)

if (length(missing_vars) > 0) {
  stop(
    paste(
      "Missing variables:",
      paste(missing_vars, collapse = ", ")
    )
  )
}




community_sizes <- model_blocked %>%
  dplyr::filter(
    !is.na(.data$community)
  ) %>%
  dplyr::count(
    .data$community,
    name = "n"
  ) %>%
  dplyr::arrange(
    dplyr::desc(.data$n)
  )


fold_sizes <- rep(0L, K)

community_fold <- integer(
  nrow(community_sizes)
)


for (i in seq_len(nrow(community_sizes))) {
  
  target_fold <- which.min(
    fold_sizes
  )
  
  community_fold[i] <- target_fold
  
  fold_sizes[target_fold] <-
    fold_sizes[target_fold] +
    community_sizes$n[i]
}


community_assignment <- data.frame(
  community = community_sizes$community,
  fold = community_fold
)



model_blocked <- model_blocked %>%
  dplyr::left_join(
    community_assignment,
    by = "community"
  )


# ------------------------------------------------------------
#  Assign genes WITHOUT STRING communities
# ------------------------------------------------------------

outside_idx <- which(
  is.na(model_blocked$community)
)

outside_positive <- outside_idx[
  model_blocked$autoantigen[outside_idx] == 1
]

outside_unlabelled <- outside_idx[
  model_blocked$autoantigen[outside_idx] == 0
]


outside_positive <- sample(
  outside_positive
)

outside_unlabelled <- sample(
  outside_unlabelled
)


model_blocked$fold[outside_positive] <-
  rep(
    1:K,
    length.out = length(outside_positive)
  )

model_blocked$fold[outside_unlabelled] <-
  rep(
    1:K,
    length.out = length(outside_unlabelled)
  )


model_blocked$fold <- as.integer(
  model_blocked$fold
)




fold_diagnostic <- model_blocked %>%
  dplyr::group_by(
    .data$fold
  ) %>%
  dplyr::summarise(
    n = dplyr::n(),
    positives = sum(.data$autoantigen),
    prevalence = mean(.data$autoantigen),
    community_genes =
      sum(!is.na(.data$community)),
    .groups = "drop"
  )



if (any(is.na(model_blocked$fold))) {
  stop("Some genes were not assigned to a fold.")
}



community_fold_check <- model_blocked %>%
  dplyr::filter(
    !is.na(.data$community)
  ) %>%
  dplyr::group_by(
    .data$community
  ) %>%
  dplyr::summarise(
    n_folds =
      dplyr::n_distinct(.data$fold),
    .groups = "drop"
  )


n_leaking <- sum(
  community_fold_check$n_folds > 1
)



if (n_leaking > 0) {
  stop(
    "STRING community leakage detected."
  )
}


# ------------------------------------------------------------
#  Metric functions
# ------------------------------------------------------------

calc_auprc <- function(y, p) {
  
  ord <- order(
    p,
    decreasing = TRUE
  )
  
  y <- y[ord]
  
  tp <- cumsum(y == 1)
  fp <- cumsum(y == 0)
  
  precision <- tp / (tp + fp)
  
  sum(
    precision[y == 1]
  ) / sum(y == 1)
}


calc_auroc <- function(y, p) {
  
  n_pos <- sum(y == 1)
  n_neg <- sum(y == 0)
  
  ranks <- rank(
    p,
    ties.method = "average"
  )
  
  (
    sum(ranks[y == 1]) -
      n_pos * (n_pos + 1) / 2
  ) /
    (n_pos * n_neg)
}



model_names <- names(model_features)

for (model_name in model_names) {
  prediction_column <- paste0("pred_", model_name)
  model_blocked[, prediction_column] <- NA_real_
}

blocked_fold_results <- list()
result_counter <- 1L


# ------------------------------------------------------------
#  Run blocked CV
# ------------------------------------------------------------


for (fold_i in 1:K) {
  
    
  train_idx <- which(model_blocked$fold != fold_i)
  test_idx  <- which(model_blocked$fold == fold_i)
  
  train_df <- model_blocked[
    train_idx,
    ,
    drop = FALSE
  ]
  
  test_df <- model_blocked[
    test_idx,
    ,
    drop = FALSE
  ]
  
  for (model_name in model_names) {
    
    features_i <- model_features[[model_name]]
    
    formula_i <- as.formula(
      paste(
        "autoantigen ~",
        paste(features_i, collapse = " + ")
      )
    )
    
    fit_i <- glm(
      formula_i,
      data = train_df,
      family = binomial()
    )
    
    pred_i <- predict(
      fit_i,
      newdata = test_df,
      type = "response"
    )
    
    prediction_column <- paste0(
      "pred_",
      model_name
    )
    
    model_blocked[
      test_idx,
      prediction_column
    ] <- pred_i
    
    fold_auprc <- calc_auprc(
      test_df$autoantigen,
      pred_i
    )
    
    fold_auroc <- calc_auroc(
      test_df$autoantigen,
      pred_i
    )
    
    fold_prevalence <- mean(
      test_df$autoantigen
    )
    
    result_row <- data.frame(
      fold = fold_i,
      model = model_name,
      n_test = nrow(test_df),
      n_positive = sum(test_df$autoantigen),
      prevalence = fold_prevalence,
      AUPRC = fold_auprc,
      AUROC = fold_auroc,
      AUPRC_over_prevalence =
        fold_auprc / fold_prevalence,
      stringsAsFactors = FALSE
    )
    
    blocked_fold_results[[result_counter]] <- result_row
    
    result_counter <- result_counter + 1L
    
      }
}


blocked_fold_results_df <- do.call(
  rbind,
  blocked_fold_results
)




for (model_name in model_names) {
  
  prediction_column <- paste0(
    "pred_",
    model_name
  )
  
  }


# ------------------------------------------------------------
# Pooled network-blocked performance
# ------------------------------------------------------------

blocked_prevalence <- mean(
  model_blocked$autoantigen
)

blocked_pooled_results <- list()

for (i in seq_along(model_names)) {
  
  model_name <- model_names[i]
  
  prediction_column <- paste0(
    "pred_",
    model_name
  )
  
  p <- model_blocked[, prediction_column]
  
  model_auprc <- calc_auprc(
    model_blocked$autoantigen,
    p
  )
  
  model_auroc <- calc_auroc(
    model_blocked$autoantigen,
    p
  )
  
  blocked_pooled_results[[i]] <- data.frame(
    model = model_name,
    AUPRC = model_auprc,
    AUROC = model_auroc,
    prevalence = blocked_prevalence,
    AUPRC_over_prevalence =
      model_auprc / blocked_prevalence,
    stringsAsFactors = FALSE
  )
}

blocked_pooled_results <- do.call(
  rbind,
  blocked_pooled_results
)


# ------------------------------------------------------------
#  comparison: M3 vs M0+M1
# ------------------------------------------------------------

m0m1 <- blocked_pooled_results[
  blocked_pooled_results$model == "M0_M1",
  ,
  drop = FALSE
]

m3 <- blocked_pooled_results[
  blocked_pooled_results$model == "M3",
  ,
  drop = FALSE
]

delta_auprc <-
  m3$AUPRC - m0m1$AUPRC

delta_auroc <-
  m3$AUROC - m0m1$AUROC


# ------------------------------------------------------------
#  Fold-specific M3 minus M0+M1
# ------------------------------------------------------------

paired_blocked <- merge(
  blocked_fold_results_df[
    blocked_fold_results_df$model == "M0_M1",
    ,
    drop = FALSE
  ],
  blocked_fold_results_df[
    blocked_fold_results_df$model == "M3",
    ,
    drop = FALSE
  ],
  by = "fold",
  suffixes = c("_M0_M1", "_M3")
)

paired_blocked$delta_AUPRC <-
  paired_blocked$AUPRC_M3 -
  paired_blocked$AUPRC_M0_M1

paired_blocked$delta_AUROC <-
  paired_blocked$AUROC_M3 -
  paired_blocked$AUROC_M0_M1




library(dplyr)
library(tidyr)
library(ggplot2)



required_objects <- c(
  "ebv_sequences",
  "human_sequences",
  "prediction_universe"
)

missing_objects <- required_objects[
  !vapply(required_objects, exists, logical(1))
]

if (length(missing_objects) > 0) {
  stop(
    "Missing objects: ",
    paste(missing_objects, collapse = ", ")
  )
}



ebv_seq <- as.character(ebv_sequences)
human_seq <- as.character(human_sequences)

names(ebv_seq) <- names(ebv_sequences)
names(human_seq) <- names(human_sequences)



human_status <- human_accession_map %>%
  
  dplyr::transmute(
    gene = as.character(gene),
    human_id = as.character(human_uniprot)
  ) %>%
  
  dplyr::left_join(
    
    prediction_universe %>%
      dplyr::transmute(
        gene = as.character(gene),
        autoantigen = as.integer(autoantigen)
      ),
    
    by = "gene"
  ) %>%
  
  dplyr::filter(
    !is.na(human_id),
    human_id != "",
    !is.na(autoantigen)
  ) %>%
  
  dplyr::distinct(
    human_id,
    .keep_all = TRUE
  )


# ============================================================
#  MATCH HUMAN ACCESSIONS TO SEQUENCE IDs
# ============================================================

human_status <- human_status %>%
  dplyr::filter(
    human_id %in% names(human_sequences)
  )



if (nrow(human_status) == 0) {
  
  stop(
    "human_uniprot IDs do not match names(human_sequences). ",
    "Paste the two ID examples printed immediately above."
  )
}

# ============================================================
# FUNCTION: UNIQUE EXACT 8-MERS FROM ONE EBV PROTEIN
# ============================================================

get_kmers <- function(sequence, k = 8) {
  
  sequence <- toupper(sequence)
  
  L <- nchar(sequence)
  
  if (
    is.na(L) ||
    L < k
  ) {
    return(character(0))
  }
  
  starts <- seq_len(L - k + 1)
  
  kmers <- substring(
    sequence,
    starts,
    starts + k - 1
  )
  
  unique(kmers)
}


# ============================================================
#  PREPARE HUMAN SEQUENCES
# ============================================================

human_ids_use <- human_status$human_id

human_seq_use <- human_seq[
  human_ids_use
]

human_seq_use <- toupper(
  as.character(human_seq_use)
)

names(human_seq_use) <- human_ids_use


auto_lookup <- human_status$autoantigen
names(auto_lookup) <- human_status$human_id



ebv_results <- vector(
  "list",
  length(ebv_seq)
)

names(ebv_results) <- names(ebv_seq)




for (i in seq_along(ebv_seq)) {
  
  viral_id <- names(ebv_seq)[i]
  
  viral_sequence <- toupper(
    ebv_seq[[i]]
  )
  
  viral_length <- nchar(
    viral_sequence
  )
  
  viral_8mers <- get_kmers(
    viral_sequence,
    k = 8
  )
  
  
 # ----------------------------------------------------------
 # Determine whether each human protein contains ANY
 # exact 8-mer from this particular EBV protein.
 # ----------------------------------------------------------
  
  matched <- vapply(
    human_seq_use,
    function(hseq) {
      
      if (
        is.na(hseq) ||
        nchar(hseq) < 8 ||
        length(viral_8mers) == 0
      ) {
        return(FALSE)
      }
      
      any(
        vapply(
          viral_8mers,
          function(kmer) {
            grepl(
              kmer,
              hseq,
              fixed = TRUE
            )
          },
          logical(1)
        )
      )
    },
    logical(1)
  )
  
  
  matched_ids <- names(
    matched
  )[matched]
  
  
  matched_status <- auto_lookup[
    matched_ids
  ]
  
  
  n_human_matches <- length(
    matched_ids
  )
  
  n_autoantigen_matches <- sum(
    matched_status == 1,
    na.rm = TRUE
  )
  
  n_unlabelled_matches <- sum(
    matched_status == 0,
    na.rm = TRUE
  )
  
  
  autoantigen_fraction <- if (
    n_human_matches > 0
  ) {
    
    n_autoantigen_matches /
      n_human_matches
    
  } else {
    
    NA_real_
    
  }
  
  
  ebv_results[[i]] <- data.frame(
    
    ebv_protein = viral_id,
    
    ebv_length = viral_length,
    
    n_unique_8mers =
      length(viral_8mers),
    
    n_matched_human_proteins =
      n_human_matches,
    
    n_matched_autoantigens =
      n_autoantigen_matches,
    
    n_matched_unlabelled =
      n_unlabelled_matches,
    
    autoantigen_fraction =
      autoantigen_fraction,
    
    stringsAsFactors = FALSE
  )
  
  
  }


ebv_size_results <- dplyr::bind_rows(
  ebv_results
)



ebv_size_results <- tibble::as_tibble(ebv_size_results)

ebv_size_results_print <- ebv_size_results %>%
  dplyr::arrange(
    dplyr::desc(n_matched_autoantigens)
  )




test_auto <- cor.test(
  ebv_size_results$ebv_length,
  ebv_size_results$n_matched_autoantigens,
  method = "spearman",
  exact = FALSE
)


# ============================================================
# CONTROL TEST
# ============================================================



test_all <- cor.test(
  ebv_size_results$ebv_length,
  ebv_size_results$n_matched_human_proteins,
  method = "spearman",
  exact = FALSE
)

# ============================================================
# SPECIFICITY TEST
# ============================================================



fraction_df <- ebv_size_results %>%
  dplyr::filter(
    n_matched_human_proteins > 0,
    !is.na(autoantigen_fraction)
  )


if (nrow(fraction_df) >= 3) {
  
  test_fraction <- cor.test(
    fraction_df$ebv_length,
    fraction_df$autoantigen_fraction,
    method = "spearman",
    exact = FALSE
  )
     
} else {
  
  }



# ============================================================
# AUTOANTIGEN MATCHES PER 100 AA
# ============================================================

ebv_size_results <- ebv_size_results %>%
  dplyr::mutate(
    
    autoantigen_matches_per_100aa =
      100 *
      n_matched_autoantigens /
      ebv_length,
    
    human_matches_per_100aa =
      100 *
      n_matched_human_proteins /
      ebv_length
  )


library(ggplot2)
library(dplyr)

# Plot 1:
# EBV length vs number of matched autoantigens
test_auto <- cor.test(
  ebv_size_results$ebv_length,
  ebv_size_results$n_matched_autoantigens,
  method = "spearman",
  exact = FALSE
)

rho_auto <- unname(test_auto$estimate)
p_auto   <- test_auto$p.value
n_auto   <- sum(
  complete.cases(
    ebv_size_results$ebv_length,
    ebv_size_results$n_matched_autoantigens
  )
)


# Plot 2:
# EBV length vs all matched human proteins
test_all <- cor.test(
  ebv_size_results$ebv_length,
  ebv_size_results$n_matched_human_proteins,
  method = "spearman",
  exact = FALSE
)

rho_all <- unname(test_all$estimate)
p_all   <- test_all$p.value
n_all   <- sum(
  complete.cases(
    ebv_size_results$ebv_length,
    ebv_size_results$n_matched_human_proteins
  )
)


# Plot 3:
# Only EBV proteins with at least one human match
fraction_df <- ebv_size_results %>%
  dplyr::filter(
    n_matched_human_proteins > 0,
    !is.na(autoantigen_fraction)
  )

test_fraction <- cor.test(
  fraction_df$ebv_length,
  fraction_df$autoantigen_fraction,
  method = "spearman",
  exact = FALSE
)

rho_fraction <- unname(test_fraction$estimate)
p_fraction   <- test_fraction$p.value
n_fraction   <- nrow(fraction_df)



format_plot_p <- function(p) {
  
  if (p < 2.22e-16) {
    return("P < 2.22 × 10^-16")
  }
  
  if (p < 0.001) {
    return(
      paste0(
        "P = ",
        format(
          p,
          scientific = TRUE,
          digits = 3
        )
      )
    )
  }
  
  paste0(
    "P = ",
    formatC(
      p,
      format = "f",
      digits = 3
    )
  )
}


library(ggplot2)
library(dplyr)


# ============================================================
# PLOT 1
# ============================================================

p1 <- ggplot(
  ebv_size_results,
  aes(
    x = ebv_length,
    y = n_matched_autoantigens
  )
) +
  geom_point(
    size = 3,
    alpha = 0.8
  ) +
  geom_smooth(
    method = "lm",
    se = TRUE
  ) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = paste0(
      "Spearman rho = ",
      sprintf("%.3f", unname(test_auto$estimate)),
      "\nP = ",
      format(
        test_auto$p.value,
        scientific = TRUE,
        digits = 3
      ),
      "\nn = ",
      nrow(ebv_size_results)
    ),
    hjust = 1.1,
    vjust = 1.2,
    size = 4.5
  ) +
  labs(
    x = "EBV protein length (amino acids)",
    y = "Number of peptide-matched human autoantigens",
    title = "EBV protein size and matched autoantigen burden"
  ) +
  theme_classic(
    base_size = 13
  ) +
  theme(
    plot.title = element_text(
      face = "bold"
    )
  )

# ============================================================
# PLOT 2
# ============================================================

p2 <- ggplot(
  ebv_size_results,
  aes(
    x = ebv_length,
    y = n_matched_human_proteins
  )
) +
  geom_point(
    size = 3,
    alpha = 0.8
  ) +
  geom_smooth(
    method = "lm",
    se = TRUE
  ) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = paste0(
      "Spearman rho = ",
      sprintf("%.3f", unname(test_all$estimate)),
      "\nP < 2.22e-16",
      "\nn = ",
      nrow(ebv_size_results)
    ),
    hjust = 1.1,
    vjust = 1.2,
    size = 4.5
  ) +
  labs(
    x = "EBV protein length (amino acids)",
    y = "Number of peptide-matched human proteins",
    title = "EBV protein size and total human peptide sharing"
  ) +
  theme_classic(
    base_size = 13
  ) +
  theme(
    plot.title = element_text(
      face = "bold"
    )
  )

# ============================================================
# PLOT 3
# ============================================================

fraction_df <- ebv_size_results %>%
  dplyr::filter(
    n_matched_human_proteins > 0,
    !is.na(autoantigen_fraction)
  )

p3 <- ggplot(
  fraction_df,
  aes(
    x = ebv_length,
    y = autoantigen_fraction
  )
) +
  geom_point(
    size = 3,
    alpha = 0.8
  ) +
  geom_smooth(
    method = "lm",
    se = TRUE
  ) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = paste0(
      "Spearman rho = ",
      sprintf("%.3f", unname(test_fraction$estimate)),
      "\nP = ",
      sprintf("%.3f", test_fraction$p.value),
      "\nn = ",
      nrow(fraction_df)
    ),
    hjust = 1.1,
    vjust = 1.2,
    size = 4.5
  ) +
  labs(
    x = "EBV protein length (amino acids)",
    y = "Fraction of matched human proteins that are autoantigens",
    title = "EBV protein size and autoantigen enrichment among matches"
  ) +
  theme_classic(
    base_size = 13
  ) +
  theme(
    plot.title = element_text(
      face = "bold"
    )
  )


# ============================================================
# EBV HOST-INTERACTION BURDEN vs PEPTIDE-MATCHED
#     HUMAN AUTOANTIGENS
# ============================================================

library(dplyr)
library(ggplot2)
library(tibble)


ppi_ebv <- ppis_std %>%
  dplyr::transmute(
    
    ebv_name = sub(
      "_.*$",
      "",
      trimws(
        as.character(Viral_Label)
      )
    ),
    
    human_gene = trimws(
      as.character(Gene)
    )
  ) %>%
  dplyr::filter(
    !is.na(ebv_name),
    ebv_name != "",
    !is.na(human_gene),
    human_gene != ""
  )


ebv_host_burden_name <- ppi_ebv %>%
  dplyr::distinct(
    ebv_name,
    human_gene
  ) %>%
  dplyr::count(
    ebv_name,
    name = "n_human_targets"
  )


map_df <- as.data.frame(
  ebv_accession_map,
  stringsAsFactors = FALSE
)

map_names <- names(map_df)


accession_candidates <- map_names[
  grepl(
    "uniprot|accession|protein.*id|^id$",
    map_names,
    ignore.case = TRUE
  )
]


name_candidates <- map_names[
  grepl(
    "name|label|viral|protein|gene|symbol",
    map_names,
    ignore.case = TRUE
  )
]


sequence_ids <- unique(
  trimws(
    as.character(
      ebv_size_results$ebv_protein
    )
  )
)


accession_overlap <- sapply(
  map_df,
  function(x) {
    
    values <- trimws(
      as.character(x)
    )
    
    sum(
      sequence_ids %in% values,
      na.rm = TRUE
    )
  }
)


accession_col <- names(
  accession_overlap
)[
  which.max(
    accession_overlap
  )
]


if (
  max(
    accession_overlap,
    na.rm = TRUE
  ) == 0
) {
  
  stop(
    paste0(
      "No column in ebv_accession_map matches ",
      "ebv_size_results$ebv_protein. ",
      "Inspect the printed accession map."
    )
  )
}


ppi_names <- unique(
  toupper(
    trimws(
      ebv_host_burden_name$ebv_name
    )
  )
)


name_overlap <- sapply(
  map_df,
  function(x) {
    
    values <- toupper(
      trimws(
        as.character(x)
      )
    )
    
    sum(
      values %in% ppi_names,
      na.rm = TRUE
    )
  }
)


name_overlap[
  accession_col
] <- -1




name_col <- names(
  name_overlap
)[
  which.max(
    name_overlap
  )
]


if (
  max(
    name_overlap,
    na.rm = TRUE
  ) <= 0
) {
  
  stop(
    paste0(
      "Could identify the UniProt accession column, ",
      "but could not identify an EBV protein-name column ",
      "matching ppis_std$Viral_Label. ",
      "Inspect the printed map."
    )
  )
}


ebv_map_clean <- data.frame(
  
  ebv_protein = trimws(
    as.character(
      map_df[[accession_col]]
    )
  ),
  
  ebv_name = trimws(
    as.character(
      map_df[[name_col]]
    )
  ),
  
  stringsAsFactors = FALSE
  
) %>%
  
  dplyr::mutate(
    
    ebv_name = sub(
      "_.*$",
      "",
      ebv_name
    ),
    
    ebv_name_upper = toupper(
      ebv_name
    )
    
  ) %>%
  
  dplyr::filter(
    !is.na(ebv_protein),
    ebv_protein != "",
    !is.na(ebv_name),
    ebv_name != ""
  ) %>%
  
  dplyr::distinct(
    ebv_protein,
    ebv_name_upper
  )


# Standardise PPI names to uppercase
ebv_host_burden_name2 <- ebv_host_burden_name %>%
  
  dplyr::mutate(
    
    ebv_name_upper = toupper(
      trimws(
        ebv_name
      )
    )
    
  )


ebv_host_burden_accession <- ebv_map_clean %>%
  
  dplyr::left_join(
    ebv_host_burden_name2 %>%
      dplyr::select(
        ebv_name_upper,
        n_human_targets
      ),
    by = "ebv_name_upper"
  )


mapping_check <- ebv_host_burden_accession %>%
  
  dplyr::group_by(
    ebv_protein
  ) %>%
  
  dplyr::summarise(
    
    n_names = dplyr::n_distinct(
      ebv_name_upper
    ),
    
    n_target_values = dplyr::n_distinct(
      n_human_targets[
        !is.na(n_human_targets)
      ]
    ),
    
    .groups = "drop"
    
  )


ebv_host_burden_accession <- ebv_host_burden_accession %>%
  
  dplyr::group_by(
    ebv_protein
  ) %>%
  
  dplyr::summarise(
    
    mapped_ebv_names = paste(
      sort(
        unique(
          ebv_name_upper
        )
      ),
      collapse = ";"
    ),
    
    n_human_targets = if (
      all(
        is.na(
          n_human_targets
        )
      )
    ) {
      
      NA_integer_
      
    } else {
      
      max(
        n_human_targets,
        na.rm = TRUE
      )
      
    },
    
    .groups = "drop"
  )



ebv_biology_results <- ebv_size_results %>%
  
  dplyr::left_join(
    ebv_host_burden_accession,
    by = "ebv_protein"
  )



ebv_target_analysis <- ebv_biology_results %>%
  
  dplyr::filter(
    !is.na(
      n_human_targets
    )
  )



if (
  nrow(
    ebv_target_analysis
  ) < 3
) {
  
  stop(
    paste0(
      "Too few EBV proteins were successfully mapped ",
      "for correlation analysis. Do not interpret the correlations."
    )
  )
}



# ============================================================
# PRIMARY TEST: Human PPI target burden
# vs peptide-matched human autoantigens
# ============================================================

test_target_auto <- cor.test(
  
  ebv_target_analysis$n_human_targets,
  
  ebv_target_analysis$n_matched_autoantigens,
  
  method = "spearman",
  
  exact = FALSE
)
# ============================================================
# CONTROL TEST
# Human PPI target burden
# vs ALL peptide-matched human proteins
# ============================================================

test_target_all <- cor.test(
  
  ebv_target_analysis$n_human_targets,
  
  ebv_target_analysis$n_matched_human_proteins,
  
  method = "spearman",
  
  exact = FALSE
)
# ============================================================
# SPECIFICITY TEST
# Does greater host-interaction burden correspond to a
# greater FRACTION of peptide matches being autoantigens?
# ============================================================

target_fraction_df <- ebv_target_analysis %>%
  
  dplyr::filter(
    n_matched_human_proteins > 0,
    !is.na(
      autoantigen_fraction
    )
  )


if (
  nrow(
    target_fraction_df
  ) >= 3
) {
  
  test_target_fraction <- cor.test(
    
    target_fraction_df$n_human_targets,
    
    target_fraction_df$autoantigen_fraction,
    
    method = "spearman",
    
    exact = FALSE
  )
  
  
          
    
    
    
} else {
  
  stop(
    "Too few mapped EBV proteins with peptide matches."
  )
}



# ============================================================
#  TEST
# ============================================================

test_target_length <- cor.test(
  
  ebv_target_analysis$n_human_targets,
  
  ebv_target_analysis$ebv_length,
  
  method = "spearman",
  
  exact = FALSE
)
# ============================================================
# PLOT 1
# HOST TARGET BURDEN vs MATCHED AUTOANTIGENS
# ============================================================

p4 <- ggplot(
  ebv_target_analysis,
  aes(
    x = n_human_targets,
    y = n_matched_autoantigens
  )
) +
  
  geom_point(
    size = 3,
    alpha = 0.8
  ) +
  
  geom_smooth(
    method = "lm",
    se = TRUE
  ) +
  
  annotate(
    "text",
    x = Inf,
    y = Inf,
    
    label = paste0(
      
      "Spearman rho = ",
      sprintf(
        "%.3f",
        unname(
          test_target_auto$estimate
        )
      ),
      
      "\nP = ",
      format.pval(
        test_target_auto$p.value,
        digits = 3
      ),
      
      "\nn = ",
      nrow(
        ebv_target_analysis
      )
    ),
    
    hjust = 1.1,
    vjust = 1.2,
    size = 4.5
  ) +
  
  labs(
    x = "Number of human proteins directly targeted by EBV protein",
    y = "Number of peptide-matched human autoantigens",
    title = "EBV host-interaction burden and matched autoantigens"
  ) +
  
  theme_classic(
    base_size = 13
  ) +
  
  theme(
    plot.title = element_text(
      face = "bold"
    )
  )




# ============================================================
# PLOT 2
# HOST TARGET BURDEN vs ALL HUMAN PEPTIDE MATCHES
# ============================================================

p5 <- ggplot(
  ebv_target_analysis,
  aes(
    x = n_human_targets,
    y = n_matched_human_proteins
  )
) +
  
  geom_point(
    size = 3,
    alpha = 0.8
  ) +
  
  geom_smooth(
    method = "lm",
    se = TRUE
  ) +
  
  annotate(
    "text",
    x = Inf,
    y = Inf,
    
    label = paste0(
      
      "Spearman rho = ",
      sprintf(
        "%.3f",
        unname(
          test_target_all$estimate
        )
      ),
      
      "\nP = ",
      format.pval(
        test_target_all$p.value,
        digits = 3
      ),
      
      "\nn = ",
      nrow(
        ebv_target_analysis
      )
    ),
    
    hjust = 1.1,
    vjust = 1.2,
    size = 4.5
  ) +
  
  labs(
    x = "Number of human proteins directly targeted by EBV protein",
    y = "Number of peptide-matched human proteins",
    title = "EBV host-interaction burden and human peptide sharing"
  ) +
  
  theme_classic(
    base_size = 13
  ) +
  
  theme(
    plot.title = element_text(
      face = "bold"
    )
  )

# ============================================================
# PLOT 3
# HOST TARGET BURDEN vs AUTOANTIGEN FRACTION
# ============================================================

p6 <- ggplot(
  target_fraction_df,
  aes(
    x = n_human_targets,
    y = autoantigen_fraction
  )
) +
  
  geom_point(
    size = 3,
    alpha = 0.8
  ) +
  
  geom_smooth(
    method = "lm",
    se = TRUE
  ) +
  
  annotate(
    "text",
    x = Inf,
    y = Inf,
    
    label = paste0(
      
      "Spearman rho = ",
      sprintf(
        "%.3f",
        unname(
          test_target_fraction$estimate
        )
      ),
      
      "\nP = ",
      format.pval(
        test_target_fraction$p.value,
        digits = 3
      ),
      
      "\nn = ",
      nrow(
        target_fraction_df
      )
    ),
    
    hjust = 1.1,
    vjust = 1.2,
    size = 4.5
  ) +
  
  labs(
    x = "Number of human proteins directly targeted by EBV protein",
    y = "Fraction of peptide matches that are autoantigens",
    title = "EBV host-interaction burden and autoantigen enrichment"
  ) +
  
  theme_classic(
    base_size = 13
  ) +
  
  theme(
    plot.title = element_text(
      face = "bold"
    )
  )
# ============================================================
# PLOT 4
# HOST TARGET BURDEN vs EBV PROTEIN LENGTH

# ============================================================

p7 <- ggplot(
  ebv_target_analysis,
  aes(
    x = ebv_length,
    y = n_human_targets
  )
) +
  
  geom_point(
    size = 3,
    alpha = 0.8
  ) +
  
  geom_smooth(
    method = "lm",
    se = TRUE
  ) +
  
  annotate(
    "text",
    x = Inf,
    y = Inf,
    
    label = paste0(
      
      "Spearman rho = ",
      sprintf(
        "%.3f",
        unname(
          test_target_length$estimate
        )
      ),
      
      "\nP = ",
      format.pval(
        test_target_length$p.value,
        digits = 3
      ),
      
      "\nn = ",
      nrow(
        ebv_target_analysis
      )
    ),
    
    hjust = 1.1,
    vjust = 1.2,
    size = 4.5
  ) +
  
  labs(
    x = "EBV protein length (amino acids)",
    y = "Number of directly targeted human proteins",
    title = "EBV protein size and host-interaction burden"
  ) +
  
  theme_classic(
    base_size = 13
  ) +
  
  theme(
    plot.title = element_text(
      face = "bold"
    )
  )



# ============================================================
# 15. UNIQUE EBV PROTEIN-LEVEL ANALYSIS
# ============================================================

library(dplyr)
library(ggplot2)
library(tibble)

ebv_variant_df <- ebv_biology_results %>%
  
  dplyr::mutate(
    
    ebv_identity = trimws(
      as.character(
        mapped_ebv_names
      )
    )
    
  ) %>%
  
  dplyr::filter(
    !is.na(ebv_identity),
    ebv_identity != "",
    !is.na(n_human_targets)
  )


variant_counts <- ebv_variant_df %>%
  
  dplyr::count(
    ebv_identity,
    name = "n_sequence_variants"
  ) %>%
  
  dplyr::arrange(
    dplyr::desc(n_sequence_variants)
  )


ebv_protein_level <- ebv_variant_df %>%
  
  dplyr::group_by(
    ebv_identity
  ) %>%
  
  dplyr::summarise(
    
    n_sequence_variants = dplyr::n(),
    
    accessions = paste(
      sort(
        unique(
          ebv_protein
        )
      ),
      collapse = ";"
    ),
    
    # --------------------------------------------------------
    # Protein length
    # --------------------------------------------------------
    
    ebv_length_median = median(
      ebv_length,
      na.rm = TRUE
    ),
    
    ebv_length_min = min(
      ebv_length,
      na.rm = TRUE
    ),
    
    ebv_length_max = max(
      ebv_length,
      na.rm = TRUE
    ),
    

    
    n_human_targets = max(
      n_human_targets,
      na.rm = TRUE
    ),
    
  
    
    n_matched_human_max = max(
      n_matched_human_proteins,
      na.rm = TRUE
    ),
    
    n_matched_human_median = median(
      n_matched_human_proteins,
      na.rm = TRUE
    ),
    
  
    
    n_matched_autoantigens_max = max(
      n_matched_autoantigens,
      na.rm = TRUE
    ),
    
    n_matched_autoantigens_median = median(
      n_matched_autoantigens,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  ) %>%
  
  dplyr::mutate(
    
   
    autoantigen_fraction_max = dplyr::if_else(
      
      n_matched_human_max > 0,
      
      n_matched_autoantigens_max /
        n_matched_human_max,
      
      NA_real_
    ),
    
    
    
    autoantigen_fraction_median = dplyr::if_else(
      
      n_matched_human_median > 0,
      
      n_matched_autoantigens_median /
        n_matched_human_median,
      
      NA_real_
    )
  )




ppi_consistency <- ebv_variant_df %>%
  
  dplyr::group_by(
    ebv_identity
  ) %>%
  
  dplyr::summarise(
    
    n_ppi_values = dplyr::n_distinct(
      n_human_targets
    ),
    
    min_targets = min(
      n_human_targets,
      na.rm = TRUE
    ),
    
    max_targets = max(
      n_human_targets,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )




if (
  any(
    ppi_consistency$n_ppi_values > 1
  )
) {
  
    
}



# ============================================================
# PRIMARY CORRELATION
# HOST-INTERACTION BURDEN
# vs
# PEPTIDE-MATCHED AUTOANTIGENS
# ============================================================

test_collapsed_auto <- cor.test(
  
  ebv_protein_level$n_human_targets,
  
  ebv_protein_level$n_matched_autoantigens_max,
  
  method = "spearman",
  
  exact = FALSE
)

# ============================================================
# CONTROL CORRELATION
# HOST-INTERACTION BURDEN
# vs
# ALL PEPTIDE-MATCHED HUMAN PROTEINS
# ============================================================

test_collapsed_all <- cor.test(
  
  ebv_protein_level$n_human_targets,
  
  ebv_protein_level$n_matched_human_max,
  
  method = "spearman",
  
  exact = FALSE
)

# ============================================================
# 15G. SPECIFICITY CORRELATION
# HOST-INTERACTION BURDEN
# vs
# FRACTION OF MATCHES THAT ARE AUTOANTIGENS
# ============================================================

ebv_fraction_level <- ebv_protein_level %>%
  
  dplyr::filter(
    n_matched_human_max > 0,
    !is.na(autoantigen_fraction_max)
  )


test_collapsed_fraction <- cor.test(
  
  ebv_fraction_level$n_human_targets,
  
  ebv_fraction_level$autoantigen_fraction_max,
  
  method = "spearman",
  
  exact = FALSE
)

# ============================================================
# LENGTH vs MATCHED AUTOANTIGENS
# ============================================================

test_collapsed_length_auto <- cor.test(
  
  ebv_protein_level$ebv_length_median,
  
  ebv_protein_level$n_matched_autoantigens_max,
  
  method = "spearman",
  
  exact = FALSE
)

# ============================================================
# LENGTH vs ALL HUMAN MATCHES
# ============================================================

test_collapsed_length_all <- cor.test(
  
  ebv_protein_level$ebv_length_median,
  
  ebv_protein_level$n_matched_human_max,
  
  method = "spearman",
  
  exact = FALSE
)

# ============================================================
# LENGTH vs AUTOANTIGEN FRACTION
# ============================================================

test_collapsed_length_fraction <- cor.test(
  
  ebv_fraction_level$ebv_length_median,
  
  ebv_fraction_level$autoantigen_fraction_max,
  
  method = "spearman",
  
  exact = FALSE
)

# ============================================================
# LENGTH vs HOST-INTERACTION BURDEN
# ============================================================

test_collapsed_length_targets <- cor.test(
  
  ebv_protein_level$ebv_length_median,
  
  ebv_protein_level$n_human_targets,
  
  method = "spearman",
  
  exact = FALSE
)

# ============================================================
# PARTIAL SPEARMAN CORRELATION
# ============================================================

partial_auto_df <- ebv_protein_level %>%
  
  dplyr::filter(
    !is.na(n_human_targets),
    !is.na(n_matched_autoantigens_max),
    !is.na(ebv_length_median)
  ) %>%
  
  dplyr::mutate(
    
    rank_targets = rank(
      n_human_targets,
      ties.method = "average"
    ),
    
    rank_auto = rank(
      n_matched_autoantigens_max,
      ties.method = "average"
    ),
    
    rank_length = rank(
      ebv_length_median,
      ties.method = "average"
    )
  )


model_targets_length <- lm(
  rank_targets ~ rank_length,
  data = partial_auto_df
)

model_auto_length <- lm(
  rank_auto ~ rank_length,
  data = partial_auto_df
)


partial_auto_df <- partial_auto_df %>%
  
  dplyr::mutate(
    
    residual_targets = residuals(
      model_targets_length
    ),
    
    residual_auto = residuals(
      model_auto_length
    )
  )


partial_auto_test <- cor.test(
  
  partial_auto_df$residual_targets,
  
  partial_auto_df$residual_auto,
  
  method = "pearson"
)


# ============================================================
#  EBNA2 / BNLF2A 
# ============================================================

contrast_df <- ebv_protein_level %>%
  
  dplyr::filter(
    ebv_identity %in% c(
      "EBNA2",
      "BNLF2A"
    )
  )



plot_label <- paste0(
  
  "Spearman rho = ",
  sprintf(
    "%.3f",
    unname(
      test_collapsed_auto$estimate
    )
  ),
  
  "\nP = ",
  format.pval(
    test_collapsed_auto$p.value,
    digits = 3
  ),
  
  "\nn = ",
  nrow(
    ebv_protein_level
  )
)


p_host_auto <- ggplot(
  
  ebv_protein_level,
  
  aes(
    x = n_human_targets,
    y = n_matched_autoantigens_max
  )
  
) +
  
  geom_point(
    size = 3,
    alpha = 0.75
  ) +
  
  geom_smooth(
    method = "lm",
    se = TRUE
  ) +
  
  geom_point(
    data = contrast_df,
    size = 4
  ) +
  
  geom_text(
    data = contrast_df,
    aes(
      label = ebv_identity
    ),
    nudge_x = 0.15,
    nudge_y = 0.15,
    check_overlap = FALSE,
    fontface = "bold"
  ) +
  
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = plot_label,
    hjust = 1.1,
    vjust = 1.2,
    size = 4.3
  ) +
  
  scale_x_continuous(
    trans = "log1p"
  ) +
  
  scale_y_continuous(
    trans = "log1p"
  ) +
  
  labs(
    x = "Number of directly targeted human proteins",
    y = "Number of peptide-matched human autoantigens",
    title = "EBV host-interaction burden and autoantigen peptide sharing"
  ) +
  
  theme_classic(
    base_size = 13
  ) +
  
  theme(
    plot.title = element_text(
      face = "bold"
    )
  )



