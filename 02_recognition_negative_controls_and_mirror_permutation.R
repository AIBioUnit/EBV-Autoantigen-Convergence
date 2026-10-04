################################################################################
## RECOGNIZED vs NON-RECOGNIZED EBV PROTEINS
################################################################################

library(dplyr)
library(tidyr)
library(stringr)

# ============================================================================
# 1. COLLAPSE EBV PROTEINS 
# ============================================================================

ppi_core <- ppis_std %>%
  dplyr::transmute(
    Gene = toupper(trimws(as.character(Gene))),
    Viral_Label = trimws(as.character(Viral_Label)),
    Viral_Core = sub("_.*$", "", trimws(as.character(Viral_Label)))
  ) %>%
  dplyr::filter(
    !is.na(Gene), Gene != "",
    !is.na(Viral_Label), Viral_Label != "",
    !is.na(Viral_Core), Viral_Core != ""
  ) %>%
  dplyr::distinct()


all_ebv_core <- sort(unique(ppi_core$Viral_Core))




# ============================================================================
# 2. GET RECOGNIZED EBV PROTEINS 
# ============================================================================

recognized_raw <- viral_sigpath_summary_all %>%
  dplyr::filter(
    !is.na(recognized_viral_proteins),
    recognized_viral_proteins != ""
  ) %>%
  dplyr::pull(recognized_viral_proteins) %>%
  paste(collapse = "; ") %>%
  stringr::str_split(";\\s*") %>%
  unlist() %>%
  trimws()

recognized_raw <- unique(
  recognized_raw[
    !is.na(recognized_raw) &
      recognized_raw != ""
  ]
)

# Apply EXACT SAME collapse
recognized_core <- sort(unique(
  sub("_.*$", "", recognized_raw)
))

# Keep only proteins actually represented in PPI universe
recognized_core <- intersect(
  recognized_core,
  all_ebv_core
)

nonrecognized_core <- setdiff(
  all_ebv_core,
  recognized_core
)



# ============================================================================
# 3. ASSIGN RECOGNITION STATUS
# ============================================================================

ppi_core <- ppi_core %>%
  dplyr::mutate(
    recognition_status = ifelse(
      Viral_Core %in% recognized_core,
      "Recognized",
      "Non_recognized"
    )
  )



print(
  data.frame(
    Status = c("Recognized", "Non-recognized"),
    N = c(
      length(recognized_core),
      length(nonrecognized_core)
    )
  )
)


# ============================================================================
# 4. MAP COLLAPSED EBV PROTEINS INTO THE 87 PATHWAYS
# ============================================================================

pathway_protein_status <- path_genes_long %>%
  dplyr::transmute(
    pathway_id = pathway_id,
    pathway_name = pathway_name,
    Gene = toupper(trimws(as.character(Gene)))
  ) %>%
  dplyr::distinct() %>%
  dplyr::inner_join(
    ppi_core %>%
      dplyr::select(
        Gene,
        Viral_Core,
        recognition_status
      ) %>%
      dplyr::distinct(),
    by = "Gene",
    relationship = "many-to-many"
  ) %>%
  dplyr::distinct(
    pathway_id,
    pathway_name,
    Viral_Core,
    recognition_status
  )


# ============================================================================
# 5. COUNT RECOGNIZED AND NON-RECOGNIZED PROTEINS PER PATHWAY
# ============================================================================

pathway_status <- pathway_protein_status %>%
  dplyr::group_by(
    pathway_id,
    pathway_name,
    recognition_status
  ) %>%
  dplyr::summarise(
    n_proteins = dplyr::n_distinct(Viral_Core),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(
    names_from = recognition_status,
    values_from = n_proteins,
    values_fill = 0
  )


# ============================================================================
# 6. MERGE WITH 87-PATHWAY CONVERGENCE TABLE
# ============================================================================

negative_control_core <- integrated_pathway_table_v2 %>%
  dplyr::select(
    pathway_id,
    pathway_name,
    n_shared_in_path,
    n_ebv_in_path,
    odds_ratio
  ) %>%
  dplyr::left_join(
    pathway_status,
    by = c("pathway_id", "pathway_name")
  ) %>%
  dplyr::mutate(
    Recognized =
      tidyr::replace_na(Recognized, 0),
    
    Non_recognized =
      tidyr::replace_na(Non_recognized, 0)
  )



# ============================================================================
# 7. NEGATIVE-CONTROL CORRELATIONS
# ============================================================================

cor_shared_rec <- cor.test(
  negative_control_core$n_shared_in_path,
  negative_control_core$Recognized,
  method = "spearman",
  exact = FALSE
)

cor_shared_nonrec <- cor.test(
  negative_control_core$n_shared_in_path,
  negative_control_core$Non_recognized,
  method = "spearman",
  exact = FALSE
)

cor_or_rec <- cor.test(
  negative_control_core$odds_ratio,
  negative_control_core$Recognized,
  method = "spearman",
  exact = FALSE
)

cor_or_nonrec <- cor.test(
  negative_control_core$odds_ratio,
  negative_control_core$Non_recognized,
  method = "spearman",
  exact = FALSE
)



################################################################################
## MERGE WITH PATHWAY TABLE
################################################################################

negative_control_core <- integrated_pathway_table %>%
  dplyr::select(
    pathway_id,
    pathway_name,
    n_shared_in_path,
    n_ebv_in_path,
    odds_ratio
  ) %>%
  dplyr::left_join(
    pathway_status,
    by = c("pathway_id", "pathway_name")
  ) %>%
  dplyr::mutate(
    Recognized = tidyr::replace_na(Recognized, 0),
    Non_recognized = tidyr::replace_na(Non_recognized, 0)
  )


################################################################################
## CORRELATIONS
################################################################################

cor_shared_rec <- cor.test(
  negative_control_core$n_shared_in_path,
  negative_control_core$Recognized,
  method = "spearman",
  exact = FALSE
)

cor_shared_nonrec <- cor.test(
  negative_control_core$n_shared_in_path,
  negative_control_core$Non_recognized,
  method = "spearman",
  exact = FALSE
)

cor_or_rec <- cor.test(
  negative_control_core$odds_ratio,
  negative_control_core$Recognized,
  method = "spearman",
  exact = FALSE
)

cor_or_nonrec <- cor.test(
  negative_control_core$odds_ratio,
  negative_control_core$Non_recognized,
  method = "spearman",
  exact = FALSE
)



################################################################################
## MIRROR NULL AUTOANTIGEN PERMUTATION
################################################################################

library(dplyr)
library(tidyr)
library(purrr)



baseline_ok <- (
  nrow(real_conv) == 87 &&
    real_summary$n_sig == 5 &&
    abs(real_summary$max_or - 5.940241) < 0.01 &&
    abs(real_summary$mean_shared - 8.770115) < 0.01 &&
    real_summary$max_shared == 39
)

if (!baseline_ok) {
  stop("Validated baseline has changed. Do not run permutations.")
}

print(real_summary, row.names = FALSE)




gene_pathway_counts <- path_genes_long_clean %>%
  dplyr::distinct(
    Gene,
    pathway_id
  ) %>%
  dplyr::count(
    Gene,
    name = "n_pathways"
  )


gene_strata_mirror <- data.frame(
  Gene = universe_genes,
  stringsAsFactors = FALSE
) %>%
  dplyr::left_join(
    gene_pathway_counts,
    by = "Gene"
  ) %>%
  dplyr::mutate(
    n_pathways = tidyr::replace_na(
      n_pathways,
      0L
    )
  )



gene_strata_mirror <- gene_strata_mirror %>%
  dplyr::mutate(
    stratum = as.character(
      dplyr::ntile(
        n_pathways,
        5
      )
    )
  )

print(
  table(
    gene_strata_mirror$stratum
  )
)


obs_strata_epi <- table(
  gene_strata_mirror$stratum[
    gene_strata_mirror$Gene %in% epi_genes
  ]
)

print(obs_strata_epi)




if (sum(obs_strata_epi) != length(epi_genes)) {
  stop("Not all autoantigens are represented in the matching strata.")
}



sample_matched_autoantigens <- function() {
  
  sampled <- lapply(
    names(obs_strata_epi),
    
    function(s) {
      
      pool <- gene_strata_mirror$Gene[
        gene_strata_mirror$stratum == s
      ]
      
      k <- as.integer(
        obs_strata_epi[[s]]
      )
      
      if (k == 0) {
        return(character(0))
      }
      
      if (length(pool) < k) {
        stop(
          paste(
            "Insufficient genes in stratum",
            s
          )
        )
      }
      
      sample(
        pool,
        size = k,
        replace = FALSE
      )
    }
  )
  
  unlist(
    sampled,
    use.names = FALSE
  )
}



set.seed(2027)

test_draw <- sample_matched_autoantigens()






if (
  length(test_draw) != length(epi_genes) ||
  length(unique(test_draw)) != length(epi_genes)
) {
  stop("Matched autoantigen sampler failed.")
}



empirical_p_mirror <- function(
    perm_values,
    observed_value) {
  
  perm_values <- perm_values[
    !is.na(perm_values) &
      is.finite(perm_values)
  ]
  
  r <- sum(
    perm_values >= observed_value
  )
  
  (r + 1) /
    (length(perm_values) + 1)
}



set.seed(2027)

n_perm_mirror <- 1000



perm_results_mirror <- purrr::map_dfr(
  
  seq_len(n_perm_mirror),
  
  function(i) {
    
    random_epi <- sample_matched_autoantigens()
    
    perm_conv <- run_pathway_convergence(
      path_genes_long_clean = path_genes_long_clean,
      epi_genes             = random_epi,
      ebv_genes             = ebv_genes,
      universe_genes        = universe_genes
    )
    
    finite_or <- perm_conv$odds_ratio[
      is.finite(
        perm_conv$odds_ratio
      )
    ]
    
    data.frame(
      
      iteration = i,
      
      n_sig = sum(
        perm_conv$fisher_fdr < 0.05,
        na.rm = TRUE
      ),
      
      max_or = max(
        finite_or,
        na.rm = TRUE
      ),
      
      mean_shared = mean(
        perm_conv$n_shared_in_path,
        na.rm = TRUE
      ),
      
      max_shared = max(
        perm_conv$n_shared_in_path,
        na.rm = TRUE
      )
    )
  }
)



mirror_results <- data.frame(
  
  metric = c(
    "Number significant pathways",
    "Maximum odds ratio",
    "Mean shared genes",
    "Maximum shared genes"
  ),
  
  observed = c(
    real_summary$n_sig,
    real_summary$max_or,
    real_summary$mean_shared,
    real_summary$max_shared
  ),
  
  null_mean = c(
    
    mean(
      perm_results_mirror$n_sig,
      na.rm = TRUE
    ),
    
    mean(
      perm_results_mirror$max_or,
      na.rm = TRUE
    ),
    
    mean(
      perm_results_mirror$mean_shared,
      na.rm = TRUE
    ),
    
    mean(
      perm_results_mirror$max_shared,
      na.rm = TRUE
    )
  ),
  
  empirical_p = c(
    
    empirical_p_mirror(
      perm_results_mirror$n_sig,
      real_summary$n_sig
    ),
    
    empirical_p_mirror(
      perm_results_mirror$max_or,
      real_summary$max_or
    ),
    
    empirical_p_mirror(
      perm_results_mirror$mean_shared,
      real_summary$mean_shared
    ),
    
    empirical_p_mirror(
      perm_results_mirror$max_shared,
      real_summary$max_shared
    )
  )
)




print(
  mirror_results,
  row.names = FALSE
)



print(
  real_summary,
  row.names = FALSE
)






