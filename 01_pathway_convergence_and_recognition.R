################################################################################
## EBV–autoantigen pathway convergence and immune-recognition analyses
################################################################################

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)

integrated_pathway_table_v2 <- integrated_pathway_table %>%
  tidyr::replace_na(list(
    n_recognized_ebv_proteins = 0,
    recognized_epitope_links  = 0,
    n_auto_epitopes           = 0,
    n_auto_genes              = 0,
    auto_gene_coverage_pct    = 0,
    n_ebv_host_genes          = 0,
    n_total_ebv_proteins      = 0
  ))




size_cols_supp <- grep("^pathway_size", names(integrated_pathway_table_v2), value = TRUE)

if (length(size_cols_supp) == 0) {
  ## fall back to recomputing from the pathway gene table
  integrated_pathway_table_v2 <- integrated_pathway_table_v2 %>%
    left_join(
      path_genes_long_clean %>%
        group_by(pathway_id) %>%
        summarise(pathway_size = n_distinct(Gene), .groups = "drop"),
      by = "pathway_id"
    )
} else if (!"pathway_size" %in% names(integrated_pathway_table_v2)) {
  integrated_pathway_table_v2$pathway_size <-
    dplyr::coalesce(!!!integrated_pathway_table_v2[size_cols_supp])
}

## sanity check: the two versions should agree
if (length(size_cols_supp) > 1) {
  disagree <- sum(integrated_pathway_table_v2[[size_cols_supp[1]]] !=
                    integrated_pathway_table_v2[[size_cols_supp[2]]], na.rm = TRUE)
  }

stopifnot(!any(is.na(integrated_pathway_table_v2$pathway_size)))


## helper: run one Spearman test and return a tidy row
spearman_row_supp <- function(x, y, label, data) {
  ok <- complete.cases(data[[x]], data[[y]])
  ct <- suppressWarnings(cor.test(data[[x]][ok], data[[y]][ok],
                                  method = "spearman", exact = FALSE))
  data.frame(
    comparison = label,
    n          = sum(ok),
    rho        = unname(ct$estimate),
    p          = ct$p.value,
    stringsAsFactors = FALSE
  )
}

## shared_pct is computed in plot_df in your script; recompute here safely
integrated_pathway_table_v2 <- integrated_pathway_table_v2 %>%
  mutate(shared_pct_v2 = 100 * n_shared_in_path / pathway_size)

correlations_v2_supp <- bind_rows(
  spearman_row_supp(
    "n_shared_in_path", "n_recognized_ebv_proteins",
    "Shared EBV-autoantigen genes vs EBV recognition",
    integrated_pathway_table_v2
  ),
  spearman_row_supp(
    "shared_pct_v2", "n_recognized_ebv_proteins",
    "Shared genes (% of pathway) vs EBV recognition",
    integrated_pathway_table_v2
  )
) %>%
  mutate(p_BH = p.adjust(p, method = "BH"))

################################################################################
## ANNOTATION-DEPTH MATCHED EBV-TARGET PERMUTATION
################################################################################

gene_pathway_count_supp <- path_genes_long_clean %>%
  dplyr::distinct(Gene, pathway_id) %>%
  dplyr::count(Gene, name = "n_pathways")

gene_strata_supp <- data.frame(
  Gene = universe_genes,
  stringsAsFactors = FALSE
) %>%
  dplyr::left_join(gene_pathway_count_supp, by = "Gene") %>%
  dplyr::mutate(
    n_pathways = tidyr::replace_na(n_pathways, 0L),
    stratum = as.character(dplyr::ntile(n_pathways, 5))
  )

obs_strata_supp <- table(
  gene_strata_supp$stratum[gene_strata_supp$Gene %in% ebv_genes]
)

sample_matched_supp <- function() {
  unlist(
    lapply(names(obs_strata_supp), function(s) {
      pool <- gene_strata_supp$Gene[gene_strata_supp$stratum == s]
      k <- as.integer(obs_strata_supp[[s]])
      sample(pool, size = k, replace = FALSE)
    }),
    use.names = FALSE
  )
}

emp_p_supp <- function(perm_vals, obs) {
  perm_vals <- perm_vals[!is.na(perm_vals) & is.finite(perm_vals)]
  (sum(perm_vals >= obs) + 1) / (length(perm_vals) + 1)
}

set.seed(2026)
n_perm_supp <- 1000

perm_results_matched_supp <- purrr::map_dfr(
  seq_len(n_perm_supp),
  function(i) {
    rand_ebv_matched <- sample_matched_supp()

    perm_conv <- run_pathway_convergence(
      path_genes_long_clean = path_genes_long_clean,
      epi_genes = epi_genes,
      ebv_genes = rand_ebv_matched,
      universe_genes = universe_genes
    )

    finite_or <- perm_conv$odds_ratio[is.finite(perm_conv$odds_ratio)]

    data.frame(
      iter = i,
      n_sig = sum(perm_conv$fisher_fdr < 0.05, na.rm = TRUE),
      max_or = max(finite_or, na.rm = TRUE),
      mean_shared = mean(perm_conv$n_shared_in_path, na.rm = TRUE),
      max_shared = max(perm_conv$n_shared_in_path, na.rm = TRUE)
    )
  }
)

empirical_matched_supp <- data.frame(
  metric = c("n_sig", "max_or", "mean_shared", "max_shared"),
  observed = c(
    real_summary$n_sig,
    real_summary$max_or,
    real_summary$mean_shared,
    real_summary$max_shared
  ),
  null_mean_matched = c(
    mean(perm_results_matched_supp$n_sig),
    mean(perm_results_matched_supp$max_or),
    mean(perm_results_matched_supp$mean_shared),
    mean(perm_results_matched_supp$max_shared)
  ),
  p_matched = c(
    emp_p_supp(perm_results_matched_supp$n_sig, real_summary$n_sig),
    emp_p_supp(perm_results_matched_supp$max_or, real_summary$max_or),
    emp_p_supp(perm_results_matched_supp$mean_shared, real_summary$mean_shared),
    emp_p_supp(perm_results_matched_supp$max_shared, real_summary$max_shared)
  )
)

sig_pathways_named_supp <- hypothesis_table2_stats %>%
  filter(!is.na(fisher_fdr), fisher_fdr < 0.05) %>%
  dplyr::select(pathway_id, pathway_name, network_group,
                pathway_size, n_shared_in_path, n_epi_in_path, n_ebv_in_path,
                shared_gene_coverage_pct, odds_ratio, fisher_p, fisher_fdr,
                combined_score, shared_genes_in_path) %>%
  arrange(fisher_fdr, desc(odds_ratio))
write.csv(sig_pathways_named_supp,
          file.path(out_dir_supp, "TableS_significant_convergent_pathways.csv"),
          row.names = FALSE)

## all EBV proteins with experimentally reported epitopes, any disease
recognized_core <- bind_rows(
  MS_viral_recognition,
  SLE_viral_recognition,
  RA_viral_recognition
) %>%
  filter(!is.na(Viral_Label), Viral_Label != "") %>%
  mutate(Viral_Core = sub("_.*$", "", Viral_Label)) %>%
  distinct(Viral_Core) %>%
  pull(Viral_Core)

## every EBV protein in the interaction data, collapsed the same way
ppis_core <- ppis_std2 %>%
  mutate(
    Viral_Core = sub("_.*$", "", Viral_Label),
    Gene       = toupper(trimws(Gene))
  ) %>%
  filter(!is.na(Gene), Gene != "", Gene %in% universe_genes) %>%
  distinct(Viral_Core, Gene)

all_core <- unique(ppis_core$Viral_Core)
nonrecognized_core <- setdiff(all_core, recognized_core)



################################################################################
## STEP 2. Host target gene sets for each arm
################################################################################

genes_recog <- ppis_core %>%
  filter(Viral_Core %in% recognized_core) %>%
  pull(Gene) %>% unique()

genes_nonrecog <- ppis_core %>%
  filter(Viral_Core %in% nonrecognized_core) %>%
  pull(Gene) %>% unique() %>%
  setdiff(genes_recog)


if (length(genes_nonrecog) < 20) {
      }


################################################################################
## STEP 3. THE TEST
################################################################################

a_rec <- sum(genes_recog    %in% epi_genes)          # recognized, autoantigen
b_rec <- length(genes_recog) - a_rec                 # recognized, not
a_non <- sum(genes_nonrecog %in% epi_genes)          # non-recognized, autoantigen
b_non <- length(genes_nonrecog) - a_non              # non-recognized, not

contingency_rec <- matrix(
  c(a_rec, b_rec, a_non, b_non),
  nrow = 2, byrow = TRUE,
  dimnames = list(
    c("Recognized EBV protein targets", "Non-recognized EBV protein targets"),
    c("Autoantigen", "Not autoantigen")
  )
)

ft_rec <- fisher.test(contingency_rec, alternative = "greater")

main_result_rec <- data.frame(
  comparison        = "Recognized vs non-recognized EBV protein host targets",
  n_recognized      = length(genes_recog),
  n_nonrecognized   = length(genes_nonrecog),
  pct_auto_recog    = round(100 * a_rec / length(genes_recog), 2),
  pct_auto_nonrecog = round(100 * a_non / max(length(genes_nonrecog), 1), 2),
  odds_ratio        = unname(ft_rec$estimate),
  ci_low            = ft_rec$conf.int[1],
  p_value           = ft_rec$p.value
)

write.csv(main_result_rec,
          file.path(out_dir_rec, "Table_recognized_vs_nonrecognized.csv"),
          row.names = FALSE)


################################################################################
## STEP 4. Per-disease 
################################################################################


recognition_by_disease <- list(
  MS  = MS_viral_recognition,
  SLE = SLE_viral_recognition,
  RA  = RA_viral_recognition
)

per_disease_rec <- bind_rows(lapply(names(recognition_by_disease), function(dis) {
  
  rec_d <- recognition_by_disease[[dis]] %>%
    filter(!is.na(Viral_Label), Viral_Label != "") %>%
    mutate(Viral_Core = sub("_.*$", "", Viral_Label)) %>%
    distinct(Viral_Core) %>% pull(Viral_Core)
  
  epi_d <- all_epi_long %>%
    filter(Disease == dis) %>%
    mutate(Gene = toupper(trimws(Gene))) %>%
    filter(Gene %in% universe_genes) %>%
    distinct(Gene) %>% pull(Gene)
  
  g_rec <- ppis_core %>% filter(Viral_Core %in% rec_d) %>% pull(Gene) %>% unique()
  g_non <- ppis_core %>% filter(!Viral_Core %in% rec_d) %>% pull(Gene) %>%
    unique() %>% setdiff(g_rec)
  
  if (length(g_rec) < 5 || length(g_non) < 5) {
    return(data.frame(Disease = dis, n_recognized_proteins = length(rec_d),
                      n_genes_recog = length(g_rec), n_genes_nonrecog = length(g_non),
                      pct_auto_recog = NA, pct_auto_nonrecog = NA,
                      odds_ratio = NA, p_value = NA))
  }
  
  tab_d <- matrix(c(
    sum(g_rec %in% epi_d), length(g_rec) - sum(g_rec %in% epi_d),
    sum(g_non %in% epi_d), length(g_non) - sum(g_non %in% epi_d)
  ), nrow = 2, byrow = TRUE)
  
  ft_d <- fisher.test(tab_d, alternative = "greater")
  
  data.frame(
    Disease               = dis,
    n_recognized_proteins = length(rec_d),
    n_genes_recog         = length(g_rec),
    n_genes_nonrecog      = length(g_non),
    pct_auto_recog        = round(100 * sum(g_rec %in% epi_d) / length(g_rec), 2),
    pct_auto_nonrecog     = round(100 * sum(g_non %in% epi_d) / length(g_non), 2),
    odds_ratio            = unname(ft_d$estimate),
    p_value               = ft_d$p.value
  )
}))

write.csv(per_disease_rec,
          file.path(out_dir_rec, "Table_recognized_per_disease.csv"),
          row.names = FALSE)


################################################################################
## STEP 5. Where do recognized-protein targets sit?
################################################################################


pathway_recog_tbl <- path_genes_long_clean %>%
  group_by(pathway_id, pathway_name) %>%
  summarise(genes = list(unique(Gene)), pathway_size = n_distinct(Gene),
            .groups = "drop") %>%
  rowwise() %>%
  mutate(
    n_recog_targets      = sum(genes_recog %in% genes),
    n_nonrecog_targets   = sum(genes_nonrecog %in% genes),
    n_recog_and_auto     = sum(intersect(genes_recog, epi_genes) %in% genes),
    n_nonrecog_and_auto  = sum(intersect(genes_nonrecog, epi_genes) %in% genes),
    pct_recog_auto       = ifelse(n_recog_targets > 0,
                                  round(100 * n_recog_and_auto / n_recog_targets, 1), NA),
    recog_auto_genes     = paste(sort(intersect(intersect(genes_recog, epi_genes),
                                                unlist(genes))), collapse = "; ")
  ) %>%
  ungroup() %>%
  dplyr::select(-genes) %>%
  arrange(desc(n_recog_and_auto))
write.csv(pathway_recog_tbl,
          file.path(out_dir_rec, "Table_pathway_recognized_targets.csv"),
          row.names = FALSE)

## antigen processing and presentation specifically
ap_row <- pathway_recog_tbl %>% filter(pathway_id == "hsa04612")
if (nrow(ap_row) > 0) {
        }


## NEGATIVE CONTROL:
## RECOGNIZED vs NON-RECOGNIZED EBV PROTEINS
################################################################################

library(dplyr)
library(tidyr)
library(stringr)

# ---------------------------------------------------------------------------
# 1. All EBV proteins in the PPI dataset
# ---------------------------------------------------------------------------

all_ebv <- ppis_std %>%
  filter(!is.na(Viral_Label), Viral_Label != "") %>%
  pull(Viral_Label) %>%
  as.character() %>%
  trimws() %>%
  unique()



# ---------------------------------------------------------------------------
# 2. Recognized EBV proteins from your EXISTING recognition analysis
# ---------------------------------------------------------------------------

recognized_ebv <- viral_sigpath_summary_all %>%
  filter(
    !is.na(recognized_viral_proteins),
    recognized_viral_proteins != ""
  ) %>%
  pull(recognized_viral_proteins) %>%
  paste(collapse = "; ") %>%
  str_split(";\\s*") %>%
  unlist() %>%
  trimws() %>%
  unique()

# Keep only proteins actually present in the PPI dataset
recognized_ebv <- intersect(recognized_ebv, all_ebv)

nonrecognized_ebv <- setdiff(all_ebv, recognized_ebv)


# ---------------------------------------------------------------------------
# 3. Assign recognition status to each EBV protein
# ---------------------------------------------------------------------------

ebv_status <- tibble(
  Viral_Label = all_ebv,
  recognition_status = ifelse(
    all_ebv %in% recognized_ebv,
    "Recognized",
    "Non_recognized"
  )
)

print(table(ebv_status$recognition_status))


# ---------------------------------------------------------------------------
# 4. EBV protein -> host gene interactions
# ---------------------------------------------------------------------------

ppi_status <- ppis_std %>%
  transmute(
    Viral_Label = trimws(as.character(Viral_Label)),
    Gene = toupper(trimws(as.character(Gene)))
  ) %>%
  filter(
    !is.na(Viral_Label),
    Viral_Label != "",
    !is.na(Gene),
    Gene != ""
  ) %>%
  distinct() %>%
  inner_join(
    ebv_status,
    by = "Viral_Label"
  )


# ---------------------------------------------------------------------------
# 5. Map EBV proteins to the 87 retained pathways through their host targets
# ---------------------------------------------------------------------------

pathway_status <- path_genes_long %>%
  transmute(
    pathway_id = pathway_id,
    pathway_name = pathway_name,
    Gene = toupper(trimws(as.character(Gene)))
  ) %>%
  distinct() %>%
  inner_join(
    ppi_status,
    by = "Gene",
    relationship = "many-to-many"
  ) %>%
  distinct(
    pathway_id,
    pathway_name,
    Viral_Label,
    recognition_status
  ) %>%
  count(
    pathway_id,
    pathway_name,
    recognition_status,
    name = "n_proteins"
  ) %>%
  pivot_wider(
    names_from = recognition_status,
    values_from = n_proteins,
    values_fill = 0
  )


# ---------------------------------------------------------------------------
# 6. Merge with pathway table
# ---------------------------------------------------------------------------

negative_control_df <- integrated_pathway_table_v2 %>%
  select(
    pathway_id,
    pathway_name,
    n_shared_in_path,
    n_ebv_in_path,
    odds_ratio,
    n_recognized_ebv_proteins
  ) %>%
  left_join(
    pathway_status,
    by = c("pathway_id", "pathway_name")
  ) %>%
  mutate(
    Recognized = replace_na(Recognized, 0),
    Non_recognized = replace_na(Non_recognized, 0)
  )


# ---------------------------------------------------------------------------
# 7. shared EBV-autoantigen genes vs recognized proteins
# ---------------------------------------------------------------------------

cor_rec_shared <- cor.test(
  negative_control_df$n_shared_in_path,
  negative_control_df$Recognized,
  method = "spearman",
  exact = FALSE
)


# ---------------------------------------------------------------------------
# 8. NEGATIVE CONTROL:
# shared EBV-autoantigen genes vs NON-recognized proteins
# ---------------------------------------------------------------------------



negative_control_df <- integrated_pathway_table_v2 %>%
  dplyr::select(
    pathway_id,
    pathway_name,
    n_shared_in_path,
    n_ebv_in_path,
    odds_ratio,
    n_recognized_ebv_proteins
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

cor_rec_shared <- cor.test(
  negative_control_df$n_shared_in_path,
  negative_control_df$Recognized,
  method = "spearman",
  exact = FALSE
)

cor_nonrec_shared <- cor.test(
  negative_control_df$n_shared_in_path,
  negative_control_df$Non_recognized,
  method = "spearman",
  exact = FALSE
)

cor_rec_or <- cor.test(
  negative_control_df$odds_ratio,
  negative_control_df$Recognized,
  method = "spearman",
  exact = FALSE
)

cor_nonrec_or <- cor.test(
  negative_control_df$odds_ratio,
  negative_control_df$Non_recognized,
  method = "spearman",
  exact = FALSE
)






