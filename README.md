# EBV-Autoantigen-Convergence
Code for investigating Epstein–Barr virus–autoantigen relationships across sequence, pathway, and protein scales in multiple sclerosis, systemic lupus erythematosus, and rheumatoid arthritis.

“Pathway-level Epstein–Barr virus perturbation, rather than exact
sequence sharing, tracks autoantigen targeting in multiple sclerosis.”

The analyses examine EBV–autoantigen relationships at sequence, pathway,
immune-recognition, and protein levels in multiple sclerosis (MS),
systemic lupus erythematosus (SLE), and rheumatoid arthritis (RA).

REPOSITORY FILES

01_pathway_convergence_and_recognition.R

Corresponds primarily to the Methods sections:

-   Pathway convergence and enrichment analysis
-   Immune recognition and protein-level comparison
-   Correlation analysis
-   Permutation analysis

This script contains the main pathway- and recognition-level analyses,
including:

-   correlations between shared EBV–autoantigen pathway burden and the
    number of immune-recognized EBV proteins;
-   pathway-size-normalized recognition correlations;
-   annotation-depth-matched EBV-target permutation testing;
-   identification/export of significantly convergent pathways;
-   comparison of host targets of immune-recognized versus
    non-recognized EBV proteins;
-   disease-specific recognition analyses for MS, SLE, and RA;
-   pathway localization of recognition-associated host targets,
    including antigen processing and presentation; and
-   recognized versus non-recognized pathway-level negative-control
    correlations.

02_recognition_negative_controls_and_mirror_permutation.R

Corresponds primarily to the Methods sections:

-   Correlation analysis
-   Permutation analysis

This script contains the dedicated recognition negative-control and
mirror-null analyses, including:

-   collapse of EBV protein labels to biological protein identities;
-   recognized versus non-recognized EBV protein counts within the 87
    pathways;
-   correlations of recognized and non-recognized EBV proteins with
    shared-gene burden and pathway enrichment odds ratios; and
-   the annotation-matched mirror null in which autoantigen genes are
    randomized while EBV-target genes are held fixed, with empirical
    statistics for significant-pathway count, maximum odds ratio, mean
    shared-gene count, and maximum shared-gene overlap.

03_sequence_and_protein_level_prediction_analyses.R

Corresponds primarily to the Methods sections:

-   EBV protein characteristics and peptide-sharing burden
-   Broad protein-level predictive modelling

This script contains:

-   construction and evaluation of the broad M0, M1, M2, M0+M1, and M3
    protein-level models using five-fold network-blocked
    cross-validation;
-   pooled AUPRC and AUROC evaluation of the broad models;
-   reconstruction of exact EBV–human 8-amino-acid peptide sharing;
-   EBV protein length, human peptide-match burden, matched-autoantigen
    burden, and autoantigen-fraction analyses;
-   integration of EBV host-interaction burden with the peptide-sharing
    analysis;
-   collapse of sequence accessions/variants to biological EBV protein
    identities for the final protein-level tests;
-   the Spearman and length-adjusted partial-Spearman tests reported for
    EBV protein characteristics; and
-   the EBNA2 and BNLF2A descriptive comparison used in the
    sequence-level results.

The hypothesis-driven MS F1–F4/pathway-perturbation analysis is
maintained separately and is not contained in this script.

04_broad_protein_predictive_modelling.R

Corresponds to the Methods section:

-   Broad protein-level predictive modelling

This script is the compact final rerun of the broad protein-level
predictive modelling analysis. It uses the prepared modelling dataset
and network-blocked fold assignments and includes:

-   M0 background protein/network characteristics;
-   M1 exact EBV–human peptide-sharing features;
-   M2 direct/local EBV interaction-network features;
-   the combined M0 + M1 model;
-   the full M3 model;
-   five-fold network-blocked cross-validation; and
-   pooled out-of-fold AUPRC and AUROC evaluation.

05_MS_F1_F4_pathway_perturbation.R

Corresponds primarily to the Methods section:
- Hypothesis-driven MS analysis

This script contains the final MS-specific F1–F4 analysis, including:
- F1: whether a human protein is targeted by an MS immune-recognized EBV protein;
- F2: length of the targeting MS immune-recognized EBV protein;
- F3: host-interaction burden of the targeting MS immune-recognized EBV protein;
- F4: pathway-level EBV perturbation;
- background-adjusted logistic regression across F1–F4 model combinations;
- estimation of odds ratios and 95% confidence intervals;
- the F1–F4 forest-plot analysis; and
- the full-human-KEGG F4 sensitivity analysis using leave-one-out pathway perturbation and the maximum value across pathways.

ANALYSIS STRUCTURE

Across the repository, the code follows the biological scales used in
the manuscript:

exact sequence sharing → pathway convergence → immune-recognition
comparison → protein-level modelling

The pathway-level analyses test whether EBV-targeted host genes and
documented autoantigens converge non-randomly within biological
pathways. Recognition analyses test whether pathway-level organization
differs according to experimentally reported EBV immune recognition and
whether recognition-associated host-target enrichment is
disease-specific. Sequence-level analyses examine exact EBV–human
peptide sharing and EBV protein characteristics. Protein-level modelling
evaluates whether EBV-related features add information beyond general
protein and network characteristics. The hypothesis-driven MS analysis
then evaluates F1–F4, with particular emphasis on pathway-level EBV
perturbation (F4), including its sensitivity across the full human KEGG
pathway universe.

DATA SOURCES

The analyses use publicly available data from DISEASES, KEGG, PHISTO,
VirHostNet, and IEDB, as described in the manuscript Methods and Data
and code availability sections. Users should consult the original
resources for their respective access conditions and licenses.

CODE USE

Copyright © 2026 Anna Onisiforou. All rights reserved.

This repository is made publicly available for transparency and
reproducibility of the analyses reported in the associated publication.
No license is granted for reuse, modification, redistribution, or
incorporation of the source code into other software or projects without
prior written permission from the author. Third-party datasets and
resources remain subject to their respective terms and licenses.
