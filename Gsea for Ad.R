title: "New thoughts"
author: "HU QIANG"
date: "2025-10-14"
output:
  html_document: default
  pdf_document: default
---

```{r setup, include=FALSE}
# Use knitr options to set the global root directory
knitr::opts_knit$set(root.dir = "C:/Users/胡强/Desktop")

```

```{r}
# Install BiocManager if not installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install required Bioconductor and dependency packages
bioc_packages <- c("tweeDEseq", "DSS", "NOISeq", "dearseq", "polyester", "edgeR", "DESeq2", "limma", "ALDEx2")
lapply(bioc_packages, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  }
})

# Install devtools and SIEVE (assuming SIEVE is stored locally)
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
library(devtools)

# Please make sure this path points to your local SIEVE package directory
sieve_path <- "C:/Users/胡强/Desktop/SIEVE"
if (dir.exists(sieve_path)) {
  devtools::install(sieve_path, upgrade = "never", quiet = TRUE)
} else {
  warning("SIEVE package path not found. Please ensure 'sieve_path' is correct.")
}

```

```{r}
# Load required libraries
# Ensure GSEA and annotation packages are loaded
if (!requireNamespace("fgsea", quietly = TRUE)) {
  BiocManager::install("fgsea", update = FALSE, ask = FALSE)
}
library(fgsea)
library(clusterProfiler)
library(org.Hs.eg.db)

library(readr); library(compositions); library(sn)
library(edgeR); library(DESeq2); library(limma)
library(tweeDEseq); library(ALDEx2); library(DSS)
library(vioplot); library(VennDiagram); library(gridExtra)
library(dplyr); library(httr); library(jsonlite)
library(SIEVE)

library(gamlss)
library(MASS); library(ggplot2)

```

```{r}
# Read input data
ad_counts <- read.csv("MayoRNAseq_RNAseq_TCX_geneCounts.csv", header = TRUE, stringsAsFactors = FALSE, row.names = 1)
ad_meta <- read.csv("MayoRNAseq_individual_metadata_031422.csv", header = TRUE, stringsAsFactors = FALSE)

# 1. Clean metadata and define groups
ad_meta <- ad_meta[-which(is.na(ad_meta$diagnosis)), ]
control_id <- (ad_meta$individualID)[ad_meta$diagnosis == "control"]
ad_id <- (ad_meta$individualID)[ad_meta$diagnosis == "Alzheimer Disease"]

# 2. Extract and merge count matrices
ad_counts_id <- colnames(ad_counts)
control_vector <- sapply(ad_counts_id, function(k) k %in% control_id)
control_counts <- ad_counts[, control_vector]
ad_vector <- sapply(ad_counts_id, function(k) k %in% ad_id)
ad_counts1 <- ad_counts[, ad_vector]
N_ad <- length(control_counts) + length(ad_counts1)
mayo_counts1 <- as.matrix(cbind(control_counts, ad_counts1),
                          nrows = dim(ad_counts)[1], ncol = N_ad)

# 3. Define group vector (0 = Control, 1 = AD)
group2 = c(rep(0, length(control_counts)), rep(1, length(ad_counts1)))

```

```{r}
# Gene filtering: CPM > 0.5 and zero proportion < 85%
CPM2 <- cpm(mayo_counts1)
keep <- rowMeans(CPM2[, 1:length(control_counts)]) > 0.5 &
  rowMeans(CPM2[, (length(control_counts)+1):N_ad]) > 0.5 &
  apply(mayo_counts1[, 1:length(control_counts)], 1, function(k) length(k[k == 0])/length(k)) < 0.85 &
  apply(mayo_counts1[, (length(control_counts)+1):N_ad], 1, function(k) length(k[k == 0])/length(k)) < 0.85

mayo_counts_filter <- mayo_counts1[keep, ]
cat(paste("Genes retained after filtering:", dim(mayo_counts_filter)[1], "\n"))

```

```{r}
# CLR transformation function
clr.transform <- function(data = NULL){
  data[data == 0] <- 1/2
  clr.count <- t(clr(t(data)))
  clr.count <- matrix(as.numeric(clr.count),
                      nrow = dim(data)[1],
                      ncol = dim(data)[2])
  row.names(clr.count) <- row.names(data)
  return(clr.count)
}

clr_counts <- clr.transform(data = mayo_counts_filter)

```

```{r}
# Plot raw vs CLR-transformed distributions
par(mfrow = c(1, 2))

raw_values <- as.vector(as.matrix(mayo_counts_filter))
raw_values <- raw_values[!is.na(raw_values) & is.finite(raw_values)]
hist(raw_values, breaks = 30, probability = TRUE, col = "gray80", border = "white",
     main = "(a) Raw Counts", xlab = "Raw count")
lines(density(raw_values), lwd = 2)
rug(raw_values)

clr_values <- as.vector(as.matrix(clr_counts))
clr_values <- clr_values[!is.na(clr_values) & is.finite(clr_values)]
hist(clr_values, breaks = 30, probability = TRUE, col = "gray90", border = "white",
     main = "(b) CLR-transformed count", xlab = "CLR-transformed count")
lines(density(clr_values), lwd = 2)
rug(clr_values)

par(mfrow = c(1, 1))

```

```{r}
# SIEVE model fitting and DE/DV/DS testing
t_sieve <- proc.time()

clrSeq_result <- clrSeq(clr_counts, group = group2)
clrSIEVE_result <- clrSIEVE(clrSeq_result = clrSeq_result,
                            alpha_level = 0.05,
                            order_DE = F, order_LFC = F, order_DS = F, order_sieve = F)

cat(paste("SIEVE run time:", as.numeric(proc.time() - t_sieve)[3], "seconds\n"))

# Extract test results and construction of U statistic for GSEA
# Extract original SIEVE results
res_de <- clrSIEVE_result$clrDE_test
res_dv <- clrSIEVE_result$clrDV_test
res_ds <- clrSIEVE_result$clrDS_test

# Ensure genes in the three tables are perfectly aligned
common_genes <- Reduce(intersect, list(rownames(res_de), rownames(res_dv), rownames(res_ds)))

df_stat <- data.frame(
  DE = res_de[common_genes, "DE"],
  DV = res_dv[common_genes, "LFC"],
  DS = res_ds[common_genes, "DS"],
  row.names = common_genes,
  check.names = FALSE
)

```

```{r}
# --- Core: Regenerate U Statistics ---
# Rule: U = the component among (DE, DV, DS) with the largest absolute value, preserving the original sign
# Tie rule: If absolute values are tied, take the first one (DE > DV > DS)
df_stat$U1 <- apply(df_stat[, c("DE","DV","DS")], 1, function(x) {
  m <- max(abs(x), na.rm = TRUE)
  idx <- which(abs(x) == m)[1]
  x[idx]
})

df_stat$U2 <- with(df_stat, sign(DE) * sqrt(DE^2 + (3/4)*DV^2 + (1/2)*DS^2))

df_stat$U3 <- apply(df_stat[, c("DE","DV","DS")], 1, function(x) {
  k <- which.max(abs(x))                               # Dominant component (max unweighted abs)
  M <- sqrt(x[1]^2 + (3/4)*x[2]^2 + (1/2)*x[3]^2)      # Weighted magnitude
  sign(x[k]) * M                                       # U3
})

# ---- Critical: ID Unification (Handling ENSG version numbers like ".12") ----
strip_ens_version <- function(x) sub("\\..*$", "", x)
rownames(df_stat) <- strip_ens_version(rownames(df_stat))

# Construct ranked lists (stats must be a named numeric vector)
de_ranks_sorted <- sort(na.omit(setNames(df_stat$DE, rownames(df_stat))), decreasing = TRUE)
dv_ranks_sorted <- sort(na.omit(setNames(df_stat$DV, rownames(df_stat))), decreasing = TRUE)
ds_ranks_sorted <- sort(na.omit(setNames(df_stat$DS, rownames(df_stat))), decreasing = TRUE)
u1_ranks_sorted  <- sort(na.omit(setNames(df_stat$U1,  rownames(df_stat))), decreasing = TRUE)
u2_ranks_sorted  <- sort(na.omit(setNames(df_stat$U2,  rownames(df_stat))), decreasing = TRUE)
u3_ranks_sorted  <- sort(na.omit(setNames(df_stat$U3,  rownames(df_stat))), decreasing = TRUE)

message(paste("DE ranked list length:", length(de_ranks_sorted)))
message(paste("DV ranked list length:", length(dv_ranks_sorted)))
message(paste("DS ranked list length:", length(ds_ranks_sorted)))
message(paste("U1 ranked list length:", length(u1_ranks_sorted)))
message(paste("U2 ranked list length:", length(u2_ranks_sorted)))
message(paste("U3 ranked list length:", length(u3_ranks_sorted)))

```

```{r}
# U1 = sgn(W_k,i) * max( |W_DE,i|, |W_DV,i|, |W_DS,i| )
check_table1 <- df_stat %>%
  as.data.frame() %>%
  # Add absolute value helper columns for quick comparison
  mutate(
    abs_DE = abs(DE),
    abs_DV = abs(DV),
    abs_DS = abs(DS)
  ) %>%
  # Record the source of U1 for each row (which column contributed the max absolute value)
  mutate(
    Source = apply(.[, c("DE", "DV", "DS")], 1, function(x) {
      names(x)[which.max(abs(x))]
    }),
    # ===== Add U1 =====
    U1 = apply(.[, c("DE", "DV", "DS")], 1, function(x) {
      k <- which.max(abs(x))          # Find dominant component
      sign(x[k]) * abs(x[k])          # sgn(W_k,i) * |W_k,i|
    })
  ) %>%
  # Take only the first 20 rows for display
  head(20)

```

```{r}
# U2 = sgn(W_DE,i) * sqrt( W_DE,i^2 + (3/4)*W_DV,i^2 + (1/2)*W_DS,i^2 )
check_table2 <- df_stat %>%
  as.data.frame() %>%
  mutate(
    contrib_DE = DE^2,
    contrib_DV = (3/4) * DV^2,
    contrib_DS = (1/2) * DS^2,
    U2 = sign(DE) * sqrt(DE^2 + (3/4) * DV^2 + (1/2) * DS^2)
  ) %>%
  rowwise() %>%
  mutate(
    Source = names(c_across(c(contrib_DE, contrib_DV, contrib_DS)))[
      which.max(c_across(c(contrib_DE, contrib_DV, contrib_DS)))
    ]
  ) %>%
  ungroup() %>%
  head(20)

```

```{r}
# Mi = sqrt( 1*W_DE,i^2 + (3/4)*W_DV,i^2 + (1/2)*W_DS,i^2 )
# Si = sgn(W_k,i),  k = arg max {|W_DE,i|, |W_DV,i|, |W_DS,i|}
# U3 = Si * Mi

check_table3 <- df_stat %>%
  as.data.frame() %>%
  mutate(
    # Weighted square contributions inside the radical
    w_DE = 1 * (DE^2),
    w_DV = (3/4) * (DV^2),
    w_DS = (1/2) * (DS^2),
    # Weighted magnitude M_i
    M_w = sqrt(w_DE + w_DV + w_DS),
    # Dominant component (compare by unweighted abs: |DE|, |DV|, |DS|)
    Source_dom = apply(.[, c("DE", "DV", "DS")], 1, function(x) {
      names(x)[which.max(abs(x))]
    }),
    # Sign S_i: take the sign of the dominant component
    S_dom = case_when(
      Source_dom == "DE" ~ sign(DE),
      Source_dom == "DV" ~ sign(DV),
      TRUE               ~ sign(DS)
    ),
    # Statistic for plotting: U_i = S_i * M_i
    U3 = S_dom * M_w
  ) %>%
  head(20)

```

```{r}
# --- Ranked List Preparation for Updated GSEA ---

# Define a convenient function to handle named vectors
get_sorted_ranks <- function(val_vec, names_vec) {
  res <- setNames(val_vec, names_vec)
  return(sort(na.omit(res), decreasing = TRUE))
}

# Construct ranked lists containing DE, DV, DS and the new U indicators
de_ranks_sorted <- get_sorted_ranks(df_stat$DE, rownames(df_stat))
dv_ranks_sorted <- get_sorted_ranks(df_stat$DV, rownames(df_stat))
ds_ranks_sorted <- get_sorted_ranks(df_stat$DS, rownames(df_stat))
u1_ranks_sorted  <- get_sorted_ranks(df_stat$U1,  rownames(df_stat))
u2_ranks_sorted  <- get_sorted_ranks(df_stat$U2,  rownames(df_stat))
u3_ranks_sorted  <- get_sorted_ranks(df_stat$U3,  rownames(df_stat))

```

```{r}
# --- Target GO IDs ---
go_ids_table4 <- c(
  "GO:0034613", "GO:0072599", "GO:0090150", "GO:0045047", "GO:0006614",
  "GO:0022610", "GO:0007155", "GO:0010647",
  "GO:0043062", "GO:0030198",
  "GO:0001568", "GO:0042060",
  "GO:0006412", "GO:0006401",
  "GO:0019884", "GO:0045619",
  "GO:0019363", "GO:0019752"
)

# --- 1) Primary Method: direct select GOALL -> ENSEMBL ---
go_map1 <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = go_ids_table4,
  keytype = "GOALL",
  columns = c("GOALL", "ONTOLOGYALL", "ENSEMBL")
)

go_map1 <- go_map1 %>%
  dplyr::filter(!is.na(ENSEMBL)) %>%
  dplyr::filter(ONTOLOGYALL == "BP") %>%
  dplyr::select(go_id = GOALL, ENSEMBL) %>%
  dplyr::distinct()

gene_sets_table4 <- split(go_map1$ENSEMBL, go_map1$go_id)

missing_gos <- setdiff(go_ids_table4, names(gene_sets_table4))
message("GO gene sets obtained via primary method: ", length(gene_sets_table4))
message("Missing GO count via primary method: ", length(missing_gos))
print(missing_gos)

# --- 2) Remedial Method: use org.Hs.egGO2ALLEGS for missing GOs (GO -> Entrez -> ENSEMBL) ---
if (length(missing_gos) > 0) {
  # org.Hs.egGO2ALLEGS is an environment: key=GO, value=Entrez list
  # Some GOs can get Entrez here, then map back to ENSEMBL
  for (go_id in missing_gos) {

    entrez_vec <- tryCatch({
      obj <- mget(go_id, envir = org.Hs.egGO2ALLEGS, ifnotfound = NA)[[1]]
      if (all(is.na(obj))) character(0) else as.character(obj)
    }, error = function(e) character(0))

    if (length(entrez_vec) == 0) next

    ens_map <- AnnotationDbi::mapIds(
      org.Hs.eg.db,
      keys = entrez_vec,
      column = "ENSEMBL",
      keytype = "ENTREZID",
      multiVals = "list"
    )

    ens_vec <- unique(unlist(ens_map))
    ens_vec <- ens_vec[!is.na(ens_vec)]

    if (length(ens_vec) > 0) {
      gene_sets_table4[[go_id]] <- ens_vec
    }
  }
}

# --- 3) Final Check ---
missing_gos_final <- setdiff(go_ids_table4, names(gene_sets_table4))
message("Final successful GO gene sets: ", length(gene_sets_table4))
message("Final missing GO count: ", length(missing_gos_final))
print(missing_gos_final)

# (Optional) View the size of each GO set
print(sort(lengths(gene_sets_table4)))

```

```{r}
# GSEA Core Calculation (Permutation Test for DE, DV, DS)
N_PERM <- 10000   # Define number of permutations (1000 is standard to ensure reliability).
set.seed(421)

# Wrap the fgsea call in suppressWarnings() to remove 'package:stats' warnings
suppressWarnings({
  fgsea_DE_results <- fgsea(pathways = gene_sets_table4, stats = de_ranks_sorted, nPermSimple = N_PERM, eps = 0)
  fgsea_DV_results <- fgsea(pathways = gene_sets_table4, stats = dv_ranks_sorted, nPermSimple = N_PERM, eps = 0)
  fgsea_DS_results <- fgsea(pathways = gene_sets_table4, stats = ds_ranks_sorted, nPermSimple = N_PERM, eps = 0)
  fgsea_U1_results  <- fgsea(pathways = gene_sets_table4, stats = u1_ranks_sorted,  nPermSimple = N_PERM, eps = 0)
  fgsea_U2_results  <- fgsea(pathways = gene_sets_table4, stats = u2_ranks_sorted,  nPermSimple = N_PERM, eps = 0)
  fgsea_U3_results  <- fgsea(pathways = gene_sets_table4, stats = u3_ranks_sorted,  nPermSimple = N_PERM, eps = 0)
})
# nPermSimple = N_PERM, eps = 0 GSEA parameters. Set permutations and set eps to 0 for more accurate low P-value estimation.
# eps (Epsilon) Minimum precision for P-value estimation. Improves accuracy of low P-values.

```

```{r}
# 1. Merge Results (for Bubble Plot)
# Extract results for each group
results_DE_summary <- fgsea_DE_results %>% dplyr::mutate(Type = "DE") %>% dplyr::select(pathway, ES, NES, padj, Type, size)
results_DV_summary <- fgsea_DV_results %>% dplyr::mutate(Type = "DV") %>% dplyr::select(pathway, ES, NES, padj, Type, size)
results_DS_summary <- fgsea_DS_results %>% dplyr::mutate(Type = "DS") %>% dplyr::select(pathway, ES, NES, padj, Type, size)
results_U1_summary <- fgsea_U1_results %>% dplyr::mutate(Type = "U1") %>% dplyr::select(pathway, ES, NES, padj, Type, size)
results_U2_summary <- fgsea_U2_results %>% dplyr::mutate(Type = "U2") %>% dplyr::select(pathway, ES, NES, padj, Type, size)
results_U3_summary <- fgsea_U3_results %>% dplyr::mutate(Type = "U3") %>% dplyr::select(pathway, ES, NES, padj, Type, size)

# Combine and add -log10(padj) column
all_gsea_results <- dplyr::bind_rows(
  results_DE_summary, 
  results_DV_summary, 
  results_DS_summary,
  results_U1_summary, 
  results_U2_summary, 
  results_U3_summary
) %>%
  dplyr::mutate(neg_log10_padj = -log10(padj))

```

```{r}
# Result Merging and Visualization (Final Integrated Version)

# 1. Merge Results (for Bubble Plot)
# Use dplyr:: to avoid naming conflicts
results_DE_summary <- fgsea_DE_results %>% dplyr::mutate(Type = "DE") %>% dplyr::select(pathway, ES, NES, padj, Type, size)
results_DV_summary <- fgsea_DV_results %>% dplyr::mutate(Type = "DV") %>% dplyr::select(pathway, ES, NES, padj, Type, size)
results_DS_summary <- fgsea_DS_results %>% dplyr::mutate(Type = "DS") %>% dplyr::select(pathway, ES, NES, padj, Type, size)
results_U1_summary <- fgsea_U1_results %>% dplyr::mutate(Type = "U1") %>% dplyr::select(pathway, ES, NES, padj, Type, size)
results_U2_summary <- fgsea_U2_results %>% dplyr::mutate(Type = "U2") %>% dplyr::select(pathway, ES, NES, padj, Type, size)
results_U3_summary <- fgsea_U3_results %>% dplyr::mutate(Type = "U3") %>% dplyr::select(pathway, ES, NES, padj, Type, size)

all_gsea_results <- dplyr::bind_rows(results_DE_summary, results_DV_summary, results_DS_summary, results_U1_summary, results_U2_summary, results_U3_summary)

# 2. Add GO Term descriptions (for visualization)
go_descriptions <- c(
    "GO:0034613" = "Cellular protein localization", "GO:0072599" = "ER localization", "GO:0090150" = "Protein localization to membrane", 
    "GO:0045047" = "Protein targeting to ER", "GO:0006614" = "SRP-dependent cotranslational protein targeting", 
    "GO:0022610" = "Biological adhesion", "GO:0007155" = "Cell adhesion", "GO:0010647" = "Pos. reg. of cell communication", 
    "GO:0043062" = "Extracellular structure organization", "GO:0030198" = "Extracellular matrix organization", 
    "GO:0001568" = "Blood vessel development", "GO:0042060" = "Wound healing", "GO:0006412" = "Translation", 
    "GO:0006401" = "RNA catabolic process", "GO:0019884" = "Antigen processing and presentation", 
    "GO:0045619" = "Reg. of lymphocyte differentiation", "GO:0019363" = "Pyridine nucleotide biosynthetic process", 
    "GO:0019752" = "Carboxylic acid metabolic process"
)

# Filter and sort data for the bubble plot
plot_data <- all_gsea_results %>%
dplyr::mutate(Term_Desc = go_descriptions[pathway]) %>%
dplyr::filter(padj < 0.2) %>%
dplyr::arrange(desc(NES)) %>%
dplyr::mutate(Term_Desc = factor(Term_Desc, levels = unique(Term_Desc)))

```

```{r}
# --- 8.1 GSEA Bubble Plot (Corrected Syntax) ---
if (nrow(plot_data) > 0) {
  message(paste("Plotting Bubble Chart: Displaying", nrow(plot_data), "enrichment results (FDR < 0.2)."))
  
  gsea_dot_plot <- ggplot(plot_data, aes(x = Type, y = Term_Desc)) +
    geom_point(aes(size = size, color = NES)) + 
    scale_color_gradient2(low = "green", mid = "white", high = "yellow", midpoint = 0, name = "NES") +
    # This geom_point is used to mark significant results
    geom_point(data = dplyr::filter(plot_data, padj < 0.05), shape = 8, size = 1.5, color = "black") + 
    
    theme_minimal(base_size = 12) +
    labs(
      title = "Functional Enrichment (GSEA) of Selected GO Terms (DE vs. DV vs. DS)",
      x = "Differential Analysis Type", 
      y = "GO Term Description", 
      size = "Gene Set Size",
      caption = "NES: Normalized Enrichment Score. Star (*) marks FDR < 0.05."
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
  print(gsea_dot_plot)
} else {
  message("No significant GSEA enrichment results found (FDR < 0.2).")
}

```

```{r}
# =========================
# Print all significant pathways and their FDR
# =========================

sig_table <- plot_data %>%
  dplyr::filter(padj < 0.05) %>%
  dplyr::select(Type, pathway, Term_Desc, NES, padj, size) %>%
  dplyr::arrange(padj)

cat("\n================ Significant Enriched Pathways (FDR < 0.05) ================\n")
if (nrow(sig_table) == 0) {
  cat("⚠️ No pathways satisfy FDR < 0.05.\n")
} else {
  print(sig_table, n = 200)
}
cat("==========================================================================\n\n")

```

```{r}
# =========================
# 8.2 Batch plotting: enrichment curve for the same GO term under 6 statistics
# =========================

selected_go_id <- "GO:0006412"
analysis_types <- c("DE","DV","DS","U1","U2","U3")

if (selected_go_id %in% names(gene_sets_table4)) {

  for (analysis_type in analysis_types) {

    # --- Dynamically select GSEA results and ranked list ---
    if (analysis_type == "DE") {
      fgsea_results  <- fgsea_DE_results
      gene_rank_list <- de_ranks_sorted

    } else if (analysis_type == "DV") {
      fgsea_results  <- fgsea_DV_results
      gene_rank_list <- dv_ranks_sorted

    } else if (analysis_type == "DS") {
      fgsea_results  <- fgsea_DS_results
      gene_rank_list <- ds_ranks_sorted

    } else if (analysis_type == "U1") {
      fgsea_results  <- fgsea_U1_results
      gene_rank_list <- u1_ranks_sorted

    } else if (analysis_type == "U2") {
      fgsea_results  <- fgsea_U2_results
      gene_rank_list <- u2_ranks_sorted

    } else if (analysis_type == "U3") {
      fgsea_results  <- fgsea_U3_results
      gene_rank_list <- u3_ranks_sorted

    } else {
      stop("Analysis type not recognized. Use one of: DE, DV, DS, U1, U2, U3.")
    }

    # Ensure the selected ranked list and result variables exist
    if (exists("gene_rank_list") && exists("fgsea_results")) {

      # Extract and format GSEA result information
      pathway_info <- fgsea_results %>%
        dplyr::filter(pathway == selected_go_id) %>%
        dplyr::mutate(
          formatted_padj = format(padj, digits = 2, scientific = TRUE),
          NES_FDR = paste0("NES=", round(NES, 2), ", FDR=", formatted_padj)
        )

      # If there are no results for this pathway under the current analysis_type, skip
      if (nrow(pathway_info) == 0) {
        message("⚠️ Skip: ", selected_go_id, " has no fgsea result row under ", analysis_type, ".")
        next
      }

      # Construct plot title
      plot_title <- paste0(
        "GSEA Enrichment Plot for ", "\n",
        go_descriptions[selected_go_id], "\n",
        " (", selected_go_id, ") in ", analysis_type, "\n",
        pathway_info$NES_FDR
      )

      # Plot
      enrichment_plot <- fgsea::plotEnrichment(
        pathway = gene_sets_table4[[selected_go_id]],
        stats = gene_rank_list
      ) +
        labs(title = plot_title,
             x = paste0(analysis_type, " Gene List Rank"),
             y = "Enrichment Score (ES)") +
        theme_minimal(base_size = 14) +
        theme(plot.title = element_text(hjust = 0.5))

      print(enrichment_plot)

      # (Optional) Print NES/padj to console for your records
      cat("\n", analysis_type, ": ",
          "NES=", round(pathway_info$NES[1], 3),
          ", padj(FDR)=", format(pathway_info$padj[1], scientific = TRUE, digits = 3),
          "\n", sep = "")

    } else {
      message("Skipping GSEA curve: ranked list or GSEA result variable is missing.")
    }
  }

} else {
  message(paste0("Skipping GSEA curve: Selected GO Term (", selected_go_id, ") gene set is missing in the database."))
}

```

```{r}
# --- 8.2 GSEA Enrichment Curve Plotting (Dynamic Selection Mode) ---

# 1. Select the pathway and analysis mode to visualize
selected_go_id <- "GO:0006401" # <-- Keep pathway unchanged
analysis_type <- "U2"          # <-- Correction: switch to U2 mode

if (selected_go_id %in% names(gene_sets_table4)) {

 # --- Dynamically select GSEA results and ranked list ---
if (analysis_type == "DE") {
  fgsea_results  <- fgsea_DE_results
  gene_rank_list <- de_ranks_sorted

} else if (analysis_type == "DV") {
  fgsea_results  <- fgsea_DV_results
  gene_rank_list <- dv_ranks_sorted

} else if (analysis_type == "DS") {
  fgsea_results  <- fgsea_DS_results
  gene_rank_list <- ds_ranks_sorted

} else if (analysis_type == "U1") {
  fgsea_results  <- fgsea_U1_results
  gene_rank_list <- u1_ranks_sorted
  
} else if (analysis_type == "U2") {
  fgsea_results  <- fgsea_U2_results
  gene_rank_list <- u2_ranks_sorted  

} else if (analysis_type == "U3") {
  fgsea_results  <- fgsea_U3_results
  gene_rank_list <- u3_ranks_sorted   
  
} else {
  stop("Analysis type not recognized. Use one of: DE, DV, DS, U1, U2, U3.")
}

# Ensure the selected ranked list and result variables exist
 if (exists("gene_rank_list") && exists("fgsea_results")) {

 # Extract and format GSEA result information
 pathway_info <- fgsea_results %>%
 dplyr::filter(pathway == selected_go_id) %>%
 dplyr::mutate(
 formatted_padj = format(padj, digits = 2, scientific = TRUE),
 NES_FDR = paste0("NES=", round(NES, 2), ", FDR=", formatted_padj)
)

 # Construct plot title
 plot_title <- paste0("GSEA Enrichment Plot for ", "\n",go_descriptions[selected_go_id],"\n",
 " (", selected_go_id, ") in ", analysis_type, "\n", pathway_info$NES_FDR)

 # Plot enrichment curve
 enrichment_plot <- fgsea::plotEnrichment(
 pathway = gene_sets_table4[[selected_go_id]],
 stats = gene_rank_list # <-- Key: use dynamic variable
 ) +
 labs(title = plot_title, 
 x = paste0(analysis_type, " Gene List Rank"), 
 y = "Enrichment Score (ES)") +
 theme_minimal(base_size = 14) +
 theme(plot.title = element_text(hjust = 0.5))

 print(enrichment_plot)
} else {
 message("Skipping GSEA curve: ranked list or GSEA result variable is missing.")
}
} else {
 message(paste0("Skipping GSEA curve: Selected GO Term (", selected_go_id, ") gene set is missing in the database."))
}

```

```{r}
# Check the size of each GO Term gene set

# Use the lengths() function to get the number of genes for each GO Term in the list
go_term_sizes <- lengths(gene_sets_table4)

# Convert results to a data frame for clear display (optional)
go_sizes_df <- data.frame(
    GO_ID = names(go_term_sizes),
    Gene_Count = as.numeric(go_term_sizes)
)

# Add GO descriptions for readability (if the go_descriptions variable is available)
if (exists("go_descriptions")) {
    go_sizes_df$GO_Description <- go_descriptions[go_sizes_df$GO_ID]
}

# Print results
print(go_sizes_df)
