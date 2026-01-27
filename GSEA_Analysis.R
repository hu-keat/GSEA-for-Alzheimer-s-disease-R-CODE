---
title: "GSEA_Analysis.R"
author: "HU QIANG"
date: "2025-10-14"
output:
  html_document: default
  pdf_document: default
---

```{r setup, include=FALSE}
# 使用 knitr 选项设置全局根目录
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
## Optional package installation commands

# install.packages("readr")

# install.packages("devtools")

# install.packages("compositions")

# BiocManager::install("polyester")

# BiocManager::install("edgeR", force = T)

# devtools::install_github("Divo-Lee/SIEVE")

# install.packages("gamlss")

# BiocManager::install("cqn")

# devtools::install_github("zjdaye/MDSeq")

# install.packages("vioplot")

# install.packages("VennDiagram")

# install.packages("gamlss")

# install.packages("gridExtra")

# install.packages("httr")

# install.packages("jsonlite")

# install.packages("MASS")

# install.packages("ggplot2")

```

```{r}
# Load required libraries
# 确保加载了 GSEA 和注释包
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

# library(MDSeq)

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
# 取 SIEVE 原始结果
res_de <- clrSIEVE_result$clrDE_test
res_dv <- clrSIEVE_result$clrDV_test
res_ds <- clrSIEVE_result$clrDS_test

# 确保三个表的基因完全对齐
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
# --- 核心：重新生成 U 统计量 ---
# 规则：U = 在 (DE, DV, DS) 中绝对值最大的那个，并保留原符号
# tie 规则：如果出现并列最大，取第一个（DE > DV > DS）
df_stat$U1 <- apply(df_stat[, c("DE","DV","DS")], 1, function(x) {
  m <- max(abs(x), na.rm = TRUE)
  idx <- which(abs(x) == m)[1]
  x[idx]
})

df_stat$U2 <- with(df_stat, sign(DE) * sqrt(DE^2 + (3/4)*DV^2 + (1/2)*DS^2))


df_stat$U3 <- apply(df_stat[, c("DE","DV","DS")], 1, function(x) {
  k <- which.max(abs(x))                       # 主导分量（未加权 abs 最大）
  M <- sqrt(x[1]^2 + (3/4)*x[2]^2 + (1/2)*x[3]^2)  # 加权幅度
  sign(x[k]) * M                               # U3
})


# ---- 非常重要：ID 统一（处理 ENSG...“.12” 这种版本号） ----
strip_ens_version <- function(x) sub("\\..*$", "", x)

rownames(df_stat) <- strip_ens_version(rownames(df_stat))

# 构造排序列表（stats 必须是 named numeric vector）
de_ranks_sorted <- sort(na.omit(setNames(df_stat$DE, rownames(df_stat))), decreasing = TRUE)
dv_ranks_sorted <- sort(na.omit(setNames(df_stat$DV, rownames(df_stat))), decreasing = TRUE)
ds_ranks_sorted <- sort(na.omit(setNames(df_stat$DS, rownames(df_stat))), decreasing = TRUE)
u1_ranks_sorted  <- sort(na.omit(setNames(df_stat$U1,  rownames(df_stat))), decreasing = TRUE)
u2_ranks_sorted  <- sort(na.omit(setNames(df_stat$U2,  rownames(df_stat))), decreasing = TRUE)
u3_ranks_sorted  <- sort(na.omit(setNames(df_stat$U3,  rownames(df_stat))), decreasing = TRUE)

message(paste("DE 排序列表长度:", length(de_ranks_sorted)))
message(paste("DV 排序列表长度:", length(dv_ranks_sorted)))
message(paste("DS 排序列表长度:", length(ds_ranks_sorted)))
message(paste("U1  排序列表长度:", length(u1_ranks_sorted)))
message(paste("U2  排序列表长度:", length(u2_ranks_sorted)))
message(paste("U3  排序列表长度:", length(u3_ranks_sorted)))
```

```{r}
# U1 = sgn(W_k,i) * max( |W_DE,i|, |W_DV,i|, |W_DS,i| )
check_table1 <- df_stat %>%
  as.data.frame() %>%
  # 增加绝对值辅助列，方便一眼看出对比
  mutate(
    abs_DE = abs(DE),
    abs_DV = abs(DV),
    abs_DS = abs(DS)
  ) %>%
  # 记录每一行 U1 的来源（是哪一列贡献了最大绝对值）
  mutate(
    Source = apply(.[, c("DE", "DV", "DS")], 1, function(x) {
      names(x)[which.max(abs(x))]
    }),
    # ===== 新增：U1 =====
    U1 = apply(.[, c("DE", "DV", "DS")], 1, function(x) {
      k <- which.max(abs(x))          # 找到主导分量
      sign(x[k]) * abs(x[k])          # sgn(W_k,i) * |W_k,i|
    })
  ) %>%
  # 只取前 20 行进行展示
  head(20)

```


```{r}
#U2 = sgn(W_DE,i) * sqrt( W_DE,i^2 + (3/4)*W_DV,i^2 + (1/2)*W_DS,i^2 )
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
#Mi = sqrt( 1·W_DE,i^2 + (3/4)·W_DV,i^2 + (1/2)·W_DS,i^2 )
#Si = sgn(W_k,i),  k = arg max {|W_DE,i|, |W_DV,i|, |W_DS,i|}
#U3 = Si · Mi

check_table3 <- df_stat %>%
  as.data.frame() %>%
  mutate(
    # 根号内的加权平方贡献
    w_DE = 1 * (DE^2),
    w_DV = (3/4) * (DV^2),
    w_DS = (1/2) * (DS^2),
# 加权幅度 M_i
    M_w = sqrt(w_DE + w_DV + w_DS),
# 主导分量（按未加权 abs 比较：|DE|,|DV|,|DS|）
    Source_dom = apply(.[, c("DE", "DV", "DS")], 1, function(x) {
      names(x)[which.max(abs(x))]
    }),
# 符号 S_i：取主导分量的符号
    S_dom = case_when(
      Source_dom == "DE" ~ sign(DE),
      Source_dom == "DV" ~ sign(DV),
      TRUE               ~ sign(DS)
    ),
# 图中统计量：U_i = S_i * M_i
    U3 = S_dom * M_w
  ) %>%
  head(20)

```

```{r}
# --- 更新后的 GSEA 排序列表准备 ---

# 定义一个便捷函数处理命名向量
get_sorted_ranks <- function(val_vec, names_vec) {
  res <- setNames(val_vec, names_vec)
  return(sort(na.omit(res), decreasing = TRUE))
}

# 构造包含 DE, DV, DS 和新指标 U 的排序列表
de_ranks_sorted <- get_sorted_ranks(df_stat$DE, rownames(df_stat))
dv_ranks_sorted <- get_sorted_ranks(df_stat$DV, rownames(df_stat))
ds_ranks_sorted <- get_sorted_ranks(df_stat$DS, rownames(df_stat))
u1_ranks_sorted  <- get_sorted_ranks(df_stat$U1,  rownames(df_stat))
u2_ranks_sorted  <- get_sorted_ranks(df_stat$U2,  rownames(df_stat))
u3_ranks_sorted  <- get_sorted_ranks(df_stat$U3,  rownames(df_stat))
```

```{r}
# --- 目标 GO IDs ---
go_ids_table4 <- c(
  "GO:0034613", "GO:0072599", "GO:0090150", "GO:0045047", "GO:0006614",
  "GO:0022610", "GO:0007155", "GO:0010647",
  "GO:0043062", "GO:0030198",
  "GO:0001568", "GO:0042060",
  "GO:0006412", "GO:0006401",
  "GO:0019884", "GO:0045619",
  "GO:0019363", "GO:0019752"
)

# --- 1) 主方法：select 直接 GOALL -> ENSEMBL ---
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
message("主方法拿到的 GO 基因集数: ", length(gene_sets_table4))
message("主方法缺失 GO 数: ", length(missing_gos))
print(missing_gos)

# --- 2) 补救方法：对缺失 GO 用 org.Hs.egGO2ALLEGS (GO -> Entrez -> ENSEMBL) ---
if (length(missing_gos) > 0) {
  # org.Hs.egGO2ALLEGS 是一个 environment：key=GO，value=Entrez 列表
  # 有些 GO 在这里能拿到 Entrez，再映射回 ENSEMBL
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

# --- 3) 最终检查 ---
missing_gos_final <- setdiff(go_ids_table4, names(gene_sets_table4))
message("最终成功 GO 基因集数: ", length(gene_sets_table4))
message("最终仍缺失 GO 数: ", length(missing_gos_final))
print(missing_gos_final)

# （可选）查看每个 GO 集合大小
print(sort(lengths(gene_sets_table4)))

```





```{r}
# GSEA 核心计算 (DE, DV, DS 的 Permutation Test)
N_PERM <- 10000   #定义置换次数。 (1000 次是标准，保证了统计显著性的可靠性)。
set.seed(421)
# 使用 suppressWarnings() 包装整个 fgsea 调用，以消除 'package:stats' 警告
suppressWarnings({
  fgsea_DE_results <- fgsea(pathways = gene_sets_table4, stats = de_ranks_sorted, nPermSimple = N_PERM, eps = 0)
  fgsea_DV_results <- fgsea(pathways = gene_sets_table4, stats = dv_ranks_sorted, nPermSimple = N_PERM, eps = 0)
  fgsea_DS_results <- fgsea(pathways = gene_sets_table4, stats = ds_ranks_sorted, nPermSimple = N_PERM, eps = 0)
  fgsea_U1_results  <- fgsea(pathways = gene_sets_table4, stats = u1_ranks_sorted,  nPermSimple = N_PERM, eps = 0)
  fgsea_U2_results  <- fgsea(pathways = gene_sets_table4, stats = u2_ranks_sorted,  nPermSimple = N_PERM, eps = 0)
  fgsea_U3_results  <- fgsea(pathways = gene_sets_table4, stats = u3_ranks_sorted,  nPermSimple = N_PERM, eps = 0)
})
#nPermSimple = N_PERM, eps = 0 GSEA 参数。 设定置换次数，并将 eps 设为 0 以获得更精确的低 P 值估计。
#eps (Epsilon)P值估计的最小精度 10^-10提高低 P 值的准确性
```

```{r}
results_DE_summary <- fgsea_DE_results %>% dplyr::mutate(Type = "DE") %>% dplyr::select(pathway, ES, NES, padj, Type, size)
results_DV_summary <- fgsea_DV_results %>% dplyr::mutate(Type = "DV") %>% dplyr::select(pathway, ES, NES, padj, Type, size)
results_DS_summary <- fgsea_DS_results %>% dplyr::mutate(Type = "DS") %>% dplyr::select(pathway, ES, NES, padj, Type, size)
results_U1_summary <- fgsea_U1_results %>% dplyr::mutate(Type = "U1") %>% dplyr::select(pathway, ES, NES, padj, Type, size)
results_U2_summary <- fgsea_U2_results %>% dplyr::mutate(Type = "U2") %>% dplyr::select(pathway, ES, NES, padj, Type, size)
results_U3_summary <- fgsea_U3_results %>% dplyr::mutate(Type = "U3") %>% dplyr::select(pathway, ES, NES, padj, Type, size)

all_gsea_results <- dplyr::bind_rows(
  results_DE_summary, 
  results_DV_summary, 
  results_DS_summary,
  results_U1_summary, 
  results_U2_summary, 
  results_U3_summary
) %>%
  # 在这里添加 -log10(padj)
  dplyr::mutate(neg_log10_padj = -log10(padj))

# (可选) 如果 padj 中存在 0，导致 -log10 变成无穷大，可以进行微调：
# all_gsea_results <- all_gsea_results %>%
#   dplyr::mutate(neg_log10_padj = -log10(ifelse(padj == 0, 1e-10, padj)))
```



```{r}

results_DE_summary <- fgsea_DE_results %>% dplyr::mutate(Type = "DE") %>% dplyr::select(pathway, ES, NES, padj, Type, size)
results_DV_summary <- fgsea_DV_results %>% dplyr::mutate(Type = "DV") %>% dplyr::select(pathway, ES, NES, padj, Type, size)
results_DS_summary <- fgsea_DS_results %>% dplyr::mutate(Type = "DS") %>% dplyr::select(pathway, ES, NES, padj, Type, size)
results_U1_summary <- fgsea_U1_results %>% dplyr::mutate(Type = "U1") %>% dplyr::select(pathway, ES, NES, padj, Type, size)
results_U2_summary <- fgsea_U2_results %>% dplyr::mutate(Type = "U2") %>% dplyr::select(pathway, ES, NES, padj, Type, size)
results_U3_summary <- fgsea_U3_results %>% dplyr::mutate(Type = "U3") %>% dplyr::select(pathway, ES, NES, padj, Type, size)

all_gsea_results <- dplyr::bind_rows(results_DE_summary, results_DV_summary, results_DS_summary,results_U1_summary, results_U2_summary, results_U3_summary)

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

plot_data <- all_gsea_results %>%
dplyr::mutate(Term_Desc = go_descriptions[pathway]) %>%
dplyr::filter(padj < 0.2) %>%
dplyr::arrange(desc(NES)) %>%
dplyr::mutate(Term_Desc = factor(Term_Desc, levels = unique(Term_Desc)))
```


```{r}
if (nrow(plot_data) > 0) {
  message(paste("绘制气泡图: 共显示", nrow(plot_data), "个富集结果 (FDR < 0.2)。"))
  
  gsea_dot_plot <- ggplot(plot_data, aes(x = Type, y = Term_Desc)) +
    geom_point(aes(size = size, color = NES)) + 
    scale_color_gradient2(low = "green", mid = "white", high = "yellow", midpoint = 0, name = "NES") +
    # 这个 geom_point 用于标记显著结果
    geom_point(data = dplyr::filter(plot_data, padj < 0.05), shape = 8, size = 1.5, color = "black") + 
    
    # labs 和 theme 应该连接在最后
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
  message("没有显著的 GSEA 富集结果 (FDR < 0.2)。")
}
```

```{r}
# =========================
# 打印所有显著通路及其 FDR
# =========================

sig_table <- plot_data %>%
  dplyr::filter(padj < 0.05) %>%
  dplyr::select(Type, pathway, Term_Desc, NES, padj, size) %>%
  dplyr::arrange(padj)

cat("\n================ 显著富集通路 (FDR < 0.05) ================\n")
if (nrow(sig_table) == 0) {
  cat("⚠️ 没有任何通路满足 FDR < 0.05。\n")
} else {
  print(sig_table, n = 200)
}
cat("===========================================================\n\n")

```


```{r}

selected_go_id <- "GO:0006412"
analysis_types <- c("DE","DV","DS","U1","U2","U3")

if (selected_go_id %in% names(gene_sets_table4)) {

  for (analysis_type in analysis_types) {

    # --- 动态选择 GSEA 结果和排序列表 ---
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
      stop("Analysis type not recognized. Use one of: DE, DV, DS, U1,U2,U3.")
    }

    # 确保选中的排序列表和结果变量存在
    if (exists("gene_rank_list") && exists("fgsea_results")) {

      # 提取和格式化 GSEA 结果信息
      pathway_info <- fgsea_results %>%
        dplyr::filter(pathway == selected_go_id) %>%
        dplyr::mutate(
          formatted_padj = format(padj, digits = 2, scientific = TRUE),
          NES_FDR = paste0("NES=", round(NES, 2), ", FDR=", formatted_padj)
        )

      # 如果当前 analysis_type 下没有这个通路结果，跳过
      if (nrow(pathway_info) == 0) {
        message("⚠️ 跳过：", selected_go_id, " 在 ", analysis_type, " 下没有 fgsea 结果行。")
        next
      }

      # 构造图标题
      plot_title <- paste0(
        "GSEA Enrichment Plot for ", "\n",
        go_descriptions[selected_go_id], "\n",
        " (", selected_go_id, ") in ", analysis_type, "\n",
        pathway_info$NES_FDR
      )

      # 画图
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

      # （可选）在控制台打印 NES/padj，方便你记录
      cat("\n", analysis_type, ": ",
          "NES=", round(pathway_info$NES[1], 3),
          ", padj(FDR)=", format(pathway_info$padj[1], scientific = TRUE, digits = 3),
          "\n", sep = "")

    } else {
      message("GSEA 曲线图跳过: 排序列表或 GSEA 结果变量缺失。")
    }
  }

} else {
  message(paste0("GSEA 曲线图跳过: 选定的 GO Term (", selected_go_id, ") 基因集在数据库中缺失。"))
}

```


```{r}

selected_go_id <- "GO:0006401" # <-- 保持通路不变
analysis_type <- "U2"          # <-- **修正：切换到 DV 模式**

if (selected_go_id %in% names(gene_sets_table4)) {

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
  stop("Analysis type not recognized. Use one of: DE, DV, DS, U1,U2,U3.")
}

 if (exists("gene_rank_list") && exists("fgsea_results")) {

 pathway_info <- fgsea_results %>%
 dplyr::filter(pathway == selected_go_id) %>%
 dplyr::mutate(
 formatted_padj = format(padj, digits = 2, scientific = TRUE),
 NES_FDR = paste0("NES=", round(NES, 2), ", FDR=", formatted_padj)
)

 plot_title <- paste0("GSEA Enrichment Plot for ", "\n",go_descriptions[selected_go_id],"\n",
 " (", selected_go_id, ") in ", analysis_type, "\n", pathway_info$NES_FDR)

 enrichment_plot <- fgsea::plotEnrichment(
 pathway = gene_sets_table4[[selected_go_id]],
stats = gene_rank_list # <-- 关键：使用动态变量
 ) +
labs(title = plot_title, 
x = paste0(analysis_type, " Gene List Rank"), 
y = "Enrichment Score (ES)") +
theme_minimal(base_size = 14) +
theme(plot.title = element_text(hjust = 0.5))

print(enrichment_plot)
} else {
 message("GSEA 曲线图跳过: 排序列表或 GSEA 结果变量缺失。")
}
} else {
 message(paste0("GSEA 曲线图跳过: 选定的 GO Term (", selected_go_id, ") 基因集在数据库中缺失。"))
}
```

```{r}

go_term_sizes <- lengths(gene_sets_table4)

go_sizes_df <- data.frame(
    GO_ID = names(go_term_sizes),
    Gene_Count = as.numeric(go_term_sizes)
)

if (exists("go_descriptions")) {
    go_sizes_df$GO_Description <- go_descriptions[go_sizes_df$GO_ID]
}

print(go_sizes_df)
```

