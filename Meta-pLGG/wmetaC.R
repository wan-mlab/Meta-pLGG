suppressMessages({
  library(ggplot2)
  library(cluster)
  library(clusterCrit)
  library(fpc)
  library(RColorBrewer)
  library(Rtsne)
  library(SHARP)
  library(umap)# ✅ wMetaC 在这里  :contentReference[oaicite:1]{index=1}
})

set.seed(123)

suppressMessages({
  library(readr)
})

# 0) 输入：三张 label 表
HCLUST_LABEL <- "/work/wanlab/bldx/glioma/0304gliomaAACR/manuscript/0601result/GMMNMFSpectral/0709results/5RP/hclustresults/all_clustering_labels.csv"
KMEANS_LABEL <- "/work/wanlab/bldx/glioma/0304gliomaAACR/manuscript/0601result/GMMNMFSpectral/0709results/5RP/kmeansresults/all_clustering_labels.csv"
SOM_LABEL    <- "/work/wanlab/bldx/glioma/0304gliomaAACR/manuscript/0601result/GMMNMFSpectral/0709results/5RP/somresults/all_clustering_labels.csv"
ALL_LABEL    <- "/work/wanlab/bldx/glioma/0304gliomaAACR/manuscript/0601result/GMMNMFSpectral/0709results/5RP/clusterresults_no11/all_clustering_labels_GMM_NMF_Spectral_K5.csv"
# 输出：给 wMetaC 用的输入文件
OUT_WMETAC_IN <- "/work/wanlab/bldx/glioma/0304gliomaAACR/manuscript/0601result/GMMNMFSpectral/0709results/5RP/wMetaC_input_6waycluster_no11.csv"

# 1) 读取 + 统一成：rownames = Sample
read_label_tbl <- function(path) {
  df0 <- read.csv(path, check.names = FALSE)
  
  # 去掉 write.csv(row.names=TRUE) 产生的空列名索引列（如果存在）
  if ("" %in% colnames(df0)) {
    df0 <- df0[, setdiff(colnames(df0), ""), drop = FALSE]
  }
  
  stopifnot("Sample" %in% colnames(df0))
  
  # ✅ 先保存真实 sample id
  ids <- as.character(df0$Sample)
  
  # 去掉 Sample，只保留 clustering 列
  df0 <- df0[, setdiff(colnames(df0), "Sample"), drop = FALSE]
  
  # 转整数标签
  df <- as.data.frame(lapply(df0, function(x) as.integer(as.factor(x))),
                      check.names = FALSE)
  
  # ✅ 关键：把 rownames 设回 BS_...
  rownames(df) <- ids
  
  return(df)
}

h_df <- read_label_tbl(HCLUST_LABEL)
k_df <- read_label_tbl(KMEANS_LABEL)
s_df <- read_label_tbl(SOM_LABEL)
a_df <- read_label_tbl(ALL_LABEL)

# 2) 对齐样本（取三者交集，保证顺序完全一致）
#common_samples <- Reduce(intersect, list(rownames(h_df), rownames(k_df)))
common_samples <- Reduce(intersect, list(rownames(h_df), rownames(k_df), rownames(s_df), rownames(a_df)))
common_samples <- sort(common_samples)

if (length(common_samples) == 0) stop("❌ 三个文件没有共同样本（rowname不一致）")

h_df <- h_df[common_samples, , drop = FALSE]
k_df <- k_df[common_samples, , drop = FALSE]
s_df <- s_df[common_samples, , drop = FALSE]
a_df <- a_df[common_samples, , drop = FALSE]

# # 3) 给列名加前缀，避免重名
# colnames(h_df) <- paste0("HCLUST__", colnames(h_df))
# colnames(k_df) <- paste0("KMEANS__", colnames(k_df))
# colnames(s_df) <- paste0("SOM__",    colnames(s_df))

# 4) 合并 & 导出（这就是 wMetaC 的 label_tbl）
#label_tbl_wmetac <- cbind(h_df, k_df)
label_tbl_wmetac <- cbind(h_df, k_df, s_df, a_df)
#label_tbl_wmetac <- k_df
write.csv(label_tbl_wmetac, OUT_WMETAC_IN, row.names = TRUE)
cat("📁 已生成 wMetaC 输入文件：\n", OUT_WMETAC_IN, "\n")
cat("📌 合并后维度：", dim(label_tbl_wmetac), " (samples x clusterings)\n")








# A) 路径 & 输入
dir <- "/work/wanlab/bldx/glioma/0304gliomaAACR/manuscript/0601result/GMMNMFSpectral/0709results/5RP"

# base clustering labels (每列一个方法/一次RP的label)
label_tbl <- read.csv(
  "/work/wanlab/bldx/glioma/0304gliomaAACR/manuscript/0601result/GMMNMFSpectral/0709results/5RP/wMetaC_input_6waycluster_no11.csv",
  row.names = 1, check.names = FALSE
)


TPM_FILE <- "/work/wanlab/bldx/glioma/0304gliomaAACR/manuscript/0601result/GMMNMFSpectral/0709results/remove11RP/tpm_filtered_50_removed11.csv"

OGdata <- as.matrix(read.csv(
  TPM_FILE,
  row.names = 1,
  check.names = FALSE
))

storage.mode(OGdata) <- "numeric"

cat("Raw TPM dimension:", dim(OGdata), "\n")


# 2) log2(TPM + 1)

OGdata <- log2(OGdata + 1)

# remove NA / Inf / NaN
OGdata[!is.finite(OGdata)] <- 0

cat("After log2(TPM+1):", dim(OGdata), "\n")
cat("NA:", sum(is.na(OGdata)), "\n")
cat("Inf:", sum(is.infinite(OGdata)), "\n")

# 3) Transpose: sample × gene

OGdata <- t(OGdata)

cat("After transpose:", dim(OGdata), "\n")

# 4) Align samples with label_tbl

common_samples <- intersect(rownames(label_tbl), rownames(OGdata))

cat("Samples in label_tbl:", nrow(label_tbl), "\n")
cat("Samples in OGdata:", nrow(OGdata), "\n")
cat("Matched samples:", length(common_samples), "\n")

label_tbl_aligned <- label_tbl[common_samples, , drop = FALSE]
OGdata_aligned <- OGdata[common_samples, , drop = FALSE]

stopifnot(all(rownames(label_tbl_aligned) == rownames(OGdata_aligned)))


# 5) Filter genes

# A. remove zero-variance genes
sds <- apply(OGdata_aligned, 2, sd, na.rm = TRUE)
keep_nonzero <- is.finite(sds) & sds > 0

OGdata_filtered <- OGdata_aligned[, keep_nonzero, drop = FALSE]

cat("Genes before zero-variance filter:", ncol(OGdata_aligned), "\n")
cat("Genes after zero-variance filter:", ncol(OGdata_filtered), "\n")

# B. select top variable genes
top_n_genes <- 2000

vars <- apply(OGdata_filtered, 2, var, na.rm = TRUE)
vars <- vars[is.finite(vars)]

top_genes <- names(sort(vars, decreasing = TRUE))[1:min(top_n_genes, length(vars))]

OGdata_filtered <- OGdata_filtered[, top_genes, drop = FALSE]

cat("Genes after top variable gene filter:", ncol(OGdata_filtered), "\n")


# 6) Scale genes

OG_scaled <- scale(OGdata_filtered)

# safety check
OG_scaled[!is.finite(OG_scaled)] <- 0

cat("Scaled data dimension:", dim(OG_scaled), "\n")
cat("Final NA:", sum(is.na(OG_scaled)), "\n")
cat("Final NaN:", sum(is.nan(OG_scaled)), "\n")
cat("Final Inf:", sum(is.infinite(OG_scaled)), "\n")

# 7) PCA

pca_eval <- prcomp(
  OG_scaled,
  center = FALSE,
  scale. = FALSE
)

# 8) Use top 20–50 PCs for silhouette / metrics

n_pc <- 10

n_pc <- min(n_pc, ncol(pca_eval$x))

OG_scaled <- pca_eval$x[, 1:n_pc, drop = FALSE]

cat("PCA evaluation space:", dim(OG_scaled), "\n")
cat("Using top PCs:", n_pc, "\n")

# 后面 tSNE / UMAP / PCA 都用 OG_scaled
n <- nrow(OG_scaled)
perp <- min(30, floor((n - 1) / 3) - 1)
if (is.na(perp) || perp < 5) perp <- max(2, floor((n - 1) / 3) - 1)
set.seed(123)
tsne_rp <- Rtsne(
  OG_scaled,
  dims = 2,
  perplexity = perp,
  pca = FALSE,
  check_duplicates = FALSE,
  verbose = FALSE
)
dataD <- data.frame(tsne1 = tsne_rp$Y[, 1], tsne2 = tsne_rp$Y[, 2])

# PCA 2D
pca_rp <- prcomp(OG_scaled, center = FALSE, scale. = FALSE)
dataPCA <- data.frame(
  pca1 = OG_scaled[, 1],
  pca2 = OG_scaled[, 2]
)

# UMAP 2D

set.seed(123)

# umap_rp <- umap(
#   OG_scaled,
#   n_neighbors = 15,
#   min_dist = 0.03,
#   spread = 1.5,
#   metric = "euclidean",
#   n_components = 2,
#   random_state = 123
# )
umap_rp <- umap(
  OG_scaled,
  n_neighbors = 30,
  min_dist = 0.1,
  spread = 2,
  metric = "euclidean",
  n_components = 2,
  random_state = 123
)

dataUMAP <- data.frame(
  umap1 = umap_rp$layout[, 1],
  umap2 = umap_rp$layout[, 2]
)
# B) 画图函数（用 tsne1/tsne2）
# plot_clusters <- function(compdata, filename, title) {
#   compdata$cluster <- as.character(compdata$cluster)
#   n_clusters <- length(unique(compdata$cluster[compdata$cluster != "0"]))
#   colors <- colorRampPalette(RColorBrewer::brewer.pal(12, "Set1"))(n_clusters)
#   
#   cluster_colors <- c(
#     "0" = "gray80",
#     setNames(colors, sort(unique(compdata$cluster[compdata$cluster != "0"])))
#   )
#   
#   p <- ggplot(compdata, aes(x = tsne1, y = tsne2, color = as.factor(cluster))) +
#     geom_point(size = 2) +
#     labs(title = title, color = "Cluster") +
#     scale_color_manual(values = cluster_colors) +
#     theme(
#       axis.text.x = element_text(size = 14, face = "bold", color = "black"),
#       axis.text.y = element_text(size = 14, face = "bold", color = "black"),
#       axis.title.x = element_text(size = 14, face = "bold", color = "black"),
#       axis.title.y = element_text(size = 14, face = "bold", color = "black"),
#       plot.title = element_text(size = 18, face = "bold", color = "black"),
#       legend.title = element_text(size = 14, face = "bold", color = "black"),
#       legend.text = element_text(size = 14, face = "bold", color = "black"),
#       panel.grid.major = element_blank(),
#       panel.grid.minor = element_blank(),
#       panel.background = element_rect(fill = "white"),
#       plot.background = element_rect(fill = "white")
#     )
#   
#   ggsave(filename = filename, plot = p, width = 8, height = 6, dpi = 300)
# }
plot_clusters <- function(compdata, filename, title, xcol, ycol) {
  compdata$cluster <- factor(compdata$cluster)
  
  p <- ggplot(compdata, aes_string(x = xcol, y = ycol, color = "cluster")) +
    geom_point(size = 2) +
    labs(title = title, color = "Cluster") +
    theme(
      axis.text.x = element_text(size = 14, face = "bold", color = "black"),
      axis.text.y = element_text(size = 14, face = "bold", color = "black"),
      axis.title.x = element_text(size = 14, face = "bold", color = "black"),
      axis.title.y = element_text(size = 14, face = "bold", color = "black"),
      plot.title = element_text(size = 18, face = "bold", color = "black"),
      legend.title = element_text(size = 14, face = "bold", color = "black"),
      legend.text = element_text(size = 14, face = "bold", color = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white")
    )
  
  ggsave(filename = filename, plot = p, width = 6, height = 6)
}
############################################################
# C) 指标函数（用 OGdata 距离算 Sil / Dunn / CH / I）
############################################################
compute_I_index <- function(data_for_metrics, clusters) {
  data <- as.matrix(data_for_metrics)
  clusters <- as.integer(clusters)
  u <- sort(unique(clusters))
  
  centroids <- do.call(rbind, lapply(u, function(k) colMeans(data[clusters == k, , drop = FALSE])))
  
  S <- sum(sapply(seq_along(u), function(i) {
    cl_data  <- data[clusters == u[i], , drop = FALSE]
    centroid <- centroids[i, ]
    sum(sqrt(rowSums((cl_data - centroid)^2)))
  }))
  
  D <- max(dist(centroids))
  D / S
}

eval_clustering_simple <- function(method, clusters, OGdata, sil_plot_file = NULL) {
  D <- dist(OGdata)
  
  sil <- silhouette(clusters, D)
  sil_mean <- mean(sil[, 3])
  
  if (!is.null(sil_plot_file)) {
    pdf(sil_plot_file, width = 8, height = 6)
    plot(sil, border = NA, main = paste("Silhouette -", method))
    dev.off()
  }
  
  dunn <- intCriteria(OGdata, as.integer(clusters), "Dunn")$dunn
  Iidx <- compute_I_index(OGdata, clusters)
  CH   <- calinhara(OGdata, as.integer(clusters))
  
  data.frame(method = method, Sil = sil_mean, Dunn = dunn, CH = CH, I = Iidx)
}

# 7) wMetaC meta-clustering：不同方法组合 top6/top5/top3/top2/top1

label_tbl_for_nc <- label_tbl_aligned

if ("Sample" %in% colnames(label_tbl_for_nc)) {
  label_tbl_for_nc <- label_tbl_for_nc[, setdiff(colnames(label_tbl_for_nc), "Sample"), drop = FALSE]
}

# A) 从列名中提取方法名

get_method_name <- function(x) {
  ifelse(
    grepl("__", x),
    sub("__.*$", "", x),
    sub("_RP[0-9]+$", "", x)
  )
}

method_each_col <- get_method_name(colnames(label_tbl_for_nc))

cat("Detected methods:\n")
print(table(method_each_col))


# B) 指定 6 个方法的顺序

method_order <- c(
  "KMeans",
  "SOM",
  "GMM",
  "Spectral",
  "Hierarchical5_no11_wardD2", 
#  "Leiden",
#  "Louvain",
  "NMF"
)

# 检查这些方法是否都在 label_tbl 中
missing_methods <- setdiff(method_order, unique(method_each_col))
if (length(missing_methods) > 0) {
  stop(
    "These methods are not found in label_tbl column names: ",
    paste(missing_methods, collapse = ", ")
  )
}

# C) 定义不同组合

method_combinations <- list(
#  top8 = method_order[1:8],
#  top7 = method_order[1:7],
  top6 = method_order[1:6],
  top5 = method_order[1:5],
  top4 = method_order[1:4],
  top3 = method_order[1:3],
  top2 = method_order[1:2],
  top1 = method_order[1:1]
)

wmetac_all_dir <- file.path(dir, "wMetaC_method_combinations")
if (!dir.exists(wmetac_all_dir)) dir.create(wmetac_all_dir, recursive = TRUE)

all_results <- data.frame()

# E) 循环运行不同组合

for (combo_name in names(method_combinations)) {
  
  cat("\n============================================================\n")
  cat("Running:", combo_name, "\n")
  cat("Methods:", paste(method_combinations[[combo_name]], collapse = ", "), "\n")
  cat("============================================================\n")
  
  selected_methods <- method_combinations[[combo_name]]
  
  keep_cols <- colnames(label_tbl_for_nc)[method_each_col %in% selected_methods]
  
  cat("Number of selected base clustering columns:", length(keep_cols), "\n")
  print(table(method_each_col[colnames(label_tbl_for_nc) %in% keep_cols]))
  
  if (length(keep_cols) == 0) {
    warning("No columns selected for ", combo_name)
    next
  }
  
  label_subset <- label_tbl_for_nc[, keep_cols, drop = FALSE]

  # Convert base clustering labels to integer matrix
  
  nC <- as.matrix(label_subset)
  nC <- apply(nC, 2, function(x) as.integer(as.factor(x)))
  
  cat("Base clustering matrix nC dimension:", dim(nC), "\n")
  
  # Run wMetaC
  
  res <- wMetaC(
    nC            = nC,
    hmethod       = "ward.D",
    enN.cluster   = NULL,
    minN.cluster  = 5,
    maxN.cluster  = 5,
    sil.thre      = 0,
    height.Ntimes = 1
  )
  
  final_meta_cluster <- as.integer(res$finalC)
  final_meta_cluster <- as.integer(as.factor(final_meta_cluster))
  
  cat("wMetaC cluster number:", length(unique(final_meta_cluster)), "\n")
  print(table(final_meta_cluster))
  
  # Check alignment
  
  stopifnot(nrow(label_tbl_aligned) == length(final_meta_cluster))
  stopifnot(nrow(OG_scaled) == length(final_meta_cluster))
  stopifnot(all(rownames(label_tbl_aligned) == rownames(OG_scaled)))
  
  # Save labels
  
  combo_dir <- file.path(wmetac_all_dir, combo_name)
  if (!dir.exists(combo_dir)) dir.create(combo_dir, recursive = TRUE)
  
  label_out <- label_tbl_aligned
  label_out$wMetaC <- final_meta_cluster
  
  write.csv(
    label_out,
    file.path(combo_dir, paste0("meta_clustering_labels_", combo_name, ".csv")),
    row.names = TRUE
  )
  
  # Evaluation
  
  sil_plot_file_meta <- file.path(combo_dir, paste0(combo_name, "_Silhouette.pdf"))
  
  resultdata <- eval_clustering_simple(
    method = paste0("wMetaC_", combo_name),
    clusters = final_meta_cluster,
    OGdata = OG_scaled,
    sil_plot_file = sil_plot_file_meta
  )
  
  resultdata$combination <- combo_name
  resultdata$n_methods <- length(selected_methods)
  resultdata$n_base_clusterings <- length(keep_cols)
  resultdata$methods_used <- paste(selected_methods, collapse = ";")
  resultdata$n_clusters <- length(unique(final_meta_cluster))
  
  all_results <- rbind(all_results, resultdata)
  
  # Plot t-SNE
  
  plot_clusters(
    compdata = data.frame(
      tsne1 = dataD$tsne1,
      tsne2 = dataD$tsne2,
      cluster = final_meta_cluster
    ),
    filename = file.path(combo_dir, paste0(combo_name, "_wMetaC_tsne_plot.pdf")),
    title    = paste0("wMetaC ", combo_name, " Clustering (t-SNE)"),
    xcol     = "tsne1",
    ycol     = "tsne2"
  )
  
  # Plot PCA
  
  plot_clusters(
    compdata = data.frame(
      pca1 = dataPCA$pca1,
      pca2 = dataPCA$pca2,
      cluster = final_meta_cluster
    ),
    filename = file.path(combo_dir, paste0(combo_name, "_wMetaC_pca_plot.pdf")),
    title    = paste0("wMetaC ", combo_name, " Clustering (PCA)"),
    xcol     = "pca1",
    ycol     = "pca2"
  )
  
  # Plot UMAP
  
  plot_clusters(
    compdata = data.frame(
      umap1 = dataUMAP$umap1,
      umap2 = dataUMAP$umap2,
      cluster = final_meta_cluster
    ),
    filename = file.path(combo_dir, paste0(combo_name, "_wMetaC_umap_plot.pdf")),
    title    = paste0("wMetaC ", combo_name, " Clustering (UMAP)"),
    xcol     = "umap1",
    ycol     = "umap2"
  )
  
  cat("Finished:", combo_name, "\n")
  print(resultdata)
}


# F) 保存所有组合的评估结果

all_results <- all_results[order(-all_results$Sil), ]

write.csv(
  all_results,
  file.path(wmetac_all_dir, "wMetaC_method_combination_metrics.csv"),
  row.names = FALSE
)

cat("\n🎉 All wMetaC combinations finished!\n")
cat("Summary metrics:\n")
print(all_results)