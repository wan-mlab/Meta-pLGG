suppressMessages({
  library(readr)
  library(ggplot2)
  library(cluster)
  library(clusterCrit)
  library(fpc)
  library(RColorBrewer)
  library(Rtsne)
  
  library(mclust)   # GMM
  library(NMF)      # NMF clustering
  library(kernlab)  # Spectral clustering
})

set.seed(123)

# 0) 路径配置

DIR_OUT <- "/work/wanlab/bldx/glioma/0304gliomaAACR/manuscript/0601result/GMMNMFSpectral/0709results/5RP/"
RESULTS_DIR <- file.path(DIR_OUT, "clusterresults_no11")
if (!dir.exists(RESULTS_DIR)) dir.create(RESULTS_DIR, recursive = TRUE)

# RP_DIR <- "/work/wanlab/bldx/glioma/1109cluster/abstract/immunecluster3"
RP_DIR <- "/work/wanlab/bldx/glioma/0304gliomaAACR/manuscript/0601result/GMMNMFSpectral/0709results/remove11RP/11_200/"
RP_PATTERN <- "random_projection_2000_%d.csv"   # 1..10
#RP_PATTERN <- "no11_random_projection2000_%d.csv"   # 1..10
RP_N <- 5

METHOD <- "GMM_Spectral_NMF_K5_no11"
#METHOD <- "NMF_K5_no11"
TARGET_K <- 5
METHODS <- c("GMM", "NMF", "Spectral")
#METHODS <- c("NMF")

TPM_FILE <- "/work/wanlab/bldx/glioma/0304gliomaAACR/manuscript/0601result/GMMNMFSpectral/0709results/remove11RP/tpm_filtered_50_removed11.csv"

# 先读第一个 RP 文件
rp_file_1 <- file.path(RP_DIR, sprintf(RP_PATTERN, 1))
rp_sample_ref <- as.matrix(read.csv(rp_file_1, row.names = 1, check.names = FALSE))
rp_samples <- rownames(rp_sample_ref)

cat("Reference RP samples:", length(rp_samples), "\n")

#   read raw TPM: gene × sample  
OG_raw <- as.matrix(read.csv(
  TPM_FILE,
  row.names = 1,
  check.names = FALSE
))

storage.mode(OG_raw) <- "numeric"

cat("Raw TPM dimension:", dim(OG_raw), "\n")

#   log2(TPM + 1)  
OG_raw <- log2(OG_raw + 1)
OG_raw[!is.finite(OG_raw)] <- 0

cat("After log2(TPM+1):", dim(OG_raw), "\n")
cat("NA:", sum(is.na(OG_raw)), "\n")
cat("Inf:", sum(is.infinite(OG_raw)), "\n")

#   transpose: sample × gene  
OG_raw <- t(OG_raw)

cat("After transpose:", dim(OG_raw), "\n")

#   align samples with RP samples  
common_samples <- intersect(rp_samples, rownames(OG_raw))

cat("Samples in RP:", length(rp_samples), "\n")
cat("Samples in OG_raw:", nrow(OG_raw), "\n")
cat("Matched samples:", length(common_samples), "\n")

# 按 RP 的样本顺序排列
common_samples <- rp_samples[rp_samples %in% common_samples]

OG_aligned <- OG_raw[common_samples, , drop = FALSE]

stopifnot(all(rownames(OG_aligned) == common_samples))

#   remove zero-variance genes  
sds <- apply(OG_aligned, 2, sd, na.rm = TRUE)
keep_nonzero <- is.finite(sds) & sds > 0

OG_filtered <- OG_aligned[, keep_nonzero, drop = FALSE]

cat("Genes before zero-variance filter:", ncol(OG_aligned), "\n")
cat("Genes after zero-variance filter:", ncol(OG_filtered), "\n")

#   select top variable genes  
top_n_genes <- 2000

vars <- apply(OG_filtered, 2, var, na.rm = TRUE)
vars <- vars[is.finite(vars)]

top_genes <- names(sort(vars, decreasing = TRUE))[1:min(top_n_genes, length(vars))]

OG_filtered <- OG_filtered[, top_genes, drop = FALSE]

cat("Genes after top variable gene filter:", ncol(OG_filtered), "\n")

#   scale genes  
OG_scaled <- scale(OG_filtered)
OG_scaled[!is.finite(OG_scaled)] <- 0

cat("Scaled data dimension:", dim(OG_scaled), "\n")
cat("Final NA:", sum(is.na(OG_scaled)), "\n")
cat("Final NaN:", sum(is.nan(OG_scaled)), "\n")
cat("Final Inf:", sum(is.infinite(OG_scaled)), "\n")

#   PCA  
pca_eval <- prcomp(
  OG_scaled,
  center = FALSE,
  scale. = FALSE
)

#   use top PCs for metrics  
n_pc <- 10
n_pc <- min(n_pc, ncol(pca_eval$x))

OGdata <- pca_eval$x[, 1:n_pc, drop = FALSE]

cat("PCA evaluation space:", dim(OGdata), "\n")
cat("Using top PCs:", n_pc, "\n")

# 2) 画图函数
plot_clusters <- function(compdata, filename, title) {
  compdata$cluster <- as.character(compdata$cluster)
  
  n_clusters <- length(unique(compdata$cluster[compdata$cluster != "0"]))
  colors <- colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(n_clusters)
  
  cluster_colors <- c(
    "0" = "gray80",
    setNames(colors, sort(unique(compdata$cluster[compdata$cluster != "0"])))
  )
  
  p <- ggplot(compdata, aes(x = tsne1, y = tsne2, color = as.factor(cluster))) +
    geom_point(size = 2) +
    labs(title = title, color = "Cluster") +
    scale_color_manual(values = cluster_colors) +
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
  
  ggsave(filename = filename, plot = p, width = 8, height = 6, dpi = 300)
}


# 3) 指标计算（用 OGdata）
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

eval_clustering <- function(method, clusters, OGdata, sil_plot_file = NULL) {
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

# 4) 循环跑 10 次 RP：

label_file   <- file.path(RESULTS_DIR, "all_clustering_labels_GMM_NMF_Spectral_K5.csv")
metrics_file <- file.path(RESULTS_DIR, "base10_metrics_GMM_NMF_Spectral_K5.csv")

label_tbl <- NULL
metrics_all <- list()
metric_counter <- 1

for (i in 1:RP_N) {
  
  rp_file <- file.path(RP_DIR, sprintf(RP_PATTERN, i))
  
  cat("\n==============================\n")
  cat("RP run:", i, "\n")
  cat("input:", rp_file, "\n")
  
  dataM <- as.matrix(read.csv(rp_file, row.names = 1, check.names = FALSE))
  storage.mode(dataM) <- "numeric"
  dataM[!is.finite(dataM)] <- 0
  dataM[is.na(dataM)] <- 0
  
  cat("dataM:", dim(dataM), "\n")
  
  #   align current RP data with OGdata metric space  
  common_run_samples <- intersect(rownames(dataM), rownames(OGdata))
  
  dataM <- dataM[common_run_samples, , drop = FALSE]
  OGdata_run <- OGdata[common_run_samples, , drop = FALSE]
  
  stopifnot(all(rownames(dataM) == rownames(OGdata_run)))
  
  cat("After alignment, dataM:", dim(dataM), "\n")
  cat("After alignment, OGdata_run:", dim(OGdata_run), "\n")
  
  #   t-SNE 仅用于可视化  
  dataM_scaled <- scale(dataM)
  dataM_scaled[!is.finite(dataM_scaled)] <- 0
  
  n <- nrow(dataM_scaled)
  perp <- min(30, floor((n - 1) / 3) - 1)
  if (is.na(perp) || perp < 5) {
    perp <- max(2, floor((n - 1) / 3) - 1)
  }
  
  set.seed(123)
  tsne_rp <- Rtsne(
    dataM_scaled,
    dims = 2,
    perplexity = perp,
    pca = TRUE,
    check_duplicates = FALSE,
    verbose = FALSE
  )
  
  dataD <- data.frame(
    tsne1 = tsne_rp$Y[, 1],
    tsne2 = tsne_rp$Y[, 2]
  )
  
  #   初始化 label table  
  if (is.null(label_tbl)) {
    label_tbl <- data.frame(Sample = rownames(dataM), stringsAsFactors = FALSE)
    rownames(label_tbl) <- label_tbl$Sample
  } else {
    label_tbl <- label_tbl[rownames(dataM), , drop = FALSE]
  }
  
  # 每个 RP 内部跑 3 种聚类方法
  for (method_name in METHODS) {
    
    cat("\n--- Method:", method_name, "RP", i, "---\n")
    
    set.seed(123)
    
    # 1) GMM clustering
    if (method_name == "GMM") {
      
      gmm_fit <- Mclust(
        data = dataM,
        G = TARGET_K,
        verbose = FALSE
      )
      
      clusters <- gmm_fit$classification
    }
    
    # 2) NMF clustering
    if (method_name == "NMF") {
      
      # NMF 要求输入矩阵非负
      # RP 数据可能有负值，所以整体平移到非负
      dataM_nmf <- dataM
      min_val <- min(dataM_nmf, na.rm = TRUE)
      
      if (min_val < 0) {
        dataM_nmf <- dataM_nmf - min_val
      }
      
      dataM_nmf <- dataM_nmf + 1e-8
      dataM_nmf[!is.finite(dataM_nmf)] <- 0
      
      # NMF 包通常要求：feature × sample
      # 你的 dataM 是 sample × feature，所以这里转置
      nmf_input <- t(dataM_nmf)
      
      nmf_fit <- nmf(
        nmf_input,
        rank = TARGET_K,
        method = "brunet",
        nrun = 1,
        seed = 123
      )
      
      # coef(nmf_fit): rank × sample
      # 每个 sample 分到 coefficient 最大的那个 component
      H <- coef(nmf_fit)
      clusters <- apply(H, 2, which.max)
      
      # 确保顺序和 dataM 的样本顺序一致
      names(clusters) <- colnames(H)
      clusters <- clusters[rownames(dataM)]
    }
    
    # 3) Spectral clustering
    if (method_name == "Spectral") {
      
      spectral_fit <- specc(
        x = as.matrix(dataM),
        centers = TARGET_K
      )
      
      clusters <- as.integer(spectral_fit)
    }
    
    clusters <- as.integer(clusters)
    
    cat("Cluster table:\n")
    print(table(clusters))
    
    #   输出目录：每个 RP + 每个 method 一个文件夹  
    run_dir <- file.path(
      DIR_OUT,
      sprintf("RP%02d_%s_K%d", i, method_name, TARGET_K)
    )
    
    if (!dir.exists(run_dir)) {
      dir.create(run_dir, recursive = TRUE)
    }
    
    #   画 t-SNE 图  
    plot_clusters(
      compdata = data.frame(
        tsne1 = dataD$tsne1,
        tsne2 = dataD$tsne2,
        cluster = clusters
      ),
      filename = file.path(run_dir, paste0(method_name, "_tSNE.png")),
      title = sprintf("%s clustering, K=%d, RP%02d t-SNE", method_name, TARGET_K, i)
    )
    
    #   silhouette plot  
    sil_pdf <- file.path(run_dir, paste0(method_name, "_Silhouette.pdf"))
    
    #   指标计算：用 OGdata_run  
    met <- eval_clustering(
      method = sprintf("%s_RP%02d", method_name, i),
      clusters = clusters,
      OGdata = OGdata_run,
      sil_plot_file = sil_pdf
    )
    
    metrics_all[[metric_counter]] <- met
    metric_counter <- metric_counter + 1
    
    #   保存 label  
    label_colname <- sprintf("%s_RP%02d", method_name, i)
    label_tbl[[label_colname]] <- clusters
  }
}

# 保存 labels
write.csv(label_tbl, label_file, row.names = TRUE)
cat("\n已保存 GMM/NMF/Spectral labels:\n", label_file, "\n")

# 保存 metrics
metrics_df <- do.call(rbind, metrics_all)
write.csv(metrics_df, metrics_file, row.names = FALSE)
cat("已保存 GMM/NMF/Spectral metrics:\n", metrics_file, "\n")

cat("\nDone!\n")