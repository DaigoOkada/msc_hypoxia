#Rhelixa scRNA-seq
#source("/Users/dokada/Dropbox/analysis/2025.4/RH_hypoxia0126_fin3.R") #2.8
out_path <- "/Users/dokada/Desktop/work/RH_hypoxia0126_fin3/"
if(file.exists(out_path)){
    unlink(out_path, recursive=TRUE)
    dir.create(out_path)
}else{
    dir.create(out_path)
}


#load
library(Seurat)
library(ggplot2)
library(clusterProfiler)
library(msigdbr)
library(DoubletFinder)
library(tidyr)
library(ggpubr)

#RNA
set.seed(1234)
data_path <- "/Users/dokada/Desktop/work/shizui_data/RH_hypoxia0116/Aggregation/filtered_feature_bc_matrix/"
sc.data <- Read10X(data.dir = data_path)
obj_rh <- CreateSeuratObject(counts = sc.data, project = "hypoxia")
donor <- sapply(colnames(obj_rh), function(chr) strsplit(chr, "-")[[1]][2])
#donor_id <- c("1"="DP1_P8_21", "2"="DP1_P8_3", "3"="DP1_P8_3_H2O2","4"="DP1_P26_21") #レリクサレポートの細胞数との一致を確認
donor_id <- c("1"="P8_norm", "2"="P8_hypo", "3"="P8_hypo_ROS","4"="P26_norm") #レリクサレポートの細胞数との一致を確認
donor[donor == 1] <- donor_id["1"]
donor[donor == 2] <- donor_id["2"]
donor[donor == 3] <- donor_id["3"]
donor[donor == 4] <- donor_id["4"]
obj_rh$donor <- donor
obj_rh[["percent.mt"]] <- PercentageFeatureSet(obj_rh, pattern = "^MT-")
features <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
for(i in 1:length(features)){
    png(paste0(out_path, "NO_QC_violin_", features[i], ".png"),width = 960, height = 960, res = 150)
    p <- VlnPlot(
        obj_rh,
        features = features[i],
        group.by = "donor",
        pt.size = 0
    ) +
    stat_summary(
        fun = median,
        geom = "point",
        size = 4,
        color = "black"
    )  + theme(
        #axis.title.x = element_blank(),
        axis.text.x = element_text(size = 28, margin = margin(t = 10)),
        axis.text.y = element_text(size = 28, margin = margin(r = 10)),
        axis.title.x = element_text(size = 28, margin = margin(t = 15)),
        axis.title.y = element_text(size = 28, margin = margin(r = 15)),
        plot.title = element_text(size = 28, hjust = 0.5, margin = margin(b = 20)),
        plot.margin = margin(t = 30, r = 30, b = 30, l = 30),
        legend.text  = element_text(size = 16),
        legend.title = element_text(size = 18),
        legend.position = "none"
    ) +
    xlab("")
    print(p)
    dev.off()

}
obj_rh <- subset(obj_rh, subset = nFeature_RNA > 200 & percent.mt < 10)
obj_rh <- NormalizeData(obj_rh, normalization.method = "LogNormalize", scale.factor = 10000)
obj_rh <- FindVariableFeatures(obj_rh, selection.method = "vst", nfeatures = 2000)
#cell cycle regression
s.genes   <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
s.genes.use   <- intersect(s.genes, rownames(obj_rh))
g2m.genes.use <- intersect(g2m.genes, rownames(obj_rh))
obj_rh <- CellCycleScoring(obj_rh, s.features = s.genes.use, g2m.features = g2m.genes.use, set.ident = FALSE)
obj_rh <- ScaleData(obj_rh, features = VariableFeatures(obj_rh))
#obj_rh <- ScaleData(obj_rh, features = VariableFeatures(obj_rh), vars.to.regress = c("S.Score","G2M.Score"))
obj_rh <- RunPCA(obj_rh, features = VariableFeatures(obj_rh), npcs = 50, verbose = FALSE)
obj_rh <- FindNeighbors(obj_rh, dims = 1:50)
obj_rh <- FindClusters(obj_rh, resolution = 0.4) #JC論文と同じ
obj_rh <- RunUMAP(obj_rh, dims = 1:50)
p <- DimPlot(obj_rh, reduction="umap", group.by="seurat_clusters", label=TRUE, label.size=10, pt.size=2) +
  theme(
    #plot.title = element_blank(),
    plot.title = element_text(size=40, hjust=0.5, face="bold"),
    axis.text = element_text(size=34),
    axis.title = element_text(size=38),
    legend.text = element_text(size=32),
    legend.title = element_text(size=34)
  ) +
  ggtitle("All") +
  xlab("UMAP1") +
  ylab("UMAP2")
ggsave(paste0(out_path, "umap_clusters.png"), p, width=12, height=10, dpi=300)
saveRDS(obj_rh, file = paste0(out_path, "obj_rh.rds"))

#DonorごとのUMAP
donors <- sort(unique(obj_rh$donor))
for(d in donors){
  cells_use <- colnames(obj_rh)[obj_rh$donor == d]
  if(length(cells_use) == 0) next

  obj_sub <- subset(obj_rh, cells = cells_use)

  p <- DimPlot(obj_sub, reduction="umap", group.by="seurat_clusters", label=TRUE, label.size=12, pt.size=2) +
    ggtitle(d) +
    xlab("UMAP1") +
    ylab("UMAP2") +
    theme(
      plot.title = element_text(size=40, hjust=0.5, face="bold"),
      axis.text = element_text(size=34),
      axis.title = element_text(size=38),
      legend.text = element_text(size=32),
      legend.title = element_text(size=34)
    )

  ggsave(paste0(out_path, "UMAP_", d, ".png"), plot=p, width=12, height=10, dpi=300)
}



# クロス集計（細胞数）
tab <- table(obj_rh$donor, obj_rh$seurat_clusters)
prop_tab <- prop.table(tab, margin = 1)
write.csv(tab, file = paste0(out_path, "RH_hypoxia_cluster_table.csv"))

#QC feature plots
obj_rh$cc_score <- pmax(obj_rh$S.Score, obj_rh$G2M.Score)
features <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "S.Score", "G2M.Score", "cc_score")
for(i in 1:length(features)){
    png(paste0(out_path, "QC_violin_", features[i], ".png"),width = 960, height = 960, res = 150)
    p <- VlnPlot(
        obj_rh,
        features = features[i],
        group.by = "donor",
        pt.size = 0
    ) +
    stat_summary(
        fun = median,
        geom = "point",
        size = 4,
        color = "black"
    )  + theme(
        axis.text.x = element_text(size = 28, margin = margin(t = 10)),
        axis.text.y = element_text(size = 28, margin = margin(r = 10)),
        axis.title.x = element_text(size = 28, margin = margin(t = 15)),
        axis.title.y = element_text(size = 28, margin = margin(r = 15)),
        plot.title = element_text(size = 28, hjust = 0.5, margin = margin(b = 20)),
        plot.margin = margin(t = 30, r = 30, b = 30, l = 30),
        legend.text  = element_text(size = 16),
        legend.title = element_text(size = 18),
        legend.position = "none"
    ) +
    xlab("")
    print(p)
    dev.off()

}


#analysis
use_dims   <- 1:50
n_boot     <- 1000
min_cells  <- 50
s_hyp   <- "P8_hypo"      # 低酸素
s_norm  <- "P8_norm"     # 定常酸素
s_highP <- "P26_norm"    # 定常酸素・高継代
s_hyph2o2 <- "P8_hypo_ROS" # 低酸素+H2O2
pc <- Embeddings(obj_rh, reduction = "pca")[, use_dims, drop = FALSE]

#grouping
cells_hyp   <- colnames(obj_rh)[obj_rh$donor == s_hyp]
cells_norm  <- colnames(obj_rh)[obj_rh$donor == s_norm]
cells_highP <- colnames(obj_rh)[obj_rh$donor == s_highP]
cells_hyph2o2 <- colnames(obj_rh)[obj_rh$donor == s_hyph2o2]

# bootstrap
centroid_dist_boot <- function(cells_a, cells_b, pc, n_boot = 2000, seed = 1) {
  set.seed(seed)
  na <- length(cells_a); nb <- length(cells_b)

  cen_a <- colMeans(pc[cells_a, , drop = FALSE])
  cen_b <- colMeans(pc[cells_b, , drop = FALSE])
  obs <- sqrt(sum((cen_a - cen_b)^2))

  boot <- numeric(n_boot)
  for (i in seq_len(n_boot)) {
    sa <- sample(cells_a, size = na, replace = TRUE)
    sb <- sample(cells_b, size = nb, replace = TRUE)
    ca <- colMeans(pc[sa, , drop = FALSE])
    cb <- colMeans(pc[sb, , drop = FALSE])
    boot[i] <- sqrt(sum((ca - cb)^2))
  }
  ci <- quantile(boot, c(0.025, 0.975), names = FALSE)
  list(obs = obs, boot = boot, ci_low = ci[1], ci_high = ci[2], n_a = na, n_b = nb)
}

# Calculate distance
r_hyp  <- centroid_dist_boot(cells_hyp,  cells_highP, pc, n_boot = n_boot, seed = 1)
r_norm <- centroid_dist_boot(cells_norm, cells_highP, pc, n_boot = n_boot, seed = 2)
r_hyp_h2o2 <- centroid_dist_boot(cells_hyph2o2, cells_highP, pc, n_boot = n_boot, seed = 3)

# delta = d(Hypoxia, HighP) - d(Normoxia, HighP)
delta_obs  <- r_hyp$obs - r_norm$obs
delta_boot <- r_hyp$boot - r_norm$boot
delta_obs2 <- r_hyp_h2o2$obs - r_hyp$obs
delta_boot2 <- r_hyp_h2o2$boot - r_hyp$boot
p_delta_lt0 <- mean(delta_boot < 0)

png(paste0(out_path, "delta_distance_hyp_vs_norm_bootstrap.png"), width = 960, height = 960)
par(mar = c(9, 9, 9, 4), mgp = c(2.5, 2.5, 0))
hist(delta_boot, main="", xlab="", ylab="", cex.axis=4, cex.lab=4, cex.main=4, xlim=c(min(delta_boot)*1.1,0))
xlab =  "Delta difference"
ylab = "Frequency"
main = ""
mtext(xlab, side=1, line=6, cex=4)
mtext(ylab, side=2, line=6, cex=4)
mtext(main, side=3, line=3, cex=4, adj=0)
abline(v=0, lwd=4, col="red")
dev.off()

png(paste0(out_path, "delta_distance_hyp_h2o2_vs_hyp_bootstrap.png"), width = 960, height = 960)
par(mar = c(9, 9, 9, 4), mgp = c(2.5, 2.5, 0))
hist(delta_boot2, main="", xlab="", ylab="", cex.axis=4, cex.lab=4, cex.main=4)
xlab =  "Delta difference"
ylab = "Frequency"
main = ""
mtext(xlab, side=1, line=6, cex=4)
mtext(ylab, side=2, line=6, cex=4)
mtext(main, side=3, line=3, cex=4, adj=0)
abline(v=0, lwd=4, col="red")
dev.off()

#細胞周Cell cycle
S.scores <- obj_rh$S.Score
G2M.scores <- obj_rh$G2M.Score


#Identify cluster
Idents(obj_rh) <- "seurat_clusters"
tab <- table(obj_rh$donor, obj_rh$seurat_clusters)
prop_tab <- prop.table(tab, margin = 2)
rep_clsts <- colnames(prop_tab)[prop_tab["P26_norm",] > 0.9]
hyp_clsts <- colnames(prop_tab)[(prop_tab["P8_hypo",] + prop_tab["P8_hypo_ROS",]) > 0.9]
write.table(rep_clsts, file = paste0(out_path, "replicative_aging_clusters.txt"), row.names = FALSE, col.names = FALSE)
write.table(hyp_clsts, file = paste0(out_path, "hypoxia_clusters.txt"), row.names = FALSE, col.names = FALSE)
clst_can <- list("rep_clsts" = rep_clsts, "hyp_clsts" = hyp_clsts)
for(clst_name in names(clst_can)){
    set.seed(1234)
    ident_can <- clst_can[[clst_name]]
    res_cl <- FindMarkers(obj_rh, ident.1 = ident_can, ident.2 = NULL, logfc.threshold = 0, min.pct = 0)
    res_cl$gene <- rownames(res_cl)
    write.csv(res_cl, file = paste0(out_path, clst_name, "_markers_full.csv"))
    res_cl2  <- res_cl[ (res_cl$pct.1 > 0 | res_cl$pct.2 > 0) & is.finite(res_cl$avg_log2FC), ]
    res_cl2 <- res_cl2[!is.na(res_cl2$avg_log2FC), , drop = FALSE]
    res_cl2 <- res_cl2[is.finite(res_cl2$avg_log2FC), , drop = FALSE]
    print(dim(res_cl2))
}

#CITE-seq module scores
unq_clueters <- c("1","5", "6")
data_path <- "/Users/dokada/Desktop/work/citeseq_dp1_0124_fin3/"
for(target_cluster in unq_clueters){
    fn <- paste0(data_path, "markers_cluster", target_cluster, "_full.csv")
    res_cl <- read.csv(fn, stringsAsFactors = FALSE, header=TRUE,row.names=1)
    res_cl <- res_cl[intersect(rownames(res_cl), rownames(obj_rh)), ]
    #markers2 <- markers[order(markers$avg_log2FC, decreasing=TRUE)[1:1000], "gene"]
    res_pos_top <- res_cl[res_cl$p_val_adj < 0.05 &  (res_cl$pct.1 > 0.1) & is.finite(res_cl$avg_log2FC) & res_cl$avg_log2FC > 0, ]
    top_n <- min(1000, nrow(res_pos_top))
    gs_use <- res_pos_top[order(res_pos_top$avg_log2FC, decreasing=TRUE)[1:top_n], "gene"]
    obj_rh <- AddModuleScore(obj_rh, features = list(gs_use), name=paste0("Cluster_score"))
    score_col <- paste0("Cluster_score", "1")  # AddModuleScoreは末尾に1が付く
    Module_score <- obj_rh[[score_col]][,1]
    obj_rh$cls_score <- Module_score
    png(paste0(out_path, "violin_cluster", target_cluster, "_marker_score.png"),width = 960, height = 960)
    p <- VlnPlot(
        obj_rh,
        features = "cls_score",
        group.by = "donor",
        pt.size = 0
    ) +
    ylab("Module Score") + 
    stat_summary(
        fun = median,
        geom = "point",
        size = 4,
        color = "black"
    )  + theme(
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 28, margin = margin(t = 10)),
        axis.text.y = element_text(size = 28, margin = margin(r = 10)),
        #axis.title.x = element_text(size = 22, margin = margin(t = 15)),
        axis.title.y = element_text(size = 28, margin = margin(r = 15)),
        plot.title = element_text(size = 28, hjust = 0.5, margin = margin(b = 20)),
        plot.margin = margin(t = 30, r = 30, b = 30, l = 30),
        legend.text  = element_text(size = 22),
        legend.title = element_text(size = 22),
        legend.position = "none"
    ) +
    ggtitle(paste0("CITE-seq C", target_cluster, "")) +
    xlab("")
    print(p)
    cat("Cluster ", target_cluster, ", ", length(gs_use), " genes used.\n")
    dev.off()

    dfp <- data.frame(donor=obj_rh$donor, cls_score=obj_rh$cls_score)
    dfp <- dfp[is.finite(dfp$cls_score) & !is.na(dfp$donor), ]
    dfp$donor <- factor(dfp$donor)
    tt <- pairwise.t.test(dfp$cls_score, dfp$donor, p.adjust.method="BH", pool.sd=FALSE)
    p_mat <- as.data.frame(as.table(tt$p.value))
    colnames(p_mat) <- c("group1","group2","p_value_adj")
    p_mat <- p_mat[!is.na(p_mat$p_value_adj), ]
    write.csv(p_mat, file=paste0(out_path, "cluster", target_cluster, "_donor_pairwise_ttest.csv"), row.names=FALSE)

}

#MMD
mmd_rbf <- function(X, Y, sigma = NULL) {
  Z <- rbind(X, Y)
  D2 <- as.matrix(dist(Z))^2
  
  if (is.null(sigma)) {
    sigma <- median(D2[D2 > 0])^0.5  # median heuristic
  }
  
  K <- exp(-D2 / (2 * sigma^2))
  n <- nrow(X); m <- nrow(Y)
  
  mean(K[1:n, 1:n]) +
  mean(K[(n+1):(n+m), (n+1):(n+m)]) -
  2 * mean(K[1:n, (n+1):(n+m)])
}
n_donor <- length(donor_id)
dist_mat <- matrix(0, nrow = n_donor, ncol = n_donor)
rownames(dist_mat) <- colnames(dist_mat)  <- donor_id
pc <- Embeddings(obj_rh, reduction = "pca")
for(i in 1:(n_donor-1)){
    for(j in (i+1):n_donor){
        donor_i <- donor_id[i]
        donor_j <- donor_id[j]
        cells_i <- colnames(obj_rh)[obj_rh$donor == donor_i]
        cells_j <- colnames(obj_rh)[obj_rh$donor == donor_j]
        X <- pc[cells_i, , drop = FALSE]
        Y <- pc[cells_j, , drop = FALSE]
        dist_ij <- mmd_rbf(X, Y)
        dist_mat[i, j] <- dist_ij
        cat(i, j, dist_ij, "\n")
    }
}
dist_mat2 <- dist_mat + t(dist_mat)
rownames(dist_mat2) <- colnames(dist_mat2)  <- donor_id
write.csv(dist_mat2, file = paste0(out_path, "donor_distance_matrix.csv"))

#MDS
dist_obj <- as.dist(dist_mat2)
mds_res <- cmdscale(dist_obj, k = 2)
mds_df <- as.data.frame(mds_res)
colnames(mds_df) <- c("MDS1", "MDS2")
mds_df$donor <- rownames(mds_df)
png(paste0(out_path, "MDS_donor_distance.png"), width = 960, height = 960)
par(mar = c(9, 9, 9, 2), mgp = c(2.5, 2.5, 0))
x_pad <- diff(range(mds_df$MDS1)) * 0.15
y_pad <- diff(range(mds_df$MDS2)) * 0.15
xlim_use <- c(min(mds_df$MDS1) - x_pad, max(mds_df$MDS1) + x_pad)
ylim_use <- c(min(mds_df$MDS2) - y_pad, max(mds_df$MDS2) + y_pad)
plot(mds_df$MDS1, mds_df$MDS2, pch = 19, cex = 2, xlab = "", ylab = "", cex.axis = 4, xlim = xlim_use, ylim = ylim_use)
mtext("MDS1", side = 1, line = 6, cex = 4)
mtext("MDS2", side = 2, line = 6, cex = 4)
mtext("", side = 3, line = 3, cex = 4)
text(mds_df$MDS1, mds_df$MDS2, mds_df$donor, pos = 3, cex = 3)
dev.off()



#MMD-based comparison
set.seed(1234)
cell_list <- list("P8_hypo"=cells_hyp, "P8_norm"=cells_norm, "P8_hypo_ROS"=cells_hyph2o2)
n_rep <- 100
res <- matrix(NA, nrow=n_rep, ncol=3)
colnames(res) <- names(cell_list)
n_sample <- 1000
for(i in 1:n_rep){
    for(j in 1:length(cell_list)){
        cells_j <- cell_list[[j]]
        X <- pc[cells_j[sample(length(cells_j), n_sample, replace = TRUE)], , drop = FALSE]
        Y <- pc[cells_highP[sample(length(cells_highP), n_sample, replace = TRUE)], , drop = FALSE]
        dist_ij <- mmd_rbf(X, Y)
        res[i, j] <- dist_ij
        cat(i, j, dist_ij, "\n")
    }
}
res <- as.matrix(res)
df <- as.data.frame(res)
df$replicate <- seq_len(nrow(df))
df_long <- pivot_longer(df, cols = c(P8_hypo, P8_norm, P8_hypo_ROS), names_to = "group", values_to = "value")
#描画
my_comp <- list(c("hypo","norm"), c("hypo","hypo_ROS"), c("norm","hypo_ROS"))
stat.test <- compare_means(value ~ group, data = df_long, comparisons = my_comp, paired = FALSE, method = "t.test", p.adjust.method = "holm")
stat.test$p.label <- paste0("P = ", signif(stat.test$p.adj, 3))
rng <- diff(range(df_long$value))
base <- max(df_long$value) + 0.05 * rng
step <- 0.08 * rng
stat.test$y.position <- base + (seq_len(nrow(stat.test)) - 1) * step
g <- ggplot(df_long, aes(x = group, y = value)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, linewidth = 0.7) +
  geom_point(position = position_jitter(width = 0.08), size = 3, alpha = 0.9) +
  theme_classic(base_size = 18) +
  theme(
    axis.title.y = element_text(size = 28),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 24),
    axis.text.y = element_text(size = 24),
    plot.title = element_text(size = 22),
    legend.text = element_text(size =20)
  ) +
  labs(y = "Distance to P26_norm (MMD)") 
ggsave(paste0(out_path, "MMD_distance_comparison.png"), plot = g, width = 8, height = 6, dpi = 300)
write.csv(df_long, file = paste0(out_path, "MMD_distance_comparison_values.csv"))
write.csv(stat.test, file = paste0(out_path, "MMD_distance_comparison_stats.csv"))

#save
sink(paste0(out_path, "session_info.txt")) 
print(sessionInfo()) 
sink()


