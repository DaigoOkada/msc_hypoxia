#DP1のCITEseq Analysis
#source("/Users/dokada/Dropbox/analysis/2025.4/citeseq_dp1_0124_fin3.R")
out_path <- "/Users/dokada/Desktop/work/citeseq_dp1_0124_fin3/"
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
library(pheatmap)
data_path <- "/Users/dokada/Desktop/work/shizui_data/citeseq_dp1/filtered_feature_bc_matrix/"
sc.data <- Read10X(data.dir = data_path)

#RNA
set.seed(1234)
obj_ori <- CreateSeuratObject(counts = sc.data[["Gene Expression"]], project = "CITEseq") #36601 features across 16046 samples
adt <- sc.data[["Antibody Capture"]]
common.cells <- intersect(colnames(obj_ori), colnames(adt))
obj_ori <- subset(x = obj_ori, cells = common.cells) #36601 features across 16046 samples
adt <- adt[, common.cells] #146 proteins
obj_ori[["ADT"]] <- CreateAssayObject(counts = adt)
DefaultAssay(obj_ori) <- "RNA"
obj_ori[["percent.mt"]] <- PercentageFeatureSet(obj_ori, pattern = "^MT-")
obj <- obj_ori
features <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
obj <- subset(obj, subset = nFeature_RNA > 200 & percent.mt < 10)
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
#cell cycle check
s.genes   <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
s.genes.use   <- intersect(s.genes, rownames(obj))
g2m.genes.use <- intersect(g2m.genes, rownames(obj))
obj <- CellCycleScoring(obj, s.features = s.genes.use, g2m.features = g2m.genes.use, set.ident = FALSE)
obj <- ScaleData(obj, features = VariableFeatures(obj))
obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = 50, verbose = FALSE)
obj <- FindNeighbors(obj, dims = 1:50)
obj <- FindClusters(obj, resolution = 0.4) #JC論文と同じ
obj <- RunUMAP(obj, dims = 1:50)
#doublet finder
sweep.res <- paramSweep(obj, PCs = 1:50, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
best.pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
nExp <- round(ncol(obj) * 0.06)  # 仮：6%
homotypic.prop <- modelHomotypic(obj$seurat_clusters)
nExp.adj <- round(nExp * (1 - homotypic.prop))
obj <- doubletFinder(obj, PCs = 1:50, pN = 0.25, pK = best.pK, nExp = nExp.adj, sct = FALSE)
df.col <- grep("^DF.classifications_", colnames(obj@meta.data), value = TRUE)
singlet.cells <- rownames(obj@meta.data)[obj@meta.data[[df.col[length(df.col)]]] == "Singlet"]
obj <- subset(obj, cells = singlet.cells)
#cell cycle regression
DefaultAssay(obj) <- "RNA"
#Re-run
obj@reductions <- list()
obj@graphs <- list()
obj@neighbors <- list()
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
# cell cycle scoring
s.genes.use   <- intersect(cc.genes$s.genes, rownames(obj))
g2m.genes.use <- intersect(cc.genes$g2m.genes, rownames(obj))
obj <- CellCycleScoring(obj, s.features = s.genes.use, g2m.features = g2m.genes.use, set.ident = FALSE)
obj <- ScaleData(obj, features = VariableFeatures(obj))
# DR
obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = 50, verbose = FALSE)
obj <- FindNeighbors(obj, dims = 1:50)
obj <- FindClusters(obj, resolution = 0.4)
obj <- RunUMAP(obj, dims = 1:50)
saveRDS(obj, file = paste0(out_path, "obj.rds"))

#Clueter5, 6 , Othersで分類
obj$cluster_5_6_other <- "Other"
obj$cluster_5_6_other[obj$seurat_clusters == 5] <- "C5"
obj$cluster_5_6_other[obj$seurat_clusters == 6] <- "C6"
obj$cluster_5_6_other <- factor(obj$cluster_5_6_other, levels = c("C5", "C6", "Other"))


#QCをしないときのクラスタリング結果の比較
set.seed(1234)
obj_noqc <- obj_ori
DefaultAssay(obj_noqc) <- "RNA"
obj_noqc <- NormalizeData(obj_noqc, normalization.method = "LogNormalize", scale.factor = 10000)
obj_noqc <- FindVariableFeatures(obj_noqc, selection.method = "vst", nfeatures = 2000)
obj_noqc <- ScaleData(obj_noqc, features = VariableFeatures(obj_noqc))
obj_noqc <- RunPCA(obj_noqc, features = VariableFeatures(obj_noqc), npcs = 50, verbose = FALSE)
obj_noqc <- FindNeighbors(obj_noqc, dims = 1:50)
obj_noqc <- FindClusters(obj_noqc, resolution = 0.4) #JC論文と同じ
obj_noqc <- RunUMAP(obj_noqc, dims = 1:50)
obj_noqc <- CellCycleScoring(obj_noqc, s.features = s.genes.use, g2m.features = g2m.genes.use, set.ident = FALSE)
saveRDS(obj_noqc, file = paste0(out_path, "obj_noqc.rds"))



#対応チェック
common_cells <- intersect(colnames(obj), colnames(obj_noqc))
tab <- table(obj$seurat_clusters[common_cells], obj_noqc$seurat_clusters[common_cells])
write.csv(tab, file = paste0(out_path, "qc_noqc_comparison_table.csv"))
map_noqc_to_qc <- apply(tab, 2, function(x) rownames(tab)[which.max(x)])
map_noqc_to_qc <- as.character(map_noqc_to_qc)  # names(map_noqc_to_qc) = noQCクラスタ
names(map_noqc_to_qc) <- colnames(tab)
qc_levels <- levels(obj$seurat_clusters)
qc_cols <- c("blue", "orange", "green", "red", "purple", "brown", "pink", "gold", "grey40")
names(qc_cols) <- qc_levels
tmp <- map_noqc_to_qc[as.character(obj_noqc$seurat_clusters)]
names(tmp) <- names(obj_noqc$seurat_clusters)
obj_noqc$seurat_clusters_mapped <- factor(tmp, levels = qc_levels)
#Draw for obj
p <- DimPlot(obj, reduction="umap", group.by="seurat_clusters", label=TRUE, label.size=12, pt.size=2, cols=qc_cols) +
  labs(title=NULL, x="UMAP1", y="UMAP2") +
  theme(
    plot.title = element_blank(),
    axis.title = element_text(size=40, face="bold"),
    axis.title.x = element_text(margin=margin(t=20)),
    axis.title.y = element_text(margin=margin(r=20)),
    axis.text = element_text(size=32),
    legend.title = element_text(size=34),
    legend.text = element_text(size=30),
    text = element_text(size=32),
  )
ggsave(paste0(out_path, "umap_clusters.png"), plot=p, width=12, height=10, dpi=300)
cluster_labels <- obj$seurat_clusters
cluster_counts <- table(cluster_labels)
tab <- round(cluster_counts / sum(cluster_counts) * 100, 2)
pie_cols <- qc_cols[names(cluster_counts)]
png(paste0(out_path, "pie_clusters.png"), width=960, height=960)
par(mar=c(9,9,9,4), mgp=c(2.5,2.5,0))
pie(cluster_counts, labels=paste0("C", names(cluster_counts), ": ", tab, "%"), col=pie_cols, cex=2)
mtext("Percentages", side=3, line=3, cex=4)
dev.off()
#Draw for obj_noqc
p_noqc <- DimPlot(obj_noqc, reduction="umap", group.by="seurat_clusters_mapped", label=TRUE, label.size=12, pt.size=2, cols=qc_cols) +
  labs(title=NULL, x="UMAP1", y="UMAP2") +
  theme(
    plot.title = element_blank(),
    axis.title = element_text(size=40, face="bold"),
    axis.title.x = element_text(margin=margin(t=20)),
    axis.title.y = element_text(margin=margin(r=20)),
    axis.text = element_text(size=32),
    legend.title = element_text(size=34),
    legend.text = element_text(size=30),
    text = element_text(size=32)
  )
ggsave(paste0(out_path, "umap_clusters_noqc_mapped_to_qc.png"), plot=p_noqc, width=12, height=10, dpi=300)
cluster_labels_noqc <- obj_noqc$seurat_clusters_mapped
cluster_counts_noqc <- table(cluster_labels_noqc)
tab_noqc <- round(cluster_counts_noqc / sum(cluster_counts_noqc) * 100, 2)
pie_cols_noqc <- qc_cols[names(cluster_counts_noqc)]
png(paste0(out_path, "pie_clusters_noqc_mapped_to_qc.png"), width=960, height=960)
par(mar=c(9,9,9,4), mgp=c(2.5,2.5,0))
pie(cluster_counts_noqc, labels=paste0("C", names(cluster_counts_noqc), ": ", tab_noqc, "%"), col=pie_cols_noqc, cex=2)
mtext("Percentages", side=3, line=3, cex=4)
dev.off()


#Feature plots for QC and Non-QC
obj$Cell.Cycle <- pmax(obj$S.Score, obj$G2M.Score)
obj_noqc$Cell.Cycle <- pmax(obj_noqc$S.Score, obj_noqc$G2M.Score)
out_dir <- paste0(out_path, "Feature_plots_QC_vs_NonQC/")
if(!file.exists(out_dir)) dir.create(out_dir)
features <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "Cell.Cycle")
obj_list <- list("QC" = obj, "Non_QC" = obj_noqc)
for(obj_name in names(obj_list)){
    obj_tmp <- obj_list[[obj_name]]
    used_group <- ifelse(obj_name == "QC", "seurat_clusters", "seurat_clusters_mapped")
    for(i in 1:length(features)){
        png(paste0(out_dir, obj_name, "_violin_", features[i], ".png"),width = 960, height = 960)
        p <- VlnPlot(obj_tmp, features = features[i], ncol = 1, group.by = used_group) +
        labs(x="Cluster") +
        theme(
            axis.text.x = element_text(size = 32, margin = margin(t = 10)),
            axis.text.y = element_text(size = 32, margin = margin(r = 10)),
            axis.title.x = element_text(size = 40, margin = margin(t = 15)),
            axis.title.y = element_text(size = 40, margin = margin(r = 15)),
            plot.title = element_text(size = 40, hjust = 0.5, margin = margin(b = 20)),
            plot.margin = margin(t = 30, r = 30, b = 30, l = 30),
            legend.text  = element_text(size = 30),
            legend.title = element_text(size = 34),
            legend.position = "none"
        )
        print(p)
        dev.off()
    }
}



#Clueter1, Cluster5、Cluster6のRNAによるマーカー同定
DefaultAssay(obj) <- "RNA"
ident_can <- c(1, 5, 6)
msig_h <- msigdbr(species = "Homo sapiens", category = "H")
term2gene_h <- msig_h[, c("gs_name", "gene_symbol")]
for(tmp_cluster in ident_can){
    #Detect merkers
    res_cl <- FindMarkers(obj, ident.1 = tmp_cluster, ident.2 = NULL, assay = "RNA", slot = "data", min.pct = 0, logfc.threshold = 0)
    res_cl$gene <- rownames(res_cl)
    write.csv(res_cl, file = paste0(out_path, paste0("markers_cluster", tmp_cluster, "_full.csv")))
    res_cl2  <- res_cl[ (res_cl$pct.1 > 0 | res_cl$pct.2 > 0) & is.finite(res_cl$avg_log2FC), ]
    res_cl2 <- res_cl2[!is.na(res_cl2$avg_log2FC), , drop = FALSE]
    res_cl2 <- res_cl2[is.finite(res_cl2$avg_log2FC), , drop = FALSE]
    print(dim(res_cl2))
}



#protein_marker_names
protein_names <- read.csv("/Users/dokada/Desktop/work/citeseq_marker_names/protein_names.csv", header=TRUE, stringsAsFactors=FALSE, row.names=1)
#rownames(protein_names) <- protein_names[,1]
DefaultAssay(obj) <- "ADT"
obj <- NormalizeData(obj, normalization.method = "CLR", margin = 2, verbose = FALSE)
prot_norm <- t(obj[["ADT"]]@data)
stopifnot(all(rownames(protein_names) == colnames(prot_norm)))
target_cols <- grep("^anti-human-", colnames(prot_norm), value = TRUE)
prot_norm2 <- prot_norm[, target_cols]
new_names <- protein_names[colnames(prot_norm2), "x_name"]
colnames(prot_norm2) <- new_names

#mean-variance plot
mean_expr <- colMeans(prot_norm2)
var_expr <- apply(prot_norm2, 2, var)
png(paste0(out_path, "mean_variance.png"), width=960, height=960)
par(mar = c(9, 9, 9, 2), mgp = c(2.5, 2.5, 0)) 
xlab = "Mean Expression"; ylab = "Variance"; main = ""
plot(mean_expr, var_expr, xlab="", ylab="", cex=2, cex.axis=4, cex.lab=4, pch=19, xlim=c(0,10))
mtext(xlab, side=1, line=6, cex=4)
mtext(ylab, side=2, line=6, cex=4)
mtext(main, side=3, line=3, cex=4)
var_thres <- 0.1
mean_thres <- 2
idx <- which(mean_expr >= mean_thres | var_expr >= var_thres)
used_proteins <- names(mean_expr)[idx]
non_used_proteins <- names(mean_expr)[-idx]
abline(h=var_thres,col="red", lwd=4)
abline(v=mean_thres,col="red", lwd=4)
text(mean_expr[idx], var_expr[idx], labels = names(mean_expr)[idx], pos = 4, cex = 2)
dev.off()
mv_res <- data.frame(Mean=mean_expr, Variance=var_expr)
write.csv(mv_res, file=paste0(out_path, "prot_mean_variance.csv"))

#negative marker check
negative_marker <- c("CD45", "CD34", "CD14", "CD11b", "CD79a", "CD19", "HLA-DR")
stopifnot(any(!(negative_marker %in% used_proteins)))

#Detect Protein Markers for Cluster5 and Cluster6
DefaultAssay(obj) <- "ADT"
res_cl5 <- FindMarkers(obj, ident.1 = 5, ident.2 = c(0,1,2,3,4))
write.csv(res_cl5, file=paste0(out_path, "prot_markers_cluster5.csv"))
res_cl6 <- FindMarkers(obj, ident.1 = 6, ident.2 = c(0,1,2,3,4))
write.csv(res_cl6, file=paste0(out_path, "prot_markers_cluster6.csv"))

#Cluster6のポジティブマーカを同定
used_proteins_id <- protein_names[protein_names[, "x_name"] %in% used_proteins,1 ]
posi6 <- intersect(rownames(res_cl6[res_cl6[,"p_val_adj"]<0.05 & res_cl6[,"avg_log2FC"] > 0,]), used_proteins_id) #CD24
posi5 <- intersect(rownames(res_cl5[res_cl5[,"p_val_adj"]<0.05 & res_cl5[,"avg_log2FC"] > 0,]), used_proteins_id) #CD105

#識別するタンパクマーカーの探索 (CD73, CD90, CD105, CD24)'s Violine Plot
marker_spe <- c("CD73", "CD90", "CD105", "CD24")
for(i in 1:length(marker_spe)){
    DefaultAssay(obj) <- "RNA"
    feat_name <- marker_spe[i]
    feat <- protein_names[protein_names$x_name == feat_name, 1]
    #Violin plot
    png(paste0(out_path, "Vln_ADT_", feat_name, ".png"), width = 2000, height = 1600, res = 300)
    p <- VlnPlot(obj, features = feat, group.by = "cluster_5_6_other",
                pt.size = 0, assay = "ADT") +
        theme(axis.text.x = element_text(size = 28, angle = 45, hjust = 1),
                axis.text.y = element_text(size = 28),
                axis.title.x = element_text(size = 28),
                axis.title.y = element_text(size = 28),
                plot.title = element_text(size = 36, hjust = 0.5),
                legend.position = "none") +
         ggtitle(feat_name) +
        xlab("")
    print(p)
    dev.off()

    #Feature plot
    DefaultAssay(obj) <- "ADT"
    png(paste0(out_path, "FeaturePlot_ADT_", feat_name, ".png"), width=2000, height=1600, res=300)
    p <- FeaturePlot(obj, features=feat, reduction="umap", cols=c("lightgrey","red")) +
        theme(axis.text.x=element_text(size=28), axis.text.y=element_text(size=28),
              axis.title.x=element_text(size=28), axis.title.y=element_text(size=28),
              plot.title=element_text(size=36, hjust=0.5),
              legend.position = "none") +
        ggtitle(feat_name) +
        xlab("UMAP1") +
        ylab("UMAP2")
    print(p)
    dev.off()
}




#Output Seurat output for RNA velocity analysis
library(SeuratDisk)
umap <- Embeddings(obj, reduction = "umap")
write.csv(umap, file = file.path(out_path, "umap.csv"), quote = FALSE)
clu <- data.frame(seurat_clusters = as.character(obj$seurat_clusters))
rownames(clu) <- colnames(obj)
write.csv(clu, file = file.path(out_path, "clusters.csv"), quote = FALSE)


#Cluster5 , 6を除外して、再度クラスタリングを行う
set.seed(1234)
target_clusters <- c(0, 1, 2, 3, 4)
obj_sub <- subset(obj, subset = seurat_clusters %in% target_clusters)
DefaultAssay(obj_sub) <- "RNA"
obj_sub <- NormalizeData(obj_sub, normalization.method = "LogNormalize", scale.factor = 10000)
obj_sub <- FindVariableFeatures(obj_sub, selection.method = "vst", nfeatures = 2000)
s.genes   <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
s.genes.use   <- intersect(s.genes, rownames(obj_sub))
g2m.genes.use <- intersect(g2m.genes, rownames(obj_sub))
#obj_sub <- CellCycleScoring(obj_sub, s.features = s.genes.use, g2m.features = g2m.genes.use, set.ident = FALSE)
obj_sub <- ScaleData(obj_sub, features = VariableFeatures(obj_sub))
obj_sub <- RunPCA(obj_sub, features = VariableFeatures(obj_sub), npcs = 50, verbose = FALSE)
obj_sub <- FindNeighbors(obj_sub, dims = 1:50)
tmp_resolution <- 0.4
obj_sub <- FindClusters(obj_sub, resolution = tmp_resolution) #JC論文と同じ
obj_sub <- RunUMAP(obj_sub, dims = 1:50)
saveRDS(obj_sub, file = paste0(out_path, "obj_sub.rds"))


#対応チェック
common_cells <- intersect(colnames(obj), colnames(obj_sub))
tab <- table(obj$seurat_clusters[common_cells], obj_sub$seurat_clusters[common_cells])
write.csv(tab, file = paste0(out_path, "qc_noqc_comparison_table.csv"))
map_sub_to_qc <- apply(tab, 2, function(x) rownames(tab)[which.max(x)])
map_sub_to_qc <- as.character(map_sub_to_qc)  # names(map_sub_to_qc) = noQCクラスタ
names(map_sub_to_qc) <- colnames(tab)
qc_levels <- levels(obj$seurat_clusters)
names(qc_cols) <- qc_levels
tmp <- map_sub_to_qc[as.character(obj_sub$seurat_clusters)]
names(tmp) <- names(obj_sub$seurat_clusters)
obj_sub$seurat_clusters_mapped <- factor(tmp, levels = qc_levels)
p_sub <- DimPlot(obj_sub, reduction="umap", group.by="seurat_clusters_mapped", label=TRUE, label.size=10, pt.size=2, cols=qc_cols) + theme(text=element_text(size=28), axis.text=element_text(size=24), axis.title=element_text(size=26), legend.text=element_text(size=22), legend.title=element_text(size=24))
ggsave(paste0(out_path, "sub_mapped_to_qc.png"), p_sub, width=12, height=10, dpi=300)

# UMAPを出力
umap <- Embeddings(obj_sub, reduction = "umap")
write.csv(umap, file = file.path(out_path, "umap_sub.csv"), quote = FALSE)
clu <- data.frame(seurat_clusters = as.character(obj_sub$seurat_clusters_mapped))
rownames(clu) <- colnames(obj_sub)
write.csv(clu, file = file.path(out_path, "clusters_sub.csv"), quote = FALSE)



#QC feature plots
features <- c("S.Score", "G2M.Score", "Cell.Cycle")
for(i in 1:length(features)){
    png(paste0(out_path, "Objsub_QC_violin_", features[i], ".png"), width=960, height=960, res=150)
    png(paste0(out_path, "Objsub_QC_violin_", features[i], ".png"), width=960, height=960)
    p <- VlnPlot(obj_sub, features=features[i], ncol=1, group.by="seurat_clusters_mapped") +
    xlab("Cluster") +
    theme(
        axis.text.x=element_text(size=32, margin=margin(t=10)),
        axis.text.y=element_text(size=32, margin=margin(r=10)),
        axis.title.x=element_text(size=32, margin=margin(t=15)),
        axis.title.y=element_text(size=32, margin=margin(r=15)),
        plot.title=element_text(size=32, hjust=0.5, margin=margin(b=20)),
        plot.margin=margin(t=30, r=30, b=30, l=30),
        legend.text=element_text(size=32),
        legend.title=element_text(size=32),
        legend.position="none"
    )
    print(p)
    dev.off()
}

#clueter1 vs cluster0,2,3,4
Idents(obj_sub) <- obj_sub$seurat_clusters_mapped
res_cl1 <- FindMarkers(obj_sub, ident.1 = 1, ident.2 = NULL, assay = "RNA", slot = "data", min.pct = 0, logfc.threshold = 0)
res_cl1$gene <- rownames(res_cl1)
write.csv(res_cl1, file = paste0(out_path, "objsub_cluster1_markers_full.csv"))



#Protein
DefaultAssay(obj_sub) <- "ADT"
obj_sub <- NormalizeData(obj_sub, normalization.method = "CLR", margin = 2, verbose = FALSE)
prot_norm <- t(obj_sub[["ADT"]]@data)
stopifnot(all(rownames(protein_names) == colnames(prot_norm)))
target_cols <- grep("^anti-human-", colnames(prot_norm), value = TRUE)
prot_norm2 <- prot_norm[, target_cols]
new_names <- protein_names[colnames(prot_norm2), "x_name"]
colnames(prot_norm2) <- new_names


#Cluster1のタンパク質マーカーを調べる
#色分けを変える
library(car)
res <- data.frame(Protein=used_proteins, Cls.P=NA, Cls.Effects=NA, Cycle.P=NA, Cycle.Effects=NA, stringsAsFactors=FALSE)
cluster_labels <- obj_sub$seurat_clusters_mapped
cycle <- factor(obj_sub$Phase, levels=c("G1","S","G2M"))
for(i in seq_along(used_proteins)){
    tmp_prot <- used_proteins[i]
    y <- prot_norm2[, tmp_prot]
    cluster <- ifelse(cluster_labels == 1, 1, 0) #Cluster1 is aged cluster
    df <- data.frame(y=y, cycle=cycle, cluster=cluster)

    fit <- lm(y ~ cluster + cycle, data=df)
    sm <- summary(fit)
    aov_tab <- Anova(fit, type=2)

    res$Cls.P[i] <- sm$coefficients["cluster","Pr(>|t|)"]
    res$Cls.Effects[i] <- sm$coefficients["cluster","Estimate"]
    res$Cycle.P[i] <- aov_tab["cycle","Pr(>F)"]
    res$Cycle.Effects[i] <- max(abs(sm$coefficients[grep("^cycle", rownames(sm$coefficients)),"Estimate"]))
}
write.csv(res, file=paste0(out_path, "cluster2_prot_marker_ccregression.csv")) #Typo for Cluster1
#idx <- which(Cls.P.adj < 0.01 & abs(res$S.Effects) < 0.05 & abs(res$G2M.Effects) < 0.05) ##CD47
Cls.P.adj <- p.adjust(res$Cls.P, method = "bonferroni")
Cycle.P.adj <- p.adjust(res$Cycle.P, method = "bonferroni")
cls_sig <- Cls.P.adj < 0.05
cc_sig <- Cycle.P.adj < 0.05
# 色設定
cols <- rep("black", nrow(res))
cols[cls_sig & !cc_sig] <- "red"
cols[!cls_sig & cc_sig] <- "blue"
cols[cls_sig & cc_sig] <- "purple"
# 効果量
cls_eff <- res$Cls.Effects
cc_eff <- -log10(pmax(res$Cycle.P, .Machine$double.xmin))
# プロット
png(paste0(out_path, "cluster1_cellcycle_effects.png"), width = 960, height = 960)
par(mar = c(9, 9, 9, 2), mgp = c(2.5, 2.5, 0))
plot(cls_eff, cc_eff, pch = 19, cex = 2, col = cols, xlab = "", ylab = "", cex.axis = 4,xlim=c(min(cls_eff)*1.2, max(cls_eff)*1.2), ylim=c(0, max(cc_eff)*1.2))
mtext("Cluster Effects", side = 1, line = 6, cex = 4)
mtext("-log10 (Cell Cycle ANOVA P value)", side = 2, line = 6, cex = 4)
top_x <- order(abs(cls_eff), decreasing = TRUE)[1:5]
top_y <- order(abs(cc_eff),  decreasing = TRUE)[1:5]
only_cls <- which(cls_sig & !cc_sig)
label_idx <- union(union(top_x, top_y), only_cls)
text(x = cls_eff[label_idx], y = cc_eff[label_idx], labels = res$Protein[label_idx], pos = 3, cex = 2)
dev.off()

#All protein violin
DefaultAssay(obj_sub) <- "ADT"
out_dir <- paste0(out_path, "Objsub_Protein_Feature_Violin/")
if(!file.exists(out_dir)) dir.create(out_dir)
for(feat in used_proteins_id){

    png(paste0(out_dir, "VioPlot_Cluster_", feat, ".png"), width=2000, height=1600, res=300)
    p <- VlnPlot(obj_sub, features=feat, group.by="seurat_clusters_mapped", pt.size=0, assay="ADT") +
        geom_boxplot(width=0.25, outlier.shape=NA, fill="white", color="black", linewidth=0.8) +
        theme_classic() +
        theme(axis.text.x=element_text(size=18, angle=45, hjust=1), axis.text.y=element_text(size=18), axis.title.x=element_text(size=20), axis.title.y=element_text(size=20), plot.title=element_text(size=24, hjust=0.5), legend.position="none") +
        xlab("Cluster")
    print(p)
    dev.off()
    png(paste0(out_dir, "VioPlot_Phase_", feat, ".png"), width=2000, height=1600, res=300)
    p <- VlnPlot(obj_sub, features=feat, group.by="Phase", pt.size=0, assay="ADT") +
        geom_boxplot(width=0.25, outlier.shape=NA, fill="white", color="black", linewidth=0.8) +
        theme_classic() +
        theme(axis.text.x=element_text(size=18, angle=45, hjust=1), axis.text.y=element_text(size=18), axis.title.x=element_text(size=20), axis.title.y=element_text(size=20), plot.title=element_text(size=24, hjust=0.5), legend.position="none") +
        xlab("Phase")
    print(p)
    dev.off()
}


# 保存するファイル名を指定
sink(paste0(out_path, "session_info.txt")) 
print(sessionInfo())
sink() 

#Run RNA velocity analysis (python script)
#conda activate velocity
#cd /Users/dokada/Desktop/work/citeseq_dp1_0124_fin3/ && python /Users/dokada/Dropbox/analysis/2025.4/rna_velo_fin.py --loom /Users/dokada/Desktop/work/citeseq_dp1_1212/velocyto_out/possorted_genome_bam_7JA39.loom --umap_csv umap_sub.csv --clusters_csv clusters_sub.csv --outdir . --plot velo1_png1_raw.png --save_h5ad adata_velo_merged.h5ad