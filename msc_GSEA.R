#GSEAの修正バージョン(2026.2.7)
#source("/Users/dokada/Dropbox/analysis/2025.4/msc_GSEA0208.R") #Run2026.2.12
out_path <- "/Users/dokada/Desktop/work/msc_GSEA0208/"
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

#analtsis
#f_vec <- c("/Users/dokada/Desktop/work/RH_hypoxia0126_fin3/hyp_clsts_gsea_results.csv","/Users/dokada/Desktop/work/RH_hypoxia0126_fin3/rep_clsts_gsea_results.csv","/Users/dokada/Desktop/work/citeseq_dp1_0124_fin3/objsub_cluster1_gsea_results.csv","/Users/dokada/Desktop/work/citeseq_dp1_0124_fin3/cluster5_gsea_hallmark.csv","/Users/dokada/Desktop/work/citeseq_dp1_0124_fin3/cluster6_gsea_hallmark.csv")
f_vec <- c("/Users/dokada/Desktop/work/RH_hypoxia0126_fin3/hyp_clsts_markers_full.csv",
           "/Users/dokada/Desktop/work/RH_hypoxia0126_fin3/rep_clsts_markers_full.csv",
           "/Users/dokada/Desktop/work/citeseq_dp1_0124_fin3/objsub_cluster1_markers_full.csv",
           "/Users/dokada/Desktop/work/citeseq_dp1_0124_fin3/markers_cluster5_full.csv",
           "/Users/dokada/Desktop/work/citeseq_dp1_0124_fin3/markers_cluster6_full.csv")
for(i in 1:length(f_vec)){
    f <- f_vec[i]
    res_cl <- read.csv(f, row.names = 1,header = TRUE)
    res_cl2  <- res_cl[ (res_cl$pct.1 > 0.01 | res_cl$pct.2 > 0.01) & is.finite(res_cl$avg_log2FC) & !is.na(res_cl$avg_log2FC), ]
    #GSEA
    gene_ranks <- res_cl2$avg_log2FC
    names(gene_ranks) <- rownames(res_cl2)
    gene_ranks <- sort(gene_ranks, decreasing = TRUE)
    msig_h <- msigdbr(species = "Homo sapiens", category = "H")
    set.seed(1000)
    gsea_h <- clusterProfiler::GSEA(geneList = gene_ranks, TERM2GENE = msig_h[, c("gs_name", "gene_symbol")], pvalueCutoff = 0.05, nPermSimple=10000, seed=TRUE)
    df <- gsea_h@result
    print(dim(df))
    clst_name <- gsub("^.*/|_gsea_results.csv$|_gsea_hallmark.csv$", "", f)
    tops <- min(10, nrow(df))
    df_top10 <- df[order(abs(df$NES), decreasing=TRUE)[tops:1], ]
    df_top10$Description <- gsub("^HALLMARK_", "", df_top10$Description)
    png(paste0(out_path, clst_name, "_gsea_top10_bar_NES.png"), width=1400, height=900)
    par(mar=c(10,55,10,2), mgp=c(5,2.2,0), tck=-0.02)
    bp <- barplot(df_top10$NES, horiz=TRUE, names.arg=df_top10$Description, las=1, xlab="", cex.axis=3.5, cex.names=3.0)
    #mtext(paste0(""), side=3, line=3.5, cex=5, font=2)
    mtext("NES", side=1, line=6.5, cex=4)
    dev.off()
}