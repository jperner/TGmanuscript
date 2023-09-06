
# PCA on top N genes with median absolute deviance
runPCA_topN <- function(vsd_sub, N=1000){
  # median absolute deviance 
  top1000 <- tail(order(rowMads(vsd_sub)), n=N)
  top1000 <- row.names(vsd_sub)[top1000]
  
  # pca 
  return(prcomp(t(vsd_sub[top1000, ])))
}

# barplot of variance explained by each PC
plot_ExplainedVar <- function(pca, main_label){
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
  names(percentVar) <- paste0("PC", 1:length(percentVar))
  barplot(head(percentVar, n=10), 
          las=2, 
          ylab="Explained variation", 
          main=main_label)
  return(percentVar)
}

# scatter plot of first and second PC with additional aesthetics 
plot_PC_1to2 <- function(toplot, 
                         setcols, 
                         main_title, 
                         vars, 
                         col_var="numeric_group", 
                         col_var_name="Group", 
                         shape_var="F2",
                         pchs=c(21, 10, 22, 7)){
  
  p <-ggplot(toplot, aes_string(x="PC1", y="PC2", color=col_var, shape=shape_var)) +
    geom_point(size=5, stroke=1.5) +
    scale_shape_manual(values = pchs, name="Cohort") +
    theme_bw() +
    ggtitle(main_title) +
    xlab(paste0("PC1 (", round(vars["PC1"]*100), "% variance)")) +
    ylab(paste0("PC2 (", round(vars["PC2"]*100), "% variance)")) 
  
  if(!is.null(setcols)) {
    p <- p + 
      scale_color_manual(values=setcols, name=col_var_name) 
  }
  
  return(p)
}

# wrapper for volcano plot with fixed parameters
plot_Volcano <- function(res, 
                         main_title, 
                         gene_labels = "", 
                         g_list=NULL, 
                         xmin=-10, 
                         xmax=10, 
                         ymax=30, 
                         fdr_cut=0.1, 
                         FCcut=2){
  
  volc <-EnhancedVolcano(toptable = res, 
                         lab = gene_labels,
                         selectLab = g_list,
                         x="log2FoldChange",
                         y= "pvalue",
                         title = main_title,
                         subtitle = NULL, 
                         caption = NULL,
                         legendPosition = "bottom", 
                         FCcutoff = FCcut,
                         pCutoffCol = "padj",
                         pCutoff = fdr_cut,
                         xlim = c(xmin, xmax),
                         ylim =c(0, ymax))
  print(volc)
} 

# select differentially expressed genes from results table
selectSignGenes <- function(res, FDRcut=0.1, absFCcut=1){
  require(tidyverse)
  as.data.frame(res) %>% 
    filter(!is.na(padj)) %>% 
    filter(padj<FDRcut & abs(log2FoldChange)>absFCcut)
}
