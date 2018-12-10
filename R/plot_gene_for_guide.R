#' Function to extract expression values for cells for a specific guide RNA and scrambles for a given gene from a CROP-seq Seurat object.
#'
#' This functions takes a Seurat object with annotation for sgRNA as input. Based on sgRNA name specification and the gene of interest
#' will extract the normalized expression values for all cells containing the guide and all scrambles. If plot_violin is true
#' (default = TRUE), will also print a violin plot of the candidate gene expression for both sgRNA and scrambles to compare
#' expression.
#'
#' @param seurat_object S4 Seurat object containing a column labeled sgRNA and normalized expression values in the data slot.
#' @param guide_ID ID of the sgRNA of interest. Format : chr_start_end_nr (example = sgRNA_19_12993205_12993305_16r_gene)
#' @param gene_name Gene of interest to extract expression values for.
#' @param plot_violin Should a violin plot for the gene of interested be drawn or not?
#' @export
#' @examples
#' run_deseq2_on_guides()
plot_gene_for_guide <- function(seurat_object,
                                guide_ID = FALSE,
                                gene_name = "HBG2",
                                plot_violin = TRUE){

  ## Subset Seurat object for guide and gene of interest and scrambled guides
  seurat_object_exp <- as.data.frame(as.matrix(seurat_object@data[gene_name,]))
  seurat_object_exp$cell <- rownames(seurat_object_exp)
  colnames(seurat_object_exp) <- c(gene_name,"cell")

  seurat_object_meta_guide <- subset(seurat_object@meta.data,
                                       sgRNA == guide_ID)

  seurat_object_meta_scramble <- subset(seurat_object@meta.data,
                                        grepl("scramble",sgRNA))


  seurat_object_exp_long <- seurat_object_exp %>%
    gather(gene,norm_exp,-cell) %>%
    subset(cell %in% rownames(seurat_object_meta_guide) |
             cell %in% rownames(seurat_object_meta_scramble))

  seurat_object_exp_long <- seurat_object_exp_long %>%
    mutate("sgRNA" = if_else(cell %in% rownames(seurat_object_meta_guide),guide_ID,"scramble"))

  violin_plot_sgRNA_gene <- ggplot(seurat_object_exp_long,aes(sgRNA,norm_exp, fill= sgRNA)) +
    geom_violin(alpha =0.6) +
    geom_jitter(alpha = 0.5) +
    labs(x = "sgRNAs",
         y = "normalized epxression",
         title = paste("Guide:",guide_ID,sep = ""),
         subtitle = paste("Gene tested:",gene_name,sep= " ")) +
    theme(legend.position = "none")

  if(plot_violin == TRUE){
    print(violin_plot_sgRNA_gene)
  }

  return(seurat_object_exp_long)

}
