#' A function to calculate DESeq2 log-fold changes for sgRNA counts from poole CRISPR screens
#'
#' This function takes as input a table containing the gRNA sequences and will output a .fasta and a .gtf file \
#' containing the merged plasmid + gRNA sequences and the gtf description to add them to a reference genome.
#' @param sgRNA_table The table containing the sgRNA sequences and annotation.
#' @export
#' @examples
#' run_deseq2_on_guides()

run_deseq2_guides <- function(design_matrix,
                                 sgRNA_count_table,
                                 design_column = "distribution")
  {

  ## Specify required libraries
  require(DESeq2)
  require(data.table)

  ## Load design matrix
  design <- read.table(design_matrix,
                       sep="\t",
                       header=T,
                       row.names= 1,
                       quote = "")

  design <- subset(design,distribution != "library")

  ## Load count tables calculated using MAGECK count
  counts <- fread(sgRNA_count_table,
                  sep="\t",
                  header=T,
                  quote = "")

  counts_filtered <- counts

  rownames(counts_filtered) <- counts_filtered$sgRNA
  counts_filtered <- counts_filtered[,-c(1,2,5)]

  ## Create DEseq2 object
  dds <- DESeqDataSetFromMatrix(countData = counts_filtered,
                                colData = design,
                                design = ~ design_column)

  dds$distribution <- factor(dds$distribution, levels = c("BOTTOM","TOP"))

  dds <- DESeq(dds)
  res <- results(dds)

  ## Add sgRNA and locus information to the DESeq2 results
  res$sgRNA <- counts$sgRNA
  res$locus <- counts$Gene

  ## Find all the guides that have NA in their logFold change
  res_na <- res %>%
    subset(is.na(log2FoldChange))

  ## Filter out the sgRNA that have NA
  res <- res %>%
    subset(!is.na(log2FoldChange))

  ## Order results by p-value
  resOrdered <- res[order(res$pvalue),]

  return(resOrdered)

}
