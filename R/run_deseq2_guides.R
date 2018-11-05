#' A function to calculate DESeq2 log-fold changes for sgRNA counts from pooled CRISPR screens. Returns a table of DEseq2
#' log-fold changes.
#'
#' This function takes as input a table containing the gRNA sequences and will output a .fasta and a .gtf file \
#' containing the merged plasmid + gRNA sequences and the gtf description to add them to a reference genome.
#' @param design_matrix Matrix containing the experimental design that is used to perform estimation of DE genes
#' @param sgRNA_count_table Table containing non-normalized sgRNA counts. Can be produced with MAGECK count.
#' @export
#' @examples
#' run_deseq2_on_guides()

run_deseq2_guides <- function(design_matrix,
                                 sgRNA_count_table)
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
                                design = ~ distribution)


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
