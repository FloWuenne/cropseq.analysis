#' A function to create a reference fasta and gtf for plasmid + gRNA fusion for CROP-seq screens
#'
#' This function takes as input a table containing the gRNA sequences and will output a .fasta and a .gtf file \
#' containing the merged plasmid + gRNA sequences and the gtf description to add them to a reference genome.
#' @param sgRNA_table The table containing the sgRNA sequences and annotation.
#' @param name The name of the column in sgRNA_table that should be used to label the fasta sequences
#' @param sequence The actual sgRNA sequence
#' @param assay The type of CRISPR assay that is being performed. If set to "activation", the tracr backbone is replaced by a tracr MS2 backbone used in VP64 screens.
#' @param outdir The directory where the .fasta and .gtf files will be saved. Default = current working directory (".")
#' @keywords Crop-seq, sgRNA, reference
#' @export
#' @examples
#' create_sgRNA_reference()

create_sgRNA_reference <- function(sgRNA_table,
                                   name = "name",
                                   sequence = "sgRNA_sequence",
                                   assay = "original",
                                   outdir = ".")
  {
  ## Create outfile names
  fasta_out = paste(outdir,"/cropseq_guides.",assay,".fasta",sep="")
  gtf_out = paste(outdir,"/cropseq_guides.",assay,".gtf" ,sep="")

  ## Store the plasmid sequences that will be assembled with the guide
  upstream_sequence <- toupper("gagggcctatttcccatgattccttcatatttgcatatacgatacaaggctgttagagagataattagaattaatttgactgtaaacacaaagatattagtacaaaatacgtgacgtagaaagtaataatttcttgggtagtttgcagttttaaaattatgttttaaaatggactatcatatgcttaccgtaacttgaaagtatttcgatttcttggctttatatatcttgtggaaaggacgaaacaccg")
  downstream_sequence_original <- toupper("gttttagagctagaaatagcaagttaaaataaggctagtccgttatcaacttgaaaaagtggcaccgagtcggtgcttttttaagcttggcgtaactagatcttgagacactgctttttgcttgtactgggtctctctggttagaccagatctgagcctgggagctctctggctaactagggaacccactgcttaagcctcaataaagcttgccttgagtgcttcaagtagtgtgtgcccgtctgttgtgtgactctggt")
  downstream_sequence_MS2 <- toupper("GTTTTAGAGCTAGGCCAACATGAGGATCACCCATGTCTGCAGGGCCTAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGGCCAACATGAGGATCACCCATGTCTGCAGGGCCAAGTGGCACCGAGTCGGTGCTTTTTTTAAGCTTGGCGTAACTAGATCTTGAGACACTGCTTTTTGCTTGTACTGGGTCTCTCTGGTTAGACCAGATCTGAGCCTGGGAGCTCTCTGGCTAACTAGGGAACCCACTGCTTAAGCCTCAATAAAGCTTGCCTTGAGTGCTTCAAGTAGTGTGTGCCCGTCTGTTGTGTGACTCTGGTAACTAGAGATCCCTCAGACCCTTTTAGTCAGTGTGGAAAATCTCTAGCAGTACGTATAGT")

  ## Check if fasta and gtf exist and if yes, delete them
  if(file.exists(fasta_out)){
    file.remove(fasta_out)
  }

  if(file.exists(gtf_out)){
    file.remove(gtf_out)
  }

  # Iterate over all entries in the sgRNA table
  for(entry in 1:nrow(sgRNA_table)){
    ## Print the number of guide working on
    print(entry)

    ## Create a fasta record
    chrom <- paste(sgRNA_table[entry,]$name,"_chrom",sep="")

    ## Check if assay has been set to activation
    if(assay == "activation"){
      assembled_gRNA <- paste(upstream_sequence,
                     sgRNA_table[entry,]$gRNA_sequence,
                     downstream_sequence_MS2,
                     sep="")
      header <- paste(">",chrom," length=",nchar(assembled_gRNA)," assay=activation",sep="")
    } else if (assay == "inhibition") {
      assembled_gRNA <- paste(upstream_sequence,
                     sgRNA_table[entry,]$gRNA_sequence,
                     downstream_sequence_original,
                     sep="")
      header <- paste(">",chrom," length=",nchar(assembled_gRNA),sep="")
    }

    ## Check if fasta file exists, if it does, make a new one
    write(header,file=fasta_out,append=TRUE)
    write(assembled_gRNA,file=fasta_out,append=TRUE)

    ## Create gtf record

    source <- "havana"
    length <- nchar(assembled_gRNA)
    if(grepl("f",sgRNA_table[entry,]$name)){
      strand <- "+"
    } else if(grepl("f",sgRNA_table[entry,]$name)){
      strand <- "-"
    }

    gene_annotation <- paste("gene_id ",
                        "\"",paste(sgRNA_table[entry,]$name,"_gene",sep=""), "\"; ",
                        "gene_name ",
                        "\"",paste(sgRNA_table[entry,]$name,"_gene",sep=""), "\"; ",
                        "gene_source ",
                        "\"","ensembl_havana", "\"; ",
                        "gene_biotype ",
                        "\"","lincRNA", "\";",
                        sep="")

    transcript_annotation <- paste("gene_id ",
                             "\"",paste(sgRNA_table[entry,]$name,"_gene",sep=""), "\"; ",
                             "transcript_id ",
                             "\"",paste(sgRNA_table[entry,]$name,"_transcript",sep=""), "\"; ",
                             "gene_name ",
                             "\"",paste(sgRNA_table[entry,]$name,"_gene",sep=""), "\"; ",
                             "gene_source ",
                             "\"","ensembl_havana", "\"; ",
                             "gene_biotype ",
                             "\"","lincRNA", "\"; ",
                             "transcript_name ",
                             "\"",paste(sgRNA_table[entry,]$name,"_transcript",sep=""), "\"; ",
                             "transcript_source ",
                             "\"","ensembl_havana", "\";",
                             sep="")

    exon_annotation <- paste("gene_id ",
                             "\"",paste(sgRNA_table[entry,]$name,"_gene",sep=""), "\"; ",
                             "transcript_id ",
                             "\"",paste(sgRNA_table[entry,]$name,"_transcript",sep=""), "\"; ",
                             "exon_number ",
                             "\"","1", "\"; ",
                             "gene_name ",
                             "\"",paste(sgRNA_table[entry,]$name,"_gene",sep=""), "\"; ",
                             "gene_source ",
                             "\"","ensembl_havana", "\"; ",
                             "gene_biotype ",
                             "\"","lincRNA", "\"; ",
                             "transcript_name ",
                             "\"",paste(sgRNA_table[entry,]$name,"_transcript",sep=""), "\"; ",
                             "transcript_source ",
                             "\"","ensembl_havana", "\"; ",
                             "exon_id ",
                             "\"",paste(sgRNA_table[entry,]$name,"_exon",sep=""), "\";",
                             sep="")

    entry_gene <- paste(chrom,source,"gene","1",length,".",strand,".",gene_annotation,sep="\t")
    entry_transcript <- paste(chrom,source,"transcript","1",length,".",strand,".",transcript_annotation,sep="\t")
    entry_exon <- paste(chrom,source,"exon","1",length,".",strand,".",exon_annotation,sep="\t")

    ## Check if gtf file exists, if it does, make a new one
    write(entry_gene,file=gtf_out,append=TRUE)
    write(entry_transcript,file=gtf_out,append=TRUE)
    write(entry_exon,file=gtf_out,append=TRUE)

    }

}
