# cropseq.analysis 

## A collection of functions to prepare and analyze single cell CRISPR perturbation screens using the CROP-seq approach or pooled CRISPR screens

## Installation
Install the package using the devtools install_github function as shown below:

library(devtools)

install_github("FloWuenne/cropseq.analysis")

## Functions 

### create_sgRNA_reference()
Uses a data frame as input and will create a fasta and gtf file containing all the guide RNA sequences merged with the 250bp upstream hU6 promoter sequences as well as the downstream gRNA backbone + 3'LTR sequence.
