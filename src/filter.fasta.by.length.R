#!/cluster/apps/r/3.6.0_openblas/x86_64/bin/Rscript

## Usage: filter.fasta.by.length.R <fasta> <minimum lenght> <maximum length>

## Value: filters a multifasta file by keeping only sequences of length X, where minlen <= X <= maxlen

## Author: simon.crameri@env.ethz.ch, Apr 2019

library(ape)

## Get arguments
args = commandArgs(trailingOnly = TRUE)

## Check arguments
if (length(args) != 3) {
  stop("3 arguments needed: 
       1) path to .fasta file to be filtered ; 
       2) minimum sequence length ; 
       3) maximum sequence length ", call.=FALSE)
}

## Set arguments
fasta = as.character(args[1])
minlen = as.numeric(as.character(args[2]))
maxlen = as.numeric(as.character(args[3]))

## Set arguments (for debugging)
# fasta = "consDalbergia11_2628.fasta"
# minlen = 1000
# maxlen = 3000

## Check arguments
stopifnot(file.exists(fasta),
          minlen >= 0,
          maxlen >= minlen)

## Read fasta
fas <- read.FASTA(fasta)

## Filter fasta
keep <- which(lengths(fas) >= minlen & lengths(fas) <= maxlen)
fas.filtered <- fas[keep]

## Write output
oname <- paste0(gsub(".fa$", "", gsub(".fas$", "", gsub(".fasta$", "", fasta))), "_l", minlen, "-", maxlen, ".fasta")
write.FASTA(x = fas.filtered, file = oname)

