#!/cluster/apps/r/3.6.0_openblas/x86_64/bin/Rscript

## Usage: replace.fasta.char.R <.fasta file> <char to replace> <char to be used as replacement>

## Value: replaces all <char to replace> characters with <char to be used as replacement> in each sequence (not headers)

## Author: simon.crameri@env.ethz.ch, Apr 2019

library(ape)

## Get arguments
args = commandArgs(trailingOnly = TRUE)

## Check arguments
if (length(args) != 3) {
  stop("3 arguments needed:
       1) <fasta|CHR>:       path to .fasta file needing character replacement ; 
       2) <toreplace|CHR>:   character to be replaced [single character] ; 
       3) <replacewith|CHR>: replacement character [single character]", call.=FALSE)
}

## Define arguments
fasta = args[1]
toreplace = args[2]
replacewith = args[3]

## Define arguments (for debugging)
# fasta = "NC_033804.1_LG_01_9556614_9557434_ID_3606_NC_033804.1_LG_01_9558269_9558494_ID_4610.all.mafft.align.fasta"
# toreplace="n"
# replacewith="-"

## Read fasta alignment
fas <- read.FASTA(fasta)

## Check arguments
stopifnot(length(unique(lengths(fas))) == 1,
          nchar(toreplace) == 1, 
          nchar(replacewith) == 1, length(replacewith) == 1)

## Replace <toreplace> with <replacewith> character
res <- apply(as.character(as.matrix(fas)), 2, function(x) {new <- x ; new[new %in% toreplace] <- replacewith ; return(new)})
# image.DNAbin(as.DNAbin(res))

## Write fasta (replaces original file!)
write.FASTA(as.DNAbin(res), file = fasta)

