#!/cluster/apps/r/3.6.0_openblas/x86_64/bin/Rscript
library("ape")

## Usage: get.fasta.lenght.R <.fasta>

## Value: writes the sequence lengths (if multifasta) or the alignment length (if alignment)

## Author: Simon Crameri, ETHZ, Apr 2019 

## Get arguments
args = commandArgs(trailingOnly=TRUE)

## Check arguments
if (!length(args) %in% c(1,2)) {
  stop("One argument needed: 
       REQUIRED
       1) <fasta|CHR>: path to .fasta ; 
       
       OPTIONAL
       2) <oname|CHR>: output file name [default: fbase.length]", call.=FALSE)
}

## Set arguments
fasta = as.character(args[1])
oname = as.character(args[2])

## Set arguments (for debugging)
# fasta = "woodmouse.fasta"
#oname = "test.length"

## Additional arguments
NA.char = "-"

## Define output name
fbase <- tools::file_path_sans_ext(fasta)
if (is.na(oname)) oname <- paste0(fbase, ".lengths")

## Read input
fas <- read.FASTA(fasta)

## Determine whether it is all empty, a multifasta or an alignment
if (all(unlist(as.character(fas)) %in% NA.char)) {
  type = "empty"
} else if (length(unique(lengths(fas))) == 1) {
  type = "aln"
} else {
  type = "mfas"
}

## Get lengths
switch(type,
       aln = {
         dd <- data.frame(name = fbase, length = unique(lengths(fas)))
         },
       mfas = {
         dd <- data.frame(name = names(fas), length = lengths(fas))
         empty <- which(sapply(as.character(fas), function(x) all(x %in% NA.char)))
         dd[empty,"length"] <- 0
       },
       empty = {
         dd <- data.frame(name = names(fas), length = rep(0, length(fas)))
       }
)

## Write output
#if (file.exists(oname)) first <- FALSE else first <- TRUE
write.table(x = dd, file = oname, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t", append = FALSE)


