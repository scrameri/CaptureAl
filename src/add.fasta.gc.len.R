#!/cluster/apps/r/3.6.0_openblas/x86_64/bin/Rscript

## Usage: add.fasta.gc.len.R

## Value: adds GC content and sequence length to fasta headers

## Author: simon.crameri@env.ethz.ch, Apr 2019

## Load libraries
library(ape)

## Get arguments
args = commandArgs(trailingOnly = TRUE)

## Check arguments
if (!length(args) %in% c(1, 2, 3, 4)) {
  stop("3 arguments taken (1 needed):
        REQUIRED
        1) <fasta|CHR>: path to .fasta file needing new header (required) ; 

	OPTIONAL
        2) <replace|BOOLEAN>: whether to replace old sequene name with <NUMER_GC_gc_LEN_length>
        3) string for header prefix ('' by default) ; 
        4) string for header suffix ('' by default)",
        call.=FALSE)
}

## Set arguments
fasta = as.character(args[1])

replace = as.logical(as.character(args[2]))
if (is.na(replace)) replace <- FALSE
prefix = as.character(args[3])
if (is.na(prefix)) prefix <- ""
suffix = as.character(args[4])
if (is.na(suffix)) suffix <- ""

## Set arguments (for debugging)
# fasta = "test.fasta"
# replace = FALSE
# prefix = "Contig"
# suffix = "seed_rbcL"

## Check parameters
stopifnot(file.exists(fasta))

## Additional parameters
verbose = TRUE

## Read fasta
fas <- read.FASTA(fasta)
onames <- names(fas)

## Modify headers
for (i in 1:length(onames)) {

  oname <- onames[i]
  seq <- fas[oname]
  
  #gc <- round(100*GC.content(seq))
  gc <- round(GC.content(seq), 2)
  len <- length(seq[[1]])
  
  if (replace) {
    nname <- paste0(prefix, i, suffix, "_", gc, "_gc_", len, "_length")
  } else {
    nname <- paste0(prefix, oname, suffix, "_", gc, "_gc_", len, "_length")
  }

  if (verbose) cat(oname, "->", nname, "\n")
  names(fas)[i] <- nname
  
}

## Safe fasta
write.FASTA(fas, file = fasta)
