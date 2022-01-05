#!/cluster/apps/r/3.6.0_openblas/x86_64/bin/Rscript

#################################
## EXTRACT ALL FASTA SEQUENCES ##
#################################

## Usage
# run script without arguments to see the arguments it takes

## Author
# simon.crameri@env.ethz.ch, Apr 2019

## Load required library
suppressPackageStartupMessages(library(ape))

## Get arguments
args = commandArgs(trailingOnly = TRUE)

## Check arguments
# options(warning.length = 2000L)
# options(width = 1000)
if (! length(args) %in% c(1,2)) {
  stop("2 argument taken (1 required):
       REQUIRED
       1) <fasta|CHR>: path to .fasta file ; 
       
       OPTIONAL 
       2) <odir|CHR>: name of output directory",
       call.=FALSE)
}

## Set arguments
fasta <- as.character(args[1])
odir <- as.character(args[2])
if (is.na(odir)) odir <- paste0(gsub(paste0(paste0(".", rev(unlist(strsplit(fasta, split = "[.]")))[1]), "$"), "", fasta), ".seqs")
  
## Set arguments (for debugging)
#fasta <- "orig.fragments.merged100_6555.fasta"
#odir <- paste0(gsub(paste0(paste0(".", rev(unlist(strsplit(fasta, split = "[.]")))[1]), "$"), "", fasta), ".seqs") # name of output directory

## Additional arguments
ext <- paste0(".", rev(unlist(strsplit(fasta, split = "[.]")))[1]) # extension of input fasta
verbose <- TRUE

## Check arguments
stopifnot(file.exists(fasta))

## Read FASTA
fas <- read.FASTA(fasta)

## Check headers (all must be unique)
stopifnot(!any(duplicated(names(fas))))

## Write output
if (!dir.exists(odir)) dir.create(odir)
for (i in seq(length(fas))) {
  fname <- names(fas)[i]
  ofile <- file.path(odir, paste0(fname, ext))
  if (!file.exists(ofile)) {
    if (verbose) cat(fname, "\n")
    write.FASTA(fas[i], file = ofile)
  } else {
    if (verbose) cat("skipping ", fname, " (file exists)\n")
  }
}
