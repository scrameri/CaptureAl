#!/cluster/apps/r/3.6.0_openblas/x86_64/bin/Rscript

## Usage: fix.fasta.headers.R

## Value: replaces ' ' and '|' in fasta headers with '_' by default

## Author: simon.crameri@env.ethz.ch, Dec 2021

library(ape)

## Get arguments
args = commandArgs(trailingOnly = TRUE)

## Check arguments
if (length(args) < 1) {
  stop("1 argument needed (2 more optional):
        REQUIRED
        1) path to .fasta file needing header fix (required)
        
        OPTIONAL
        2) string of characters to replace (' |\\|' by default) 
        3) replacement character ('_' by default)
        4) verbose (TRUE by default)", call.=FALSE)
}

## Set parameters
fasta = args[1]
toreplace = as.character(args[2])
replacewith = as.character(args[3])
verbose = as.logical(args[4])

## Check parameters
stopifnot(file.exists(fasta))
if (is.na(toreplace)) toreplace <- " |\\|"
if (is.na(replacewith)) replacewith <- "_"
if (is.na(verbose)) verbose <- TRUE

## Fix parameters (' |\\|' will be interpreted as " |\\\\|" by commandArgs())
if (toreplace == " |\\\\|") toreplace <- " |\\|"

## Read fasta
fas <- read.FASTA(fasta)
fnames <- names(fas)

## Fix headers
nnames <- gsub(toreplace, replacewith, fnames)
repeat{
  rw <- paste0(rep(replacewith, 2), collapse = "")
  nnames <- gsub(rw, replacewith, nnames)
  test <- !any(grepl(rw, nnames))
  if (test) break()
}
if (identical(fnames, nnames)) {
  if (verbose) print(paste0("No characters <", toreplace, "> found - no fixes needed!"))
} else {
  names(fas) <- nnames
  if (verbose) print(paste0("Fixed ", length(which(nnames != fnames)), "/", length(fnames), " headers"))
  
  ## Safe fasta
  nfile <- paste0(tools::file_path_sans_ext(fasta), ".badheaders.", tools::file_ext(fasta))
  invisible(file.rename(from = fasta, to = nfile))
  write.FASTA(fas, file = fasta)
}

