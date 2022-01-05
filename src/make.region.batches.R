#!/cluster/home/crameris/bin/R/bin/Rscript

## Generate region subsets

## Description
# generates files with <n> region subsets of a reference FASTA file

## Author
# simon.crameri@env.ethz.ch, Apr 2020

## Load required library

## Get arguments
args = commandArgs(trailingOnly = TRUE)

## Usage
options(warning.length = 2000L)
options(width = 1000)
if (! length(args) %in% c(2)) {
  stop("2 arguments needed:
       REQUIRED
       1) <ref|CHR>:     path to reference FASTA ;
       2) <n|NUM>:  number of regions per batch",
       call.=FALSE)
}

## Get arguments
ref=args[1]
n=as.numeric(args[2])

## Get argumetns (for debugging)
# ref="test.fasta"
# n=20

## Additional arguments
odir="batches"
if (!dir.exists(odir)) dir.create(odir)

## Check arguments
stopifnot(file.exists(ref), n > 0)

## Get region names
l <- readLines(ref)
regions <- gsub("^>", "", l[grep("^>", l)]) ; rm(l)

## Make region subsets
starts <- seq(1, length(regions), by = n)
for (i in starts) {
  r <- na.omit(regions[i:(i+n-1)])
  writeLines(text = r, con = file.path(odir, paste0("batch_", which(starts==i),".txt")))
}

