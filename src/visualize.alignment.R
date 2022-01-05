#!/cluster/apps/r/3.6.0_openblas/x86_64/bin/Rscript

### Visualize (short) FASTA ALIGNMENTS

## Usage: Rscript visualize.alignment.R <alignment.fasta file>
## Author: Simon Crameri, ETHZ, 7.1.2019 

## Get arguments
args = commandArgs(trailingOnly=TRUE)

if (length(args) < 1 | length(args) > 3) {
  stop("At least 1 argument needed: 
        
        required:
        1) path to first alignment.fasta file ; 
	
	optional:
	2) path to second alignment.fasta file ;
        3) path to third alignment.fasta file ;
        ", call.=FALSE)
}

fname1 <- as.character(args[1])
if (!file.exists(fname1)) stop("file ", fname1, "not found!")

fname2 <- as.character(args[2])
if (!is.na(fname2)) {if (!file.exists(fname2)) stop("file ", fname2, "not found!")}

fname3 <- as.character(args[3])
if (!is.na(fname3)) {if (!file.exists(fname3)) stop("file ", fname3, "not found!")}

## Display alignment
library(ape)

x1 <- read.FASTA(fname1)
if (!is.na(fname2)) x2 <- read.FASTA(fname2)
if (!is.na(fname3)) x3 <- read.FASTA(fname3) 

oname <- paste0(gsub(".fasta$", "", basename(fname1)), ".pdf")
pdf(oname, width = 15, height = 7)
image.DNAbin(x1, sub = paste(fname1), cex.lab = 0.3)
if (!is.na(fname2)) {image.DNAbin(x2, sub = paste(fname2), cex.lab = 0.3)}
if (!is.na(fname3)) {image.DNAbin(x3, sub = paste(fname3), cex.lab = 0.3)}
graphics.off()
