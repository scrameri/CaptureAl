#!/cluster/apps/r/3.6.0_openblas/x86_64/bin/Rscript

###############################################
## REMOVE EMPTY TAXA FROM A FASTA ALIGNMENT ###
###############################################

## Usage
# run script without arguments to see the arguments it takes

## Value: returns a a FASTA alignment without taxa lacking any sequence data

## Author
# simon.crameri@env.ethz.ch, May 2019

## Load required library
library(ape)

## Get arguments
args = commandArgs(trailingOnly = TRUE)

## Check arguments
options(warning.length = 2000L)
if (!length(args) %in% c(1:2)) {
  stop("1 arguments required (2 taken): 
       REQUIRED
       1) <aln|CHR>: path to .fasta alignment ; 
       
       OPTIONAL
       2) <oname|CHR>: name of output file [DEFAULT: <alnbase>.sub.<alnformat> ]",
       call.=FALSE)
}

## Additional arguments
ext <- paste0(".", sapply(strsplit(basename(args[1]), split = "[.]"), function(x) {rev(x)[1]})) # alignment file extension
NA.char <- "-" # gap character

## Set arguments
aln <- as.character(args[1])

oname <- as.character(args[2])
if (is.na(oname)) oname <- paste0(gsub(paste0(ext, "$"), "", basename(aln)), ".sub", ext)

## Set arguments (for debugging)
# aln <- "test.fasta"
# oname <- "test.sub.fasta"

## Check arguments
stopifnot(file.exists(aln))

## Define helperfunctions
# quiet stop
stopQuietly <- function(...) {
  blankMsg <- sprintf("\r%s\r", paste(rep(" ", getOption("width")-1L), collapse=" "))
  stop(simpleError(blankMsg))
}

##########################################################################################

## Read fasta
fas.all <- read.FASTA(aln)

## Determine empty taxa
seqlen <- apply(as.character(as.matrix(fas.all)), 1, function(x) {length(x[! x %in% NA.char])})

## Subset fasta
nonzero <- names(seqlen[seqlen > 0])
mat <- as.character(as.matrix(fas.all[nonzero]))
fas <- as.DNAbin(mat)
              
## Move any existing output file
if (file.exists(oname)) {
  system(paste("mv", oname, paste0(oname, ".bak")))
}

## Write subsetted FASTA
#if (!dir.exists(odir)) suppressWarnings(dir.create(odir))
write.FASTA(x = fas, file = oname)
