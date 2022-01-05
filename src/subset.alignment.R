#!/cluster/apps/r/3.6.0_openblas/x86_64/bin/Rscript

###############################
## SUBSET A FASTA ALIGNMENT ###
###############################

## Usage
# run script without arguments to see the arguments it takes

## Value: returns a subset of a FASTA alignment containing a specified taxon set

## Author
# simon.crameri@env.ethz.ch, May 2019

## Load required library
library(ape)

## Get arguments
args = commandArgs(trailingOnly = TRUE)

## Check arguments
options(warning.length = 2000L)
if (!length(args) %in% c(2:5)) {
  stop("1 argument required (5 taken): 
       REQUIRED
       1) <aln|CHR>: path to .fasta alignment ; 
       2) <sfile|CHR>: sample file (taxon subset). If set to 'complete', will subset all taxa with at least one base of data ;
       
       OPTIONAL
       3) <ignore.gaps|BOOLEAN>: if TRUE, any column consisting of gaps only will be removed from the output [DEFAULT: FALSE] ;
       4) <oname|CHR>: output file name [DEFAULT: <aln>.sub.<ext> ]
       5) <odir|CHR>: output directory [DEFAULT: current directory ]",
       call.=FALSE)
}

## Additional arguments
falsestrings <- c("F","FALSE","f","false","False") # accepted false strings for BOOLEAN argument
truestrings <- c("T","TRUE","t","true","True") # accepted true strings for BOOLEAN argument
ext <- paste0(".", sapply(strsplit(basename(args[1]), split = "[.]"), function(x) {rev(x)[1]})) # alignment file extension
NA.char <- "-" # gap character

## Set arguments
aln <- as.character(args[1])
sfile <- as.character(args[2])

ignore.gaps <- as.character(args[3])
if (is.na(ignore.gaps)) {
  ignore.gaps <- FALSE
} else {
  if (ignore.gaps %in% truestrings) ignore.gaps <- TRUE
  if (ignore.gaps %in% falsestrings) ignore.gaps <- FALSE
}
oname <- as.character(args[4])
if (is.na(oname)) oname <- paste0(gsub(paste0(ext, "$"), "", basename(aln)), ".sub", ext)
odir <- as.character(args[5])
if (is.na(odir)) odir <- getwd()

## Set arguments (for debugging)
# aln <- "test.fasta"
# sfile <- "samples.txt"
# ignore.gaps <- FALSE
# oname <- NA
# odir <- NA

## Check arguments
stopifnot(file.exists(aln),
          is.logical(ignore.gaps))
if (! sfile %in% c("complete")) {
  stopifnot(file.exists(sfile))
}

## Define helperfunctions
# quiet stop
stopQuietly <- function(...) {
  blankMsg <- sprintf("\r%s\r", paste(rep(" ", getOption("width")-1L), collapse=" "))
  stop(simpleError(blankMsg))
}

##########################################################################################

## Read fasta
fas.all <- read.FASTA(aln)

## Read sample file
if (! sfile %in% c("complete")) {
  samples <- readLines(sfile)
  mat <- as.character(as.matrix(fas.all[samples]))
} else {
  mat.all <- as.character(as.matrix(fas.all))
  samples <- rownames(mat.all[apply(mat.all, 1, function(x) {!all(x %in% NA.char)}), , drop = FALSE])
  mat <- mat.all[samples,]
  rm(mat.all)
}

## Check conformity of sample file and alignment file
if (!all(samples %in% names(fas.all))) {
  notfound <- samples[!samples %in% names(fas.all)]
  cat(paste0("Found ", length(notfound), " sample(s) in ", sfile, " that are missing in ", basename(aln), ":\n"))
  print(notfound)
  stopQuietly()
}


## Subset columns if ignore.gaps = TRUE
if (ignore.gaps) mat <- mat[,apply(mat, 2, function(x) {!all(x %in% NA.char)}), drop = FALSE]
fas <- as.DNAbin(mat)
              
## Move any existing output file
ofile <- file.path(odir, oname)
if (file.exists(ofile)) {
  system(paste("mv", ofile, paste0(ofile, ".bak")))
}

## Write consensus
if (!dir.exists(odir)) suppressWarnings(dir.create(odir))
write.FASTA(x = fas, file = ofile)


