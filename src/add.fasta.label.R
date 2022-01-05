#!/cluster/apps/r/3.6.0_openblas/x86_64/bin/Rscript

############################################
## ADD LABEL (METADATA) TO FASTA HEADERS ###
############################################

## Usage
# run script without arguments to see the arguments it takes

## Value: returns a FASTA file with updated headers (oldheader_correspondingnewlabel)

## Author
# simon.crameri@env.ethz.ch, Apr 2019

## Load required library
library(ape)

## Get arguments
args = commandArgs(trailingOnly = TRUE)

## Check arguments
options(warning.length = 2000L)
if (! length(args) %in% c(2:4)) {
  stop("2 arguments required: 
       REQUIRED
       1) <fas|CHR>: path to FASTA file ; 
       2) <meta|CHR>: path to metadata file mapping taxon ID (1st column, must exactly match header in <fas>) to metadata (2nd column). Expects header and tab-delimiter. Can contain additional taxa and different order with respect to <fas> ;
       
       OPTIONAL
       3) <sortit|BOOLEAN>: if TRUE, will sort the sequences according to 1) metadata, 2) old header [DEFAULT: FALSE] ; 
       4) <odir|CHR>: name of output directory. Replaces <fas> if <odir> is the input directory [DEFAULT: meta]",
       call.=FALSE)
}

## Set arguments
fas <- as.character(args[1])
meta <- as.character(args[2])

sortit <- as.character(args[3])

odir <- as.character(args[4])
if (is.na(odir)) odir <- "meta"

## Set arguments (for debugging)
# fas <- "test.fasta"
# meta <- "mapfile.txt"
# sortit <- FALSE
# odir <- getwd()

## Additional arguments
# input reading
falsestrings <- c("F","FALSE","f","false","False")
truestrings <- c("T","TRUE","t","true","True")

# output name
# ext <- paste0(".", sapply(strsplit(basename(fas), split = "[.]"), function(x) {rev(x)[1]})) # FASTA file extension
# oname <- paste0(gsub(paste0(ext, "$"), "", basename(fas)), ".cons", ext) # name of output FASTA (renamed)

## Handle booleans
if (is.na(sortit)) {
  sortit <- FALSE
} else {
  if (sortit %in% falsestrings) sortit <- FALSE
  if (sortit %in% truestrings) sortit <- TRUE
  if (! sortit %in% c(truestrings, falsestrings)) stop("\n<sortit> must be logical ", paste(c(truestrings, falsestrings), collapse = ", "), "!\n")
}

## Check arguments
stopifnot(file.exists(fas),
          file.exists(meta),
          is.logical(sortit))

## Define helperfuncitions
# quiet stop
stopQuietly <- function(...) {
  blankMsg <- sprintf("\r%s\r", paste(rep(" ", getOption("width")-1L), collapse=" "))
  stop(simpleError(blankMsg))
}

##########################################################################################

## Read fasta
fas.old <- read.FASTA(fas)

## Read sample file
dmeta <- read.delim(meta)
names(dmeta)[c(1,2)] <- c("ID", "GROUP")
rownames(dmeta) <- dmeta$ID

## Check conformity of sample file and alignment file
if (!all(names(fas.old) %in% dmeta$ID)) {
  notfound <- dmeta$ID[!dmeta$ID %in% names(fas.old)]
  cat(paste0("Found ", length(notfound), " sample(s) in ", meta, " that are missing in ", basename(fas), ":\n"))
  print(notfound)
  stopQuietly()
}
dmeta <- dmeta[names(fas.old),] # bring to same order

## Rename fasta headers
dmeta$LABEL <-  paste(dmeta$ID, dmeta$GROUP, sep = "_")
fas.new <- fas.old
names(fas.new) <- dmeta$LABEL

## Sort if specified
if (sortit) {
  fas.new <- fas.new[order(dmeta$GROUP, dmeta$ID)]
}

## Write output
if (!dir.exists(odir)) suppressWarnings(dir.create(odir))
write.FASTA(x = fas.new, file = file.path(odir, basename(fas)))
