#!/cluster/apps/r/3.6.0_openblas/x86_64/bin/Rscript

### Rename fasta headers in batch mode

## Usage: Rscript rename.fasta.headers.R <path to file> <pattern1> <pattern2> <splitstring pattern> 
## Author: Simon Crameri, ETHZ, Apr 2018 

## Get arguments
args = commandArgs(trailingOnly=TRUE)

## Check arguments
if (!length(args) %in% c(4)) {
  stop("4 arguments needed:
        REQUIRED:
        1) <file|CHR>:         path to .fasta file ; 
        2) <pattern1|CHR|LOG>: regex pattern to replace with ''. If FALSE ('FALSE', 'F', 'false', 'f', 'False'), no replacement is done ;
        3) <pattern2|CHR|LOG>: regex pattern to replace with ''. If FALSE ('FALSE', 'F', 'false', 'f', 'False'), no replacement is done ;
        4) <splitchar|CHR|LOG>: strsplit character (only left part will be retained). If FALSE ('FALSE', 'F', 'false', 'f', 'False'), no string splitting is done", call.=FALSE)
}

## Set arguments
file <- as.character(args[1])
pattern1 <- as.character(args[2]) # this character will be replaced by ""
pattern2 <- as.character(args[3]) # this character will be replaced by ""
splitchar <- as.character(args[4]) # this character will be interpreted as a strsplit split pattern, and only the left side of the splitted string will be kept

## Additional arguments
replacement1 <- ""  # replace pattern1 with replacement1
replacement2 <- ""  # replace pattern2 with replacement2
splitselect <- 1    # keep left side of split string 

## Check parameters
stopifnot(file.exists(file))

## Additional parameters
falsestrings <- c("FALSE", "F", "f", "false", "False")

## Handle FALSE strings
if (pattern1 %in% falsestrings) pattern1 <- FALSE
if (pattern2 %in% falsestrings) pattern2 <- FALSE
if (splitchar %in% falsestrings) splitchar <- FALSE
if (pattern1 == F & pattern2 == F  & splitchar == F) rename <- FALSE else rename <- TRUE

if (rename) {

  ## Rename file
  cat(basename(file), "\n")
  d <- read.table(file, stringsAsFactors = FALSE, sep = "\t")
  idhead <- which(apply(d, 1, grep, pattern = "^>")==1)

  for (j in 1:length(idhead)) {
    # pattern1
    if (pattern1 != FALSE) d[idhead[j],] <- gsub("^>>", ">", paste0(">", gsub(pattern1, replacement1, d[idhead[j],])))

    # pattern2
    if (pattern1 != FALSE) d[idhead[j],] <- gsub("^>>", ">", paste0(">", gsub(pattern2, replacement2, d[idhead[j],])))

    # split pattern
    if (splitchar != FALSE) d[idhead[j],] <- sapply(strsplit(d[idhead[j],], split = splitchar), "[", splitselect)
  }

  write.table(d, file, quote = FALSE, row.names = FALSE, col.names = FALSE)
} else {
  cat("skipping", basename(file), "\n")
}
