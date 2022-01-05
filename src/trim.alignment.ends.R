#!/cluster/apps/r/3.6.0_openblas/x86_64/bin/Rscript

######################################################
## TRIM FASTA ALIGNMENT ENDS (MISSINGNESS / NUCDIV) ##
######################################################

## Usage
# run script without arguments to see the arguments it takes

## Author
# simon.crameri@env.ethz.ch, May 2019

## Load required library
suppressPackageStartupMessages(library(ape))

## Get arguments
args = commandArgs(trailingOnly = TRUE)

## Check arguments
options(warning.length = 2000L)
if (! length(args) %in% c(3:9)) {
  stop("9 arguments taken (3 required):
       REQUIRED
       1) <sfile|CHR>:          path to sample mapping file (header and tab-separation expected). Can contain group membership in second column, which will be used to sort the aligned sequences.
       2) <aln|CHR>:            path to .fasta alignment ; 
       3) <odir|CHR>:           base name of output folder ; 
       
       OPTIONAL (if any is provided, it must be provided in this order)
       4) <completeness|NUM>:   trimming is performed from alignment ends until >= <completeness> bases are encountered. Can be an integer > 1 (number of sequences with data) or a fraction (0 keeps all sites, 1 keeps only complete sites). Individuals with no sequence data are not considered for completeness. [DEFAULT: 0.5]
       5) <maxnucdiv|NUM>:      sites with <= maximum nucleotide diversity are kept [DEFAULT: 0.25]
       6) <NA.char|CHR>:        character(s) that denote missing or unknown data [DEFAULT: -] 
       7) <visualize|BOOLEAN>:  if TRUE, the trimming procedure will be visualized step by step (saved in .pdf) [DEFAULT: TRUE]
       8) <plot.width|NUM>:     width of saved plots [DEFAULT: 15]
       9) <plot.height|NUM>:    height of saved plots [DEFAULT: 7]", 
       call.=FALSE)
}

## Set arguments
sfile <- as.character(args[1])
aln <- as.character(args[2])
odir <- as.character(args[3])

completeness <- as.numeric(args[4])
if (is.na(completeness)) completeness <- 0.5
maxnucdiv <- as.numeric(args[5])
if (is.na(maxnucdiv)) maxnucdiv <- 0.25
NA.char <- as.character(args[6])
if (is.na(NA.char)) NA.char <- "-"
visualize <- as.character(args[7])
if (is.na(visualize)) visualize <- TRUE
plot.width <- as.numeric(args[8])
if (is.na(plot.width)) plot.width <- 15
plot.height <- as.numeric(args[9])
if (is.na(plot.height)) plot.height <- 7

# booleans
get.boolean <- function(x, truestrings = c("T","TRUE", "t","true", "True"), 
                        falsestrings = c("F","FALSE","f","false","False")) {
  if (x %in% truestrings) {
    return(TRUE)
  } else {
    if (x %in% falsestrings) {
      return(FALSE)
    } else stop("<", x, "> cannot be interpreted. Please specify one of ", 
                paste(c(truestrings, falsestrings), collapse = ", "), ".")
  }
}
visualize <- get.boolean(visualize)

# ## Set arguments (for debugging)
# sfile <- "mapfile.txt"
# aln <- "test.fasta"
# odir <- "trimmed"
# completeness = 0.5
# maxnucdiv = 0.25
# NA.char = "-"
# visualize = TRUE
# plot.width  = 15
# plot.height = 7

## Additional arguments
# site trimming
ignore.empty <- TRUE               # if TRUE, <completeness> and <maxnucdiv> are interpreted based only on individuals with sequence data

# visualization
cex.indlab <- 0.2 # size of individual labels in image.DNAbin

# nucleotide diversity assessment
levs <- c("a","c","g","t","n","-") # characters expected in input (read.FASTA reads A as a etc.)

# output file
ext <- paste0(".", rev(unlist(strsplit(aln, split = "[.]")))[1]) # file extension
infix <- ".etr" # infix of trimmed output (before file extension)

## Check arguments
stopifnot(file.exists(sfile), file.exists(aln),
          completeness >= 0,
          maxnucdiv >= 0, maxnucdiv <= 1,
          is.logical(visualize),
          plot.width > 0,
          plot.height > 0,
          cex.indlab >= 0)
          
## Define helperfunction
# nucleotide diversity
assess <- function(x, NA.char = "-", levs = c("a","c","g","t","n","-")) {
  tab <- table(factor(x, levels = levs))
  nseq <- length(x)
  ngap <- length(x[x %in% NA.char])
  if (ngap < nseq-1) {
    nucdiv <- (tab["a"] * tab["c"] + tab["a"] * tab["g"] + tab["a"] * tab["t"] + tab["c"] * tab["g"] + tab["c"] * tab["t"] + tab["g"] * tab["t"]) / ((nseq - ngap) * (nseq - ngap -1) / 2)
  } else {
    nucdiv <- 0
  }
  res <- c(ngap, 1 - (ngap / nseq), nucdiv)
  names(res) <- c("ngap", "completeness", "nucdiv")
  return(res)
}

#################################################################

## Read sfile
ds <- read.delim(sfile, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
if (ncol(ds) == 1) ds[,grpvar1] <- rep("Undefined", nrow(ds))
names(ds) <- c("ID","GROUP")
samples <- ds[,"ID"]
ds[ds[,"GROUP"] %in% c(NA, "NA", ""), "GROUP"] <- "NA" # sets <NA> or <> group values to "NA"

## Read fasta alignment
fas <- read.FASTA(aln)
# image.DNAbin(fas)

## Order fasta alignment according to ds$GROUP
if (length(unique(ds$GROUP)) > 1) {
  
  if (!any(ds$ID %in% names(fas))) {
    stop(paste0("no sample in <", sfile, "> found in <", aln, ">\n"))
  }
  
  if (!all(ds$ID %in% names(fas))) {
    mis.id <- ds$ID[which(!ds$ID %in% names(fas))]
    warning(paste0(length(mis.id), " samples in <", sfile, "> not found in <", aln, ">:\n", paste(mis.id, collapse = ", "), "\n"))
  }
  
  if (!all(names(fas) %in% ds$ID)) {
    mis.fas <- names(fas)[!names(fas) %in% ds$ID]
    warning(paste0(length(mis.fas), " samples in <", aln, "> not found in <", sfile, ">:\n", paste(mis.fas, collapse = ", ")), "\n")
  }
  
  ds <- ds[ds$ID %in% names(fas),]
  ds <- ds[order(ds$GROUP, ds$ID),]
  fas <- fas[ds$ID]
}

## Check alignment
flen <- unique(lengths(fas))
if (length(flen) != 1) stop(paste0("Detected unequal sequence lengths, did you align sequences? [", aln, "]"))

## Walk through sequence from both ends until <completeness> AND <maxnucdiv> are met
dd <- as.character(as.matrix(fas))
nseq <- nrow(dd) # number of sequences
nseq.z <- nrow(dd[!apply(dd, 1, FUN = function(x) all(x %in% NA.char)),,drop=F]) # number of non-empty sequences
nsite <- ncol(dd) # number of sites
# turn <completeness> to fraction
if (completeness > 1) {
  if (ignore.empty) {
    completeness <- completeness/nseq.z
  } else {
    completeness <- completeness/nseq
  }
}
if (completeness > 1) completeness <- 1 # if completeness fraction is higher than 1, interpret as full completeness

# assess
if (ignore.empty) {
  dd.z <- dd[!apply(dd, 1, FUN = function(x) all(x %in% NA.char)),,drop=F]
  di <- apply(dd.z, MARGIN = 2, FUN = assess, NA.char = NA.char, levs = levs)
} else {
  di <- apply(dd, MARGIN = 2, FUN = assess, NA.char = NA.char, levs = levs)
}

# left
for (i in 1:flen) {
  # test
  test.completeness <- di["completeness", i] >= completeness
  test.nucdiv <- di["nucdiv", i] <= maxnucdiv
  keep.start <- all(c(test.completeness, test.nucdiv))
  if (keep.start) {
    idx.start <- i
    break
  }
}

# right
for (i in flen:1) {
  # test
  test.completeness <- di["completeness", i] >= completeness
  test.nucdiv <- di["nucdiv", i] <= maxnucdiv
  keep.end <- sum(c(test.completeness, test.nucdiv)) == 2
  if (keep.end) {
    idx.end <- i
    break
  }
}

## Handle sequences that were completely trimmed
if (!exists("idx.start")) idx.start <- 0 #nsite + 1
if (!exists("idx.end")) idx.end <- 0
  
## Trim ends
fas.end <- as.DNAbin(dd[,idx.start:idx.end, drop = FALSE])
# image.DNAbin(fas.end)

## Be verbose
if (idx.start > 0 & idx.end > 0) {
  nkept <- length(idx.start:idx.end)
  cat("trimmed first", idx.start -1, "bases\n")
  cat("trimmed last", nsite - idx.end, "bases\n")
} else {
  nkept <- 0
  cat("trimmed first", nsite, "bases\n")
  cat("trimmed last", nsite, "bases\n")
}
nbase <- sum(base.freq(fas, freq = TRUE))
nbase.end <- sum(base.freq(fas.end, freq = TRUE))
nbase.trimmed <- nbase - nbase.end
nsite.trimmed <- nsite - nkept

cat(paste0("new trimmed alignment coordinates: ", idx.start, ", ", idx.end, "\n"))
cat(paste0("trimmed ", nbase.trimmed, " / ", nbase, " (", round(100*nbase.trimmed/nbase, 2), "%) bases at ", nsite.trimmed, " / ", nsite, " (", round(100*nsite.trimmed/nsite, 2), "%) sites\n"))
cat(paste0("kept ", nkept, " / ", nsite, " (", round(100*nkept/nsite,2), "%) sites\n"))

## Write trimmed alignment (only if it was not completely trimmed)
oname <- paste0(gsub(ext, "", basename(aln)), infix, ext)
if (!dir.exists(odir)) suppressWarnings(dir.create(odir))
if (idx.start > 0 & idx.end > 0) {
  write.FASTA(x = fas.end, file = file.path(odir, oname))
}

## Visualize trimming step by step
if (visualize) {
  pname <- paste0(gsub(ext, "", basename(aln)), infix, ".pdf")
  pdir <- paste0(odir, ".viz")
  if (!dir.exists(pdir)) suppressWarnings(dir.create(pdir))
  pdf(file = file.path(pdir, pname), width = plot.width, height = plot.height)
  image.DNAbin(fas, cex.lab = cex.indlab)
  title(sub = paste("RAW:", basename(aln)))
  if (idx.start > 0 & idx.end > 0) {
    image.DNAbin(fas.end, cex.lab = cex.indlab)
  } else {
    plot(1:10,1:10,type="n",axes=F,ylab="",xlab="")
    text(5,5, "All aligned sequences were completely trimmed!")
  }
  title(sub = paste0("ENDTRIMMED (", idx.start, ":", idx.end, ") ", basename(aln)))
  graphics.off()
}



