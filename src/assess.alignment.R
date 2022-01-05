#!/cluster/apps/r/3.6.0_openblas/x86_64/bin/Rscript
library("ape")
#library("pegas") # used for theta.s()

# author: Simon Crameri, ETHZ, Apr 2019 

## Get arguments
args = commandArgs(trailingOnly=TRUE)

## Check arguments
if (!length(args) %in% c(1,2)) {
  stop("1 argument needed (2 taken): 
       REQUIRED
       1) <aln|CHR>: path to alignment .fasta ; 
       
       OPTIONAL
       2) <oname|>: output file name [default: alnbase.assess.txt]", call.=FALSE)
}

## Set arguments
aln = as.character(args[1])
oname = as.character(args[2])

## Set arguments (for debugging)
# aln = "file.fasta"
# oname = "test.assess.txt"

## Define output name
alnbase <- tools::file_path_sans_ext(aln)
if (is.na(oname)) oname <- paste0(alnbase, ".assess.txt")

## Additional arguments
na.char <- "-"                     # character denoting gaps
levs <- c("a","c","g","t","n","-") # characters expected in input (read.FASTA reads A as a etc.)

## Helperfunctions
# nucleotide diversity
assess <- function(x, na.char = "-", levs = c("a","c","g","t","n","-")) {
  tab <- table(factor(x, levels = levs))
  nseq <- length(x)
  ngap <- length(x[x %in% na.char])
  if (ngap < nseq) {
    nucdiv <- (tab["a"] * tab["c"] + tab["a"] * tab["g"] + tab["a"] * tab["t"] + tab["c"] * tab["g"] + tab["c"] * tab["t"] + tab["g"] * tab["t"]) / ((nseq - ngap) * (nseq - ngap -1) / 2)
  } else {
    nucdiv <- NA
  }
  return(nucdiv)
}

# parsimony informative sites
get.PIS <- function(x, NA.char = "-") {
  
  # count NA.char, ref alleles, alt alleles
  tab <- table(as.character(x))
  if (any(names(tab) %in% NA.char)) {
    tab.base <- tab[!names(tab) %in% NA.char]
    tab.na <- tab[names(tab) %in% NA.char]
  } else {
    tab.base <- tab
    tab.na <- 0 ; names(tab.na) <- NA.char[1]
  }
  
  # determine whether the site is parsimony informative
  'A site is parsimony-informative if it contains at least two types 
  of nucleotides (or amino acids), and at least two of them occur with a minimum frequency of two.'
  ispoly <- ifelse(length(tab.base) >= 2, TRUE, FALSE)
  if (ispoly) {
    isinf <- ifelse(all(sort(tab.base, decreasing = T)[1:2] >= 2), TRUE, FALSE)
  } else {
    isinf <- FALSE
  }
  
  # return results
  isparsinf <- ispoly & isinf
  return(c(fmis = as.numeric(sum(tab.na))/length(x), poly = ispoly, pis = isparsinf))
}

## Read fasta
fas <- read.FASTA(aln)

## Get alignment length
alnlen <- unique(lengths(fas))
nseqs <- length(fas)
alndim <- nseqs*alnlen

## Check that it is aligned data
stopifnot(length(unique(alnlen)) == 1)
fas <- as.matrix(fas)

## Get alignment statistics
# Number of gaps
ngtot <- sum(as.character(fas) %in% na.char)

# Gap Ratio
# nuctot <- sum(!as.character(fas) %in% na.char)
# ngpernuctot <- ngtot / nuctot 
ngpernuctot <- ngtot / (alnlen * nseqs) # ngtot / (alnlen * (nseqs-1)) is the formula in Stefan Zoller's script

# Total Nucleotide Diversity
# library(pegas)
# nuc.div(fas, pairwise.deletion = FALSE) # gives NaN
nucdivpersite <- apply(as.character(fas), 2, assess, na.char, levs)
nuctotdiversity <- sum(nucdivpersite, na.rm = TRUE)

# Average Nucleotide Diversity per Site
avgnucdivpersite <- nuctotdiversity / alnlen # as in Stefan Zoller's script

# GC content
# base.freq(fas.rec)
gc <- GC.content(fas)

# theta.s
#thetaS <- theta.s(fas) # needs pegas

# Alignment length
length <- alnlen

# Number and proportion of segregationg sites
# nSEG <- length(seg.sites(fas)) # counts monomorphic sites with gaps '-' as polymorphic
pis <- apply(as.character(fas), 2, get.PIS, na.char)
nSEG <- sum(pis["poly",])
pSEG <- nSEG / length

# parsimony informative sites (at least two sequences deviating from the most frequently observed character)
# count number of Parsimony Informative Sites 
nPIS <- sum(pis["pis",])
pPIS <- nPIS / length

# file
file <- basename(aln)

## Compile result
dd <- data.frame(ngtot = ngtot, ngpernuctot = ngpernuctot, 
                 nuctotdiversity = nuctotdiversity, avgnucdivpersite = avgnucdivpersite,
                 gc = gc, 
                 #thetaS = thetaS,
                 length = length,
                 nSEG = nSEG, pSEG = pSEG,
                 nPIS = nPIS, pPIS = pPIS,
                 file = file)

## Write output
#if (file.exists(oname)) first <- FALSE else first <- TRUE
#write.table(x = dd, file = oname, col.names = first, row.names = FALSE, quote = FALSE, sep = "\t", append = !first)
write.table(x = dd, file = oname, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
