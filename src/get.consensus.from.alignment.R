#!/cluster/apps/r/3.6.0_openblas/x86_64/bin/Rscript

########################################################################
## GENERATE A CONSENSUS SEQUENCE FROM A (SUBSAMPLED) FASTA ALIGNMENT ###
########################################################################

## Usage
# run script without arguments to see the arguments it takes

## Value: returns a consensus sequence from a (subsampled) alignment fasta (and a plot if visualize = T and a log if logfile = T)

## Details
# only taxa in <sfile> are considered for deriving the consensus
# the consensus is as follows:
#
# - ambiguity coding follows the IUPAC rules as stated here:
#   https://droog.gs.washington.edu/parc/images/iupac.html
#
# - if ignore.empty = TRUE, all empty sequences are not considered for allele and base frequency calculation
#
# - if ignore.gaps = TRUE
#   - (minor) allele frequencies are calculated based on all base (but no gap) characters at that site
#   - sites with a gap consensus are removed from the consensus sequence (which will be shorter than the alignment)
#
# - if ignore.gaps = FALSE
#   - (minor) allele frequencies are calculated based on all base and gap characters at that site
#   - sites with a gap consensus are retained in the consensus sequence
#
# - if ignore.n = TRUE
#   - sites with a 'n' consensus (all actg above minallfreq) are removed from the consensus sequence (which will be shorter than the alignment)
#
# - if there is sufficient base frequency AND at least one allele above minallfreq, the ambiguity is calculated as follows:
#   - if only 1 allele is at >= minallfreq, that allele is the consensus
#   - if two or more alleles are at >= minallfreq <= 0.5, the respective IUPAC consensus is returned
#   - if minallfreq > 0.5, the most commonly observed base is returned irrespective of its frequency (if two bases are equally common, one is randomly sampled)
#
# - if there is insufficient base frequency OR no allele above minallfreq, a gap is returned
#
# - if ignore.gaps = TRUE, consensus positions with a '-' character are removed (i.e., sites with little base or allele frequency)
# - if ignore.n = TRUE, consensus positions with a 'n' character are removed (i.e., sites with a, c, g and t at high allele frequency)

## Author
# simon.crameri@env.ethz.ch, May 2019

## Load required library
library(ape)

## Get arguments
args = commandArgs(trailingOnly = TRUE)

## Check arguments
options(warning.length = 2000L)
if (! length(args) %in% c(2, 12)) {
  stop("2 arguments required (12 taken): 
       REQUIRED
       1) <aln|CHR>: path to .fasta alignment ; 
       2) <sfile|CHR>: sample file (consensus will be derived with regard to these taxa) ;

       OPTIONAL (if any is provided, all must be provided in this order)
       3) <odir|CHR>: base name of output folder [default: cons]; 
       4) <minallfreq|NUM>: minimum allele count (if >1) or frequency to call an alternate allele, which leads to an ambiguous consensus (based on available bases only if ignore.gaps = TRUE). If minallfreq > 0.5, the most commonly observed base is returned irrespective of its frequency (if two bases are equally common, one is randomly sampled) [default: 0.05] ;
       5) <minbasefreq|NUM>: minimum number (if > 1) or frequency of bases to call a consensus (a gap is returned otherwise). Individuals with no sequence data are not considered for completeness. [default: 0.5] ;
       6) <ignore.gaps|BOOLEAN>: whether to ignore gaps when computing a consensus [default: TRUE] ;  
       7) <ignore.n|BOOLEAN>: whether to omit all 'n' characters from the consensus sequence [default: TRUE] ;
       8) <prefix|CHR>: prefix of consensus sequence (string before basename of <aln>) [default: cons_] ; 
       9) <suffix|CHR>: suffix of consensus sequence (string after basename of <aln>) [default: ''];
       10) <visualize|BOOLEAN>: whether to visualize the consensus [default: TRUE] ;
       11) <plot.width|NUM>: plot width [default: 15] ; 
       12) <plot.height|NUM>: plot height [default: 7]", 
       call.=FALSE)
}

## Set arguments
aln <- as.character(args[1])
sfile <- as.character(args[2])

odir <- as.character(args[3])
if (is.na(odir)) odir <- "cons"
minallfreq <- as.numeric(as.character(args[4]))
if (is.na(minallfreq)) minallfreq <- 0.05
minbasefreq <- as.numeric(as.character(args[5]))
if (is.na(minbasefreq)) minbasefreq <- 0.5
ignore.gaps <- as.logical(as.character(args[6]))
if (is.na(ignore.gaps)) ignore.gaps <- TRUE
ignore.n <- as.logical(as.character(args[7]))
if (is.na(ignore.n)) ignore.n <- TRUE
prefix <- as.character(args[8])
if (is.na(prefix)) prefix <- "cons_"
suffix <- as.character(args[9])
if (is.na(suffix)) suffix <- ""
visualize <- as.logical(as.character(args[10]))
if (is.na(visualize)) visualize <- TRUE
plot.width <- as.numeric(as.character(args[11]))
if (is.na(plot.width)) plot.width <- 15
plot.height <- as.numeric(as.character(args[12]))
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
ignore.gaps <- get.boolean(ignore.gaps)
ignore.n <- get.boolean(ignore.n)
visualize <- get.boolean(visualize)

## Set arguments (for debugging)
# aln <- "test.fasta"
# sfile <- "samples.txt"
# odir <- "cons"
# minallfreq <- 0.3
# minbasefreq <- 0.2
# ignore.gaps <- TRUE
# ignore.n <- TRUE
# prefix <- "cons_"
# suffix <- ""
# visualize <- TRUE
# plot.width <- 15
# plot.height <- 7

## Additional arguments
# handing empty sequences
ignore.empty <- TRUE               # if TRUE, <minallfreq> and <minbasefreq> are interpreted based only on individuals with sequence data

# expected base characters
NA.char <- c("n", "-") # character denoting gaps or unknown characters (will be ignored at the consensus calculation step if ignore.gaps = TRUE)
levs <- c("a","c","g","t","n","-") # characters expected in input (read.FASTA reads A as a etc.). A warning is issued if there are other characters present, or if ambiguous characters are detected in the input.

# plotting parameters
cex.indlab <- 0.2 # size of individual labels in image.DNAbin

# output directories and files
logfile <- TRUE # whether to log the consensus information to <.logs>
ext <- paste0(".", sapply(strsplit(basename(aln), split = "[.]"), function(x) {rev(x)[1]})) # alignment file extension
ldir <- paste0(odir, ".logs") # name of LOG directory
pdir <- paste0(odir, ".viz") # name of PDF directory
oname <- paste0(gsub(paste0(ext, "$"), "", basename(aln)), ".cons", ext) # name of output FASTA (consensus)
pname <- paste0(gsub(paste0(ext, "$"), "", basename(aln)), ".cons.pdf") # name of output PDF (consensus)
lname <- paste0(gsub(paste0(ext, "$"), "", basename(aln)), ".cons.log") # name of output LOG

## Handle empty strings (prefix and suffix)
falsestrings = c("F","FALSE","f","false","False")
if (prefix %in% falsestrings) prefix <- ""
if (suffix %in% falsestrings) suffix <- ""

# Define IUPAC ambiguity codes
# taken from https://droog.gs.washington.edu/parc/images/iupac.html
IUPACamb <-   c("a", "c", "g", "t", "m",  "r",  "w",  "s",  "y",  "k",  "v",   "h",   "d",   "b",   "n",    "-")
BASES <-      c("a", "c", "g", "t", "ac", "ag", "at", "cg", "ct", "gt", "acg", "act", "agt", "cgt", "acgt", "-")
COMPLEMENT <- c("t", "g", "c", "a", "k",  "y",  "w",  "s",  "r",  "m",  "b",   "d",   "h",   "v",   "n",    "-")
dIUPACamb <- data.frame(IUPACamb = IUPACamb, BASES = BASES, COMPLEMENT = COMPLEMENT, stringsAsFactors = FALSE)

## Check arguments
stopifnot(file.exists(aln),
          any(c(is.na(sfile), file.exists(sfile))),
          minallfreq >= 0,
          minbasefreq >= 0,
          is.logical(ignore.gaps),
          is.logical(visualize),
          plot.width > 0,
          plot.height > 0,
          is.character(NA.char), is.character(levs),
          is.numeric(cex.indlab), cex.indlab >= 0,
          is.logical(logfile))

##########################################################################################

## Define helperfunctions
# consensus
get.IUPACamb <- function(x, db = dIUPACamb, minallfreq = 0.05, minbasefreq = 0.5,
                         levs = c("a","c","g","t","n","-"), NA.char = c("n","-"), 
                         ignore.gaps = TRUE) {
  ## RULES
  # if ignore.gaps = TRUE, (minor) allele frequencies are calculated based on all sequences with data
  # if ignore.gaps = FALSE, (minor) allele frequencies are calculated based on all sequences (including sequences with gaps)
  
  # if there is sufficient base frequency AND at least one base above minallfreq (up to 0.5), the ambiguity is calculated
  # if there is insufficient base frequency OR no base above minallfreq (up to 0.5), a gap is returned
  
  # if there is sufficient base frequency and minallfreq is above 0.5, the most commonly observed base is returned. If two bases are equally common, one is randomly sampled.
  
  # check that the string does not include any ambiguity codes beside n
  if (any(x %in% db[nchar(db$BASES)>1 & db$IUPACamb != "n", "IUPACamb"])) {
    warning(paste0(names(table(x)), collapse = ""), " contains ambiguity characters!")
  }
  
  # turn minallfreq and minbasefreq to fractions
  if (minallfreq > 1) minallfreq <- minallfreq/length(x) 
  if (minbasefreq > 1) minbasefreq <- minbasefreq/length(x)
  if (minallfreq > 1) minallfreq <- 1 # if minallfreq was specified over the upper limit
  if (minbasefreq > 1) minbasefreq <- 1 # if minbasefreq was specified over the upper limit
  
  # get allele frequencies (only alleles above minallfreq will be considered)
  tab.base.gaps <- table(factor(x, levels = levs[!levs %in% "n"]))
  tab.base <- table(factor(x[!x %in% NA.char], levels = levs[!levs %in% NA.char])) 
  if (ignore.gaps) {
    tab <- tab.base
  } else {
    tab <- tab.base.gaps
  }
  allfreq <- tab/sum(tab) # allele frequency depends on ignore.gaps
  
  # check if site is above minbasefreq
  basefreq <- sum(tab.base)/sum(tab.base.gaps) # basefreq does not depend on ignore.gaps
  highbasefreq <- basefreq >= minbasefreq
  
  # get alleles above minallfreq (up to 0.5) if minbasefreq is met
  if (highbasefreq) {
    idx <- which(allfreq > 0 & allfreq >= minallfreq)
  } else {
    idx <- integer()
  }
  
  # paste allele names (to match them with db$BASES and return the corresponding db$IUPACamb)
  if (minallfreq <= 0.5 & length(idx) > 0) {
    # get consensus
    bases <- paste(sort(names(tab[idx])), collapse = "")
  } 
  if (minallfreq <= 0.5 & length(idx) == 0) {
    # get gap
    bases <- "-"
  }
  if (minallfreq > 0.5 & highbasefreq) {
    # get most commonly observed allele if minbasefreq is met
    tab.base.max <- tab.base[tab.base == max(tab.base)]
    bases <- sample(names(tab.base.max), size = 1)
  }
  if (minallfreq > 0.5 & !highbasefreq) {
    # get gap
    bases <- "-"
  }
  
  # match allele names with db$IUPACamb
  recode <- db[db$BASES == bases, "IUPACamb"]
  if (length(recode) != 1) {
    recode <- "n"
    warning("WARNING: no consensus could be found for ", paste(x, collapse = ""), ". Returned 'n'.")
  }  
  return(recode)
}

# quiet stop
stopQuietly <- function(...) {
  blankMsg <- sprintf("\r%s\r", paste(rep(" ", getOption("width")-1L), collapse=" "))
  stop(simpleError(blankMsg))
}

##########################################################################################

## Read fasta
fas.all <- read.FASTA(aln)

## Read sample file
samples <- readLines(sfile)

## Check conformity of sample file and alignment file
if (!all(samples %in% names(fas.all))) {
  notfound <- samples[!samples %in% names(fas.all)]
  warning(length(notfound), " samples in <", sfile, "> not found in <", basename(aln), ">:\n", 
          paste(notfound, collapse = ", "), "\n")
}

## Subset fasta
fas <- fas.all[samples]

## Determine minimum allele frequency (if minimum allele count provided)
nseq <- length(fas)

## Compute consensus base
mat <- as.character(as.matrix(fas))
if (ignore.empty) {
  matcons <- apply(mat[!apply(mat, 1, FUN = function(x) all(x %in% NA.char)),,drop=F], 2, get.IUPACamb, db = dIUPACamb, 
                   minallfreq = minallfreq, minbasefreq = minbasefreq,
                   levs = levs, NA.char = unique(c(NA.char, "n")), ignore.gaps = ignore.gaps)
} else {
  matcons <- apply(mat[,], 2, get.IUPACamb, db = dIUPACamb, 
                   minallfreq = minallfreq, minbasefreq = minbasefreq,
                   levs = levs, NA.char = unique(c(NA.char, "n")), ignore.gaps = ignore.gaps)
}

## Remove 'n' from consensus if ignore.n = TRUE
if (ignore.n) matcons <- matcons[!matcons %in% "n"]

## Revmoe '-' from consensus if ignore.gaps = TRUE
if (ignore.gaps) matcons <- matcons[!matcons %in% "-"]

## Create named DNAbin object from consensus
cons <- as.DNAbin(list(matcons))
names(cons) <- paste0(prefix, gsub(paste0(ext, "$"), "", basename(aln)), suffix)

## Create log file
if (logfile) {
  #taborig <- table(factor(unlist(as.character(fas)), levels = dIUPACamb$IUPACamb))
  consseq <- as.character(unlist(as.character(cons)))
  tabambig <- table(factor(consseq, levels = dIUPACamb$IUPACamb))
  
  ntax <- length(samples)
  nall <- length(fas.all)
  laln <- ncol(mat)
  lcons <- length(matcons)
  
  dlog <- c(
    paste0("alignment: ", aln),
    paste0("sample file: ", sfile),
    paste0(""),
    paste0("output consensus FASTA: ", file.path(odir, oname)),
    paste0("output PDF (if any): ", file.path(pdir, pname)),
    paste0("output LOG (if any): ", file.path(ldir, lname)),
    paste0(""),
    paste0("minimum allele frequency: ", minallfreq),
    paste0("minimum base frequency: ", minbasefreq),
    paste0("ignore.gaps: ", as.character(ignore.gaps)),
    paste0("ignore.n: ", as.character(ignore.n)),
    paste0(""),
    paste0("prefix: ", prefix),
    paste0("suffix: ", suffix),
    paste0("visualize: ", as.character(visualize)),
    paste0("plot.width: ", plot.width),
    paste0("plot.height: ", plot.height),
    paste0(""),
    paste0("number of considered samples: ", ntax, "/", nall, " (", round(100*ntax/nall,2), "%)"),
    paste0("alignment length: ", laln),
    paste0("consensus length: ", lcons, " (", round(100*lcons/laln,2), "%)"),
    paste0("")
  )
  
  dd <- rbind(dIUPACamb$IUPACamb,
              dIUPACamb$BASES,
              #taborig,
              tabambig,
              round(tabambig/ncol(mat), 4))
  
  if (!dir.exists(ldir)) suppressWarnings(dir.create(ldir))
  writeLines(text = dlog, con = file.path(ldir, lname))
  write.table(dd, file = file.path(ldir, lname), row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE, sep = "\t")
}

## Write consensus
if (!dir.exists(odir)) suppressWarnings(dir.create(odir))
write.FASTA(x = cons, file = file.path(odir, oname))

## Visualize consensus
if (visualize) {
  if (!dir.exists(pdir)) suppressWarnings(dir.create(pdir))
 
  pdf(file = file.path(pdir, pname), width = plot.width, height = plot.height)
  
  image.DNAbin(fas, cex.lab = cex.indlab)
  title(sub = paste("ALIGNMENT:", basename(aln)))
  
  # ape < 5.2 cannot display a single sequence
  image.DNAbin(as.matrix(rbind(as.matrix(cons),as.matrix(cons))),
               show.labels =  FALSE)
  title(sub = paste("CONSENSUS:", basename(oname)))
  
  graphics.off()
}

