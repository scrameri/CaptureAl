#!/cluster/apps/r/3.6.0_openblas/x86_64/bin/Rscript

########################################################
## TRIM FASTA ALIGNMENT (MISSINGNESS / MISASSEMBLIES) ##
########################################################

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
if (! length(args) %in% c(3:13)) {
  stop("13 arguments taken (3 required): 
       REQUIRED
       1) <sfile|CHR>:          path to sample mapping file (header and tab-separation expected). Can contain group membership in second column, which will be used to sort the aligned sequences.
       2) <aln|CHR>:            path to .fasta alignment
       3) <odir|CHR>:           base name of output folder
       
       OPTIONAL (if any is provided, it must be provided in this order)
       4) <completeness|NUM>:   only sites with >= <completeness> are kept. Can be an integer > 1 (number of sequences with data) or a fraction (0 keeps all sites, 1 keeps only complete sites). Individuals with no sequence data are not considered for completeness. [DEFAULT: 0.3]
       5) <winsize|INT>:        size of sliding window used to mask potentially misaligned / misassembled sites in each sequence [DEFAULT: 20]
       6) <wincons|NUM>:        minimum fraction of base matches between the window and consensus sequence. Sites with IUPAC ambiguities (base frequency above 1 and minor allele(s) frequency above 0.3 are not trimmed). If the fraction of base matches is < wincons, that sequence window is masked. If set to 0, no window-based trimming is performed (and the trimming becomes much faster) [DEFAULT: 0.5]
       7) <step|INT>:           step size for sliding window approach. Increase to half the window size to gain speed at the expense of a less saving trimming [DEFAULT: 1]
       8) <internal|BOOLEAN>:   if FALSE, only masked bases at sequence ends (+- a tolerance of 10 bases) are trimmed. If TRUE, all masked bases (also at internal sites) are trimmed [DEFAULT: FALSE]
       9) <NA.char|CHR>:        character(s) that denote missing or unknown data [DEFAULT: -]
       10) <visualize|BOOLEAN>:  if TRUE, the trimming procedure will be visualized step by step (saved in .pdf) [DEFAULT: TRUE]
       11) <parsimony|BOOLEAN>: if TRUE, a the number of missing as well as parsimony informative sites (PIS) in the trimmed output is saved and plotted as density plots (saved in .pdf and .txt) [DEFAULT: FALSE]
       12) <plot.width|NUM>:    width of saved plots [DEFAULT: 15]
       13) <plot.height|NUM>:   height of saved plots [DEFAULT: 7]", 
       call.=FALSE)
}

## Set arguments
sfile <- as.character(args[1])
aln <- as.character(args[2])
odir <- as.character(args[3])

completeness <- as.numeric(args[4])
if (is.na(completeness)) completeness <- 0.3
winsize <- as.numeric(args[5])
if (is.na(winsize)) winsize <- 20
wincons <- as.numeric(args[6])
if (is.na(wincons)) wincons <- 0.5
step <- as.numeric(args[7])
if (is.na(step)) step <- 1
internal <- as.character(args[8])
if (is.na(internal)) internal <- FALSE
NA.char <- as.character(args[9])
if (is.na(NA.char)) NA.char <- "-"
visualize <- as.character(args[10])
if (is.na(visualize)) visualize <- TRUE
parsimony <- as.character(args[11])
if (is.na(parsimony)) parsimony <- FALSE
plot.width <- as.numeric(args[12])
if (is.na(plot.width)) plot.width <- 15
plot.height <- as.numeric(args[13])
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
internal <- get.boolean(internal)
visualize <- get.boolean(visualize)
parsimony <- get.boolean(parsimony)

# ## Set arguments (for debugging)
# sfile <- "mapfile.txt"
# aln <-  "test.fasta"
# odir <- "trimmed"
# 
# completeness = 0.4
# winsize = 20
# wincons = 0.5
# step = 1
# internal = FALSE
# NA.char = "-"
# visualize = TRUE
# parsimony =  FALSE
# plot.width  = 15
# plot.height = 7

## Additional arguments
# site trimming
ignore.empty <- TRUE               # if TRUE, <completeness>, <minallfreq> and <minbasefreq> are interpreted based only on individuals with sequence data

# consensus calculation (in the window-based trimming, each base is compared against unambiguous, i.e. conserved, consensus bases)
levs <- c("a","c","g","t","n","-") # characters expected in input (read.FASTA reads A as a etc.). A warning is issued if there are other or ambiguous characters in the input.
minbasefreq <- 0.2                 # minimum base frequency to call a consensus (gap is returned otherwise). Comparisons between window sequence and consensus are only at sites with an unambiguous consensus base.
minallfreq <- 0.3                  # minimum minor allele frequency to call a IUPAC ambiguity (if none is above minallfreq, a gap is returned). Comparisons between window sequence and consensus are only at sites with an unambiguous consensus base.
ignore.gaps <- TRUE                # if TRUE, minallfreq is computed without considering gaps / missing data

# visualization
cex.indlab <- 0.2                  # size of sample labels in image.DNAbin rows
tol <- 20                          # basepair tolerance to define sequence ends if internal = FALSE. A masked sequence end is defined as a masked subsequence that is connected to another masked subequence by up to <tol> bases AND at some point connected to a sequence end by up to <tol> bases

# output file name
ext <- paste0(".", rev(unlist(strsplit(aln, split = "[.]")))[1]) # file extension
infix <- ".itr"                                                  # infix of trimmed output (before file extension)

## Check parameters
stopifnot(file.exists(sfile), file.exists(aln),
          completeness >= 0,
          winsize >= 1,
          wincons >= 0, wincons <= 1,
          plot.width > 0,
          plot.height > 0,
          cex.indlab >= 0)

# Define IUPAC ambiguity codes
# taken from https://droog.gs.washington.edu/parc/images/iupac.html
IUPACamb <-   c("a", "c", "g", "t", "m",  "r",  "w",  "s",  "y",  "k",  "v",   "h",   "d",   "b",   "n",    "-")
BASES <-      c("a", "c", "g", "t", "ac", "ag", "at", "cg", "ct", "gt", "acg", "act", "agt", "cgt", "acgt", "-")
COMPLEMENT <- c("t", "g", "c", "a", "k",  "y",  "w",  "s",  "r",  "m",  "b",   "d",   "h",   "v",   "n",    "-")
dIUPACamb <- data.frame(IUPACamb = IUPACamb, BASES = BASES, COMPLEMENT = COMPLEMENT, stringsAsFactors = FALSE)

## Define helperfunctions
# site trimming
get.complete.sites <- function(dd, completeness, NA.char = "-", ignore.empty = TRUE) {
  if (ignore.empty) {
    dd.z <- dd[!apply(dd, 1, FUN = function(x) all(x %in% NA.char)),,drop=F]
  } else {
    dd.z <- dd
  }
  nttab <- apply(dd.z, 2, table)
  basefreq <- lapply(nttab, function(x) sum(x[!names(x) %in% NA.char])/sum(x))
  if (completeness > 1) completeness <- completeness / nrow(dd.z)
  if (completeness > 1) completeness <- 1
  keep <- which(unlist(basefreq) >= completeness)
  return(keep)
}

# window-based trimming
mask.bad.windows <- function(sequence, fasta, ref, window, winsize, wincons, NA.char, bases = c("a","c","g","t")) {
  indices <- window:(window+winsize-1)
  
  # query bases (bases, ev. gaps and/or Ns)
  query <- factor(as.character(as.character(fasta[sequence, indices])), levels = dIUPACamb$IUPACamb)
  base.indices <- indices[!query%in% NA.char] # masking is only performed at sites with bases
  
  # consensus bases (bases, ev. with ambiguities) to compare with
  refseq <- factor(ref[indices], levels = dIUPACamb$IUPACamb)
  
  # table of query bases (without gaps or Ns) vs. consensus bases (without gaps, Ns or ambiguity codes):
  # sum of comparisons is sum of window sites with a base (no gap, no IUPAC ambiguity) in both the query and refseq
  tab <- table(refseq, query)[bases,bases]
  ncomp <- sum(tab) # number of comparisons
  nmatches <- sum(diag(tab)) # number of base matches
  # ncomp - nmatches # number of base mismatches
  
  # decide on whether to mask this window
  if (ncomp > 0) {
    mask <- nmatches / ncomp < wincons
  } else {
    mask <- FALSE # sequence windows consisting of gaps only are not masked
  }
  if (mask) dd.mask[sequence, base.indices] <<- 1
}

# consensus
get.IUPACamb <- function(x, db = dIUPACamb, minallfreq = 0.05, minbasefreq = 0.5,
                         levs = c("a","c","g","t","n","-"), NA.char = c("n","-"), 
                         ignore.gaps = TRUE) {
  ## RULES
  # if ignore.gaps = TRUE, (minor) allele frequencies are calculated based on all sequences with data
  # if ignore.gaps = FALSE, (minor) allele frequencies are calculated based on all sequences (including sequences with gaps)
  
  # if there is sufficient base frequency AND at least one base above minallfreq, the ambiguity is calculated
  # if there is insufficient base frequency OR no base above minallfreq, a gap is returned
  
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
  
  # get alleles above minallfreq if minbasefreq is met
  if (highbasefreq) {
    idx <- which(allfreq > 0 & allfreq >= minallfreq)
  } else {
    idx <- integer()
  }
  
  # paste allele names (to match them with db$BASES and return the corresponding db$IUPACamb)
  if (length(idx) > 0) {
    bases <- paste(sort(names(tab[idx])), collapse = "")
  } else {
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

##############################################################################

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

## Site trimming I: find sites with sufficient data
dd <- as.character(as.matrix(fas)) # all sequences
nseq <- nrow(dd)  # number of sequences
nsite <- ncol(dd) # number of sites

# keep only sites with >= completeness data
keep <- get.complete.sites(dd = dd, completeness = completeness, NA.char = NA.char, ignore.empty = ignore.empty)
nkeep <- length(keep)
cat(paste0("kept ", nkeep, " / ", nsite, " (", round(100*nkeep/nsite, 2), "%) sites with completeness above ", completeness, "\n"))
fas.data <- as.matrix(fas)[, keep]
# image.DNAbin(fas.data)

## Window-based trimming: sliding window approach to mask windows with excessive mismatches to the consensus
if (winsize > nkeep) stop("Window size (", winsize, ") is larger than alignment length after site trimming I (", nkeep, ")!")
if (wincons > 0) {
  
  # get consensus at each alignment site
  if (ignore.empty) {
    cons.data <- apply(dd[!apply(dd, 1, FUN = function(x) all(x %in% NA.char)),keep,drop=F], 2, get.IUPACamb, db = dIUPACamb, 
                     minallfreq = minallfreq, minbasefreq = minbasefreq,
                     levs = levs, NA.char = unique(c(NA.char, "n")), ignore.gaps = ignore.gaps)
  } else {
    cons.data <- apply(dd[,keep], 2, get.IUPACamb, db = dIUPACamb, 
                     minallfreq = minallfreq, minbasefreq = minbasefreq,
                     levs = levs, NA.char = unique(c(NA.char, "n")), ignore.gaps = ignore.gaps)
  }
  
  # define windows
  last <- (ncol(fas.data)-(winsize-1))
  allwindows <- 1:last
  windows <- allwindows[seq(1,length(allwindows), by = step)]
  windows <- unique(c(windows, last))
  
  # mask globally (internal and at sequence ends, takes some time)
  dd.mask <- array(data = 0, dim = dim(fas.data), dimnames = list(rownames(fas.data), colnames(fas.data)))
  for (window in windows) {
    # cat(which(windows %in% window), "\n") 
    sapply(seq(nrow(fas.data)), FUN = function(x) {
      mask.bad.windows(sequence = x, fasta = fas.data, ref = cons.data, 
                       window = window, winsize = winsize, wincons = wincons, 
                       NA.char = NA.char)
    })
  }
  
  # mask only at sequence ends +- <tol> (i.e. revert internal masking)
  if (!internal) {
    dd.endmask <- array(data = 0, dim = dim(fas.data), dimnames = list(rownames(fas.data), colnames(fas.data)))
    for (x in seq(nrow(fas.data))) {
      bases <- which(!as.character(fas.data)[x,] %in% NA.char) # indices of bases relative to aln
      masked <- which(dd.mask[x,] == 1)                        # indices of masked bases relative to aln
      
      if (length(masked) > 0) {
        seqlen <- length(bases)
        ismasked <- rep(FALSE, seqlen)                           # has length of bases
        ismasked[bases %in% masked] <- TRUE                      # sum is length of masked bases
        
        # walk from sequence ends towards other end and record indices of masked bases 
        # walk is breaked if there is no masked base at a sequence end +- <tol>
        leftwalk <- numeric()
        for (i in seq(length(bases))) {
          win <- i:(i+tol) ; win <- win[!win > seqlen]
          if (any(ismasked[win])) leftwalk <- c(leftwalk, i) else break
        }
        rightwalk <- numeric()
        for (i in rev(seq(length(bases)))) {
          win <- (i-tol):i ; win <- win[!win < 1]
          if (any(ismasked[win])) rightwalk <- c(rightwalk, i) else break
        }
        idx.end.tol <- unique(sort(c(bases[leftwalk], bases[rightwalk]))) # indices of masked bases at sequence ends (including <tol> unmasked sites)
        idx.end <- idx.end.tol[idx.end.tol %in% masked]                   # indices of masked bases at sequence ends (only masked sites)
        
        # stopifnot(all(idx.end.tol %in% bases))
        # stopifnot(all(idx.end %in% masked))
        
        # define masked bases for sequence <x>
        dd.endmask[x, idx.end] <- 1
      } 
    }
    
    # define masked bases for all sequences
    tomask <- which(dd.endmask == 1) # indices of masked bases
    imask <- sum(rowSums(dd.endmask) > 0) # number of sequences with >= 1 masked base
    smask <- sum(colSums(dd.endmask) > 0) # number of sites with >= 1 masked base
  } else {
    tomask <- which(dd.mask == 1) # indices of masked bases
    imask <- sum(rowSums(dd.mask) > 0) # number of sequences with >= 1 masked base
    smask <- sum(colSums(dd.mask) > 0) # number of sites with >= 1 masked base
  }
} else {
  tomask <- integer() # indices of masked bases
  imask <- 0 # number of sequences with >= 1 masked base
  smask <- 0 # number of sites with >= 1 masked base
}
# image.DNAbin(as.DNAbin(cons.data))

# log masking
ndata <- ncol(fas.data) # number of alignment sites with >= completeness data
nmask <- length(tomask) # number of masked bases (depends on winsize, wincons, step and internal)
sumbases <- sum(base.freq(fas.data, freq = T)) # total number of bases in fas.data
cat(paste0("trimmed ", nmask, " / ", sumbases, " (", round(100*nmask/sumbases, 2), "%) bases at ", smask, " / ", ndata, " (", round(100*smask/ndata, 2), "%) sites in ", imask, " / ", nseq, " (", round(100*imask/nseq, 2), "%) sequences\n"))

# trim masked bases
dd.masked <- as.character(fas.data)
dd.masked[tomask] <- "" # actual masking ind dd.masked / fas.masked (appears white on image.DNAbin)
fas.masked <- as.DNAbin(dd.masked)
# image.DNAbin(fas.masked)

dd.trimmed <-  dd.masked
dd.trimmed[tomask] <- NA.char[1] # actual trimming ind dd.trimmed / fas.trimmed (appears black on image.DNAbin)
fas.trimmed <- as.DNAbin(dd.trimmed)
# image.DNAbin(fas.trimmed)

## Site trimming II: find sites with sufficient data
# keep only sites with >= completeness data
keep.trimmed <- get.complete.sites(dd = dd.trimmed, completeness = completeness, NA.char = NA.char, ignore.empty = ignore.empty)
nkeep.trimmed <- length(keep.trimmed)
nsite.trimmed <- ncol(dd.trimmed)
cat(paste0("kept ", nkeep.trimmed, " / ", nsite.trimmed, " (", round(100*nkeep.trimmed/nsite.trimmed, 2), "%) sites with completeness above ", completeness, "\n"))
fas.end <- as.matrix(fas.trimmed)[, keep.trimmed]
# image.DNAbin(fas.end)

## Write trimmed alignment (only if it was not completely trimmed)
if (!dir.exists(odir)) suppressWarnings(dir.create(odir))
oname <- paste0(gsub(ext, "", basename(aln)), infix, ext)
if (ncol(fas.end) > 0) {
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
  
  image.DNAbin(fas.data, cex.lab = cex.indlab)
  title(sub = paste("REDUCED:", basename(aln)))
  
  if (wincons > 0) {
    fas.cons <- as.matrix(as.DNAbin(cons.data))
    rownames(fas.cons) <- "CONSENSUS"
    image.DNAbin(rbind(fas.cons, fas.masked), cex.lab = cex.indlab)
  } else {
    image.DNAbin(fas.masked, cex.lab = cex.indlab)
  }
  title(sub = paste("MASKED:", basename(aln)))
  
  # image.DNAbin(fas.trimmed, cex.lab = cex.indlab)
  # title(sub = paste("TRIMMED:", basename(aln)))
  
  image.DNAbin(fas.end, cex.lab = cex.indlab)
  title(sub = paste("TRIMMED+REDUCED:", basename(aln)))
  
  graphics.off()
}

## Plot number of parsimony informative sites (EXPERIMENTAL!)
if (parsimony) {
  # count number of Parsimony Informative Sites 
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
    return(c(fmis = as.numeric(sum(tab.na))/length(x), pis = isparsinf))
  }
  
  # pis <- apply(as.character(fas.trimmed), 2, get.PIS)
  # dpis <- data.frame(site = seq(ncol(fas.trimmed)), pis = pis[2,], mis = pis[1,])
  pis <- apply(as.character(fas.end), 2, get.PIS)
  dpis <- data.frame(site = seq(ncol(fas.end)), pis = pis[2,], mis = pis[1,])

  # plot density of PIS along the alignment
  library(ggplot2)
  library(ggpubr)
  plot.pis <- function(dat, bw) {
    p <- ggplot(subset(dat, pis == 1)) +
      geom_density(aes(site), bw = bw) +
      theme_bw() +
      ggtitle(paste0("Parsimony Informative Sites ; bw = ", bw))
    return(p)
  }
  
  d1 <- plot.pis(dpis, bw = 2)
  d2 <- plot.pis(dpis, bw = 10)
  d3 <- plot.pis(dpis, bw = 20)
  d4 <- plot.pis(dpis, bw = "nrd0")
 
  pisname <- paste0(gsub(ext, "", basename(aln)), infix, ".parsimony.pdf")
  pisdir <- paste0(odir, ".parsimony")
  if (!dir.exists(pisdir)) suppressWarnings(dir.create(pisdir))
  pdf(file = file.path(pisdir, pisname), width = plot.width, height = plot.height)
  d0 <- ggplot(dpis) +
    geom_point(aes(site, mis)) +
    ylab("Fraction of missing bases") +
    geom_hline(aes(yintercept = 1-completeness)) +
    theme_bw() +
    ggtitle(basename(aln))
  print(d0)
  print(ggarrange(d1, d2, d3, d4))
  graphics.off()
  
  # write output
  pistabname <- paste0(gsub(ext, "", basename(aln)), infix, ".parsimony.txt")
  write.table(dpis, file = file.path(pisdir, pistabname), row.names = F, quote = F, sep = "\t")
}


