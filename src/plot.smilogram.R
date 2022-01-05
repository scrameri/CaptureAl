#!/cluster/apps/r/3.6.0_openblas/x86_64/bin/Rscript

## Usage
# run script without arguments to see the arguments it takes

## Author: simon.crameri@env.ethz.ch, Feb 2019

## Load required libraries
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(dplyr))
#suppressPackageStartupMessages(library(adegenet)) # for funky colors
suppressPackageStartupMessages(library(ggplot2))

## Get arguments
args = commandArgs(trailingOnly = TRUE)

## Check arguments
options(warning.length = 2000L)
if (! length(args) %in% c(1:4)) {
  stop("at least 1 argument needed: 
       REQUIRED
       1) <infolder|CHR>: path to folder with alignments ;

       OPTIONAL (if 3rd is given, 2nd must be given, etc.)
       2) <NA.char|CHR>: character(s) that denote missing or unknown data [default: '-'];
       3) <perclass|BOOLEAN>: whether to create smilograms for each class [default: FALSE] ;
       4) <mapfile|CHR>: file mapping sequences to classes [default: NULL] ;
       ", 
       call.=FALSE)
}

## Set parameters
infolder <- as.character(args[1])
stopifnot(dir.exists(infolder))

NA.char <- as.character(args[2])
if (is.na(NA.char)) NA.char <- "-"
perclass <- as.logical(args[3])
if (is.na(perclass)) perclass <- FALSE
mapfile <- as.character(args[4])
if (is.na(mapfile)) mapfile <- NULL
if (perclass & !is.null(mapfile)) stopifnot(file.exists(mapfile))

## Set parameters (debugging)
# infolder <- "mafft.1e-04.0.80.100.5000.80.trim-0.4-20-10-1"
# NA.char = "-"
# perclass <- FALSE
# mapfile <- NULL

## Additional parameters
exp.char = c("a","t","g","c","A","T","G","C")

## Define helperfunctions
# read fastas
read.fastas <- function(infolder, format = "fasta", verbose = FALSE) {
  library(ape)
  fastas <- list.files(infolder, pattern = paste0(".", format), full.names = TRUE)
  flist <- list()
  for (fasta in fastas) {
    if (verbose) cat(which(fastas == fasta), "/", length(fastas), "\n")
    probeID <- sapply(strsplit(basename(fasta), split = "__"), "[", 1)
    probeID <- sapply(strsplit(probeID, split = "_"), function(x) rev(x)[1])
    aln <- read.FASTA(fasta)
    flist[[probeID]] <- aln
    rm(list=c("aln","probeID","fasta"))
  }
  return(flist)
}

# features of alignment columns
getAlnFeatures <- function(aln.list, exp.char = c("a","t","g","c"), NA.char = "-", progress = TRUE) {
  
  ## helperfunctions
  # number of sequences
  nseq <- function(x) {length(x)}
  
  # number of bases
  nbases <- function(x, NA.char) {length(x[!x %in% NA.char])}
  
  # number of alleles
  nall <- function(x, exp.char) {length(unique(x[x %in% exp.char]))}
  
  # parsimony informative
  pis <- function(x, exp.char) {
    tab <- table(x[x %in% exp.char])
    ispoly <- ifelse(length(tab) >= 2, TRUE, FALSE)
    if (ispoly) {
      isinf <- ifelse(all(sort(tab, decreasing = TRUE)[1:2] >= 2), TRUE, FALSE)
    } else {
      isinf <- FALSE
    }
    return(isinf)
  }
  
  # frequency of variant bases
  freqvar <- function(x, exp.char) {
    bases <- x %in% exp.char
    tab <- table(x[bases])/sum(bases)
    varfreq <- sum(tab[-which.max(tab)])
    return(varfreq)
  }
  
  # combine above statistics for each alignment, add position and distance from alignment center
  get.Features <- function(aln, exp.char, NA.char, verbose = FALSE) {
    if (verbose) cat(".")
    library(ape)
    stopifnot(inherits(aln, "DNAbin"))
    if (length(unique(lengths(as.list(aln)))) > 1) stop("Unequal sequence lengths detected!") 
    mat <- as.character(as.matrix(aln))
    x.position <- seq(ncol(mat))
    x.distc <- floor(scale(x.position, scale = FALSE, center = TRUE))[,1]
    x.nseq <- apply(mat, 2, nseq)
    x.nbases <- apply(mat, 2, nbases, NA.char = NA.char)
    x.nall <- apply(mat, 2, nall, exp.char = exp.char)
    x.pis <- apply(mat, 2, pis, exp.char = exp.char)
    x.freqvar <- apply(mat, 2, freqvar, exp.char = exp.char)
    res <- rbind(position = x.position, distc = x.distc, nseq = x.nseq, nbases = x.nbases, 
                 nall = x.nall, pis = x.pis, freqvar = x.freqvar)
    return(t(res))
  }
  
  ## get features for list of alignments
  lres <- lapply(aln.list, get.Features, exp.char = exp.char, NA.char = NA.char, verbose = progress)
  
  ## bind to data.frame
  library(dplyr)
  df <- bind_rows(lapply(lres, data.frame), .id = "locus")
  df$mis <- 1 - (df$nbases / df$nseq)
  attr(df, "class") <- c("data.frame", "aln.features")
  return(df)
}


# smoothing over median frequency of variant bases
get.smilogram <- function(df, x = "distc", y = "freqvar", fun = median, f = 0.05, variant.only = TRUE, pis.only = FALSE, 
                          max.mis = NULL, max.nall = NULL, max.dist = NULL) {
  
  stopifnot(inherits(df, "aln.features"))
  
  # Filter input
  if (variant.only) df <- df[df$freqvar > 0,]
  if (pis.only) df <- df[df$pis == 1,]
  if (!is.null(max.mis)) {df <- df[df$mis <= max.mis,]}
  if (!is.null(max.nall)) df <- df[df$nall <= max.nall,]
  if (!is.null(max.dist)) df <- df[abs(df$distc) <= max.dist,]
  
  # Running median (median at each distance from center of alignment)
  agg <- aggregate(df[,y], by = list(df[,x]), FUN = fun)
  
  # Smooth through running median
  dl <- lowess(agg$Group.1, agg$x, f = f)
  
  # Compile summarized data
  dx <- data.frame(agg$Group.1, agg$x, dl$y)
  names(dx) <- c(x, paste("agg", y, sep = "_"), paste("smoothed", y, sep = "_"))
  attr(dx, "x") <- x
  attr(dx, "y") <- y
  attr(dx, "fun") <- substitute(fun)
  attr(dx, "f") <- f
  attr(dx, "variant.only") <- variant.only
  attr(dx, "pis.only") <- pis.only
  attr(dx, "max.mis") <- max.mis
  attr(dx, "max.nall") <- max.nall
  attr(dx, "max.dist") <- max.dist
  return(dx)
}


# apply above functions and plot results
process <- function(aln.list, level = "all", infolder, exp.char, NA.char, progress = TRUE) {
  
  nseqs <- nrow(as.matrix(aln.list[[1]]))
  
  ## Get alignment features
  cat("computing alignment features in", length(aln.list), "alignments, each with", nseqs, "sequences...\n")
  df <- getAlnFeatures(aln.list = aln.list, exp.char = exp.char, NA.char = NA.char, progress = progress)
  save(df, file = paste(infolder, gsub(" ", "_", level), "rda", sep = "."))
  
  # Get freqvar smilogram data (smooth over per-locus-median frequency of variant bases as a function of distance from alignment center)
  dx <- get.smilogram(df, x = "distc", y = "freqvar", fun = median, variant.only = TRUE, max.dist = NULL)
  
  # Get mis smilogram data (smooth over per-locus-median frequency of missing bases as a function of distance from alignment center)
  dmis <- get.smilogram(df, x = "distc", y = "mis", fun = median, variant.only = FALSE, max.dist = NULL)
  
  ## Plots
  cat("\nplotting results...\n")
  pdf(file = paste(basename(infolder), gsub(" ", "_", level), "pdf", sep = "."), height = 10, width = 10)
  
  library(ggplot2)
  p0 <- ggplot(dmis, aes(distc, agg_mis)) +
    geom_line() +
    geom_line(aes(distc, smoothed_mis), col = "tomato", size = 2) +
    labs(x = "Position in alignment (bp)", y = "Frequency of missing bases") +
    ggtitle(level) +
    theme_bw()
  
  print(p0)
  
  p1 <- ggplot(dx, aes(distc, agg_freqvar)) +
    geom_line() +
    geom_line(aes(distc, smoothed_freqvar), col = "tomato", size = 2) +
    labs(x = "Position in alignment (bp)", y = "Frequency of variant bases") +
    ggtitle(level) +
    theme_bw()
  
  print(p1)
  
  # # per locus
  # set.seed(1040223)
  # locs <- sample(unique(df$locus), 100)
  # library(adegenet) # for funky colors
  # 
  # p2 <- ggplot(subset(df, freqvar > 0 & locus %in% locs), aes(distc, freqvar)) +
  #   geom_point(aes(col = freqvar), size = 0.2) +
  #   geom_smooth(method = "loess", formula = y~x, se = FALSE) +
  #   scale_color_continuous(low = funky(4)[1], high = funky(4)[4]) +
  #   facet_wrap(~locus) +
  #   theme_bw()
  # 
  # print(p2)
  
  dev.off()

}

## ----------------------------------------------------------------------- ##

cat("reading fastas...\n")
flist <- read.fastas(infolder)
save(flist, file = paste0(infolder, ".flist.rda"))

## Get metadata
if (perclass) {
  
  ## Read metadata
  dmap <- read.delim(mapfile, header = T)
  names(dmap) <- c("ID", "CLASS")
  
  ## Merge metadata
  dd <- merge(data.frame(ID=names(flist[[1]])), dmap, by = "ID", all.x = T, all.y = F, sort = F)
  rownames(dd) <- dd$ID
  
  ## Get Classes with >= 2 sequences
  dtab <- table(droplevels(dd$CLASS))
  dtab <- dtab[dtab > 1]
  
  ## Get Smilogram per Class
  for (level in names(dtab)) {
    cat("subsetting fastas for class level(s):", paste(level, collapse = ", "), "...\n")
    idsel <- which(dd$CLASS %in% level)
    aln.list <- lapply(flist, function(x) {as.matrix(x)[idsel,]})
    process(aln.list = aln.list, level = level, infolder = infolder, exp.char = exp.char, NA.char = NA.char)
  }
  
} else {
  aln.list <- flist
  level <- "all"
  
  ## Get Smilogram
  process(aln.list = aln.list, level = level, infolder = infolder, exp.char = exp.char, NA.char = NA.char)
}
