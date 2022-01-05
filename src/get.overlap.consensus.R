#!/cluster/apps/r/3.6.0_openblas/x86_64/bin/Rscript

###########################
## GET OVERLAP CONSENSUS ##
############################

## Usage
# run script without arguments to see the arguments it takes

## Author
# simon.crameri@env.ethz.ch, Apr 2019

## Load required library
suppressPackageStartupMessages(library(ape))

## Get arguments
args = commandArgs(trailingOnly = TRUE)

## Check arguments
if (! length(args) %in% c(1)) {
  stop("1 arguments taken (1 required):
       REQUIRED
       1) <aln|CHR>: alignment of overlapping loci",
       call.=FALSE)
}

## Set arguments
aln <- as.character(args[1])

## Set arguments (for debugging)
# aln <- "test.mafft.align.fasta"

## Additional arguments
outdir <- dirname(aln) # where output .fasta and .pdf is written to
ext <- paste0(".", rev(unlist(strsplit(aln, split = "[.]")))[1]) # file extension
obase <- gsub(paste0(".aln", ext, "$"), ".all.aln", aln) # basename of output files (obase.ext, obase.log, obase.pdf)

NA.char <- "-"         # character encoding gaps
cex.lab <- 0.2         # plot labels size
cex.sub <- 0.5	       # plot subtitle size
plot.width <- 15       # plot width 
plot.height <- 7       # plot height
verbose <- FALSE

## Check arguments
stopifnot(file.exists(aln),
          dir.exists(outdir))

## Read alignment
faln <- as.matrix(read.FASTA(aln))

## Get SFS
laln <- apply(as.character(faln), 2, function(x, na.char = NA.char) {names(sort(table(x[!x %in% na.char]), decreasing = TRUE))})

## Loop through samples
samples <- unique(rownames(faln))

caln <- as.character(as.matrix(faln)[0,])
qaln <- as.numeric(as.matrix(faln)[0,])
for (sample in samples) {
  if (verbose) cat(sample, "\n")
  
  # potentially overlapping contigs
  saln <- as.character(faln)[rownames(faln) %in% sample,]
  
  # number of base states at each position (0 = gap, 1 = 1 base, 2 = 2 bases, etc.)
  naln <- apply(saln, 2, function(x, na.char = NA.char) {length(unique(x[!x %in% na.char]))})
  qaln <- rbind(qaln, naln)
  rownames(qaln)[nrow(qaln)] <- sample
  
  # if there is at least 1 base (naln > 0), take that base (if naln == 1) or take the base that is most frequent at that alignment site (if naln > 1)
  scons <- character()
  for (j in seq(length(naln))) {
    # if there is at least 1 base
    if (j %in% which(naln > 0)) {
      uniq <- unique(saln[,j])
      uniq <- uniq[!uniq %in% NA.char] # bases

      # if there are >1 bases, the base that is most frequent at that alignment site is the consensus
      if (length(uniq) > 1) {
        for (z in seq(length(laln[[j]]))) {
          cons <- uniq[uniq == laln[[j]][z]]
          if (length(cons) == 1) break
        }
      } else {
        cons <- uniq # if there is 1 base, that base in the consensus
      }
    } else {
      cons <- NA.char # if there is no base, a gap character is the consensus
    }
    scons <- c(scons, cons)
  }
  caln <- rbind(caln, scons)
  rownames(caln)[nrow(caln)] <- sample
}
caln <- as.DNAbin(caln)

## Score it
alndim <- nrow(qaln)*ncol(qaln)

# number of bases in alignment
nbases <- sum(qaln > 0)

# number of mismatches between overlapping contigs of the same individual
npara <- sum(qaln > 1)

# alignment score
qscore <- 1- (npara / nbases)

qscores <- c(paste0("Number of bases:\t", nbases, "\t", paste0(obase, ext)),
             paste0("Fraction of bases\t", nbases/alndim, "\t", paste0(obase, ext)),
             paste0("Number of paralogs:\t", npara, "\t", paste0(obase, ext)),
             paste0("Fraction of paralogs:\t", npara/alndim, "\t", paste0(obase, ext)),
             paste0("Overall alignment score:\t", qscore, "\t", paste0(obase, ext)))

# qmean <- apply(qaln, 2, mean)
# qposmean <- apply(qaln, 2, function(x) {mean(x[x > 0])})
# qsum <- sum(qmean > 1)
# qscore <- qsum/length(qmean)
# 
# qpossum <- sum(qposmean > 1)
# qposscore <- qpossum/length(qposmean)
# 
# qscores <- c(paste0("Mean number of non-homolog sites\t", qsum),
#              paste0("Mean number of non-homolog bases\t", qpossum),
#              paste0("Fraction of non-homolog sites\t", qscore),
#              paste0("Fraction of non-homolog bases\t", qposscore))

# GOOD summary(qmean)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.01786 0.05357 0.25000 0.41568 0.82143 0.92857 

# BAD summary(qmean)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.01786 0.01786 0.01786 0.25105 0.26786 1.82143 

# Visualize
# image(seq(ncol(qaln)),seq(nrow(qaln)),t(qaln))
# qmean2 <- rbind(qmean,qmean)
# image(seq(length(qmean)),1:2,t(qmean2))

## Write merged alignment
write.FASTA(caln, file = paste0(obase, ext))

## Write alignment score
writeLines(qscores, con = paste0(obase, ".log"))

## Visualize
pdf(paste0(obase, ".pdf"), width = plot.width, height = plot.height)
image.DNAbin(faln, sub = paste0("OVERLAPPING (SCORE = ", round(qscore, 2), "): ", basename(aln)), cex.sub = cex.sub, cex.lab = cex.lab)
image.DNAbin(caln, sub = paste0("MERGED / CONSENSUS: ", basename(paste0(obase, ext))), cex.sub = cex.sub, cex.lab = cex.lab)
graphics.off()

