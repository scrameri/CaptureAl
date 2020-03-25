#!/gdc_home4/scrameri/bin/Rscript

#########################
## PLOT CONTIG LENGTHS ##
#########################

## Usage
# run script without arguments to see the arguments it takes

## Author
# simon.crameri@env.ethz.ch, Apr 2019

## Load required library
suppressPackageStartupMessages(library(gplots))

## Get arguments
args = commandArgs(trailingOnly = TRUE)

## Check arguments
options(warning.length = 2000L)
if (! length(args) %in% c(1, 2, 3)) {
  stop("2 arguments taken (1 required): 
       REQUIRED
       1) <lentab|CHR>: path to table with contig ID (1st column) and contigs lengths (2nd to nth columns, 1 column for each locus). Expects no header and tab-delimiter ; 

       OPTIONAL
       2) <filter.loci|BOOLEAN>: whether to filter out loci where all individuals have empty sequences [default: TRUE] ;
       3) <meta|CHR>: path to metadata file mapping taxon ID (1st column) to metadata (2nd column). Expects header and tab-delimiter",
       call.=FALSE)
}

## Set arguments
lentab <- as.character(args[1])
filter.loci <- as.logical(as.character(args[2]))
if (is.na(filter.loci)) filter.loci <- TRUE
meta <- as.character(args[3])

## Set arguments (for debugging)
# lentab <- "multifasta.11.2605.locus.lengths"
# filter.loci <- TRUE
# meta <- "mapfile.orig.cpgroup.txt"

## Additional arguments
hclustmethod <- "ward.D2" # locus clustering method in heatmap
nulllength <- 20 # length of a null sequence, e.g. 20 for (--------------------)
ofile <- paste0(lentab, ".pdf") # name of output file
plot.width <- 15 # output plot width
plot.height <- 15 # output plot height
find.loci <- TRUE # if TRUE and filter.loci = TRUE, attempts to filter multifastas with only zero-contigs in all individuals
suffix1 <- ".locus"
suffix2 <- ".lengths"
rmdirname <- "loci_rm_allzerolength"

## Check arguments
stopifnot(file.exists(lentab))
if (!is.na(meta)) stopifnot(file.exists(meta))

## Define helperfunctions
# create a matrix for heatmaps
do.mat <- function(df, var, NA.replacement = 0, fac = "locus", ids = "ID") {
  mat <- matrix(as.numeric(df[,var]), ncol = nlevels(df[,fac]))
  rownames(mat) <- unique(df[,ids])
  mat[is.na(mat)] <- NA.replacement
  return(mat)
}

# extrafun to heatmap.2
myfun <- function() {
  plot(0,0,type="n",xlab="",ylab="",axes=F,bty="n")
  legend("topleft", legend = levels(dmerge$GROUP), 
         pch = 15, 
         pt.cex = 3, 
         cex = 1, 
         col = unique(RowSideColors), 
         ncol = ceiling(nlevels(dmerge$GROUP)/2),
         bty = "n", xpd = TRUE)
}

# wrapper around heatmap.2
do.heatmap <- function(mat, max.nbreaks = 200, labRow, Colv = TRUE, reverse.colors = FALSE, highquant = 1, lmat, lwid, lhei, hclustmethod, RowSideColors, extrafun, key.title) {
  
  stopifnot(is.matrix(mat), highquant > 0)
  
  # set heatmap legend color breaks
  nstates <- length(unique(as.numeric(mat)))
  matmax <- max(as.numeric(mat, na.rm = TRUE))
  if (nstates > max.nbreaks) {
    nbreaks <- max.nbreaks
  } else {
    nbreaks <- max(c(nstates, matmax))
  }
  nhigh <- floor(nbreaks/10)
  if (nhigh < 1) nhigh = 1
  
  if (highquant == 1) {
    breaks <- c(0, seq(1, matmax, length.out = nbreaks)-0.001)
  } else {
    if (highquant > 1) {
      if (highquant > matmax) highquant <- matmax
      breaks <- unique(c(seq(0, highquant, length.out = nbreaks - nhigh), 
                         seq(highquant, matmax, length.out = nhigh + 2))) 
    } else {
      breaks <- unique(c(seq(0, quantile(mat, probs = highquant, na.rm = TRUE), length.out = nbreaks - nhigh), 
                         seq(quantile(mat, probs = highquant, na.rm = TRUE), matmax, length.out = nhigh + 2)))
    }
  }
  
  # set heatmap legend color breaks (color does not change much after the highquant percentile)
  
  # pie(rep(1/nbreaks, nbreaks), col = rainbow(nbreaks, start = 0/6, end = 6/6))
  if (reverse.colors) colfun <- rev else colfun <- identity
  col <- c("white", do.call(what = colfun, args = list(rainbow(nbreaks-1, start = 0/6, end = 4/6)))) # 0/6 is red, 4/6 is blue
  if (length(breaks) < length(col) + 1) breaks <- c(0, seq(1, matmax, length.out = length(col))-0.001)
  
  # kex.xtickfun
  if (length(breaks) <= 25) {
    key.xtickfun <- function() {
      breaks <- parent.frame()$breaks
      return(list(
        at=parent.frame()$scale01(round(breaks, 1)),
        labels=c(as.character(round(breaks, 1)))
      ))
    }
  } else {
    key.xtickfun = NULL
  }
  
  # plot heatmap
  heatmap.2(x = mat,
            Rowv = FALSE, # no reordering of rows (individuals)
            Colv = Colv, # would reorder columns (loci) according to decreasing mean value (but separately for each cluster)
            labRow = labRow, labCol = NA,
            add.expr = box(),
            #hclustfun = function(x) hclust(x, method = hclustmethod), # here, x is going to be distfun(mat)
            #distfun = dist,                                           # 
            col = col,
            scale = "none",
            margins = c(3,0), # bottom margin will be 3
            trace = "none", 
            symkey = FALSE, 
            symbreaks = FALSE, 
            dendrogram = "none",
            density.info = "histogram",
            denscol = "black",
            breaks = breaks,
            key.title = key.title,
            key.xlab = "",
            key.ylab = "",
            key.xtickfun = key.xtickfun,
            keysize = 1,
            extrafun = extrafun,
            RowSideColors = RowSideColors,
            # RowSideColors=RowSideColors,
            key.par = list(mar = c(2,0,2,0)), # margins around key (bottom, left, top, right)
            lmat = lmat, lhei = lhei, lwid = lwid # to understand this, see show.layout(layout(lmat,lhei,lwid))
  )
}

# ggplot2 colors
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

######################################################################################################

## Read lenfile
dlen <- read.delim(lentab, header = FALSE)
names(dlen)[1] <- "ID"

## Read metadata
if (!is.na(meta)) {
  dmeta <- read.delim(meta)
  names(dmeta)[1] <- "ID"
  names(dmeta)[2] <- "GROUP"
}

## Merge it
if (!is.na(meta)) {
  dmerge <- merge(dlen, dmeta, by = "ID", all.x = TRUE, all.y = FALSE, sort = FALSE)
  dmerge$GROUP <- as.character(dmerge$GROUP)
  dmerge$GROUP[is.na(dmerge$GROUP)] <- "NA"
  dmerge$GROUP <- factor(dmerge$GROUP)
} else {
  dmerge <- dlen
  dmerge$GROUP <- factor("NA")
}
stopifnot(all(as.character(dmerge$ID) %in% as.character(dlen$ID)))
nind <- nrow(dlen)
nloci <- ncol(dlen)-1

## Order by group
dmerge <- dmerge[order(dmerge$GROUP, dmerge$ID),]
rownames(dlen) <- dlen$ID
dlen <- dlen[as.character(dmerge$ID),]

## Get matrix
mat <- as.matrix(dlen[,-1])
rownames(mat) <- dlen[,1]

## Handle NAs
mat[is.na(mat)] <- 0
mat[mat == nulllength] <- 0

## Plot
# cat("\ncreating plots...\n")

## Set heatmap labels and group colors
# labRow <- rownames(mat)
labRow <- paste(rownames(mat), dmerge$GROUP)

grcols <- gg_color_hue(nlevels(dmerge$GROUP))
RowSideColors <- grcols[dmerge$GROUP]

## Set heatmap layout (read the ?heatmap.2 documentation, or my own explanations below)
# heatmap.2 creates at least 4 plot components that must be placed at the desired position on the figure
# - if ColSideColors or RowSideColors are added, there are 6 plot components
# - if an extrafun is passed that plots something, there are 7 plot components
# - in our case, we plot 6 components (just RowSideColors [for individual groups] and 1 extrafun, no ColSideColors [for loci])
#
# Order of plot components by heatmap.2:
# - RowSideColors = 1 
# - heatmap matrix = 2
# - rowdendrogram = 3 (not plotted if dendrogram="none")
# - coldendrogram = 4 (not plotted if dendrogram="none")
# - heatmap color legend = 5
# - extrafun = 6, etc. (legend for RowSideColors)
#
# <lmat> codes the LOCATION of each plot component on the figure
# - here, the matrix [2] will end up in the bottom layout row, second from the left
# - while the legend [5] will end up in the second row, second from the left
lmat = rbind(c(0, 6, 0, 0),
             c(4, 5, 0, 0), 
             c(3, 2, 1, 0))

# <lhei> and <lwid> code the SIZE (height, width) of each plot component on the figure
# - these are the RELATIVE heights of each layout row (e.g. 3rd row 5/1.5 times as high as 2nd column)
lhei =            c(0.4,
                    1.5,
                    5)

# - these are the RELATIVE widths of each layout column (e.g. 2nd column 10 times as wide as 4th column)
lwid =      c(0.4, 10, 0.2, 1)

# show plot component layout and sizes
# layout.show(layout(lmat, heights = lhei, widths = lwid)) # matrix [2] gets the large space at the bottom, 2nd from left

pdf(file = ofile, width = plot.width, height = plot.height)
do.heatmap(mat, max.nbreaks = 200, labRow, Colv = TRUE, reverse.colors = FALSE, highquant = 0.975, lmat, lwid, lhei, hclustmethod, RowSideColors, myfun, paste0("Contig length (", nloci, " loci)"))
graphics.off()

## Look for all-zero-length loci or individuals
# loci
lzero <- which(apply(mat, 2, function(x) {all(x==0)}))
cat("found", length(lzero), "loci with zero length in all individuals!\n")

if (length(lzero) > 0) {
  writeLines(names(lzero), con = paste0(gsub(".lengths$", "", lentab), ".zeroloci"))
  
  if (find.loci) {
    ldir <- gsub(suffix1, "", lentab)
    adir <- gsub(paste0(suffix1, suffix2, "$"), "", lentab)
    if (dir.exists(ldir)) {
      zeroloci <- list.files(ldir)[lzero]
      tofilter <- file.path(adir, gsub(paste0(suffix2, "$"), ".fasta", zeroloci))
      print(basename(tofilter), max = 10, quote = FALSE)
      writeLines(tofilter, con = paste0(gsub(paste0(suffix2, "$"), "", lentab), ".zeroloci"))
    }
  }
  if (filter.loci) {
    if (dir.exists(adir)) {
      rmdir <- file.path(adir, rmdirname)
      if (!dir.exists(rmdir)) dir.create(path = rmdir)
      for (aln in tofilter) {
        if (file.exists(aln)) {
          system(command = paste("mv", aln, rmdir))
        }
      }
      cat("filtered", length(lzero), "loci with zero length in all individuals!\n")
    }
  }
}

# individuals
izero <- which(apply(mat, 1, function(x) {all(x==0)}))
cat("found", length(izero), "individuals with zero length in all loci!\n")

if (length(izero) > 0) {
  zeroinds <- names(izero)
  print(zeroinds, max = 10, quote = FALSE)
  writeLines(zeroinds, con = paste0(gsub(paste0(suffix2, "$"), "", lentab), ".zeroinds"))
}

