#!/cluster/apps/r/3.6.0_openblas/x86_64/bin/Rscript

## Usage: plot.coverage.stats.R <samples.txt> <mapQ> <OPT: mapfile> <OPT: min.len> <OPT: min.cov> <OPT: max.cov> <OPT: min.ratio> <OPT: max.perc>

## Load libraries
suppressPackageStartupMessages(library(ggplot2)) # ggplot
suppressPackageStartupMessages(library(tidyr))  # complete
suppressPackageStartupMessages(library(gplots)) # heatmap.2

# t1 <- Sys.time()

## Author: simon.crameri@env.ethz.ch, Feb 2020

## Get arguments
# print(Sys.time())
args <- commandArgs(trailingOnly=TRUE)

## Check arguments
if (!length(args) >= 2) {
  stop("8 arguments taken (2 needed):
       REQUIRED
       1) <sfile>     CHR sample file (path to sample directories with coverage stats in SAMPLE.bwa-mem.mapped.QMAPQ.sorted.bam.coverage.txt; 
       2) <mapQ>      CHR mapping quality threshold;
       
       OPTIONAL
       3) <meta>      CHR mapping file (w/ header, sample in 1st field, group in 2nd field, used to filter and sort coverage heatmap accordingly) [DEFAULT: NA];
       4) <min.len>   INT minimum length in .bam [DEFAULT: 500];
       5) <min.cov>   INT minimum coverage in .bam [DEFAULT: 10];
       6) <max.cov>   INT maximum coverage in .bam [DEFAULT: 1000];
       7) <min.ratio> NUM minimum alignment fraction [DEFAULT: 0.5];
       8) <max.perc>  NUM maximum percentage of non-conforming sample per group [DEFAULT: 0.1]",
       call.=FALSE)
}

## Set arguments
sfile <- as.character(args[1])
mapQ <- as.character(args[2])
meta <- as.character(args[3])

min.len <- as.numeric(args[4])
if (is.na(min.len)) min.len <- 500
min.cov <- as.numeric(args[5])
if (is.na(min.cov)) min.cov <- 10
max.cov <- as.numeric(args[6])
if (is.na(max.cov)) max.cov <- 1000
min.ratio <- as.numeric(args[7])
if (is.na(min.ratio)) min.ratio <- 0.5
max.perc <- as.numeric(args[8])
if (is.na(max.perc)) max.perc <- 0.1


## Set arguments (for debugging)
# sfile <- "samples.txt"
# mapQ <- "10"
# meta <- "mapfile.dalbergia.63.txt"
# 
# min.len <- 500
# min.cov <- 8
# max.cov <- 1000
# min.ratio <- 0.5
# max.perc <- 0.1


## Additional arguments
# paths
outpdf <- "coverage_stats.pdf" # name of output file
outtxt1 <- "coverage_stats.txt"
outtxt2 <- "loci_kept.txt"
outtxt3 <- "loci_rm.txt"
logtxt <- "coverage_stats.log"
pathinsubdir <- "SAMPLE.bwa-mem.mapped.QMAPQ.sorted.bam.coverage.txt"


# what to plot
pvar1 <- "length.in.bam"
pvar2 <- "length.in.refseq"
pvar3 <- "avg.coverage.bam" # summarized
pvar4 <- "avg.coverage.refseq"


# plotting
hclustmethod <- "ward.D2" # locus clustering method in heatmap
sortmat = TRUE # if TRUE, will sort loci based on <hclustmethod> in heatmap plots
plot.width <- 15 # output plot width
plot.height <- 15 # output plot height


## Check arguments
stopifnot(file.exists(sfile))
if (!is.na(meta)) stopifnot(file.exists(meta))


## Define helperfunctions
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# create a matrix for heatmaps
do.mat <- function(df, var, NA.replacement = 0, fac = "locus", ids = "ID") {
  mat <- matrix(as.numeric(df[,var]), ncol = nlevels(factor(df[,fac])))
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

###################################################################################

## Read sample file
samples <- as.character(read.delim(sfile, check.names = FALSE, header = FALSE)$V1)


# check subdirectories for each sample
scheck <- sapply(samples, FUN = dir.exists)
if (any(scheck == FALSE)) {
  stop(paste0("the following directories could not be found:\n",
              paste(names(scheck[scheck == FALSE]), collapse = ", ")))
}


## Read coverage data
cat("reading coverage data...\n")
covdata <- data.frame(array(NA, dim = c(0, 6), dimnames = list(NULL, c("ID", "region","length.in.bam","length.in.refseq","avg.coverage.bam","avg.coverage.refseq"))))
for (i in samples) {
  cat(which(samples == i), " / ", length(samples), "\r")
  p <- file.path(i, paste0(gsub("MAPQ", mapQ, gsub("SAMPLE", i, pathinsubdir))))
  if (!file.exists(p)) warning("file ", p, " not found,\ndid you specify the correct path (line 47 in script)?")
  dd <- read.delim(p, header = T, stringsAsFactors = FALSE)
  dd <- data.frame(ID = rep(i, nrow(dd)), dd, stringsAsFactors = FALSE)
  covdata <- rbind(covdata, dd)
}


# fill up regions with no coverage
covdata$ID <- factor(covdata$ID)
covdata$region <- factor(covdata$region)
dd <- data.frame(tidyr::complete(covdata, ID, region))


## Read metadata
if (!is.na(meta)) {
  dmeta <- read.delim(meta)
  names(dmeta)[1] <- "ID"
  names(dmeta)[2] <- "GROUP"
}


## Merge data
cat("\n\nmerging data...\n")
if (!is.na(meta)) {
  dmerge <- merge(dd, dmeta, by = "ID", all.x = TRUE, all.y = FALSE, sort = FALSE)
  dmerge$GROUP <- as.character(dmerge$GROUP)
} else {
  dmerge <- dd
  dmerge$GROUP <- "NA"
}

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


## Get number of individuals and groups, number of target loci
dmerge[,"GROUP"][is.na(dmerge[,"GROUP"])] <- "NA"
dmerge[,"GROUP"] <- factor(dmerge[,"GROUP"])
ntaxa <- nlevels(dmerge$ID) # all taxa
ngroups <- length(unique(dmerge$GROUP)) # all mapped groups
nloci <- nlevels(dmerge$region) # all target loci
dmerge.notconsidered <- subset(dmerge[!duplicated(dmerge$ID),c("ID","GROUP")], GROUP == "NA")
ntaxa.notconsidered <- nrow(dmerge.notconsidered)


## Print grouping
grtab <- table(dmerge[!duplicated(dmerge$ID),"GROUP"])
cat(paste0("\n", ngroups, " taxon grouping(s) found:\n"))
print(grtab)
if (!is.na(meta) & ntaxa.notconsidered > 0) {
  cat(paste0("\nWARNING: ", ntaxa.notconsidered, " taxa are either missing or <NA> in ", meta, " and will not be considered for locus filtering:\n"))
  print(as.character(dmerge.notconsidered$ID))
}


## Sort according to groups (needed to create matrix)
dmerge <- dmerge[order(dmerge["region"], dmerge[,"GROUP"], dmerge[,"ID"]),]


## Add number of individuals per group
dmerge$LABEL <- as.character(dmerge$GROUP)
for (group in unique(dmerge$GROUP)) {
  dmerge$LABEL[dmerge$LABEL == group] <- paste0(names(grtab[group]), " (n = ", grtab[group], ")")
}
dmerge$LABEL <- as.factor(dmerge$LABEL)


## Add ratios
covmax <- quantile(dmerge[,pvar3], 0.975, na.rm = TRUE)
dmerge$avg.coverage.bam_log <- log(dmerge$avg.coverage.bam)
dmerge$length_ratio <- dmerge$length.in.bam / dmerge$length.in.refseq
dmerge$coverage_per_length <- dmerge$avg.coverage.bam / dmerge$length.in.bam

## Get stats
if (all(dmerge[,"GROUP"] == "NA")) {
  groups <- "NA"
} else {
  groups <- levels(dmerge[,"GROUP"])
  groups <- groups[!groups %in% "NA"]
}

lowlen <- lowratio <- lowcov <- highcov <- lowmap <- character()
for (group in groups) {
  dsub <- subset(dmerge, GROUP == group)
  
  s1 <- aggregate(dsub[, pvar1], by = list(dsub[,"region"]), FUN = function(x) {
    length(x[x<min.len & !is.na(x)])/length(x[!is.na(x)])})
  torm1 <- as.character(subset(s1, x > max.perc)[,1])
  lowlen <- sort(unique(c(lowlen, torm1)))

  s2 <- aggregate(dsub[, "length_ratio"], by = list(dsub[,"region"]), FUN = function(x) {
    length(x[x<min.ratio & !is.na(x)])/length(x[!is.na(x)])})
  torm2 <- as.character(subset(s2, x > max.perc)[,1])
  lowratio <- sort(unique(c(lowratio, torm2)))
  
  
  s3 <- aggregate(dsub[, pvar3], by = list(dsub[,"region"]), FUN = function(x) {
    length(x[x<min.cov & !is.na(x)])/length(x[!is.na(x)])})
  torm3 <- as.character(subset(s3, x > max.perc)[,1])
  lowcov <- sort(unique(c(lowcov, torm3)))
  
  s4 <- aggregate(dsub[, pvar3], by = list(dsub[,"region"]), FUN = function(x) {
    length(x[x>max.cov & !is.na(x)])/length(x[!is.na(x)])})
  torm4 <- as.character(subset(s4, x > max.perc)[,1])
  highcov <- sort(unique(c(highcov, torm4)))
  
  s5 <- aggregate(dsub[, pvar3], by = list(dsub[,"region"]), FUN = function(x) {
    length(x[is.na(x)])/length(x)})
  torm5 <- as.character(subset(s5, x > max.perc)[,1])
  lowmap <- sort(unique(c(lowmap, torm5)))

}
torm <- sort(unique(c(lowlen, lowratio, lowcov, highcov, lowmap)))
tokeep <- unique(dmerge$region)[!unique(dmerge$region) %in% torm]
lfailed <- length(torm)
lpassed <- length(tokeep)


## Filtering summary
cat("\n")
paste1 <- paste0("Loci with length <", min.len, " in >", max.perc, " of individuals in each group:")
paste2 <- paste0("Loci with coverage <", min.cov, " in >", max.perc, " of individuals in each group:")
paste3 <- paste0("Loci with coverage >", max.cov, " in >", max.perc, " of individuals in each group:")
paste4 <- paste0("Loci with alignment fraction <", min.ratio, " in >", max.perc, " of individuals in each group:")
paste5 <- paste0("Loci without mapped reads in >", max.perc, " of individuals in each group:")

fillto <- max(c(nchar(paste1), nchar(paste2), nchar(paste3), nchar(paste4), nchar(paste5))) + 2
pasteA <- paste0("Total number of target loci:", paste(rep(" ", fillto-28), collapse = ""), nloci, " (100%)\n")

paste1 <- paste0(paste1, paste(rep(" ", fillto-nchar(paste1)), collapse = ""), length(lowlen),  " (", round(100*length(lowlen)/nloci,2),"%)\n")
paste2 <- paste0(paste2, paste(rep(" ", fillto-nchar(paste2)), collapse = ""), length(lowcov),  " (", round(100*length(lowcov)/nloci,2),"%)\n")
paste3 <- paste0(paste3, paste(rep(" ", fillto-nchar(paste3)), collapse = ""), length(highcov)," (", round(100*length(highcov)/nloci,2),"%)\n")
paste4 <- paste0(paste4, paste(rep(" ", fillto-nchar(paste4)), collapse = ""), length(lowratio),  " (", round(100*length(lowratio)/nloci,2),"%)\n")
paste5 <- paste0(paste5, paste(rep(" ", fillto-nchar(paste5)), collapse = ""), length(lowmap),  " (", round(100*length(lowmap)/nloci,2),"%)\n")
paste6 <- paste0("Cumulative loci to remove:", paste(rep(" ", fillto-26), collapse = ""), lfailed, " (", round(100*lfailed/nloci,2),"%)\n")
paste7 <- paste0("Cumulative loci to keep:", paste(rep(" ", fillto-24), collapse = ""), lpassed, " (", round(100*lpassed/nloci,2),"%)\n")

cat(pasteA)
cat(paste1)
cat(paste2)
cat(paste3)
cat(paste4)
cat(paste5)
cat(paste6)
cat(paste7)


## Filter loci
dmerge.rm <- subset(dmerge, region %in% torm)
dmerge.rm$region <- droplevels(dmerge.rm$region)
dmerge.keep <- subset(dmerge, ! region %in% torm)
dmerge.keep$region <- droplevels(dmerge.keep$region)


## Get colors
grcols <- gg_color_hue(nlevels(dmerge[,"GROUP"]))
RowSideColors <- grcols[dmerge[!duplicated(dmerge$ID),][,"GROUP"]]

## Create matrices
dd1 <- do.mat(df = dmerge, var = pvar1, NA.replacement = 0, fac = "region", ids = "ID")
# dd2 <- do.mat(df = dmerge, var = pvar2, NA.replacement = 0, fac = "region", ids = "ID")
dd3 <- do.mat(df = dmerge, var = pvar3, NA.replacement = 0, fac = "region", ids = "ID")
# dd4 <- do.mat(df = dmerge, var = pvar4, NA.replacement = 0, fac = "region", ids = "ID")
labRow <- rownames(dd1)

# to rm
dd1.rm <- do.mat(df = dmerge.rm, var = pvar1, NA.replacement = 0, fac = "region", ids = "ID")
# dd2.rm <- do.mat(df = dmerge.rm, var = pvar2, NA.replacement = 0, fac = "region", ids = "ID")
dd3.rm <- do.mat(df = dmerge.rm, var = pvar3, NA.replacement = 0, fac = "region", ids = "ID")
# dd4.rm <- do.mat(df = dmerge.rm, var = pvar4, NA.replacement = 0, fac = "region", ids = "ID")

# to keep
dd1.keep <- do.mat(df = dmerge.keep, var = pvar1, NA.replacement = 0, fac = "region", ids = "ID")
# dd2.keep <- do.mat(df = dmerge.keep, var = pvar2, NA.replacement = 0, fac = "region", ids = "ID")
dd3.keep <- do.mat(df = dmerge.keep, var = pvar3, NA.replacement = 0, fac = "region", ids = "ID")
# dd4.keep <- do.mat(df = dmerge.keep, var = pvar4, NA.replacement = 0, fac = "region", ids = "ID")


## Write LOG
suffix <- paste(paste0("-q", mapQ),min.len,min.cov,max.cov,min.ratio,paste0(max.perc,".txt"), sep = "-")
outpdf <- gsub(".pdf$", gsub(".txt$", ".pdf", suffix), outpdf)
outtxt1 <- outtxt1
outtxt2 <- gsub(".txt$", suffix, outtxt2)
outtxt3 <- gsub(".txt$", suffix, outtxt3)
logtxt <- gsub(".log", gsub(".txt$", ".log", suffix), logtxt)

dlog <- c(
  paste0("=== Filter loci based on visualization of coverage statistics ==="),
  paste0(""),
  paste0("Input data:"),
  paste0("-----------"),
  paste0("Sample file:     ", sfile),
  paste0("Group file:      ", meta),
  paste0("Mapping quality: ", mapQ),
  paste0(""),
  paste0("Thresholds for locus filtering:"),
  paste0("-------------------------------"),
  paste0("Min. length:                                         ", min.len),
  paste0("Min. average coverage:                               ", min.cov),
  paste0("Min. alignment fraction:                             ", min.ratio),
  paste0("Max. percentage of non-conforming samples per group: ", max.perc),
  paste0(""),
  paste0("Other (internal) parameters:"),
  paste0("----------------------------"),
  paste0("Path to coverage stats: ", pathinsubdir),
  paste0("Hierarchical clustering of loci in heatmaps: ", as.character(sortmat)),
  paste0("Hierarchical clustering method: ", hclustmethod),
  paste0("PDF height (inches): ", plot.height),
  paste0("PDF width (inches): ", plot.width),
  paste0(""),
  paste0("Locus filtering results:"),
  paste0("------------------------"),
  # paste0("Total number of target loci: ", nloci),
  substring(pasteA, 1, nchar(pasteA)-1),
  substring(paste1, 1, nchar(paste1)-1),
  substring(paste2, 1, nchar(paste2)-1),
  substring(paste3, 1, nchar(paste3)-1),
  substring(paste4, 1, nchar(paste4)-1),
  substring(paste5, 1, nchar(paste5)-1),
  substring(paste6, 1, nchar(paste6)-1),
  substring(paste7, 1, nchar(paste7)-1),
  paste(""),
  paste0("Output files:"),
  paste0("-------------"),
  paste0("Log file:                               ", logtxt),
  paste0("Visualization of coverage stats:        ", outpdf),
  paste0("Coverage stats:                         ", outtxt1),
  paste0("Kept loci:                              ", outtxt2),
  paste0("Removed loci:                           ", outtxt3),
  paste0("")
)


## Write TXT
cat("\nwriting results .txt and .log...\n")
write.table(dmerge, file = outtxt1, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
write.table(tokeep, file = outtxt2, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
write.table(torm,   file = outtxt3, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
writeLines(text = dlog, con = logtxt)


## Plot stats
p1 <- ggplot(dmerge, aes_string(x = "LABEL", y = pvar1, fill = "GROUP")) +
  geom_boxplot(alpha = 0.5, na.rm = TRUE) +
  geom_violin(alpha = 0.5, na.rm = TRUE) +
  labs(x = "", y = pvar1) +
  geom_hline(yintercept = min.len, col = "tomato") +
  scale_fill_discrete(guide = FALSE) +
  coord_flip() +
  theme_bw()

p2 <- ggplot(dmerge, aes_string(x = "LABEL", y = pvar3, fill = "GROUP")) +
  geom_boxplot(alpha = 0.5, na.rm = TRUE) +
  geom_violin(alpha = 0.5, na.rm = TRUE) +
  labs(x = "", y = paste0(pvar3, " (>97.5% quantile removed)")) +
  ylim(c(0, covmax)) +
  geom_hline(yintercept = min.cov, col = "tomato") +
  geom_hline(yintercept = covmax, col = "black") +
  scale_fill_discrete(guide = FALSE) +
  coord_flip() +
  theme_bw()

p3 <- ggplot(dmerge, aes_string(x = "LABEL", y = paste0(pvar3, "_log"), fill = "GROUP")) +
  geom_boxplot(alpha = 0.5, na.rm = TRUE) +
  geom_violin(alpha = 0.5, na.rm = TRUE) +
  labs(x = "", y = paste0("log(", pvar3, ")")) +
  geom_hline(yintercept = log(min.cov), col = "tomato") +
  scale_fill_discrete(guide = FALSE) +
  coord_flip() +
  theme_bw()

p4 <- ggplot(dmerge, aes_string(x = "LABEL", y = "length_ratio", fill = "GROUP")) +
  geom_boxplot(alpha = 0.5, na.rm = TRUE) +
  geom_violin(alpha = 0.5, na.rm = TRUE) +
  labs(x = "", y = "alignment fraction") +
  geom_hline(yintercept = min.ratio, col = "tomato") +  
  scale_fill_discrete(guide = FALSE) +
  coord_flip() +
  theme_bw()

p5 <- ggplot(dmerge, aes_string(x = "length.in.bam", y = "avg.coverage.bam", col = "GROUP")) +
  geom_point(alpha = 0.1, na.rm = TRUE) +
  geom_density_2d(aes_string(colour = "GROUP"), alpha = 1, na.rm = TRUE) +
  labs(x = "length in .bam", y = "average coverage in .bam") +
  geom_hline(yintercept = max.cov, col = "tomato") +  
  scale_fill_discrete(guide = FALSE) +
  coord_flip() +
  scale_colour_discrete(guide = FALSE) +
  facet_wrap(~GROUP) +
  theme_bw()


## Plot heatmaps
cat("\nplotting heatmaps...")
pdf(outpdf, width = plot.width, height = plot.height)

print(p1)
print(p2)
print(p3)
print(p4)
print(p5)

do.heatmap(dd1, max.nbreaks = 200, labRow, Colv = sortmat, reverse.colors = TRUE, highquant = 0.975, lmat, lwid, lhei, hclustmethod, RowSideColors, myfun, paste0(pvar1, ", ALL, nloc = ", ncol(dd1)))
# do.heatmap(dd1.rm, max.nbreaks = 200, labRow, Colv = sortmat, reverse.colors = TRUE, highquant = 0.975, lmat, lwid, lhei, hclustmethod, RowSideColors, myfun, paste0(pvar1, ", REMOVED, nloc = ", ncol(dd1.rm)))
do.heatmap(dd1.keep, max.nbreaks = 200, labRow, Colv = sortmat, reverse.colors = TRUE, highquant = 0.975, lmat, lwid, lhei, hclustmethod, RowSideColors, myfun, paste0(pvar1, ", FILTERED, nloc = ", ncol(dd1.keep)))

do.heatmap(dd3, max.nbreaks = 200, labRow, Colv = sortmat, reverse.colors = TRUE, highquant = 0.975, lmat, lwid, lhei, hclustmethod, RowSideColors, myfun, paste0(pvar3, ", ALL, nloc = ", ncol(dd3)))
# do.heatmap(dd3.rm, max.nbreaks = 200, labRow, Colv = sortmat, reverse.colors = TRUE, highquant = 0.975, lmat, lwid, lhei, hclustmethod, RowSideColors, myfun, paste0(pvar3, ", REMOVED, nloc = ", ncol(dd3.rm)))
do.heatmap(dd3.keep, max.nbreaks = 200, labRow, Colv = sortmat, reverse.colors = TRUE, highquant = 0.975, lmat, lwid, lhei, hclustmethod, RowSideColors, myfun, paste0(pvar3, ", FILTERED, nloc = ", ncol(dd3.keep)))

graphics.off()

cat("done\n")


