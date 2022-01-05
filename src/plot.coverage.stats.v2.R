#!/cluster/apps/r/3.6.0_openblas/x86_64/bin/Rscript

## Usage: plot.coverage.stats.R <samples.txt> <mapQ> <OPT: min.len> <OPT: min.cov> <OPT: max.cov> <OPT: min.ratio> <OPT: max.frac> <OPT: mapfile.txt>

## Load libraries
suppressPackageStartupMessages(library(ggplot2)) # ggplot
suppressPackageStartupMessages(library(tidyr))  # complete
suppressPackageStartupMessages(library(dplyr))  # sym

## Author: simon.crameri@env.ethz.ch, Mar 2020

## Get arguments
args <- commandArgs(trailingOnly=TRUE)

## Check arguments
if (!length(args) >= 2) {
  stop("8 arguments taken (2 needed):
       REQUIRED
       1) <sfile>     CHR sample file (path to sample directories with coverage stats in SAMPLE.bwa-mem.mapped.QMAPQ.sorted.bam.coverage.txt; 
       2) <mapQ>      CHR mapping quality threshold;
       
       OPTIONAL
       3) <min.len>   INT minimum length in .bam [DEFAULT: 500];
       4) <min.cov>   INT minimum coverage in .bam [DEFAULT: 10];
       5) <max.cov>   INT maximum coverage in .bam [DEFAULT: 1000];
       6) <min.ratio> NUM minimum alignment fraction [DEFAULT: 0.5];
       7) <max.frac>  NUM maximum fraction of non-conforming samples per group [DEFAULT: 0.1];
       8) <meta>      CHR mapping file (w/ header, sample in 1st field, group in 2nd field, used to filter regions and sort heatmaps accordingly) [DEFAULT: NA]",
       call.=FALSE)
}

## Set arguments
sfile <- as.character(args[1])
mapQ <- as.character(args[2])

min.len <- as.numeric(args[3])
if (is.na(min.len)) min.len <- 500
min.cov <- as.numeric(args[4])
if (is.na(min.cov)) min.cov <- 10
max.cov <- as.numeric(args[5])
if (is.na(max.cov)) max.cov <- 1000
min.ratio <- as.numeric(args[6])
if (is.na(min.ratio)) min.ratio <- 0.5
max.frac <- as.numeric(args[7])
if (is.na(max.frac)) max.frac <- 0.1

meta <- as.character(args[8])

t1 <- Sys.time()
paste0("Starting time: ", t1)


## Set arguments (for debugging)
# sfile <- "samples.txt"
# mapQ <- "10"
# 
# min.len <- 1
# min.cov <- 1
# max.cov <- 10000
# min.ratio <- 0
# max.frac <- 0.1
#
# meta <- "mapfile.txt"


## Additional arguments
# input paths (path to coverage stats file)
# will look for <covstatfile> in each sample subdirectory. 
# <SAMPLE> and <MAPQ> can be part of the path and will be handled using regular expressions.
covstatfile <- "SAMPLE.bwa-mem.mapped.QMAPQ.sorted.bam.coverage.txt" 

# ouput paths
suffix <- paste(paste0("-q", mapQ),min.len,min.cov,max.cov,min.ratio,max.frac, sep = "-")
outpdf <- paste0("coverage_stats", suffix, ".pdf")  # output .pdf (visualization of coverage stats)
outtxt1 <- paste0("coverage_stats", ".txt")         # output .txt (merged coverage stats)
outtxt2 <- paste0("regions_kept", suffix, ".txt")   # passed regions 
outtxt3 <- paste0("regions_rm", suffix, ".txt")     # failed regions
logtxt <- paste0("coverage_stats", suffix, ".log")  # output .log
 
# variables in <covstatfile>
regvar <- "region"              # variable name in <covstatfile> denoting <region id>
lvar1 <- "length.in.refseq"     # variable name in <covstatfile> denoting <length of reference sequence>
pvar1 <- "length.in.bam"        # variable name in <covstatfile> denoting <length of mapped region>
cvar2 <- "avg.coverage.refseq"  # variable name in <covstatfile> denoting <average coverage of reference sequence>
pvar2 <- "avg.coverage.bam"     # variable name in <covstatfile> denoting <average coverage of mapped region>

# variables created in merged data.frame
idvar <- "ID"                   # variable name denoting <sample ID>
grpvar1 <- "GROUP"              # variable name denoting <group ID>
grpvar2 <- "LABEL"              # variable name denoting <group ID (n = N)> 
pvar3 <- "alignment_fraction"   # computed as pvar1 / lvar1
navar <- "NA"                   # string used to denote membership to a not considered group

# plotting
qlow <- 0.1                     # lower quantile (shown on violin plots and used to cap heatmaps of FAILED loci)
qhigh <- 0.9                    # upper quantile (shown on violin plots and used to cap heatmaps of ALL and PASSED loci)
hclustmethod <- "ward.D2"       # region clustering method in heatmap
sortx = TRUE                    # if TRUE, will sort regions based on <hclustmethod> in heatmap plots
ycols = c("#A6CEE3","#99CD91","#B89B74","#F06C45","#ED8F47","#825D99","#B15928") # group colors, will be recycled
low = "#A50026"                 # color used for the low end of the heatmap gradient
mid = "#F0E442"                 # color used for the mid point of the heatmap gradient
high = "#081D58"                # color used for the high end of the heatmap gradient
plot.width <- 15                # output plot width
plot.height <- 15               # output plot height

##########################################################################################

## Check arguments
stopifnot(file.exists(sfile), min.len >= 1, min.cov >= 1, is.numeric(max.cov), min.ratio >= 0, min.ratio <= 1, max.frac >= 0, max.frac <= 1)
if (!is.na(meta)) stopifnot(file.exists(meta))


## Define helperfunctions
# create a ggplot2 heatmap (y sorted for groups, z sorted for clusters)
do.heatmap <- function(dat, xfac, yfac, znum, ygr, ycols = NULL, limit = 100, title = znum,
                       sortx = TRUE, hclustmethod = "ward.D2", low = "#A50026", mid = "#F0E442", high = "#081D58") {
  
  # Usage
  # do.heatmap(dat, xfac, yfac, znum, ygr, ycols = NULL, limit = 100, title = znum,
  #            hclustmethod = "ward.D2", low = "#A50026", mid = "#F0E442", high = "#081D58")
  #
  # Arguments
  # dat	         data.frame with factor variables for x and y axes of the heatmap, 
  #              a numeric variable for the z axis of the heatmap, as well as a y grouping factor
  # xfac         character denoting variable name of x factor
  # yfac         character denoting variable name of y factor
  # znum         character denoting variable name of z numeric
  # ygr          character denoting variable name of y grouping factor
  # ycols        color vector for y grouping. Will be recycled.
  # limit        numeric denoting the limit the heatmap color scale to c(0, limit). 
  #              Values higher than the limit will be reset to the limit.
  # title        character used to inform the heatmap title
  # sortx        boolean denoting whether or not to sort x columns in the heatmap using hierarchical
  #              clustering, or character vector of region names used for sorting
  # chlustmethod method used for the hierarchical clustering. See ?hclust()
  # low          color used for the low end of the heatmap gradient
  # mid          color used for the mid point of the heatmap gradient
  # high         color used for the high end of the heatmap gradient
  
  # Author simon.crameri@env.ethz.ch, Mar 2020
  
  # helperfunctions
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  do.mat <- function(df, var, NA.replacement = 0, fac, ids) {
    mat <- matrix(as.numeric(df[,var]), ncol = nlevels(factor(df[,fac])))
    rownames(mat) <- unique(df[,ids])
    mat[is.na(mat)] <- NA.replacement
    return(mat)
  }
  
  # check input
  nygr <- nlevels(dat[,ygr])
  if (is.null(ycols)) {
    ycols <- gg_color_hue(nygr)
  } else {
    ycols <- rep(ycols, ceiling(nygr/length(ycols)))[1:nygr]
  }
  stopifnot(is.data.frame(dat), is.factor(dat[,xfac]), is.factor(dat[,yfac]), is.numeric(dat[,znum]),
            is.factor(dat[,ygr]), length(ycols) == nygr, is.numeric(limit),
            any(is.logical(sortx), is.character(sortx)))
  
  if (nrow(dat) > 0) {
    # sort for xfac, ygr and yfac
    dat <- dat[order(dat[,xfac], dat[,ygr], dat[,yfac]),]
    
    # cluster xfac
    if (is.logical(sortx)) {
      if (sortx) {
        dd.mat <- do.mat(df = dat, var = znum, NA.replacement = 0, fac = xfac, ids = yfac)
        dd.dendro <- as.dendrogram(hclust(d = dist(x = t(dd.mat)), method = hclustmethod))
        dd.order <- order.dendrogram(dd.dendro)
        dat[,xfac] <- factor(dat[,xfac], levels = levels(dat[,xfac])[dd.order], ordered = T)
      }
    } else {
      if (any(!sortx %in% levels(dat[,xfac]))) warning("some regions in <sortx> not found <dat>!")
      if (any(!levels(dat[,xfac]) %in% sortx)) warning("some regions in <dat> not found in <sortx>!")
      dat[,xfac] <- factor(dat[,xfac], levels = sortx[sortx %in% levels(dat[,xfac])], ordered = T)
    }
    
    # limit znum
    zmax <- round(max(dat[,znum], na.rm = TRUE), 2)
    dat[,znum][which(dat[,znum] > limit)] <- limit
    dat[,yfac] <- factor(as.character(dat[,yfac]), levels = rev(unique(dat[,yfac])))
    
    # plot
    dummy <- data.frame(x = (1:nygr), cols = levels(dat[,ygr]))
    legrows <- function(x) {ifelse(x>21, 4, ifelse(x>14, 3, ifelse(x>7, 2, 1)))}
    ygrfac <- rev(dat[!duplicated(dat[,yfac]),ygr])
    ybreaks <- sapply(1:nygr, FUN = function(x) {max(which(ygrfac == rev(levels(ygrfac))[x]))})
    p <- ggplot(dat, aes_string(xfac, yfac, fill = znum)) + 
      geom_tile() +
      geom_hline(yintercept = c(0, ybreaks)+0.5, size = 0.5) +
      geom_point(data = dummy, aes(x = x, y = x, colour = cols), alpha = 0, inherit.aes = FALSE) +
      scale_fill_gradient2(low = low, mid = mid, high = high, na.value = "white", name = "",
                           midpoint = median(dat[,znum], na.rm = TRUE), limits = c(0, limit)) +
      scale_color_manual(name = "", values = ycols) +
      guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1, name = ""), nrow = legrows(nygr), byrow = TRUE)) +
      labs(x = NULL, y = NULL) +
      theme_bw() +
      theme(axis.text.y = element_text(size = 4, colour = ycols[ygrfac]),
            axis.title.x = element_blank(), axis.text.x = element_blank(),
            axis.ticks.x = element_blank(), legend.position = "bottom",
            legend.margin = margin()) +
      ggtitle(paste0("Heatmap of ", title, " (nind = ", nlevels(dat[,yfac]), 
                     "; nloc= ", nlevels(dat[,xfac]), "), capped at ", limit, " (maximum = ", zmax, ")"))
    return(p)
  } else {
    warning("No regions found!\n")
  }
}

##########################################################################################

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
covdata <- data.frame(array(NA, dim = c(0, 6), dimnames = list(NULL, c(idvar,regvar,pvar1,lvar1,pvar2,cvar2))))
for (i in samples) {
  cat(which(samples == i), " / ", length(samples), "\r")
  p <- file.path(i, paste0(gsub("MAPQ", mapQ, gsub("SAMPLE", i, covstatfile))))
  if (!file.exists(p)) warning("file <", p, "> not found,\ndid you specify the correct path (line 69 in script)?")
  dd <- read.delim(p, header = TRUE, stringsAsFactors = FALSE)
  dd <- data.frame(ID = rep(i, nrow(dd)), dd, stringsAsFactors = FALSE)
  if (!regvar %in% names(dd)) {stop("variable <", regvar, "> not found in <", p, ">, please check line 79 in script")}
  if (!pvar1 %in% names(dd))  {stop("variable <", pvar1,  "> not found in <", p, ">, please check line 80 in script")}
  if (!lvar1 %in% names(dd))  {stop("variable <", lvar1,  "> not found in <", p, ">, please check line 81 in script")}
  if (!pvar2 %in% names(dd))  {stop("variable <", pvar2,  "> not found in <", p, ">, please check line 82 in script")}
  if (!cvar2 %in% names(dd))  {stop("variable <", cvar2,  "> not found in <", p, ">, please check line 83 in script")}
  covdata <- rbind(covdata, dd)
}
names(covdata)[1] <- c(idvar)

# fill up regions with no coverage
covdata[,idvar] <- factor(covdata[,idvar])
covdata[,regvar] <- factor(covdata[,regvar])
dd <- data.frame(tidyr::complete(covdata, !! sym(idvar), !! sym(regvar)))


## Read metadata
if (!is.na(meta)) {
  dmeta <- read.delim(meta)
  names(dmeta)[1] <- idvar
  names(dmeta)[2] <- grpvar1
}


## Merge data
cat("\n\nmerging data...\n")
if (!is.na(meta)) {
  dmerge <- merge(dd, dmeta, by = idvar, all.x = TRUE, all.y = FALSE, sort = FALSE)
  dmerge[,grpvar1] <- as.character(dmerge[,grpvar1])
} else {
  dmerge <- dd
  dmerge[,grpvar1] <- navar
}


## Get number of individuals and groups, number of regions
dmerge[,grpvar1][is.na(dmerge[,grpvar1])] <- navar
dmerge[,grpvar1] <- factor(dmerge[,grpvar1])
ntaxa <- nlevels(dmerge[,idvar]) # all taxa
ngroups <- length(unique(dmerge[,grpvar1])) # all mapped groups
nreg <- nlevels(dmerge[,regvar]) # all regions
dmerge.notconsidered <- dmerge[!duplicated(dmerge[,idvar]) & dmerge[,grpvar1] == navar,c(idvar,grpvar1)]
ntaxa.notconsidered <- nrow(dmerge.notconsidered)


## Print grouping
grtab <- table(dmerge[!duplicated(dmerge[,idvar]),grpvar1])
cat(paste0("\n", ngroups, " taxon grouping(s) found:\n"))
print(grtab)
if (!is.na(meta) & ntaxa.notconsidered > 0) {
  cat(paste0("\nWARNING: ", ntaxa.notconsidered, " taxa are either missing or <NA> in ", meta, " and will not be considered for region filtering:\n"))
  print(as.character(dmerge.notconsidered[,idvar]))
}


## Sort according to groups (needed to create matrix)
dmerge <- dmerge[order(dmerge[,regvar], dmerge[,grpvar1], dmerge[,idvar]),]


## Add number of individuals per group
dmerge[,grpvar2] <- as.character(dmerge[,grpvar1])
for (group in unique(dmerge[,grpvar1])) {
  dmerge[,grpvar2][dmerge[,grpvar2] == group] <- paste0(names(grtab[group]), " (n = ", grtab[group], ")")
}
dmerge[,grpvar2] <- as.factor(dmerge[,grpvar2])


## Add alignmnet fraction
dmerge[,pvar3] <- dmerge[,pvar1] / dmerge[,lvar1]


## Filter regions
if (all(dmerge[,grpvar1] == navar)) {
  groups <- navar
} else {
  groups <- levels(dmerge[,grpvar1])
  groups <- groups[!groups %in% navar]
}
lowlen <- lowratio <- lowcov <- highcov <- lowmap <- character()

for (group in groups) {
  dsub <- dmerge[dmerge[,grpvar1] == group,]
  
  s1 <- aggregate(dsub[,pvar1], by = list(dsub[,regvar]), FUN = function(x) {
    length(x[x<min.len & !is.na(x)])/length(x[!is.na(x)])})
  torm1 <- as.character(subset(s1, x > max.frac)[,1])
  lowlen <- sort(unique(c(lowlen, torm1)))
  
  s2 <- aggregate(dsub[,pvar3], by = list(dsub[,regvar]), FUN = function(x) {
    length(x[x<min.ratio & !is.na(x)])/length(x[!is.na(x)])})
  torm2 <- as.character(subset(s2, x > max.frac)[,1])
  lowratio <- sort(unique(c(lowratio, torm2)))
  
  s3 <- aggregate(dsub[,pvar2], by = list(dsub[,regvar]), FUN = function(x) {
    length(x[x<min.cov & !is.na(x)])/length(x[!is.na(x)])})
  torm3 <- as.character(subset(s3, x > max.frac)[,1])
  lowcov <- sort(unique(c(lowcov, torm3)))
  
  s4 <- aggregate(dsub[,pvar2], by = list(dsub[,regvar]), FUN = function(x) {
    length(x[x>max.cov & !is.na(x)])/length(x[!is.na(x)])})
  torm4 <- as.character(subset(s4, x > max.frac)[,1])
  highcov <- sort(unique(c(highcov, torm4)))
  
  s5 <- aggregate(dsub[,pvar2], by = list(dsub[,regvar]), FUN = function(x) {
    length(x[is.na(x)])/length(x)})
  torm5 <- as.character(subset(s5, x > max.frac)[,1])
  lowmap <- sort(unique(c(lowmap, torm5)))
  
}
torm <- sort(unique(c(lowlen, lowratio, lowcov, highcov, lowmap)))
tokeep <- unique(dmerge[,regvar])[!unique(dmerge[,regvar]) %in% torm]
lfailed <- length(torm)
lpassed <- length(tokeep)


## Filtering summary
cat("\n")
paste1 <- paste0("Regions with length <", min.len, " in >", 100*max.frac, "% of individuals in each group:")
paste2 <- paste0("Regions with average coverage <", min.cov, " in >", 100*max.frac, "% of individuals in each group:")
paste3 <- paste0("Regions with average coverage >", max.cov, " in >", 100*max.frac, "% of individuals in each group:")
paste4 <- paste0("Regions with alignment fraction <", min.ratio, " in >", 100*max.frac, "% of individuals in each group:")
paste5 <- paste0("Regions without mapped reads in >", 100*max.frac, "% of individuals in each group:")

fillto <- max(c(nchar(paste1), nchar(paste2), nchar(paste3), nchar(paste4), nchar(paste5))) + 2
pasteA <- paste0("Number of ALL regions:", paste(rep(" ", fillto-22), collapse = ""), nreg, " (100%)\n")

paste1 <- paste0(paste1, paste(rep(" ", fillto-nchar(paste1)), collapse = ""), length(lowlen),  " (", round(100*length(lowlen)/nreg,2),"%)\n")
paste2 <- paste0(paste2, paste(rep(" ", fillto-nchar(paste2)), collapse = ""), length(lowcov),  " (", round(100*length(lowcov)/nreg,2),"%)\n")
paste3 <- paste0(paste3, paste(rep(" ", fillto-nchar(paste3)), collapse = ""), length(highcov)," (", round(100*length(highcov)/nreg,2),"%)\n")
paste4 <- paste0(paste4, paste(rep(" ", fillto-nchar(paste4)), collapse = ""), length(lowratio),  " (", round(100*length(lowratio)/nreg,2),"%)\n")
paste5 <- paste0(paste5, paste(rep(" ", fillto-nchar(paste5)), collapse = ""), length(lowmap),  " (", round(100*length(lowmap)/nreg,2),"%)\n")
paste6 <- paste0("Cumulative regions FAILED:", paste(rep(" ", fillto-26), collapse = ""), lfailed, " (", round(100*lfailed/nreg,2),"%)\n")
paste7 <- paste0("Cumulative regions PASSED:", paste(rep(" ", fillto-26), collapse = ""), lpassed, " (", round(100*lpassed/nreg,2),"%)\n")

cat(pasteA)
cat(paste1)
cat(paste2)
cat(paste3)
cat(paste4)
cat(paste5)
cat(paste6)
cat(paste7)


## Filter coverage data
dmerge.rm <- dmerge[dmerge[,regvar] %in% torm,]
dmerge.rm[,regvar] <- droplevels(dmerge.rm[,regvar])
dmerge.keep <- dmerge[!dmerge[,regvar] %in% torm,]
dmerge.keep[,regvar] <- droplevels(dmerge.keep[,regvar])


## Write LOG
dlog <- c(
  paste0("=== Filter regions based on visualization of coverage statistics ==="),
  paste0("Starting time: ", t1),
  paste0(""),
  paste0("Input data:"),
  paste0("-----------"),
  paste0("Sample file:     ", sfile),
  paste0("Group file:      ", meta),
  paste0("Mapping quality: ", mapQ),
  paste0(""),
  paste0("Thresholds for region filtering:"),
  paste0("-------------------------------"),
  paste0("Min. length:                                         ", min.len),
  paste0("Min. average coverage:                               ", min.cov),
  paste0("Max. average coverage:                               ", max.cov),
  paste0("Min. alignment fraction:                             ", min.ratio),
  paste0("Max. fraction of non-conforming samples per group:   ", max.frac),
  paste0(""),
  paste0("Other (internal) parameters:"),
  paste0("----------------------------"),
  paste0("Path to coverage stats: ", covstatfile),
  paste0("Hierarchical clustering of regions in heatmaps: ", as.character(sortx)),
  paste0("Hierarchical clustering method: ", hclustmethod),
  paste0("PDF height (inches): ", plot.height),
  paste0("PDF width (inches): ", plot.width),
  paste0(""),
  paste0("Region filtering results:"),
  paste0("------------------------"),
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
  paste0("Merged coverage stats:                  ", outtxt1),
  paste0("Passed regions:                         ", outtxt2),
  paste0("Failed regions:                         ", outtxt3),
  paste0("Visualization of coverage stats:        ", outpdf),
  paste0("Log file:                               ", logtxt),
  paste0("")
)


## Write TXT
cat("\nwriting results .txt and .log...\n")
write.table(dmerge, file = outtxt1, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
write.table(tokeep, file = outtxt2, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
write.table(torm,   file = outtxt3, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")


## Plot stats
nygr <- nlevels(dmerge[,grpvar2])
ycols <- rep(ycols, ceiling(nygr/length(ycols)))[1:nygr]

p1q <- round(quantile(dmerge[,pvar1], c(qlow, qhigh), na.rm = TRUE))
p1 <- ggplot(dmerge, aes_string(x = grpvar2, y = pvar1, fill = grpvar2)) +
  geom_boxplot(alpha = 0.5, na.rm = TRUE) +
  geom_violin(alpha = 0.5, na.rm = TRUE) +
  labs(x = "", y = pvar1) +
  geom_hline(yintercept = c(p1q[1], min.len, p1q[2]), 
             linetype = c(2, 1, 2), col = c(1, "tomato", 1)) +
  scale_fill_manual(guide = FALSE, values = ycols) +
  coord_flip() +
  theme_bw()

p2q <- round(quantile(dmerge[,pvar2], c(qlow, qhigh), na.rm = TRUE))
p2 <- ggplot(dmerge, aes_string(x = grpvar2, y = pvar2, fill = grpvar2)) +
  geom_boxplot(alpha = 0.5, na.rm = TRUE) +
  geom_violin(alpha = 0.5, na.rm = TRUE) +
  labs(x = "", y = pvar2) +
  geom_hline(yintercept = c(p2q[1], min.cov, max.cov, p2q[2]), 
             linetype = c(2, 1, 1, 2), col = c(1, "tomato", "tomato", 1)) +
  scale_fill_manual(guide = FALSE, values = ycols) +
  scale_y_log10() +
  coord_flip() +
  theme_bw()

p3q <- round(quantile(dmerge[,pvar3], c(qlow, qhigh), na.rm = TRUE), 2)
p3 <- ggplot(dmerge, aes_string(x = grpvar2, y = pvar3, fill = grpvar2)) +
  geom_boxplot(alpha = 0.5, na.rm = TRUE) +
  geom_violin(alpha = 0.5, na.rm = TRUE) +
  labs(x = "", y = "Alignment fraction (length of mapped region / length of reference sequence)") +
  geom_hline(yintercept = c(p3q[1], min.ratio, p3q[2]), 
             linetype = c(2, 1, 2), col = c(1, "tomato", 1)) +
  scale_fill_manual(guide = FALSE, values = ycols) +
  coord_flip() +
  theme_bw()

p4 <- ggplot(dmerge, aes_string(x = pvar3, y = pvar2)) +
  geom_vline(xintercept = c(p3q[1], min.ratio, p3q[2]),
             linetype = rep(c(2, 1, 2), nygr), col = rep(c(1, "tomato", 1), nygr)) +
  geom_hline(yintercept = c(p2q[1], min.cov, max.cov, p2q[2]),
             linetype = rep(c(2, 1, 1, 2), nygr), col = rep(c(1, "tomato", "tomato", 1), nygr)) +
  geom_point(aes_string(colour = grpvar2), size = .05, alpha = 0.3, na.rm = TRUE) +
  geom_density_2d(alpha = 1, colour = "black", na.rm = TRUE) +
  xlab("Alignment fraction (length of mapped region / length of reference sequence)") +
  ylab("Average coverage of mapped region") +
  scale_colour_manual(guide = FALSE, values = ycols) +
  scale_y_log10() +
  facet_wrap(as.formula(paste("~", grpvar2))) +
  theme_bw()

p5 <- ggplot(dmerge, aes_string(x = pvar1, y = pvar2)) +
  geom_vline(xintercept = c(p3q[1], min.ratio, p3q[2]),
             linetype = rep(c(2, 1, 2), nygr), col = rep(c(1, "tomato", 1), nygr)) +
  geom_hline(yintercept = c(p1q[1], min.cov, max.cov, p1q[2]),
             linetype = rep(c(2, 1, 1, 2), nygr), col = rep(c(1, "tomato", "tomato", 1), nygr)) +
  geom_point(aes_string(colour = grpvar2), size = .05, alpha = 0.3, na.rm = TRUE) +
  geom_density_2d(alpha = 1, colour = "black", na.rm = TRUE) +
  xlab("Length of mapped region") +
  ylab("Average coverage of mapped region") +
  scale_colour_manual(guide = FALSE, values = ycols) +
  scale_y_log10() +
  facet_wrap(as.formula(paste("~", grpvar2))) +
  theme_bw()

## Plot heatmaps
cat("\nplotting heatmaps...")
pdf(outpdf, width = plot.width, height = plot.height)

print(p1)
print(p2)
print(p3)
print(p4)
print(p5)

# length of mapped region
p6 <- suppressWarnings(do.heatmap(dat = dmerge,       xfac = regvar, yfac = idvar, znum = pvar1, ygr = grpvar2, ycols = ycols, limit = p1q[2], sortx = sortx,                     hclustmethod = hclustmethod, low = low, mid = mid, high = high, title = paste("ALL regions,", pvar1)))
suppressWarnings(print(p6))

p7 <- suppressWarnings(do.heatmap(dat = dmerge.keep,  xfac = regvar, yfac = idvar, znum = pvar1, ygr = grpvar2, ycols = ycols, limit = p1q[2], sortx = levels(p6$data[,regvar]),  hclustmethod = hclustmethod, low = low, mid = mid, high = high, title = paste("PASSED regions,", pvar1)))
suppressWarnings(print(p7))

p8 <- suppressWarnings(do.heatmap(dat = dmerge.rm,    xfac = regvar, yfac = idvar, znum = pvar1, ygr = grpvar2, ycols = ycols, limit = p1q[1], sortx = levels(p6$data[,regvar]),  hclustmethod = hclustmethod, low = low, mid = mid, high = high, title = paste("FAILED regions,", pvar1)))
suppressWarnings(print(p8))

# average coverage in mapped region
p9 <- suppressWarnings(do.heatmap(dat = dmerge,       xfac = regvar, yfac = idvar, znum = pvar2, ygr = grpvar2, ycols = ycols, limit = p2q[2], sortx = sortx,                     hclustmethod = hclustmethod, low = low, mid = mid, high = high, title = paste("ALL regions,", pvar2)))
suppressWarnings(print(p9))

p10 <- suppressWarnings(do.heatmap(dat = dmerge.keep, xfac = regvar, yfac = idvar, znum = pvar2, ygr = grpvar2, ycols = ycols, limit = p2q[2], sortx = levels(p9$data[,regvar]),  hclustmethod = hclustmethod, low = low, mid = mid, high = high, title = paste("PASSED regions,", pvar2)))
suppressWarnings(print(p10))

p11 <- suppressWarnings(do.heatmap(dat = dmerge.rm,   xfac = regvar, yfac = idvar, znum = pvar2, ygr = grpvar2, ycols = ycols, limit = p2q[1], sortx = levels(p9$data[,regvar]),  hclustmethod = hclustmethod, low = low, mid = mid, high = high, title = paste("FAILED regions,", pvar2)))
suppressWarnings(print(p11))

# alignment fraction
p12 <- suppressWarnings(do.heatmap(dat = dmerge,      xfac = regvar, yfac = idvar, znum = pvar3, ygr = grpvar2, ycols = ycols, limit = 1,      sortx = sortx,                     hclustmethod = hclustmethod, low = low, mid = mid, high = high, title = paste("ALL regions,", pvar3)))
suppressWarnings(print(p12))

p13 <- suppressWarnings(do.heatmap(dat = dmerge.keep, xfac = regvar, yfac = idvar, znum = pvar3, ygr = grpvar2, ycols = ycols, limit = 1,      sortx = levels(p12$data[,regvar]), hclustmethod = hclustmethod, low = low, mid = mid, high = high, title = paste("PASSED regions,", pvar3)))
suppressWarnings(print(p13))

p14 <- suppressWarnings(do.heatmap(dat = dmerge.rm,   xfac = regvar, yfac = idvar, znum = pvar3, ygr = grpvar2, ycols = ycols, limit = 1,      sortx = levels(p12$data[,regvar]), hclustmethod = hclustmethod, low = low, mid = mid, high = high, title = paste("FAILED regions,", pvar3)))
suppressWarnings(print(p14))

graphics.off()

t2 <- Sys.time()
cat("\n") ; paste0("Finish time: ", t2)

dlog <- c(dlog, paste0("Finish time: ", t2))
writeLines(text = dlog, con = logtxt)
