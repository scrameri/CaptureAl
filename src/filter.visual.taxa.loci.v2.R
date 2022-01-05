#!/cluster/apps/r/3.6.0_openblas/x86_64/bin/Rscript

## Usage: filter.visual.taxa.loci.R <samples.txt> <loci_stats.txt> <reference.fasta> <OPT: minploci> <OPT: minptaxa> <OPT: maxncontigs> <OPT: minnormbbestscore> <OPT: mintalsn> <OPT: mintfrac> <OPT: minbestscore> <OPT: minbestlength> <OPT: min.frac> <OPT: mapfile.txt>

## Load libraries
suppressPackageStartupMessages(library(ape)) # read.fasta
suppressPackageStartupMessages(library(ggplot2)) # ggplot

## Author: simon.crameri@env.ethz.ch, Mar 2020

## Get arguments
args <- commandArgs(trailingOnly=TRUE)

## Check arguments
options(warning.length = 3000L)
if (! length(args) %in% c(3:13)) {
  stop("At least 3 arguments needed (13 taken):
       REQUIRED
       1) <sfile|CHR>:              path to sample file. No header expected, sample ID in first column.
       2) <file|CHR>:               path to collected alignment stats. Header and tab separation 
                                    expected, sample ID expected in first column. Must contain further
                                    variables defined in lines 114-122 of script. Only alignment stats 
                                    of taxa in <sfile> will be read (warns or stops if there is a mismatch).
       3) <refseqs|CHR>             path to region reference sequences (.fasta). Only regions in <file>
                                    will be considered (warns or stops if there is a mismatch).
       
       OPTIONAL
       # taxon filtering criteria
       4) <min.pregion|NUM>:        minimum proportion of regions recovered in a taxon (i.e., taxon has at least 1 
                                    successfully aligned contig in <min.pregion>*nreg regions) [DEFAULT: 0.3]
       
       # locus filtering criteria
       5) <min.ptaxa|NUM>:          minimum proportion of taxa recovered in a locus (i.e., locus has at least 1 
                                    successfully aligned contig in <min.ptaxa>*ntaxa taxa) [DEFAULT: 0.3]
       6) <max.ncontigs|NUM>:       maximum number of (non-zero) contigs per region [DEFAULT: 1]
       7) <min.bestscore.norm|NUM>: minimum normalized alignment score (i.e., raw score / target alignment length) 
                                    of best alignment [DEFAULT: 1]
       8) <min.taln|NUM>:           minimum length of best alignment (in reference) [DEFAULT: 50]
       9) <min.tfrac|NUM>:          minimum fraction of best alignment (in reference) [DEFAULT: 0.2]
       10) <min.bestscore|NUM>:     minimum raw score of best alignment [DEFAULT: 1]
       11) <min.bestlength|NUM>:    minimum length of best aligned contig [DEFAULT: 1]
       12) <min.frac|NUM>:          minimum fraction of conforming taxa (i.e., meet criteria 6 to 11) [DEFAULT: 0.9];
       
       # grouping for locus filtering
       13) <meta|CHR>:              path to group file that maps taxon labels (1st column) to groups (2nd column). 
                                    Groups labelled as 'NA' will be displayed in plots but will not be considerd 
                                    during locus filtering. Header and tab separation expected. More taxa of different 
                                    order than in <file> are allowed, missing taxa will be labeled as 'NA'. If NULL,
                                    all taxa are assumed to constitute one group. If provided, the locus filtering 
                                    criteria 5 to 12 must be met for each considered group [DEFAULT: NULL]", 
       call.=FALSE)
}

## Set arguments
sfile <- as.character(args[1])
file <- as.character(args[2])
refseqs <- as.character(args[3])

min.pregion <- as.numeric(as.character(args[4]))
if (is.na(min.pregion)) min.pregion <- 0.3
min.ptaxa <- as.numeric(as.character(args[5]))
if (is.na(min.ptaxa)) min.ptaxa <- 0.3
max.ncontigs <- as.numeric(as.character(args[6]))
if (is.na(max.ncontigs)) max.ncontigs <- 1
min.bestscore.norm <- as.numeric(as.character(args[7]))
if (is.na(min.bestscore.norm)) min.bestscore.norm <- 1

min.taln <- as.numeric(as.character(args[8]))
if (is.na(min.taln)) min.taln <- 50
min.tfrac <- as.numeric(as.character(args[9]))
if (is.na(min.tfrac)) min.tfrac <- 0.2
min.bestscore <- as.numeric(as.character(args[10]))
if (is.na(min.bestscore)) min.bestscore <- 1
min.bestlength <- as.numeric(as.character(args[11]))
if (is.na(min.bestlength)) min.bestlength <- 1
min.frac <- as.numeric(args[12])
if (is.na(min.frac)) min.frac <- 0.9

meta <- as.character(args[13])

t1 <- Sys.time()
paste0("Starting time: ", t1)


## Set arguments (for debugging)
# sfile = "samples.txt"
# file = "loci_stats.txt"
# # ranges <- "loci_ranges.txt"
# refseqs = "consFabaceae_4c_1747.fasta"
# 
# min.pregion <- 0.2      # minimum proportion of loci with non-zero contigs (filters taxa)
# min.ptaxa <- 0.01       # minimum proportion of taxa with non-zero contigs (filters loci)
# max.ncontigs <- 2       # maximum number of non-zero contigs (filters loci)
# min.bestscore.norm <- 2 # minimum best normalized alignment score (filters loci)
# min.taln <- 80          # minimum alignment length (filters loci)
# min.tfrac <- 0          # minimum alignment length / target length (filters loci)
# min.bestscore <- 1      # minimum best alignment score (filters loci)
# min.bestlength <- 1     # minimum best contig length (filters loci)
# min.frac <- 0.5         # minimum fraction of conforming taxa (filters loci)
# 
# meta <- "mapfile.txt"


## Additional arguments
# ouput paths
suffix <- paste(min.pregion, min.ptaxa, max.ncontigs, min.bestscore.norm, min.taln, min.frac, sep = "-")
outpdf = paste0("loci_stats-", suffix, ".pdf")    # output .pdf (visualization of alignment stats)
outtxt1 = paste0("loci_stats-", suffix, ".txt")   # output .txt (merged alignment stats)
outtxt2 = paste0("loci_kept-", suffix, ".txt")    # filtered loci
outtxt3 = paste0("taxa_kept-", min.pregion, ".txt") # filtered taxa
logtxt = paste0("loci_stats-", suffix, ".log")    # output .log 

# variables in <file>
regvar <- "locus"               # variable name in <file> denoting <region id>
ncontigs <- "ncontigs"          # variable name in <file> denoting <number of contigs per region>
bestscore <- "bestscore"        # variable name in <file> denoting <raw alignment score of best-matching contig>
bestlength <- "bestlength"      # variable name in <file> denoting <length of best-matching contig>
query <- "query"                # variable name in <file> denoting <contig id>
qstart <- "qstart"              # variable name in <file> denoting <alignment start in query (contig)>
qend <- "qend"                  # variable name in <file> denoting <alignment end in query (contig)>
tstart <- "tstart"              # variable name in <file> denoting <alignment start in target (reference)>
tend <- "tend"                  # variable name in <file> denoting <alignment end in target (reference)>

# variables created in merged data.frame
idvar <- "ID"                   # variable name denoting <sample ID>
grpvar1 <- "GROUP"              # variable name denoting <group ID>
grpvar2 <- "LABEL"              # variable name denoting <group ID (n = N)>
ncontigsc <- "ncontigscat"      # variable name denoting <number of contigs category>
bestscoren <- "bestscore.norm"  # variable name denoting <normalized alignment score of best-matching contig>
tlen <- "tlen"                  # variable name denoting <alignment length of best-matching contig>
tgc <- "tgc"                    # variable name denoting <GC content in target region>
qaln <- "qaln"                  # variable name denoting <query (contig) alignment length> computed as qend-qstart
taln <- "taln"                  # variable name denoting <target (reference) alignment length> computed as tend-tstart
qfrac <- "qfrac"                # variable name denoting <query (contig) alignment fraction> computed as qaln / bestlength
tfrac <- "tfrac"                # variable name denoting <target (reference) alignment fraction> computed as taln / tlen
navar <- "NA"                   # string used to denote membership to a not considered group

# plotting
qlow <- 0.1                     # lower quantile (shown on violin plots and used to cap heatmaps of FAILED loci)
qhigh <- 0.9                    # upper quantile (shown on violin plots and used to cap heatmaps of ALL and PASSED loci)
XandMore <- 3                   # will create number of contigs categories for 0, 1, 2, ..., >=XandMore contigs
Ncont <- paste0(ncontigsc, "N") # used to plot heatmap of contig categories as: as.numeric(dmerge[,ncontigsc]-1)
hclustmethod <- "ward.D2"       # region clustering method in heatmap
sortx <- TRUE                   # if TRUE, will sort regions based on <hclustmethod> in heatmap plots
ycols <- c("#A6CEE3","#99CD91","#B89B74","#F06C45","#ED8F47","#825D99","#B15928") # group colors, will be recycled
low <- "#A50026"                # color used for the low end of the heatmap gradient
mid <- "#F0E442"                # color used for the mid point of the heatmap gradient
high <- "#081D58"               # color used for the high end of the heatmap gradient
plot.width <- 15                # output plot width
plot.height <- 15               # output plot height

##########################################################################################

## Check arguments
stopifnot(file.exists(sfile),
          file.exists(file),
          file.exists(refseqs),
          min.pregion >= 0, min.pregion <= 1,
          min.ptaxa >= 0, min.ptaxa <= 1,
          max.ncontigs >= 1,
          min.bestscore.norm >= 0,
          min.taln >= 1,
          min.tfrac >= 0, min.tfrac <= 1,
          min.bestscore >= 1,
          min.bestlength >= 1)
if (!is.na(meta)) stopifnot(file.exists(meta))



## Define helperfunctions
# find non-conforming regions
find.fails <- function(df, var, by, thresh, method, min.frac) {
  stopifnot(is.data.frame(df), is.character(var), var %in% names(df), 
            is.character(by), by %in% names(df), is.numeric(thresh), 
            method %in% c("min", "max"))
  switch(method,
         min = {fun <- function(x) {length(x[x<thresh & !is.na(x) & x>0])/length(x[!is.na(x) & x>0])}},
         max = {fun <- function(x) {length(x[x>thresh & !is.na(x) & x>0])/length(x[!is.na(x) & x>0])}})
  s <- aggregate(df[,var], by = list(dsub[,by]), FUN = fun)
  torm <- as.character(subset(s, x > (1-min.frac))[,1])
  return(torm)
}

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
      if (any(!sortx %in% levels(dat[,xfac]))) warning("some regions in <sortx> not found in <dat>!")
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
samples <- as.character(read.delim(sfile, header = FALSE, check.names = FALSE)$V1)


## Read alignment stats
cat("\nreading alignment stats...\n")
fnames <- c(regvar, ncontigs, bestscore, bestlength, query, qstart, qend, tstart, tend)
dc <- read.delim(file, check.names = FALSE)
names(dc)[1] <- idvar

# check alignment stats
if (any(!samples %in% dc[,idvar])) {
  warning("some taxa in <", sfile, "> not found in <", file, "> (sample ID expected in first column):\n", 
          paste(samples[!samples %in% dc[,idvar]], collapse = ", "), "\n")
}
if (any(!fnames %in% names(dc))) {
  stop("some expected variables defined in script (lines 114-122) not found in <", file, ">:\n", 
       paste(fnames[!fnames %in% names(dc)], collapse = ", "))
}

# subset alignment stats (only for taxa in sfile)
dc <- dc[dc[,idvar] %in% samples,]
for (i in (c(idvar, regvar, query))) dc[,i] <- droplevels(dc[,i])
ntaxa <- length(unique(dc[,idvar])) # number of taxa
nreg <- length(unique(dc[,regvar])) # number of regions
cat(paste0("\nalignment stats have ", nrow(dc), " rows, ", ntaxa, " taxa, ", nreg, " regions\n\n"))
if (nrow(dc) != ntaxa * nreg) {
  warning(ntaxa * nreg, " rows in alignment stats expected, but only ", nrow(dc), " found!\n")
}

## Read metadata
if (!is.na(meta)) {
  dmeta <- read.delim(meta)
} else {
  dmeta <- data.frame(samples, rep("Undefined", length(samples)))
}
names(dmeta)[1] <- idvar ; names(dmeta)[2] <- grpvar1

## Read reference fasta
ref <- ape::read.FASTA(refseqs)
dgc <- numeric()
for (i in seq(length(ref))) {
  dgc <- c(dgc, ape::GC.content(ref[i]))
}
dref <- data.frame(names(ref), lengths(ref), dgc)
names(dref) <- c(regvar, tlen, tgc)

## Merge data
cat("merging data...\n")

# merge alignment data (taxa in sfile) with mapfile (matching taxa, missing taxa will be set to NA)
dmerge <- merge(dc, dmeta, by = idvar, all.x = TRUE, all.y = FALSE, sort = FALSE)
dmerge[,grpvar1] <- as.character(dmerge[,grpvar1]) ; dmerge[,grpvar1][is.na(dmerge[,grpvar1])] <- navar
dmerge[,grpvar1] <- as.factor(dmerge[,grpvar1])
ngroups <- length(unique(dmerge[,grpvar1])) # number of groups

# merge with reference sequence data
if (any(!levels(dmerge[,regvar]) %in% levels(dref[,regvar]))) {
  stop(nreg, " region names in <", file, "> do not match ", length(ref), " region names in <", refseqs, ">:\n", 
       paste(levels(dmerge[,regvar])[!levels(dmerge[,regvar]) %in% levels(dref[,regvar])], collapse = ", "))
}
if (any(!levels(dref[,regvar]) %in% levels(dmerge[,regvar]))) {
  exreg <- which(!levels(dref[,regvar]) %in% levels(dmerge[,regvar]))
  warning(length(exreg), " regions in <", refseqs, "> not found in <", file, "> (region ID expected in column <", regvar, ">):\n", 
          paste(levels(dref[,regvar])[!levels(dref[,regvar]) %in% levels(dmerge[,regvar])], collapse = ", "), "\n")
}
dmerge <- merge(dmerge, dref, by = regvar, all.x = TRUE, all.y = FALSE, sort = FALSE)

# print grouping
tnocons <- as.character(dmerge[!duplicated(dmerge[,idvar]) & dmerge[,grpvar1] == navar,idvar])
cat(paste0("\n", ngroups, " taxon grouping(s) found:\n"))
grtab <- table(dmerge[!duplicated(dmerge[,idvar]),grpvar1])
print(grtab)
if (length(tnocons) > 0) {
  cat(paste0("\nWARNING: ", length(tnocons), " taxa are either missing or <NA> in <", meta, "> and will not be considered for region filtering:\n"))
  print(tnocons)
  cat("\n")
}


## Add variables
# add number of individuals per group (all taxa, all groups)
dmerge[,grpvar2] <- as.character(dmerge[,grpvar1])
for (group in unique(dmerge[,grpvar1])) {
  dmerge[,grpvar2][dmerge[,grpvar2] == group] <- paste0(names(grtab[group]), " (n = ", grtab[group], ")")
}
dmerge[,grpvar2] <- as.factor(dmerge[,grpvar2])

# add number of contigs category (ncontigsc)
XandMorelevs <- as.character(0:XandMore)
XandMorelevs[length(XandMorelevs)] <- paste0(">=", XandMore)
dmerge[,ncontigsc] <- as.character(dmerge[,ncontigs])
dmerge[,ncontigsc][dmerge[,ncontigs] >= XandMore ] <- paste0(">=", XandMore)
dmerge[,ncontigsc] <- factor(dmerge[,ncontigsc], levels = XandMorelevs)

# get alignment lengths and fractions (if available)
dmerge[,qaln] <- abs(dmerge[,qend] - dmerge[,qstart]) # query alignment length
dmerge[,taln] <- abs(dmerge[,tend] - dmerge[,tstart]) # target alignment length
dmerge[,qfrac] <- dmerge[,qaln] / dmerge[,bestlength] # aligned portion in query
dmerge[,tfrac] <- dmerge[,taln] / dmerge[,tlen] # aligned portion in target

# normalize bestscore for alignment length
# in the case of a Smith-Waterman alignment (affine:local, exhaustive yes), the raw scores are the sum of the substitution matrix scores and the gap penalties
# to normalize it, we divide the raw score by the target alignment length
# in the case of a Smith-Waterman alignment (affine:local, exhaustive yes), the upper limit of the normalized score is usually 5 (follows from the Smith-Waterman scoring algorithm)
#dmerge$bestscore.norm <- dmerge$bestscore / dmerge$length # if normalized for target length, normalized scores are low if the alignment length is small
dmerge$bestscore.norm <- dmerge$bestscore / dmerge$taln

# add number of contigs (capped numeric variable to plot the heatmap)
dmerge[,Ncont] <- as.numeric(dmerge[,ncontigsc])-1


## Sort according to groups (needed to create matrix)
varorder <- c(regvar, grpvar1, grpvar2, idvar, ncontigsc, Ncont, ncontigs, query, bestscore, bestscoren, bestlength, tlen, tgc, qaln, taln, qfrac, tfrac, names(dc)[!names(dc) %in% c(idvar, regvar, ncontigs, bestscore, bestlength, query)])
dmerge <- dmerge[order(dmerge[,regvar], dmerge[,grpvar1], dmerge[,idvar]),varorder]


## Filter taxa based on min.pregion (minimum proportion of loci with data)
min.nreg <-  min.pregion * nreg # minimum required number of recovered loci
tax.preg <- tapply(X = dmerge[,ncontigs], INDEX = dmerge[,idvar], FUN = function(x) {sum(x > 0, na.rm = TRUE)/length(x)})
tax.reg <- tapply(X = dmerge[,ncontigs], INDEX = dmerge[,idvar], FUN = function(x) {sum(x > 0, na.rm = TRUE)})
tpassed <- names(tax.preg[tax.preg >= min.pregion]) # passed taxa
ntpassed <- length(tpassed) # number of passed taxa
dpreg <- data.frame(dmerge[match(names(tax.preg), table = dmerge[,idvar]),c(idvar, grpvar1, grpvar2)], tax.reg, tax.preg)
dpreg$passed <- ifelse(dpreg[,idvar] %in% tpassed, 1, 0)
dpreg$considered <- ifelse(dpreg[,grpvar1] == navar, 0, 1)
gpassed <- levels(droplevels(dpreg[dpreg$passed == 1,grpvar1])) # passed groups
ngpassed <- length(gpassed) # number of passed groups
tcons <- levels(droplevels(dpreg[dpreg$passed == 1 & dpreg$considered == 1,idvar])) # considered taxa
ntcons <- length(tcons) # number of considered taxa
gcons <- levels(droplevels(dpreg[dpreg$passed == 1 & dpreg$considered == 1,grpvar1])) # considered groups
ngcons<- length(gcons) # number of considered groups
cat(paste0(ntpassed, " (", round(100*ntpassed/ntaxa,2), "%) taxa passed, ", 
           ntcons, " (", round(100*ntcons/ntaxa,2), "%) taxa from ", 
           ngcons, " (", round(100*ngcons/ngroups,2), "%) groups considered\n"))


## Filter regions based on min.ptaxa, max.ncontigs, min.bestscore.norm, min.taln, min.tfrac, min.bestscore, min.bestlength and min.frac
lowptaxa <- highncontigs <- lowbestscoren <- lowtaln <- lowtfrac <- lowbestscore <- lowbestlength <- character()
for (group in gcons) {
  dsub <- dmerge[dmerge[,grpvar1] == group & dmerge[,idvar] %in% tcons,] # only considered taxa from considered groups
  
  # regions with zero contigs in > (1-min.ptaxa) of taxa per group
  torm1 <- aggregate(dsub[,ncontigs], by = list(dsub[,regvar]), FUN = function(x) {sum(x==0)/length(x)}) # fraction missing
  lowptaxa <- as.character(subset(torm1, x > (1-min.ptaxa))[,1])
  
  # regions with ncontigs > max.ncontigs in > (1-min.frac) of taxa per group
  torm2 <- find.fails(dsub, ncontigs, regvar, max.ncontigs, "max", min.frac)
  highncontigs <- sort(unique(c(highncontigs, torm2)))
  
  # regions with bestscore.norm < min.bestscore.norm in > (1-min.frac) of taxa per group
  torm3 <- find.fails(dsub, bestscoren, regvar, min.bestscore.norm, "min", min.frac)
  lowbestscoren <- sort(unique(c(lowbestscoren, torm3)))
  
  # regions with taln < min.taln in > (1-min.frac) of taxa per group
  torm4 <- find.fails(dsub, taln, regvar, min.taln, "min", min.frac)
  lowtaln <- sort(unique(c(lowtaln, torm4)))
  
  # regions with tfrac < min.tfrac in > (1-min.frac) of taxa per group
  torm5 <- find.fails(dsub, tfrac, regvar, min.tfrac, "min", min.frac)
  lowtfrac <- sort(unique(c(lowtfrac, torm5)))
  
  # regions with bestscore < min.bestscore in > (1-min.frac) of taxa per group
  torm6 <- find.fails(dsub, bestscore, regvar, min.bestscore, "min", min.frac)
  lowbestscore <- sort(unique(c(lowbestscore, torm6)))
  
  # regions with bestlength < min.bestlength in > (1-min.frac) of taxa per group
  torm7 <- find.fails(dsub, bestlength, regvar, min.bestlength, "min", min.frac)
  lowbestlength <- sort(unique(c(lowbestlength, torm7)))
}
torm <- sort(unique(c(lowptaxa, highncontigs, lowbestscoren, lowtaln, lowtfrac, lowbestscore, lowbestlength)))
tokeep <- unique(dmerge[,regvar])[!unique(dmerge[,regvar]) %in% torm]
ptaxa.passed <- nreg - length(lowptaxa)
ncontigs.passed <- nreg - length(highncontigs)
bestscoren.passed <- nreg - length(lowbestscoren)
taln.passed <- nreg - length(lowtaln)
tfrac.passed <- nreg - length(lowtfrac)
bestscore.passed <- nreg - length(lowbestscore)
bestlength.passed <- nreg - length(lowbestlength)
lfailed <- length(torm)
lpassed <- length(tokeep)


## Region filtering summary
cat("\n")
paste1 <- paste0("Regions with contigs in <", 100*min.ptaxa, "% of considered taxa in each considered group:")
paste2 <- paste0("Regions with >", max.ncontigs, " contigs in >", 100*(1-min.frac), "% of considered + aligned taxa in each considered group:")
paste3 <- paste0("Regions with normalized best alignment score <", min.bestscore.norm, " in >", 100*(1-min.frac), "% as above:")
paste4 <- paste0("Regions with target alignment length <", min.taln, " in >", 100*(1-min.frac), "% as above:")
paste5 <- paste0("Regions with target alignment fraction <", min.tfrac, " in >", 100*(1-min.frac), "% as above:")
paste6 <- paste0("Regions with best alignment score <", min.bestscore, " in >", 100*(1-min.frac), "% as above:")
paste7 <- paste0("Regions with best contig length <", min.bestlength, " in >", 100*(1-min.frac), "% as above:")

fillto <- max(c(nchar(paste1), nchar(paste2), nchar(paste3), nchar(paste4), nchar(paste5))) + 2
pasteA <- paste0("Number of ALL regions:", paste(rep(" ", fillto-22), collapse = ""), nreg, " (100%)\n")

paste1 <- paste0(paste1, paste(rep(" ", fillto-nchar(paste1)), collapse = ""), length(lowptaxa),  " (", round(100*length(lowptaxa)/nreg,2),"%)\n")
paste2 <- paste0(paste2, paste(rep(" ", fillto-nchar(paste2)), collapse = ""), length(highncontigs),  " (", round(100*length(highncontigs)/nreg,2),"%)\n")
paste3 <- paste0(paste3, paste(rep(" ", fillto-nchar(paste3)), collapse = ""), length(lowbestscoren)," (", round(100*length(lowbestscoren)/nreg,2),"%)\n")
paste4 <- paste0(paste4, paste(rep(" ", fillto-nchar(paste4)), collapse = ""), length(lowtaln),  " (", round(100*length(lowtaln)/nreg,2),"%)\n")
paste5 <- paste0(paste5, paste(rep(" ", fillto-nchar(paste5)), collapse = ""), length(lowtfrac),  " (", round(100*length(lowtfrac)/nreg,2),"%)\n")
paste6 <- paste0(paste6, paste(rep(" ", fillto-nchar(paste6)), collapse = ""), length(lowbestscore),  " (", round(100*length(lowbestscore)/nreg,2),"%)\n")
paste7 <- paste0(paste7, paste(rep(" ", fillto-nchar(paste7)), collapse = ""), length(lowbestlength),  " (", round(100*length(lowbestlength)/nreg,2),"%)\n")
paste8 <- paste0("Cumulative regions FAILED:", paste(rep(" ", fillto-26), collapse = ""), lfailed, " (", round(100*lfailed/nreg,2),"%)\n")
paste9 <- paste0("Cumulative regions PASSED:", paste(rep(" ", fillto-26), collapse = ""), lpassed, " (", round(100*lpassed/nreg,2),"%)\n")

cat(pasteA)
cat(paste1)
cat(paste2)
cat(paste3)
cat(paste4)
cat(paste5)
cat(paste6)
cat(paste7)
cat(paste8)
cat(paste9)


## Filter region data
dmerge.rm <- dmerge[dmerge[,regvar] %in% torm,]
dmerge.rm[,regvar] <- droplevels(dmerge.rm[,regvar])
dmerge.keep <- dmerge[!dmerge[,regvar] %in% torm,]
dmerge.keep[,regvar] <- droplevels(dmerge.keep[,regvar])


## Write LOG
dlog <- c(
  paste0("=== Filter taxa and loci based on visualization of alignment statistics ==="),
  paste0("Starting time: ", t1),
  paste0(""),
  paste0("Input data:"),
  paste0("-----------"),
  paste0("Sample file:               ", sfile),
  paste0("Group file:                ", meta),
  paste0("Alignment statistics file: ", file),
  paste0("Reference .fasta:          ", refseqs),
  paste0(""),
  paste0("Thresholds for taxon filtering:"),
  paste0("-------------------------------"),
  paste0("Min. proportion of regions recovered in a taxon:    ", min.pregion),
  paste0(""),
  paste0("Thresholds for region filtering:"),
  paste0("--------------------------------"),
  paste0("Min. proportion of taxa recovered in a region:      ", min.ptaxa),
  paste0("Max. number of (non-zero) contigs per region:       ", max.ncontigs),
  paste0("Min. best normalized alignment score:               ", min.bestscore.norm),
  paste0("Min. best target alignment length:                  ", min.taln),
  paste0("Min. best target alignment fraction:                ", min.tfrac),
  paste0("Min. best raw alignment score:                      ", min.bestscore),
  paste0("Min. best contig length:                            ", min.bestlength),
  paste0("Min. fraction of conforming taxa per group:         ", min.frac),
  paste0(""),
  paste0("Other (internal) parameters:"),
  paste0("----------------------------"),
  paste0("Contig number summary categories:                   ", 0, " - >=", XandMore),
  paste0("Hierarchical clustering of regions in heatmaps:     ", as.character(sortx)),
  paste0("Hierarchical clustering method:                     ", hclustmethod),
  paste0("PDF height (inches):                                ", plot.height),
  paste0("PDF width (inches):                                 ", plot.width),
  paste0(""),
  paste0("Taxon filtering results:"),
  paste0("------------------------"),
  paste0("Number of ALL taxa:                                 ", ntaxa),
  paste0("Number of PASSED taxa:                              ", ntpassed, " (", round(100*ntpassed/ntaxa,2), "%)"),
  paste0("Number of considered taxa (PASSED on non-NA group): ", ntcons, " (", round(100*ntcons/ntaxa,2), "%)"),
  paste(""),
  paste0("Number of ALL groups:                               ", ngroups),
  paste0("Number of PASSED groups (at least 1 PASSED taxon):  ", ngpassed, " (", round(100*ngpassed/ngroups,2), "%)"),
  paste0("Number of considered groups (PASSED non-NA group):  ", ngcons, " (", round(100*ngcons/ngroups,2), "%)"),
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
  substring(paste8, 1, nchar(paste8)-1),
  substring(paste9, 1, nchar(paste9)-1),
  paste(""),
  paste0("Output files:"),
  paste0("-------------"),
  paste0("Merged region stats and metadata: ", outtxt1),
  paste0("Passed regions:                   ", outtxt2),
  paste0("Passed taxa:                      ", outtxt3),
  paste0("Visualization of alignment stats: ", outpdf),
  paste0("Log file:                         ", logtxt),
  paste0("")
)


## Write TXT
cat("\nwriting results .txt and .log...\n")
write.table(dmerge, file = outtxt1, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
write.table(tokeep, file = outtxt2, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
write.table(tpassed, file = outtxt3, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")


## Plot stats
nygr <- nlevels(dmerge[,grpvar2])
ycols <- rep(ycols, ceiling(nygr/length(ycols)))[1:nygr]

# Taxon filtering
p0 <- ggplot(dpreg, aes_string(x = grpvar2, y = "tax.reg", fill = grpvar2)) +
  geom_boxplot(alpha = 0.5) +
  # geom_violin(alpha = 0.5) +
  geom_point(aes(color = LABEL)) +
  geom_hline(aes(yintercept = min.nreg), col = "tomato") +
  geom_hline(aes(yintercept = nreg)) +
  scale_fill_manual(guide = FALSE, values = ycols) +
  scale_colour_manual(guide = FALSE, values = ycols) +
  scale_y_continuous(sec.axis = sec_axis(~./nreg, breaks = seq(0, 1, by = 0.1), labels = seq(0, 1, by = 0.1)), limits = c(0, nreg)) +
  labs(x = "", y = "Number of target regions recovered in a taxon") +
  coord_flip() +
  theme_bw() +
  ggtitle(paste0("Number of regions: ", nreg, " ; required: ", min.nreg, " (", round(100*min.pregion,2), "%)\n",
                 "Number of ALL taxa: ", ntaxa, " ; PASSED: ", ntpassed, " (", round(100*ntpassed/ntaxa,2), "%)"))

p1q <- round(quantile(dmerge[dmerge[,ncontigs] > 0, ncontigs], c(qlow, qhigh, 0.99), na.rm = TRUE))
p1 <- ggplot(dmerge[dmerge[,ncontigs] > 0,], aes_string(x = grpvar2, y = ncontigs, fill = grpvar2)) +
  geom_boxplot(alpha = 0.5, na.rm = TRUE) +
  geom_violin(alpha = 0.5, na.rm = TRUE) +
  labs(x = "", y = "Number of non-zero contigs per region") +
  geom_hline(yintercept = c(p1q[1], max.ncontigs, p1q[2]), 
             linetype = c(2, 1, 2), col = c(1, "tomato", 1)) +
  scale_fill_manual(guide = FALSE, values = ycols) +
  scale_y_log10(breaks=c(1:10,seq(20,100,by=10))[c(1:10,seq(20,100,by=10))<max(dmerge[,ncontigs],na.rm=T)]) +
  coord_flip() +
  theme_bw() +
  ggtitle(paste0("Number of ALL regions: ", nreg, " ; PASSED: ", ncontigs.passed, " (", round(100*ncontigs.passed/nreg,2), "%)"))

p2q <- round(quantile(dmerge[,bestscoren], c(qlow, qhigh), na.rm = TRUE), 2)
p2 <- ggplot(dmerge, aes_string(x = grpvar2, y = bestscoren, fill = grpvar2)) +
  geom_boxplot(alpha = 0.5, na.rm = TRUE) +
  geom_violin(alpha = 0.5, na.rm = TRUE) +
  labs(x = "", y = "Normalized alignment score of best-matching contig") +
  geom_hline(yintercept = c(p2q[1], min.bestscore.norm, p2q[2]), 
             linetype = c(2, 1, 2), col = c(1, "tomato", 1)) +
  scale_fill_manual(guide = FALSE, values = ycols) +
  # scale_y_log10() +
  coord_flip() +
  theme_bw() +
  ggtitle(paste0("Number of ALL regions: ", nreg, " ; PASSED: ", bestscoren.passed, " (", round(100*bestscoren.passed/nreg,2), "%)"))

p3q <- round(quantile(dmerge[,bestlength], c(qlow, qhigh), na.rm = TRUE))
p3 <- ggplot(dmerge, aes_string(x = grpvar2, y = bestlength, fill = grpvar2)) +
  geom_boxplot(alpha = 0.5, na.rm = TRUE) +
  geom_violin(alpha = 0.5, na.rm = TRUE) +
  labs(x = "", y = "Length of best-matching contig") +
  geom_hline(yintercept = c(p3q[1], min.bestlength, p3q[2]), 
             linetype = c(2, 1, 2), col = c(1, "tomato", 1)) +
  scale_fill_manual(guide = FALSE, values = ycols) +
  coord_flip() +
  theme_bw() +
  ggtitle(paste0("Number of ALL regions: ", nreg, " ; PASSED: ", bestlength.passed, " (", round(100*bestlength.passed/nreg,2), "%)"))

p4q <- round(quantile(dmerge[,taln], c(qlow, qhigh), na.rm = TRUE))
p4 <- ggplot(dmerge, aes_string(x = grpvar2, y = taln, fill = grpvar2)) +
  geom_boxplot(alpha = 0.5, na.rm = TRUE) +
  geom_violin(alpha = 0.5, na.rm = TRUE) +
  labs(x = "", y = "Length of best alignment (in reference)") +
  geom_hline(yintercept = c(p4q[1], min.taln, p4q[2]), 
             linetype = c(2, 1, 2), col = c(1, "tomato", 1)) +
  scale_fill_manual(guide = FALSE, values = ycols) +
  coord_flip() +
  theme_bw() +
  ggtitle(paste0("Number of ALL regions: ", nreg, " ; PASSED: ", taln.passed, " (", round(100*taln.passed/nreg,2), "%)"))

p5q <- round(quantile(dmerge[,tfrac], c(qlow, qhigh), na.rm = TRUE),2)
p5 <- ggplot(dmerge, aes_string(x = grpvar2, y = tfrac, fill = grpvar2)) +
  geom_boxplot(alpha = 0.5, na.rm = TRUE) +
  geom_violin(alpha = 0.5, na.rm = TRUE) +
  labs(x = "", y = "Fraction of best alignment (reference length / alignment length in reference)") +
  geom_hline(yintercept = c(p5q[1], min.tfrac, p5q[2]), 
             linetype = c(2, 1, 2), col = c(1, "tomato", 1)) +
  scale_fill_manual(guide = FALSE, values = ycols) +
  coord_flip() +
  theme_bw() +
  ggtitle(paste0("Number of ALL regions: ", nreg, " ; PASSED: ", tfrac.passed, " (", round(100*tfrac.passed/nreg,2), "%)"))

p6 <- ggplot(dmerge, aes_string(x = bestscoren, y = tgc)) +
  geom_vline(xintercept = c(p2q[1], min.bestscore.norm, p2q[2]),
             linetype = rep(c(2, 1, 2), nygr), col = rep(c(1, "tomato", 1), nygr)) +
  geom_point(aes_string(colour = grpvar2), size = .05, alpha = 0.3, na.rm = TRUE) +
  geom_density_2d(alpha = 1, colour = "black", na.rm = TRUE) +
  xlab("Normalized alignment score of best-matching contig") +
  ylab("GC content in reference region") +
  scale_colour_manual(guide = FALSE, values = ycols) +
  scale_y_log10() +
  facet_wrap(as.formula(paste("~", grpvar2))) +
  theme_bw() +
  ggtitle(paste0("Number of ALL regions: ", nreg, " ; PASSED: ", bestscoren.passed, " (", round(100*bestscoren.passed/nreg,2), "%)"))

dref$PASSED <- factor(ifelse(dref$locus %in% tokeep, "PASSED", "FAILED"), levels = c("PASSED","FAILED"))
dcon <- aggregate(dmerge[,Ncont], by = list(dmerge[,regvar]), FUN = mean, na.rm = TRUE)
dref[,ncontigs] <- dcon[match(dref$locus, dcon$Group.1),]$x
p7 <- ggplot(dref, aes_string("tlen", "tgc", colour = "PASSED", alpha = ncontigs, size = ncontigs)) +
  geom_point() +
  geom_density_2d(aes(group = PASSED), colour = "black") +
  scale_colour_manual(guide = FALSE, values = c(high, low)) +
  scale_alpha_continuous(guide = FALSE, range = c(0.8,0.3)) +
  scale_size_continuous(range = c(0.5,3.5), name = "Nb. of\nconotigs", labels = XandMorelevs) +
  labs(x = "Region length (in .fasta)", y = "Region GC content") +
  theme_bw() +
  facet_wrap(~PASSED) +
  ggtitle(paste0("Number of ALL regions: ", nreg, " ; PASSED: ", lpassed, " (", round(100*lpassed/nreg,2), "%)"))


## Plot heatmaps
cat("\nplotting heatmaps...")
pdf(outpdf, width = plot.width, height = plot.height)

print(p0)
print(p1)
print(p2)
print(p3)
print(p4)
print(p5)
print(p6)
print(p7)

# number of contigs (min.ptaxa, max.ncontigs)
p8 <- suppressWarnings(do.heatmap(dat = dmerge,       xfac = regvar, yfac = idvar, znum = Ncont,       ygr = grpvar2, ycols = ycols, limit = XandMore, sortx = sortx,                     hclustmethod = hclustmethod, low = low, mid = mid, high = high, title = paste("ALL regions,", "capped number of contigs")))
suppressWarnings(print(p8))

p9 <- suppressWarnings(do.heatmap(dat = dmerge.keep,  xfac = regvar, yfac = idvar, znum = Ncont,       ygr = grpvar2, ycols = ycols, limit = XandMore, sortx = levels(p8$data[,regvar]),  hclustmethod = hclustmethod, low = low, mid = mid, high = high, title = paste("PASSED regions,", "capped number of contigs")))
suppressWarnings(print(p9))

p10 <- suppressWarnings(do.heatmap(dat = dmerge.rm,   xfac = regvar, yfac = idvar, znum = Ncont,       ygr = grpvar2, ycols = ycols, limit = XandMore, sortx = levels(p8$data[,regvar]),  hclustmethod = hclustmethod, low = low, mid = mid, high = high, title = paste("FAILED regions,", "capped number of contigs")))
suppressWarnings(print(p10))

# best normalized alignment score (min.bestscore.norm)
p11 <- suppressWarnings(do.heatmap(dat = dmerge,      xfac = regvar, yfac = idvar, znum = bestscoren, ygr = grpvar2, ycols = ycols, limit = p2q[2],    sortx = sortx,                     hclustmethod = hclustmethod, low = low, mid = mid, high = high, title = paste("ALL regions,", "norm. alignment score")))
suppressWarnings(print(p11))

p12 <- suppressWarnings(do.heatmap(dat = dmerge.keep, xfac = regvar, yfac = idvar, znum = bestscoren, ygr = grpvar2, ycols = ycols, limit = p2q[2],    sortx = levels(p11$data[,regvar]), hclustmethod = hclustmethod, low = low, mid = mid, high = high, title = paste("PASSED regions,", "norm. alignment score")))
suppressWarnings(print(p12))

p13 <- suppressWarnings(do.heatmap(dat = dmerge.rm,   xfac = regvar, yfac = idvar, znum = bestscoren, ygr = grpvar2, ycols = ycols, limit = p2q[2],    sortx = levels(p11$data[,regvar]), hclustmethod = hclustmethod, low = low, mid = mid, high = high, title = paste("FAILED regions,", "norm. alignment score")))
suppressWarnings(print(p13))

# best contig length (min.bestlength)
p14 <- suppressWarnings(do.heatmap(dat = dmerge,      xfac = regvar, yfac = idvar, znum = bestlength, ygr = grpvar2, ycols = ycols, limit = p3q[2],    sortx = sortx,                     hclustmethod = hclustmethod, low = low, mid = mid, high = high, title = paste("ALL regions,", "contig length")))
suppressWarnings(print(p14))

p15 <- suppressWarnings(do.heatmap(dat = dmerge.keep, xfac = regvar, yfac = idvar, znum = bestlength, ygr = grpvar2, ycols = ycols, limit = p3q[2],    sortx = levels(p14$data[,regvar]), hclustmethod = hclustmethod, low = low, mid = mid, high = high, title = paste("PASSED regions,", "contig length")))
suppressWarnings(print(p15))

p16 <- suppressWarnings(do.heatmap(dat = dmerge.rm,   xfac = regvar, yfac = idvar, znum = bestlength, ygr = grpvar2, ycols = ycols, limit = p3q[2],    sortx = levels(p14$data[,regvar]), hclustmethod = hclustmethod, low = low, mid = mid, high = high, title = paste("FAILED regions,", "contig length")))
suppressWarnings(print(p16))

# region alignment length (min.taln)
p17 <- suppressWarnings(do.heatmap(dat = dmerge,      xfac = regvar, yfac = idvar, znum = taln,       ygr = grpvar2, ycols = ycols, limit = p4q[2],    sortx = sortx,                     hclustmethod = hclustmethod, low = low, mid = mid, high = high, title = paste("ALL regions,", "alignment length")))
suppressWarnings(print(p17))

p18 <- suppressWarnings(do.heatmap(dat = dmerge.keep, xfac = regvar, yfac = idvar, znum = taln,       ygr = grpvar2, ycols = ycols, limit = p4q[2],    sortx = levels(p17$data[,regvar]), hclustmethod = hclustmethod, low = low, mid = mid, high = high, title = paste("PASSED regions,", "alignment length")))
suppressWarnings(print(p18))

p19 <- suppressWarnings(do.heatmap(dat = dmerge.rm,   xfac = regvar, yfac = idvar, znum = taln,       ygr = grpvar2, ycols = ycols, limit = p4q[2],    sortx = levels(p17$data[,regvar]), hclustmethod = hclustmethod, low = low, mid = mid, high = high, title = paste("FAILED regions,", "alignment length")))
suppressWarnings(print(p19))

# region alignment fraction (min.tfrac)
p20 <- suppressWarnings(do.heatmap(dat = dmerge,      xfac = regvar, yfac = idvar, znum = tfrac,      ygr = grpvar2, ycols = ycols, limit = 1,         sortx = sortx,                     hclustmethod = hclustmethod, low = low, mid = mid, high = high, title = paste("ALL regions,", "alignment fraction")))
suppressWarnings(print(p20))

p21 <- suppressWarnings(do.heatmap(dat = dmerge.keep, xfac = regvar, yfac = idvar, znum = tfrac,      ygr = grpvar2, ycols = ycols, limit = 1,         sortx = levels(p20$data[,regvar]), hclustmethod = hclustmethod, low = low, mid = mid, high = high, title = paste("PASSED regions,", "alignment fraction")))
suppressWarnings(print(p21))

p22 <- suppressWarnings(do.heatmap(dat = dmerge.rm,   xfac = regvar, yfac = idvar, znum = tfrac,      ygr = grpvar2, ycols = ycols, limit = 1,         sortx = levels(p20$data[,regvar]), hclustmethod = hclustmethod, low = low, mid = mid, high = high, title = paste("FAILED regions,", "alignment fraction")))
suppressWarnings(print(p22))

graphics.off()

t2 <- Sys.time()
cat("\n") ; paste0("Finish time: ", t2)

dlog <- c(dlog, paste0("Finish time: ", t2))
writeLines(text = dlog, con = logtxt)

