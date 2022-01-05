#!/cluster/apps/r/3.6.0_openblas/x86_64/bin/Rscript

## Usage: filter.visual.coverages.R <sfile> <stats> <refseqs> <min.pregion=0.3> <min.ptaxa=0.3> <min.len=500> <min.cov=10> <max.cov=1000> <min.ratio=0.5> <min.frac=0.9>

## Needs: ggplot2, ape, grid, VennDiagramm, tidyr

## Load libraries
suppressPackageStartupMessages(library(ape)) # read.fasta, GC.content
suppressPackageStartupMessages(library(ggplot2)) # ggplot, sym

## Author: simon.crameri@env.ethz.ch, Mar 2020

## Get arguments
args <- commandArgs(trailingOnly=TRUE)

## Check arguments
options(warning.length = 5000L)
if (!length(args) %in% c(3:10)) {
  stop("3 arguments needed (10 taken):
       REQUIRED
       1) <sfile|CHR>:              path to samples file. Header and tab-separation expected.
                                    Sample IDs must be in the first column. Group IDs can be specified in the 
                                    second column (if not specified, all samples are assumed to constitute one group).
                                    The group ID is used to apply region filtering criteria 4-9 within all considered 
                                    groups, to determine regions passing the filtering criteria in all groups.
                                    Samples that do not belong to any specified group (second column empty or 'NA') 
                                    will be displayed in summary plots but will not be considerd during region filtering. 
                                    Additional columns are ignored.
       2) <stats|CHR>:              path to alignment stats. Header and tab-separation expected.
                                    Sample IDs must be in the first column. Alignment stats must be in the following
                                    columns as defined in lines 114-122 of this script. Only alignment stats 
                                    of samples in <sfile> will be read (warns or stops if there is a mismatch).
       3) <refseqs|CHR>             path to region reference sequences. Fasta format expected. Used to correlate alignment
                                    stats with reference sequence lengths and GC content. Only regions in <stats>
                                    will be considered (warns or stops if there is a mismatch).

       OPTIONAL
       # sample filtering criteria (sample quality)
       4) <min.pregion|NUM>:        minimum fraction of regions recovered in a sample (i.e., sample has at least 1 
                                    mapped read in <min.pregion>*nregions regions) [DEFAULT: 0.3]
       
       # region filtering criteria (mapping sensitivity)
       5) <min.ptaxa|NUM>:          minimum fraction of samples recovered in a region (i.e., region has at least 1 
                                    mapped read in <min.ptaxa>*nsamples samples) [DEFAULT: 0.3]
                                    
       # region filtering criteria (mapping specificity)
       6) <min.len|NUM>:            minimum length in .bam [DEFAULT: 500];
       7) <min.cov|NUM>:            minimum coverage in .bam [DEFAULT: 10];
       8) <max.cov|NUM>:            maximum coverage in .bam [DEFAULT: 1000];
       9) <min.ratio|NUM>:          minimum alignment fraction [DEFAULT: 0.5];
       10) <min.frac|NUM>:          minimum fraction of samples conforming to the absolute filtering criteria 6-9
                                    (i.e., regions must meet criteria 5-8 in (100*<min.frac>)% of considered samples,
                                    separately for each considered group) [DEFAULT: 0.9]",
       call.=FALSE)
}

## Set arguments
sfile <- as.character(args[1])
stats <- as.character(args[2])
refseqs <- as.character(args[3])

min.pregion <- as.numeric(args[4])
if (is.na(min.pregion)) min.pregion <- 0.3
min.ptaxa <- as.numeric(args[5])
if (is.na(min.ptaxa)) min.ptaxa <- 0.3
min.len <- as.numeric(args[6])
if (is.na(min.len)) min.len <- 500
min.cov <- as.numeric(args[7])
if (is.na(min.cov)) min.cov <- 10
max.cov <- as.numeric(args[8])
if (is.na(max.cov)) max.cov <- 1000
min.ratio <- as.numeric(args[9])
if (is.na(min.ratio)) min.ratio <- 0.5
min.frac <- as.numeric(args[10])
if (is.na(min.frac)) min.frac <- 0.9

t1 <- Sys.time()
paste0("Starting time: ", t1)


## Set arguments (for debugging)
# sfile <- "mapfile.txt"
# stats <- "coverage_stats.txt"
# refseqs = "test.fasta"
# 
# min.pregion <- 0.95
# min.ptaxa <- 0.95
# min.len <- 100
# min.cov <- 8
# max.cov <- 1000
# min.ratio <- 0
# min.frac <- 0.1


## Additional arguments
# ouput paths
prefix <- tools::file_path_sans_ext(basename(stats))
suffix <- paste(min.pregion,min.ptaxa,min.len,min.cov,max.cov,min.ratio,min.frac, sep = "-")
outpdf <- paste0(prefix, "-", suffix, ".pdf")  # output .pdf (visualization of coverage stats)
outtxt1 <- paste0(prefix, "-", suffix, ".txt")         # output .txt (merged coverage stats)
outtxt2 <- paste0(prefix, "-regions-", suffix, ".txt")   # filtered regions
outtxt3 = paste0(prefix, "-taxa-", min.pregion, ".txt")  # filtered samples
logtxt <- paste0(prefix, "-", suffix, ".log")  # output .log
write.outtxt1 = TRUE                                 # if TRUE, will write merged alignment stats

# variables expected in <stat>
regvar <- "region"              # variable name in <stat> denoting <region id>
linref <- "linref"              # variable name in <stat> denoting <length of reference sequence>
linbam <- "linbam"              # variable name in <stat> denoting <length of mapped region>
cinref <- "cinref"              # variable name in <stat> denoting <average coverage of reference sequence>
cinbam <- "cinbam"              # variable name in <stat> denoting <average coverage of mapped region>

# variables created in merged data.frame
idvar <- "ID"                   # variable name denoting <sample ID>
grpvar1 <- "GROUP"              # variable name denoting <group ID>
grpvar2 <- "LABEL"              # variable name denoting <group ID (n = N)> 
faln <- "alignment_fraction"    # computed as linbam / linref
tlen <- "tlen"                  # variable name denoting <alignment length of best-matching contig>
tgc <- "tgc"                    # variable name denoting <GC content in target region>
navar <- "NA"                   # string used to denote membership to a not considered group

# plotting
qlow <- 0.1                     # lower quantile (shown on violin plots)
qhigh <- 0.9                    # upper quantile (shown on violin plots and used to cap heatmaps of ALL and PASSED loci)
hclustmethod <- "ward.D2"       # region clustering method in heatmap
sortx = TRUE                    # if TRUE, will sort regions based on <hclustmethod> in heatmap plots
ycols = c("#A6CEE3","#99CD91","#B89B74","#F06C45","#ED8F47","#825D99","#B15928") # group colors, will be recycled
low = "#A50026"                 # color used for the low end of the heatmap gradient
mid = "#F0E442"                 # color used for the mid point of the heatmap gradient
high = "#081D58"                # color used for the high end of the heatmap gradient
draw.venn = TRUE                # if TRUE, will plot VENN diagram of regions passing in each group (requires grid and VennDiagram libraries)
plot.pca = FALSE                # if TRUE, will plot PCAs of coverage statistics per region
plot.width <- 15                # output plot width
plot.height <- 15               # output plot height

##########################################################################################

## Check arguments
stopifnot(file.exists(sfile), 
          file.exists(stats),
          file.exists(refseqs),
          min.pregion >= 0, min.pregion <= 1,
          min.ptaxa >= 0, min.ptaxa <= 1,
          min.len >= 1,
          min.cov >= 1, 
          min.ratio >= 0, min.ratio <= 1, 
          min.frac >= 0, min.frac <= 1)


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
  
  # Internal arguments
  low.limit <- 0 # specifies the lower variable limit (corresponds to lower capping of heatmap color gradient)
  
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
    dat <- dat[order(dat[,xfac], (dat[,ygr]), (dat[,yfac])),]
    
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
                           midpoint = mean(c(low.limit, limit)), limits = c(low.limit, limit)) +
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

# PCA with ggplot2
ggpca <- function(dat, colvar, cols, x = 1, y = 2, scale = TRUE, center = TRUE, varscale = 5, labsize = 3) {
  library(ggplot2)
  isnum <- sapply(colnames(dat), FUN = function(x) {is.numeric(dat[,x])})
  isvar <- apply(dat[,isnum], 2, function(x) {var(x) > 0})
  dpca <- prcomp(dat[,isnum][,isvar], scale. = scale, center = center)
  dvar <- round(100*dpca$sdev^2/sum(dpca$sdev^2), 2)
  dPCA <- data.frame(dat[,names(isnum[!isnum])], dpca$x) ; names(dPCA)[1] <- colvar
  drot <- data.frame(dpca$rotation)[,c(x,y)] ; drot$labs <- rownames(drot)
  drot[,paste0("exp_",x)] <- drot[,x]*varscale ; drot[,paste0("exp2_",x)] <- drot[,x]*varscale*1.25
  drot[,paste0("exp_",y)] <- drot[,y]*varscale ; drot[,paste0("exp2_",y)] <- drot[,y]*varscale*1.1
  p <- ggplot(dPCA, aes_string(paste0("PC", x), paste0("PC", y), colour = colvar)) +
    geom_point() +
    geom_density_2d(na.rm = TRUE) +
    geom_segment(data = drot, size = 1, colour = "tomato", inherit.aes = FALSE,
                 aes_string(x = 0, y = 0, 
                            xend = paste0("exp_", x), 
                            yend = paste0("exp_", y)),
                 arrow = arrow(length = unit(4, "mm"))) +
    geom_text(data = drot, size = labsize, inherit.aes = FALSE,
              aes_string(x = paste0("exp2_", x), y = paste0("exp2_", y), 
                         label = "labs")) +
    scale_colour_manual(name = "", values = cols) +
    labs(x = paste0("PC ", x, " (", dvar[x], "%)"), y = paste0("PC ", y, " (", dvar[y], "%)")) +
    theme_bw()
  invisible(p)
}

##########################################################################################

## Read sample / group file
ds <- read.delim(sfile, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
if (ncol(ds) == 1) ds[,grpvar1] <- rep("Undefined", nrow(ds))
names(ds) <- c(idvar, grpvar1)
samples <- ds[,idvar]
ds[ds[,grpvar1] %in% c(NA, "NA", ""),grpvar1] <- navar # sets <NA> or <> group values to "NA"


## Read coverage stats
cat("\nreading coverage stats...\n")
fnames <- c(regvar, linref, linbam, cinref, cinbam)
dc <- read.delim(stats, check.names = FALSE)
names(dc)[1] <- idvar ; if (any(names(dc) %in% grpvar1)) dc <- dc[,-which(names(dc) %in% grpvar1)]
ssamples <- unique(dc[,idvar])

# check conformity of samples and alignment stats
if (any(!samples %in% ssamples)) {
  warning(length(samples[!samples %in% ssamples]), " samples in <", sfile, "> not found in <", stats, "> (sample ID expected in first column):\n", 
          paste(samples[!samples %in% ssamples], collapse = ", "), "\n")
  samples <- samples[samples %in% ssamples]
  ds <- ds[ds[,idvar] %in% samples,]
}
if (any(!ssamples %in% samples)) {
  warning(length(ssamples[!ssamples %in% samples]), " samples in <", stats, "> not found in <", sfile, "> (sample ID expected in first column):\n", 
          paste(ssamples[!ssamples %in% samples], collapse = ", "), "\n")
}
if (any(!fnames %in% names(dc))) {
  stop("some expected variables defined in script (lines 78-82) not found in <", stats, ">:\n", 
       paste(fnames[!fnames %in% names(dc)], collapse = ", "))
}

# subset alignment stats (only for samples in sfile)
dc <- dc[dc[,idvar] %in% samples,] ; rownames(dc) <- seq(nrow(dc))
for (i in (c(idvar, regvar))) dc[,i] <- droplevels(dc[,i])
ntaxa <- length(unique(dc[,idvar])) # number of samples
nreg <- length(unique(dc[,regvar])) # number of regions
ngroups <- length(unique(ds[,grpvar1])) # number of groups
nnsamples <- sum(ds[,grpvar1] == navar) # number of samples without specified group
cat(paste0("\nmerged coverage stats have ", nrow(dc), " rows, ", ntaxa, " samples from ", ngroups, " groups (", ifelse(nnsamples > 0, "1", "0"), " thereof set to NA), ", nreg, " regions\n\n"))
if (nrow(dc) != ntaxa * nreg) {
  warning(ntaxa * nreg, " rows in coverage stats expected, but only ", nrow(dc), " found.\nWill expand alignment stats for missing combinations (tidyr package needed).\n")
  dc <- data.frame(tidyr::complete(dc, !! sym(idvar), !! sym(regvar)))
}


## Read reference fasta
ref <- ape::read.FASTA(refseqs)
dgc <- numeric()
for (i in seq(length(ref))) {
  dgc <- c(dgc, ape::GC.content(ref[i]))
}
dref <- data.frame(names(ref), lengths(ref), dgc)
names(dref) <- c(regvar, tlen, tgc)


## Merge data
# merge alignment data (samples in sfile) with mapfile (matching samples, missing samples will be set to NA)
dmerge <- merge(dc, ds, by = idvar, all.x = TRUE, all.y = FALSE, sort = FALSE) ; rm(dc)
dmerge[,grpvar1] <- as.factor(dmerge[,grpvar1])

# merge with reference sequence data
if (any(!levels(dmerge[,regvar]) %in% levels(dref[,regvar]))) {
  exreg <-  which(!levels(dmerge[,regvar]) %in% levels(dref[,regvar]))
  stop(length(exreg), " regions in <", stats, "> not found in <", refseqs, "> (region ID expected in column <", regvar, ">):\n", 
       paste(levels(dmerge[,regvar])[!levels(dmerge[,regvar]) %in% levels(dref[,regvar])], collapse = ", "))
}
if (any(!levels(dref[,regvar]) %in% levels(dmerge[,regvar]))) {
  exreg <- which(!levels(dref[,regvar]) %in% levels(dmerge[,regvar]))
  warning(length(exreg), " regions in <", refseqs, "> not found in <", stats, "> (region ID expected in column <", regvar, ">):\n", 
          paste(levels(dref[,regvar])[!levels(dref[,regvar]) %in% levels(dmerge[,regvar])], collapse = ", "), "\n")
  dref <- dref[dref[,regvar] %in% levels(dmerge[,regvar]),]
  dref[,regvar] <- droplevels(dref[,regvar])
}
dmerge <- merge(dmerge, dref, by = regvar, all.x = TRUE, all.y = FALSE, sort = FALSE)

# print grouping
tnocons <- as.character(dmerge[!duplicated(dmerge[,idvar]) & dmerge[,grpvar1] == navar,idvar])
cat(paste0("sample grouping(s):\n"))
grtab <- table(dmerge[!duplicated(dmerge[,idvar]),grpvar1])
print(grtab)
cat("\n")
if (length(tnocons) > 0) {
  cat(paste0("NOTE: ", length(tnocons), " sample(s) have no group membership in <", sfile, "> and will not be considered for region filtering:\n"))
  print(tnocons)
  cat("\n")
}


## Add variables
# add number of individuals per group
dmerge[,grpvar2] <- as.character(dmerge[,grpvar1])
for (group in unique(dmerge[,grpvar1])) {
  dmerge[,grpvar2][dmerge[,grpvar2] == group] <- paste0(names(grtab[group]), " (n = ", grtab[group], ")")
}
dmerge[,grpvar2] <- as.factor(dmerge[,grpvar2])

# add alignmnet fraction
dmerge$mapped <- ifelse(!is.na(dmerge[,linbam]), 1, 0)
dmerge[,faln] <- dmerge[,linbam] / dmerge[,linref]

# log representations
dmerge[,paste0("log.", cinbam)] <- log10(dmerge[,cinbam]) ; dmerge[,paste0("log.", cinref)] <- log10(dmerge[,cinref])

## Sort according to groups (needed to create matrix)
varorder <- c(idvar, grpvar1, grpvar2, regvar, linbam, linref, cinbam, paste0("log.", cinbam), cinref, paste0("log.", cinref), faln, tlen, tgc, "mapped")
dmerge <- dmerge[order(dmerge[,regvar], dmerge[,grpvar1], dmerge[,idvar]),varorder]


## Filter samples based on min.pregion (minimum fraction of loci with data)
min.nreg <-  min.pregion * nreg # minimum required number of recovered loci
tax.preg <- tapply(X = dmerge[,linbam], INDEX = dmerge[,idvar], FUN = function(x) {sum(x > 0, na.rm = TRUE)/length(x)})
tax.reg <- tapply(X = dmerge[,linbam], INDEX = dmerge[,idvar], FUN = function(x) {sum(x > 0, na.rm = TRUE)})
tpassed <- names(tax.preg[tax.preg >= min.pregion]) # passed samples
ntpassed <- length(tpassed) # number of passed samples
ntfailed <- length(samples) - length(tpassed) # number of failed samples
dpreg <- data.frame(dmerge[match(names(tax.preg), table = dmerge[,idvar]),c(idvar, grpvar1, grpvar2)], tax.reg, tax.preg)
dpreg$passed <- ifelse(dpreg[,idvar] %in% tpassed, 1, 0)
dpreg$considered <- ifelse(dpreg[,grpvar1] == navar, 0, 1)
gpassed <- levels(droplevels(dpreg[dpreg$passed == 1,grpvar1])) # passed groups
ngpassed <- length(gpassed) # number of passed groups
tcons <- levels(droplevels(dpreg[dpreg$passed == 1 & dpreg$considered == 1,idvar])) # considered samples
ntcons <- length(tcons) # number of considered samples
gcons <- levels(droplevels(dpreg[dpreg$passed == 1 & dpreg$considered == 1,grpvar1])) # considered groups
ngcons<- length(gcons) # number of considered groups
cat(paste0(ntpassed, " (", round(100*ntpassed/ntaxa,2), "%) samples passed, ", 
           ntcons, " (", round(100*ntcons/ntaxa,2), "%) samples from ", 
           ngcons, " (", round(100*ngcons/ngroups,2), "%) groups considered\n"))

## Filter regions
regs.passed <- list() ; all.regs <- levels(droplevels(dmerge[,regvar]))
lowptaxa <- lowlen <- lowcov <- highcov <- lowratio <- character()
for (group in gcons) {
  dsub <- dmerge[dmerge[,grpvar1] == group & dmerge[,idvar] %in% tcons,] # only considered samples from considered groups
  
  # regions with zero reads in > (1-min.ptaxa) of samples per group
  torm1 <- as.character(subset(aggregate(dsub[,linbam], by = list(dsub[,regvar]), FUN = function(x) {sum(is.na(x))/length(x)}), x > (1-min.ptaxa))[,1])
  lowptaxa <- sort(unique(c(lowptaxa, torm1)))
  
  # regions with length in bam < min.len in > (1-min.frac) of samples per group
  torm2 <- find.fails(dsub, linbam, regvar, min.len, "min", min.frac)
  lowlen <- sort(unique(c(lowlen, torm2)))
  
  # regions with coverage in bam < min.cov in > (1-min.frac) of samples per group
  torm3 <- find.fails(dsub, cinbam, regvar, min.cov, "min", min.frac)
  lowcov <- sort(unique(c(lowcov, torm3)))
  
  # regions with coverage in bam > max.cov in > (1-min.frac) of samples per group
  torm4 <- find.fails(dsub, cinbam, regvar, max.cov, "max", min.frac)
  highcov <- sort(unique(c(highcov, torm4)))
  
  # regions with coverage in bam < faln in > (1-min.frac) of samples per group
  torm5 <- find.fails(dsub, faln, regvar, min.ratio, "min", min.frac)
  lowratio <- sort(unique(c(lowratio, torm5)))
  
  # list of regions passed in group
  regs.passed[[group]] <- all.regs[!all.regs %in% sort(unique(c(torm1, torm2, torm3, torm4, torm5)))]
  s.cons <- length(unique(dmerge[dmerge[,grpvar1] == group & dmerge[,idvar] %in% tpassed,idvar]))
  s.passed <- length(unique(dsub[,idvar]))
  r.passed <- length(regs.passed[[group]])
  names(regs.passed)[names(regs.passed) == group] <- paste0(group, " (", "samples passed: ", s.passed, " / ", s.cons, " [", round(100*s.passed/s.cons,2), "%] ; loci passed: ", r.passed, " / ", nreg, " [", round(100*r.passed/nreg,2), "%])")
}
torm <- sort(unique(c(lowptaxa, lowlen, lowcov, highcov, lowratio)))
tokeep <- all.regs[!all.regs %in% torm]
ptaxa.passed <- nreg - length(lowptaxa)
linbam.passed <- nreg - length(lowlen)
cinbam.passed <- nreg - length(unique(c(lowcov, highcov)))
faln.passed <- nreg - length(lowratio)
cinbam.faln.passed <- nreg - length(unique(c(lowcov, highcov, lowratio)))
cinbam.linbam.passed <- nreg - length(unique(c(lowcov, highcov, lowlen)))
lfailed <- length(torm)
lpassed <- length(tokeep)


## Filtering summary
cat("\n")
paste1 <- paste0("Regions without mapped reads in >", 100*(1-min.ptaxa), "% of considered samples in at least 1 group:")
paste2 <- paste0("Regions with mapped length <", min.len, " in >", 100*(1-min.frac), "% as above:")
paste3 <- paste0("Regions with average coverage <", min.cov, " in >", 100*(1-min.frac), "% as above:")
paste4 <- paste0("Regions with average coverage >", max.cov, " in >", 100*(1-min.frac), "% as above:")
paste5 <- paste0("Regions with alignment fraction <", min.ratio, " in >", 100*(1-min.frac), "% as above:")

fillto <- max(c(nchar(paste1), nchar(paste2), nchar(paste3), nchar(paste4), nchar(paste5))) + 2
pasteA <- paste0("Number of ALL regions:", paste(rep(" ", fillto-22), collapse = ""), nreg, " (100%)\n")

paste1 <- paste0(paste1, paste(rep(" ", fillto-nchar(paste1)), collapse = ""), length(lowptaxa),  " (", round(100*length(lowptaxa)/nreg,2),"%)\n")
paste2 <- paste0(paste2, paste(rep(" ", fillto-nchar(paste2)), collapse = ""), length(lowlen),  " (", round(100*length(lowlen)/nreg,2),"%)\n")
paste3 <- paste0(paste3, paste(rep(" ", fillto-nchar(paste3)), collapse = ""), length(lowcov),  " (", round(100*length(lowcov)/nreg,2),"%)\n")
paste4 <- paste0(paste4, paste(rep(" ", fillto-nchar(paste4)), collapse = ""), length(highcov)," (", round(100*length(highcov)/nreg,2),"%)\n")
paste5 <- paste0(paste5, paste(rep(" ", fillto-nchar(paste5)), collapse = ""), length(lowratio),  " (", round(100*length(lowratio)/nreg,2),"%)\n")
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


## Filter region data
dmerge.rm <- dmerge[dmerge[,regvar] %in% torm,]
dmerge.rm[,regvar] <- droplevels(dmerge.rm[,regvar])
dmerge.keep <- dmerge[!dmerge[,regvar] %in% torm,]
dmerge.keep[,regvar] <- droplevels(dmerge.keep[,regvar])


## Write LOG
dlog <- c(
  paste0("=== Filter samples and regions based on visualization of coverage statistics ==="),
  paste0("Starting time: ", t1),
  paste0(""),
  paste0("Input data:"),
  paste0("-----------"),
  paste0("Sample / Grouping file:                               ", sfile),
  paste0("Coverage statistics file:                             ", stats),
  paste0("Reference .fasta:                                     ", refseqs),
  paste0(""),
  paste0("Thresholds for sample filtering:"),
  paste0("-------------------------------"),
  paste0("Min. fraction of regions recovered in a sample:       ", min.pregion),
  paste0(""),
  paste0("Thresholds for region filtering:"),
  paste0("--------------------------------"),
  paste0("Min. fraction of samples recovered in a region:       ", min.ptaxa),
  paste0("Min. length:                                          ", min.len),
  paste0("Min. average coverage:                                ", min.cov),
  paste0("Max. average coverage:                                ", max.cov),
  paste0("Min. alignment fraction:                              ", min.ratio),
  paste0("Min. fraction of conforming samples per group:        ", min.frac),
  paste0(""),
  paste0("Other (internal) parameters:"),
  paste0("----------------------------"),
  paste0("Hierarchical clustering of regions in heatmaps:       ", as.character(sortx)),
  paste0("Hierarchical clustering method:                       ", hclustmethod),
  paste0("PDF height (inches):                                  ", plot.height),
  paste0("PDF width (inches):                                   ", plot.width),
  paste0(""),
  paste0("Taxon filtering results:"),
  paste0("------------------------"),
  paste0("Number of ALL samples:                                ", ntaxa),
  paste0("Number of PASSED samples:                             ", ntpassed, " (", round(100*ntpassed/ntaxa,2), "%)"),
  paste0("Number of considered samples (PASSED + non-NA group): ", ntcons, " (", round(100*ntcons/ntaxa,2), "%)"),
  paste(""),
  paste0("Number of ALL groups:                                 ", ngroups),
  paste0("Number of PASSED groups (at least 1 PASSED sample):   ", ngpassed, " (", round(100*ngpassed/ngroups,2), "%)"),
  paste0("Number of considered groups (PASSED + non-NA group):  ", ngcons, " (", round(100*ngcons/ngroups,2), "%)"),
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
  paste0("Merged coverage stats and metadata:                   ", outtxt1, ifelse(write.outtxt1, "", " (not written)")),
  paste0("PASSED regions:                                       ", outtxt2),
  paste0("PASSED samples:                                       ", outtxt3),
  paste0("Visualization of coverage stats:                      ", outpdf),
  paste0("Log file:                                             ", logtxt),
  paste0("")
)


## Write TXT
if (write.outtxt1) write.table(dmerge, file = outtxt1, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
write.table(tokeep, file = outtxt2, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
write.table(tpassed, file = outtxt3, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")


## Plot stats
nygr <- nlevels(dmerge[,grpvar2]) ; grpvar3 <- "LABELREV"
ycols <- rev(rep(ycols, ceiling(nygr/length(ycols)))[1:nygr])
dmerge[,grpvar3] <- factor(dmerge[,grpvar2], levels = rev(levels(dmerge[,grpvar2]))) # reverse label order for top-to-bottom y axis
dpreg[,grpvar3] <- factor(dpreg[,grpvar2], levels = rev(levels(dpreg[,grpvar2])))
dmerge <- dmerge[order(dmerge[,regvar], dmerge[,grpvar3], dmerge[,idvar]),c(varorder,grpvar3)]

# Taxon filtering
min.nreg2 <- floor(min(c(min(dpreg[,"tax.reg"], na.rm = TRUE), min.nreg)))
min.nreg2lab <- pretty(seq(min.nreg2/nreg, 1, length.out = 10))
p0 <- ggplot(dpreg, aes_string(x = grpvar3, y = "tax.reg", fill = grpvar3)) +
  geom_boxplot(alpha = 0.5) +
  geom_point(aes_string(colour = grpvar3)) +
  geom_hline(aes(yintercept = min.nreg), col = "tomato") +
  geom_hline(aes(yintercept = nreg)) +
  scale_fill_manual(guide = "none", values = ycols) +
  scale_colour_manual(guide = "none", values = ycols) +
  scale_y_continuous(sec.axis = sec_axis(~./nreg, breaks = min.nreg2lab, labels = min.nreg2lab), limits = c(min.nreg2, nreg)) +
  labs(x = "", y = "Number of regions recovered in a sample") +
  coord_flip() +
  theme_bw() +
  ggtitle(paste0("Number of regions: ", nreg, " ; required: ", min.nreg, " (", round(100*min.pregion,2), "%)\n",
                 "Number of ALL samples: ", ntaxa, " ; PASSED: ", ntpassed, " (", round(100*ntpassed/ntaxa,2), "%) ; FAILED: ", ntfailed, " (", round(100*ntfailed/ntaxa,2), "%)"))

p1q <- round(quantile(dmerge[,linbam], c(qlow, qhigh), na.rm = TRUE))
p1 <- ggplot(dmerge, aes_string(x = grpvar3, y = linbam, fill = grpvar3)) +
  geom_boxplot(alpha = 0.5, na.rm = TRUE) +
  geom_violin(alpha = 0.5, na.rm = TRUE) +
  labs(x = "", y = linbam) +
  geom_hline(yintercept = c(p1q[1], min.len, p1q[2]), 
             linetype = c(2, 1, 2), col = c(1, "tomato", 1)) +
  scale_fill_manual(guide = "none", values = ycols) +
  scale_y_log10() +
  coord_flip() +
  theme_bw() +
  ggtitle(paste0("Number of ALL regions: ", nreg, " ; PASSED: ", linbam.passed, " (", round(100*linbam.passed/nreg,2), "%) ; FAILED: ", nreg-linbam.passed, " (", round(100*(nreg-linbam.passed)/nreg,2), "%)"))

p2q <- round(quantile(dmerge[,cinbam], c(qlow, qhigh), na.rm = TRUE))
p2 <- ggplot(dmerge, aes_string(x = grpvar3, y = cinbam, fill = grpvar3)) +
  geom_boxplot(alpha = 0.5, na.rm = TRUE) +
  geom_violin(alpha = 0.5, na.rm = TRUE) +
  labs(x = "", y = cinbam) +
  geom_hline(yintercept = c(p2q[1], min.cov, max.cov, p2q[2]), 
             linetype = c(2, 1, 1, 2), col = c(1, "tomato", "tomato", 1)) +
  scale_fill_manual(guide = "none", values = ycols) +
  scale_y_log10() +
  coord_flip() +
  theme_bw() +
  ggtitle(paste0("Number of ALL regions: ", nreg, " ; PASSED: ", cinbam.passed, " (", round(100*cinbam.passed/nreg,2), "%) ; FAILED: ", nreg-cinbam.passed, " (", round(100*(nreg-cinbam.passed)/nreg,2), "%)"))

p3q <- round(quantile(dmerge[,faln], c(qlow, qhigh), na.rm = TRUE), 2)
p3 <- ggplot(dmerge, aes_string(x = grpvar3, y = faln, fill = grpvar3)) +
  geom_boxplot(alpha = 0.5, na.rm = TRUE) +
  geom_violin(alpha = 0.5, na.rm = TRUE) +
  labs(x = "", y = "Alignment fraction (length of mapped region / length of reference sequence)") +
  geom_hline(yintercept = c(p3q[1], min.ratio, p3q[2]), 
             linetype = c(2, 1, 2), col = c(1, "tomato", 1)) +
  scale_fill_manual(guide = "none", values = ycols) +
  coord_flip() +
  theme_bw() +
  ggtitle(paste0("Number of ALL regions: ", nreg, " ; PASSED: ", faln.passed, " (", round(100*faln.passed/nreg,2), "%) ; FAILED: ", nreg-faln.passed, " (", round(100*(nreg-faln.passed)/nreg,2), "%)"))

p4 <- ggplot(dmerge, aes_string(x = faln, y = cinbam)) +
  geom_vline(xintercept = c(p3q[1], min.ratio, p3q[2]),
             linetype = rep(c(2, 1, 2), nygr), col = rep(c(1, "tomato", 1), nygr)) +
  geom_hline(yintercept = c(p2q[1], min.cov, max.cov, p2q[2]),
             linetype = rep(c(2, 1, 1, 2), nygr), col = rep(c(1, "tomato", "tomato", 1), nygr)) +
  geom_point(aes_string(colour = grpvar2), size = .05, alpha = 0.3, na.rm = TRUE) +
  geom_density_2d(alpha = 1, colour = "black", na.rm = TRUE) +
  xlab("Alignment fraction (length of mapped region / length of reference sequence)") +
  ylab("Average coverage of mapped region") +
  scale_colour_manual(guide = "none", values = rev(ycols)) +
  scale_y_log10() +
  facet_wrap(as.formula(paste("~", grpvar2))) +
  theme_bw() +
  ggtitle(paste0("Number of ALL regions: ", nreg, " ; PASSED: ", cinbam.faln.passed, " (", round(100*cinbam.faln.passed/nreg,2), "%) ; FAILED: ", nreg-cinbam.faln.passed, " (", round(100*(nreg-cinbam.faln.passed)/nreg,2), "%)"))

p5 <- ggplot(dmerge, aes_string(x = linbam, y = cinbam)) +
  geom_vline(xintercept = c(p1q[1], min.len, p1q[2]),
             linetype = rep(c(2, 1, 2), nygr), col = rep(c(1, "tomato", 1), nygr)) +
  geom_hline(yintercept = c(p1q[1], min.cov, max.cov, p1q[2]),
             linetype = rep(c(2, 1, 1, 2), nygr), col = rep(c(1, "tomato", "tomato", 1), nygr)) +
  geom_point(aes_string(colour = grpvar2), size = .05, alpha = 0.3, na.rm = TRUE) +
  geom_density_2d(alpha = 1, colour = "black", na.rm = TRUE) +
  xlab("Length of mapped region") +
  ylab("Average coverage of mapped region") +
  scale_colour_manual(guide = "none", values = rev(ycols)) +
  scale_y_log10() +
  scale_x_log10() +
  facet_wrap(as.formula(paste("~", grpvar2))) +
  theme_bw() +
  ggtitle(paste0("Number of ALL regions: ", nreg, " ; PASSED: ", cinbam.linbam.passed, " (", round(100*cinbam.linbam.passed/nreg,2), "%) ; FAILED: ", nreg-cinbam.linbam.passed, " (", round(100*(nreg-cinbam.linbam.passed)/nreg,2), "%)"))

dref$PASSED <- factor(ifelse(dref[,regvar] %in% tokeep, "PASSED", "FAILED"), levels = c("PASSED","FAILED"))
dcon <- aggregate(dmerge[,faln], by = list(dmerge[,regvar]), FUN = mean, na.rm = TRUE)
dref[,faln] <- dcon[match(dref[,regvar], dcon$Group.1),]$x
p6 <- ggplot(dref, aes_string("tlen", "tgc", colour = "PASSED", alpha = faln, size = faln)) +
  geom_point(na.rm = TRUE) +
  geom_density_2d(aes(group = PASSED), colour = "black") +
  scale_colour_manual(guide = "none", values = c(high, low)) +
  scale_alpha_continuous(guide = "none", range = c(0.8,0.3)) +
  scale_size_continuous(name = "Alignment\nfraction") +
  labs(x = "Region length (in .fasta)", y = "Region GC content") +
  scale_x_log10() +
  theme_bw() +
  facet_wrap(~PASSED) +
  ggtitle(paste0("Number of ALL regions: ", nreg, " ; PASSED: ", lpassed, " (", round(100*lpassed/nreg,2), "%) ; FAILED: ", nreg-lpassed, " (", round(100*(nreg-lpassed)/nreg,2), "%)"))

# PCA 
nvar <- c(linbam, faln, tlen, tgc, "mapped", paste0("log.", cinbam), paste0("log.", cinref))
dmedL <- sapply(nvar, FUN = function(x) {tapply(X = dmerge[,x], INDEX = dmerge[,regvar], FUN = function(y) {median(y, na.rm = T)})})
datL <- na.omit(data.frame(dmedL, PASSED = dref[match(rownames(dmedL), dref[,regvar]),c("PASSED")]))

dmedT <- sapply(nvar, FUN = function(x) {tapply(X = dmerge[,x], INDEX = dmerge[,idvar], FUN = function(y) {median(y, na.rm = T)})})
datT <- na.omit(data.frame(dmedT, GROUP = dmerge[match(rownames(dmedT), dmerge[,idvar]),grpvar2]))

if (plot.pca) {
	pPCA1 <- ggpca(datL, colvar = "PASSED", cols = c(high, low)) +
  		facet_wrap(~PASSED) +
  		guides(color = "none") +
  		ggtitle(paste0("Number of ALL regions: ", nreg, " ; PASSED: ", lpassed, " (", round(100*lpassed/nreg,2), "%) ; FAILED: ", nreg-lpassed, " (", round(100*(nreg-lpassed)/nreg,2), "%)"))

	pPCA2 <- ggpca(datT, colvar = grpvar2, cols = rev(ycols), varscale = 2) +
  		facet_wrap(as.formula(paste("~", grpvar2))) +
  		guides(color = "none") +
  		ggtitle(paste0("Number of ALL regions: ", nreg, " ; PASSED: ", lpassed, " (", round(100*lpassed/nreg,2), "%) ; FAILED: ", nreg-lpassed, " (", round(100*(nreg-lpassed)/nreg,2), "%)"))
}

## Plots
cat("\nplotting...")
pdf(outpdf, width = plot.width, height = plot.height)

# single variables
print(p0)
print(p1)
print(p2)
print(p3)
print(p4)
print(p5)
print(p6)

if (plot.pca) {
	suppressWarnings(print(pPCA1))
	suppressWarnings(print(pPCA2))
}

# VENN diagram for regions that passed filters in all considered groups
if (draw.venn) {  
  # get color vector for passed + considered groups
  grlab <- levels(dmerge[,grpvar2])[levels(dmerge[,grpvar1]) %in% gcons]
  grcol <- rev(ycols)[which(grlab %in% levels(dmerge[,grpvar2]))]
  
  if (length(regs.passed) <= 5) {
    # load libraries
    suppressPackageStartupMessages(library(grid))
    suppressPackageStartupMessages(library(VennDiagram))
    
    # set plot margins and color function
    oldmar <- par()$mar
    par(mar = rep(0,4))
    
    # plot VENN diagram for passed loci
    plot.new()
    futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
    grid.draw(venn.diagram(regs.passed, filename = NULL, 
                           height = 3000, width = 3000,
                           cat.dist = rep(0, length(grcol)),
                           cat.col = grcol,
                           col = grcol,
                           fill = grcol,
                           alpha = 0.6,
                           cex = 2, cat.cex = 0,
                           lty = rep(2, length(grcol))))
    legend("topleft", names(regs.passed), text.col = grcol, bty = "n", cex = 1)
    par(mar = oldmar)
  }
}

# Barplot for regions that passed filters in each considered group (more than 5 allowed)
dp <- data.frame(array(data = NA, dim = c(0, 3)))
for (group in names(regs.passed)) {
  d <- data.frame(rep(group),
                  rep(gcons[which(names(regs.passed) == group)], nreg), 
                  c(rep("PASSED", length(regs.passed[[group]])), 
                    rep("FAILED", nreg - length(regs.passed[[group]]))),
                  stringsAsFactors = FALSE)
  dp <- rbind(dp, d)
}
names(dp) <- c("LAB", grpvar2, "PASSED") ; dp[,"ALLREG"] <- 1 ; dp[,"PASSED"] <- factor(dp[,"PASSED"], levels = c("FAILED","PASSED"))

dp2 <- aggregate(dp$PASSED, by = list(dp$LAB, dp[,grpvar2]), FUN = function(x) {length(x[x=="PASSED"])})
dp2[,"ALLREG"] <- nreg
names(dp2) <- c("LAB", grpvar2, "PASSED", "ALLREG") ; dp2 <- dp2[order(dp2[,"PASSED"]),] ; dp2$y_pos <- rep(1, nrow(dp2)) ; dp2$x_pos <- 1:nrow(dp2)
dp[,grpvar2] <- factor(dp[,grpvar2], levels = unique(dp2[,grpvar2]))

ggplot(dp, aes_string(grpvar2, "ALLREG", fill = "PASSED")) +
  geom_bar(stat="identity") +
  geom_text(data = dp2, inherit.aes = FALSE,
            aes_string(x = "x_pos", y = "y_pos", label = "LAB"),
            color = "white", hjust = 0, size = 3.5) +
  scale_fill_manual(name = "", values = c(low, high)) +
  labs(x = "", y = "Number of PASSED / FAILED regions") +
  coord_flip() +
  theme_bw() +
  theme(axis.text.y = element_blank()) +
  ggtitle(paste0("Number of ALL regions: ", nreg, " ; PASSED: ", lpassed, " (", round(100*lpassed/nreg,2), "%) ; FAILED: ", nreg-lpassed, " (", round(100*(nreg-lpassed)/nreg,2), "%)"))


# Heatmaps
# length of mapped region
p7 <- suppressWarnings(do.heatmap(dat = dmerge,       xfac = regvar, yfac = idvar, znum = linbam, ygr = grpvar2, ycols = rev(ycols), limit = p1q[2], sortx = sortx,                     hclustmethod = hclustmethod, low = low, mid = mid, high = high, title = paste("ALL regions,", linbam)))
suppressWarnings(print(p7))

p8 <- suppressWarnings(do.heatmap(dat = dmerge.keep,  xfac = regvar, yfac = idvar, znum = linbam, ygr = grpvar2, ycols = rev(ycols), limit = p1q[2], sortx = levels(p7$data[,regvar]),  hclustmethod = hclustmethod, low = low, mid = mid, high = high, title = paste("PASSED regions,", linbam)))
suppressWarnings(print(p8))

p9 <- suppressWarnings(do.heatmap(dat = dmerge.rm,    xfac = regvar, yfac = idvar, znum = linbam, ygr = grpvar2, ycols = rev(ycols), limit = min.len,sortx = levels(p7$data[,regvar]),  hclustmethod = hclustmethod, low = low, mid = mid, high = high, title = paste("FAILED regions,", linbam)))
suppressWarnings(print(p9))

# average coverage in mapped region
p10 <- suppressWarnings(do.heatmap(dat = dmerge,      xfac = regvar, yfac = idvar, znum = cinbam, ygr = grpvar2, ycols = rev(ycols), limit = p2q[2], sortx = sortx,                     hclustmethod = hclustmethod, low = low, mid = mid, high = high, title = paste("ALL regions,", cinbam)))
suppressWarnings(print(p10))

p11 <- suppressWarnings(do.heatmap(dat = dmerge.keep, xfac = regvar, yfac = idvar, znum = cinbam, ygr = grpvar2, ycols = rev(ycols), limit = p2q[2], sortx = levels(p10$data[,regvar]), hclustmethod = hclustmethod, low = low, mid = mid, high = high, title = paste("PASSED regions,", cinbam)))
suppressWarnings(print(p11))

p12 <- suppressWarnings(do.heatmap(dat = dmerge.rm,   xfac = regvar, yfac = idvar, znum = cinbam, ygr = grpvar2, ycols = rev(ycols), limit = min.cov,sortx = levels(p10$data[,regvar]), hclustmethod = hclustmethod, low = low, mid = mid, high = high, title = paste("FAILED regions,", cinbam)))
suppressWarnings(print(p12))

# alignment fraction
p13 <- suppressWarnings(do.heatmap(dat = dmerge,      xfac = regvar, yfac = idvar, znum = faln, ygr = grpvar2, ycols = rev(ycols), limit = 1,        sortx = sortx,                     hclustmethod = hclustmethod, low = low, mid = mid, high = high, title = paste("ALL regions,", faln)))
suppressWarnings(print(p13))

p14 <- suppressWarnings(do.heatmap(dat = dmerge.keep, xfac = regvar, yfac = idvar, znum = faln, ygr = grpvar2, ycols = rev(ycols), limit = 1,        sortx = levels(p13$data[,regvar]), hclustmethod = hclustmethod, low = low, mid = mid, high = high, title = paste("PASSED regions,", faln)))
suppressWarnings(print(p14))

p15 <- suppressWarnings(do.heatmap(dat = dmerge.rm,   xfac = regvar, yfac = idvar, znum = faln, ygr = grpvar2, ycols = rev(ycols), limit = min.frac, sortx = levels(p13$data[,regvar]), hclustmethod = hclustmethod, low = low, mid = mid, high = high, title = paste("FAILED regions,", faln)))
suppressWarnings(print(p15))

graphics.off()

t2 <- Sys.time()
cat("\n") ; paste0("Finish time: ", t2)

dlog <- c(dlog, paste0("Finish time: ", t2))
writeLines(text = dlog, con = logtxt)

