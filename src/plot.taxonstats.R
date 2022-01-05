#!/cluster/apps/r/3.6.0_openblas/x86_64/bin/Rscript

## Usage: plot.taxonstats.R <taxon stats>

## Load libraries
suppressPackageStartupMessages(library(data.table)) # data.table
suppressPackageStartupMessages(library(ggplot2)) # ggplot

## Author: simon.crameri@env.ethz.ch, Mar 2020

## Get arguments
args <- commandArgs(trailingOnly=TRUE)

## Check arguments
options(warning.length = 5000L)
if (! length(args) %in% c(1,2)) {
  stop("1 arguments needed (2 taken):
       REQUIRED
       1) <stats|CHR>:  path to taxon stats
       
       OPTIONAL
       2) <sfile|CHR>: path to sample mapping file",
       call.=FALSE)
}

## Set arguments
stats <- as.character(args[1])
sfile <- as.character(args[2])

## Set arguments (for debugging)
# stats <- "test.tstats.txt"
# sfile <- "mapfile.txt"

## Additional arguments
ycols = c("#A6CEE3","#99CD91","#B89B74","#F06C45","#ED8F47","#825D99","#B15928") # group colors, will be recycled
highquant = 0.9
width = 15
height = 15
oname = gsub(".txt$", ".pdf", stats)

## Helperfunctions
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


## Check arguments
stopifnot(file.exists(stats))
if (!is.na(sfile)) stopifnot(file.exists(sfile))

## read stats
ds <- fread(stats, stringsAsFactors = TRUE)

## read sfile
if (!is.na(sfile)) {
  dmeta <- read.delim(sfile, stringsAsFactors = F)
  names(dmeta)[1] <- "ID"
  names(dmeta)[2] <- "GROUP"
} else {
  dmeta <- data.frame(ID = unique(unlist(ds[,1])), GROUP = "NA", stringsAsFactors = F)
}
dmeta$ID <- as.factor(dmeta$ID)
dmeta$GROUP[dmeta$GROUP %in% c("NA", NA)] <- "NA"
dmeta$GROUP <- factor(dmeta$GROUP, levels = sort(unique(dmeta$GROUP)))
dmeta <- dmeta[dmeta$ID %in% unique(ds$ID),]
dmeta$ID <- droplevels(dmeta$ID)
dmeta$GROUP <- droplevels(dmeta$GROUP)

## merge
dmerge <- merge(ds, dmeta, by = "ID", all.x = FALSE, all.y = TRUE, sort = FALSE)

## add number of individuals per group (all samples, all groups)
dmerge[,LABEL := as.character(GROUP)]
grtab <- table(dmerge[!duplicated(ID),GROUP])
for (group in levels(unlist(dmerge[,GROUP]))) {
  dmerge[GROUP == group, LABEL := paste0(names(grtab[group]), " (n = ", grtab[group], ")")] 
}
dmerge[,LABEL := factor(LABEL)]
dmerge[,completeness := 1-(mis/len)]

## correct numbers for sequence lengths
dmerge[,c_trans := len*(trans/(len-mis))]
dmerge[,c_transv := len*(transv/(len-mis))]
dmerge[,c_ag := len*(ag/(len-mis))]
dmerge[,c_ct := len*(ct/(len-mis))]
dmerge[,c_ac := len*(ac/(len-mis))]
dmerge[,c_at := len*(at/(len-mis))]
dmerge[,c_gc := len*(gc/(len-mis))]
dmerge[,c_gt := len*(gt/(len-mis))]

## sort
dmerge <- dmerge[order(dmerge$GROUP, dmerge$ID, dmerge$loc),]

## handle completely missing
dmerge[mis == len,c("trans","transv","ag","ct","ac","at","gc","gt")] <- NA

## Plot
pdf(file = oname, height = height, width = width)

nygr <- nlevels(dmerge$GROUP)
ycols <- rep(ycols, ceiling(nygr/length(ycols)))[1:nygr]
dpl <- data.frame(dmerge)

# Completeness
suppressWarnings(do.heatmap(dat = dpl, xfac = "loc", yfac = "ID", znum = "completeness",  ygr = "LABEL", ycols = ycols, limit = 1,  title = "Completeness"))

# Transitions
suppressWarnings(do.heatmap(dat = dpl, xfac = "loc", yfac = "ID", znum = "c_trans",  ygr = "LABEL", ycols = ycols, limit = quantile(dpl$c_trans, highquant, na.rm = TRUE),  title = "Transititons"))

suppressWarnings(do.heatmap(dat = dpl, xfac = "loc", yfac = "ID", znum = "c_ag",     ygr = "LABEL", ycols = ycols, limit = quantile(dpl$c_ag, highquant, na.rm = TRUE),     title = "A>G / T>C"))

# common transition in post-mortem damage
suppressWarnings(do.heatmap(dat = dpl, xfac = "loc", yfac = "ID", znum = "c_ct",     ygr = "LABEL", ycols = ycols, limit = quantile(dpl$c_ct, highquant, na.rm = TRUE),     title = "C>T / G>A"))

# Transversions
suppressWarnings(do.heatmap(dat = dpl, xfac = "loc", yfac = "ID", znum = "c_transv", ygr = "LABEL", ycols = ycols, limit = quantile(dpl$c_transv, highquant, na.rm = TRUE), title = "Transversions"))

suppressWarnings(do.heatmap(dat = dpl, xfac = "loc", yfac = "ID", znum = "c_ac",     ygr = "LABEL", ycols = ycols, limit = quantile(dpl$c_ac, highquant, na.rm = TRUE),     title = "A>C / T>G"))
suppressWarnings(do.heatmap(dat = dpl, xfac = "loc", yfac = "ID", znum = "c_at",     ygr = "LABEL", ycols = ycols, limit = quantile(dpl$c_at, highquant, na.rm = TRUE),     title = "A>T / T>A"))
suppressWarnings(do.heatmap(dat = dpl, xfac = "loc", yfac = "ID", znum = "c_gc",     ygr = "LABEL", ycols = ycols, limit = quantile(dpl$c_gc, highquant, na.rm = TRUE),     title = "G>C / C>G"))
suppressWarnings(do.heatmap(dat = dpl, xfac = "loc", yfac = "ID", znum = "c_gt",     ygr = "LABEL", ycols = ycols, limit = quantile(dpl$c_gt, highquant, na.rm = TRUE),     title = "G>T / C>A"))

graphics.off()


