#!/cluster/apps/r/3.6.0_openblas/x86_64/bin/Rscript

#########################
## PLOT CONTIG LENGTHS ##
#########################

## Usage
# plot.contig.lengths.R <lengthtab> <OPT: mapfile>

## Author
# simon.crameri@env.ethz.ch, Apr 2019

## Load required library
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(data.table))

## Get arguments
args = commandArgs(trailingOnly = TRUE)

## Check arguments
options(warning.length = 2000L)
if (! length(args) %in% c(1:2)) {
  stop("2 arguments taken (1 required): 
       REQUIRED
       1) <lentab|CHR>: path to table with contig ID (1st column) and contigs lengths (2nd to nth columns, 1 column for each locus). Expects no header and tab-delimiter ; 
       OPTIONAL
       2) <meta|CHR>: path to metadata file mapping taxon ID (1st column) to metadata (2nd column). Expects header and tab-delimiter",
       call.=FALSE)
}

## Set arguments
lentab <- as.character(args[1])
meta <- as.character(args[2])

## Set arguments (for debugging)
# lentab <- "multifasta.11.2605.locus.lengths"
# meta <- "mapfile.orig.cpgroup.txt"

## Additional arguments
# clusterin
hclustmethod <- "ward.D2" # locus clustering method in heatmap
nulllength <- 20 # length of a null sequence, e.g. 20 for (--------------------)
sortx <- TRUE                   # if TRUE, will sort regions based on <hclustmethod> in heatmap plots

# plotting
ycols <- c("#A6CEE3","#99CD91","#B89B74","#F06C45","#ED8F47","#825D99","#B15928") # group colors, will be recycled
low <- "#A50026"                # color used for the low end of the heatmap gradient
mid <- "#F0E442"                # color used for the mid point of the heatmap gradient
high <- "#081D58"               # color used for the high end of the heatmap gradient
qlow <- 0.1                     # lower quantile (shown on violin plots)
qhigh <- 0.9                    # upper quantile (shown on violin plots and used to cap heatmaps of ALL and PASSED regions)
plot.width <- 15 # output plot width
plot.height <- 15 # output plot height

# paths
# suffix1 <- ".locus"
suffix2 <- ".lengths"
# rmdirname <- "loci_rm_allzerolength"
ofile <- paste0(lentab, ".pdf") # name of output file


## Check arguments
stopifnot(file.exists(lentab))
if (!is.na(meta)) stopifnot(file.exists(meta))

## Define helperfunctions
# create a ggplot2 heatmap (y sorted for groups, z sorted for clusters)
do.heatmap <- function(dat, xfac, yfac, znum, ygr, ycols = NULL, limit = 100, title = znum,
                       sortx = TRUE, hclustmethod = "ward.D2", low = "#A50026", mid = "#F0E442", high = "#081D58") {
  
  library(ggplot2)
  
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
  low.limit <- min(dat[,znum], na.rm = TRUE) # specifies the lower variable limit (corresponds to lower capping of heatmap color gradient)
  
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

######################################################################################################

## Read lenfile
dlen <- fread(lentab, stringsAsFactors = TRUE)
names(dlen)[1] <- "ID"

## Melt lenfile if there are more than 10 columns (assuming this means a wide data format)
if (ncol(dlen) > 10) {
  dlen <- melt(dlen, id.vars = "ID", variable.name = "target", value.name = "length")
}

## Read metadata
if (!is.na(meta)) {
  dmeta <- read.delim(meta, stringsAsFactors = FALSE)
  names(dmeta)[1] <- "ID"
  names(dmeta)[2] <- "GROUP"
}

## Merge it
if (!is.na(meta)) {
  dmerge <- merge(dlen, dmeta, by = "ID", all.x = TRUE, all.y = FALSE, sort = FALSE)
  dmerge[is.na(GROUP),"GROUP"] <- "NA"
  dmerge$GROUP <- factor(dmerge$GROUP)
  dmerge$ID <- as.factor(dmerge$ID)
} else {
  dmerge <- dlen
  dmerge$GROUP <- factor("NA")
}

## Add number of individuals per group (all samples, all groups)
grtab <- table(dmerge[!duplicated(dmerge[,"ID"]),"GROUP"])
dmerge[,"LABEL"] <- dmerge[,"GROUP"]
for (group in levels(unlist(dmerge[,"GROUP"]))) {
  dmerge[GROUP == group, "LABEL"] <- paste0(names(grtab[group]), " (n = ", grtab[group], ")")
}
dmerge[,"LABEL"] <- droplevels(dmerge[,"LABEL"])


## Order by group
dmerge <- dmerge[order(dmerge$GROUP, dmerge$ID),]


## Get color vector
nreg <- length(unique(dmerge$target))
nygr <- nlevels(unlist(dmerge[,"LABEL"]))
ycols <- rev(rep(ycols, ceiling(nygr/length(ycols)))[1:nygr])


## Prepare Plots
# aap lengths
pq <- round(quantile(dmerge[,"length"], c(qlow, qhigh), na.rm = TRUE))


# add null lengths for heatmap
dmerge[length == nulllength,"length"] <- as.integer(NA_integer_) # as.integer(0)


# add group variable with levels in reversed order
dmerge[,"LABEL2"] <- factor(unlist(dmerge[,"LABEL"]), levels = rev(levels(unlist(dmerge[,"LABEL"])))) # reverse label order for top-to-bottom y axis


# plot length distributions per group
p1 <- ggplot(dmerge, aes_string(x = "LABEL2", y = "length", fill = "LABEL2")) +
  geom_boxplot(alpha = 0.5, na.rm = TRUE) +
  geom_violin(alpha = 0.5, na.rm = TRUE) +
  labs(x = "", y = "Contig length") +
  geom_hline(yintercept = c(pq[1], pq[2]), 
             linetype = c(2, 2), col = c(1, 1)) +
  scale_fill_manual(guide = FALSE, values = ycols) +
  # scale_y_log10() +
  coord_flip() +
  theme_bw() +
  ggtitle(paste0("Number of ALL regions: ", nreg))

p2 <- ggplot(dmerge, aes_string(x = "LABEL2", y = "length", fill = "LABEL2")) +
  geom_boxplot(alpha = 0.5, na.rm = TRUE) +
  geom_violin(alpha = 0.5, na.rm = TRUE) +
  labs(x = "", y = "Contig length") +
  geom_hline(yintercept = c(pq[1], pq[2]), 
             linetype = c(2, 2), col = c(1, 1)) +
  scale_fill_manual(guide = FALSE, values = ycols) +
  scale_y_log10() +
  coord_flip() +
  theme_bw() +
  ggtitle(paste0("Number of ALL regions: ", nreg))


## Plot
pdf(file = ofile, width = plot.width, height = plot.height)

print(p1)
print(p2)

suppressWarnings(do.heatmap(dat = data.frame(dmerge), xfac = "target", yfac = "ID", znum = "length", ygr = "LABEL",
                            ycols = rev(ycols), limit = pq[2], sortx = sortx, hclustmethod = hclustmethod, 
                            low = low, mid = mid, high = high, title = paste("ALL regions,", "capped contig length")))

graphics.off()


## Look for all-zero-length loci or individuals
# loci
lzero <- which(tapply(dmerge$length, INDEX = dmerge$target, FUN = function(x) all(x == 0)))
cat("found", length(lzero), "loci with zero length in all individuals!\n")

# individuals
izero <- which(tapply(dmerge$length, INDEX = dmerge$ID, FUN = function(x) all(x == 0)))
cat("found", length(izero), "individuals with zero length in all loci!\n")

if (length(izero) > 0) {
  zeroinds <- names(izero)
  print(zeroinds, max = 10, quote = FALSE)
  writeLines(zeroinds, con = paste0(gsub(paste0(suffix2, "$"), "", lentab), ".zeroinds"))
}

