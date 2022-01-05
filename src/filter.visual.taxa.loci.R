#!/cluster/apps/r/3.6.0_openblas/x86_64/bin/Rscript

############################################################################
## Filter taxa and loci based on visualization of (aggregated) statistics ##
############################################################################

## Usage
# run script without arguments to see the arguments it takes

## Author
# simon.crameri@env.ethz.ch, May 2019

## Load required library
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gplots))
# suppressPackageStartupMessages(library(grid))
# suppressPackageStartupMessages(library(VennDiagram))

## Get arguments
args = commandArgs(trailingOnly = TRUE)

## Check arguments
options(warning.length = 2000L)
if (! length(args) %in% c(2:11)) {
  stop("At least 2 arguments needed (11 taken):
       REQUIRED
       1) <file|CHR>: path to collected contig stats (loci_stats.txt)
       2) <refseqs|CHR> path to target locus reference sequences (.fasta)
       
       OPTIONAL
       # taxon filtering criteria
       3) <min.ploci|NUM>: minimum proportion of target loci covered by a taxon (i.e., taxon has at least 1 successfully aligned contig in <min.ploci>*nloci target loci) [DEFAULT: 0.3]
       
       # locus filtering criteria
       4) <min.ptaxa|NUM>: minimum median proportion of taxa covered by a locus (i.e., locus has at least 1 successfully aligned contig in <min.ptaxa>*ntaxa taxa) (per group) [DEFAULT: 0.3]
       5) <max.ncontigs|NUM>: maximum median number of (non-zero) contigs per target locus (per group) [DEFAULT: 1]
       6) <min.bestscore.norm|NUM>: minimum median normalized alignment score (i.e., raw score / target alignment length) of best alignment (per group) [DEFAULT: 1]
       7) <min.taln|NUM>: minimum median target alignment length (per group) [DEFAULT: 50]
       8) <min.tfrac|NUM>: minimum median target alignment fraction (per group) [DEFAULT: 0.2]
       9) <min.bestscore|NUM>: minimum median alignment score of best contig (per group) [DEFAULT: 1]
       10) <min.bestlength|NUM>: minimum median length of best contig (per group) [DEFAULT: 1]
       
       # grouping for locus filtering
       11) <meta|CHR>: path to metadata file mapping taxon labels (1st column) to groups (2nd column). Groups labelled as 'NA' will be displayed in plots but will not be considerd during locus filtering. Header and tab separation expected. More taxa of different order than in <file> are ok, missing taxa will be labeled as 'NA'. [DEFAULT: NA]
       ", call.=FALSE)
}

## Set arguments
file <- as.character(args[1])
refseqs <- as.character(args[2])

min.ploci <- as.numeric(as.character(args[3]))
if (is.na(min.ploci)) min.ploci <- 0.3
min.ptaxa <- as.numeric(as.character(args[4]))
if (is.na(min.ptaxa)) min.ptaxa <- 0.3
max.ncontigs <- as.numeric(as.character(args[5]))
if (is.na(max.ncontigs)) max.ncontigs <- 1
min.bestscore.norm <- as.numeric(as.character(args[6]))
if (is.na(min.bestscore.norm)) min.bestscore.norm <- 1

min.taln <- as.numeric(as.character(args[7]))
if (is.na(min.taln)) min.taln <- 50
min.tfrac <- as.numeric(as.character(args[8]))
if (is.na(min.tfrac)) min.tfrac <- 0.2
min.bestscore <- as.numeric(as.character(args[9]))
if (is.na(min.bestscore)) min.bestscore <- 1
min.bestlength <- as.numeric(as.character(args[10]))
if (is.na(min.bestlength)) min.bestlength <- 1

meta <- as.character(args[11])

## Sample starting time
t1 <- Sys.time()
t1str <- strsplit(as.character(t1), split = " ")
t1str <- paste0(".", sapply(t1str, "[", 1), "_", gsub(":", "-", sapply(t1str, "[", 2)))
cat("### Filter taxa and loci based on visualization of (aggregated) statistics ###\n")
print(t1)

## Set arguments (for debugging)
# file = "loci_stats.txt"
# # ranges <- "loci_ranges.txt"
# refseqs = "test.fasta"
# min.ploci <- 0.2 # minimum proportion of loci with non-zero contigs (filters taxa)
# min.ptaxa <- 0.01 # minimum proportion of taxa with non-zero contigs (filters loci)
# max.ncontigs <- 2 # maximum number of non-zero contigs (filters loci)
# min.bestscore.norm <- 2 # minimum best normalized alignment score (filters loci)
# min.taln <- 80 # minimum target alignment length (filters loci)
# min.tfrac <- 0 # minimum target alignment length / target length (filters loci)
# min.bestscore <- 1 # minimum best alignment score (filters loci)
# min.bestlength <- 1 # minimum best contig length (filters loci)
# meta <- "mapfile.txt"

## Additional arguments
# group summary statistic for locus filtering
sumfun <- median
funlab <- "Median" # function label in plots (needs to match <fun>)

# expected field names in <file> after merging with <refseqs> and <meta>
fnames <- c("locus","GROUP","LABEL","ID","ncontigs","ncontigscat","scores","lengths","bestscore","bestscore.norm","bestlength","tlen",
            "tgc","whichbest","query","qstart","qend","tstart","tend","score",
            "qaln","taln","qfrac","tfrac","tqfrac")

# definition of filtering criteria (filtervars) and extreme value filtering (filterexts)
filtervars <- c("ptaxa", "ncontigs", "bestscore.norm", "taln", "tfrac", "bestscore", "bestlength")
filterexts <- c("min",   "max",      "min",            "min",  "min",   "min",       "min")

# visualization
XandMore = 3 # will create number of contigs categories for 0, 1, 2, ..., >=XandMore contigs
sortmat = TRUE # if TRUE, will sort loci based on <hclustmethod> in heatmap plots
hclustmethod = "ward.D2" # method to cluster loci on heatmaps
plot.width = 15 # plot width in output .pdf
plot.height = 15 # plot height in output .pdf

# output files
obase <- paste(min.ploci, min.ptaxa, max.ncontigs, min.bestscore.norm, min.taln, min.tfrac, min.bestscore, min.bestlength, sep = "-")
obase <- paste(min.ploci, min.ptaxa, max.ncontigs, sep = "-")
outpdf = paste0("loci_stats-", obase, ".pdf") # output .pdf (visualization of (aggregated) stats)
outtxt = paste0("loci_stats-", obase, ".txt") # output .txt (merged stats)
aggtax = paste0("taxa_stats-", min.ploci, ".txt") # output .txt (aggregation by taxon)
aggtxt = paste0("loci_aggr-", obase, ".txt") # output .txt (aggregation by group)
loctxt = paste0("loci_kept-", obase, ".txt") # filtered loci
taxtxt = paste0("taxa_kept-", min.ploci, ".txt") # filtered taxa
logtxt = paste0("loci_stats-", obase, ".log") # output .log 

## Check arguments
stopifnot(file.exists(file),
          file.exists(refseqs),
          min.ploci >= 0, min.ploci <= 1,
          min.ptaxa >= 0, min.ptaxa <= 1,
          max.ncontigs >= 1,
          min.bestscore.norm >= 0,
          min.taln >= 1,
          min.tfrac >= 0, min.tfrac <= 1,
          min.bestscore >= 1,
          min.bestlength >= 1)
if (!is.na(meta)) stopifnot(file.exists(meta))

## Move any existing output files
if (file.exists(outpdf)) {system(paste("mv", outpdf, gsub(".pdf$", paste0(t1str, ".pdf"), outpdf)))}
if (file.exists(outtxt)) {system(paste("mv", outtxt, gsub(".txt$", paste0(t1str, ".txt"), outtxt)))}
if (file.exists(aggtax)) {system(paste("mv", aggtax, gsub(".txt$", paste0(t1str, ".txt"), aggtax)))}
if (file.exists(aggtxt)) {system(paste("mv", aggtxt, gsub(".txt$", paste0(t1str, ".txt"), aggtxt)))}
if (file.exists(loctxt)) {system(paste("mv", loctxt, gsub(".txt$", paste0(t1str, ".txt"), loctxt)))}
if (file.exists(taxtxt)) {system(paste("mv", taxtxt, gsub(".txt$", paste0(t1str, ".txt"), taxtxt)))}
if (file.exists(logtxt)) {system(paste("mv", logtxt, gsub(".log$", paste0(t1str, ".log"), logtxt)))}

## Define helperfunctions
# aggregate data
get.agg <- function(df, sumfun) {
  # check input
  stopifnot(all(fnames %in% names(df)))
  
  # proportion of loci covered by taxa
  dptaxa <- aggregate(df$ncontigs, by = list(df$LABEL, df$locus), FUN = length, drop = FALSE)
  dptaxa$x[is.na(dptaxa$x)] <- 0 # all zero contigs in a group
  names(dptaxa) <- c("LABEL", "locus", "nloci")
  df.map <- df[!duplicated(df$ID),c("ID","LABEL")]
  df.agg <- aggregate(df.map$ID, by = list(df.map$LABEL), FUN = length, drop = FALSE)
  names(df.agg) <- c("LABEL", "ntaxa")
  dptaxa <- merge(dptaxa, df.agg, by = "LABEL", sort = FALSE)
  dptaxa$ptaxa <- dptaxa$nloci / dptaxa$ntaxa
  dptaxa <- dptaxa[order(dptaxa$locus, dptaxa$LABEL),]
  rownames(dptaxa) <- seq(nrow(dptaxa))
  
  # number of distinct contigs with alignment to target locus
  dncontigs <- aggregate(df$ncontigs, by = list(df$LABEL, df$locus), FUN = function(x) {do.call(sumfun, args = list(x, na.rm = TRUE))}, drop = FALSE)
  
  # raw and normalized alignment score of best contig
  dbestscore <- aggregate(df$bestscore, by = list(df$LABEL, df$locus), FUN = function(x) {do.call(sumfun, args = list(x, na.rm = TRUE))}, drop = FALSE)
  dbestscore.norm <- aggregate(df$bestscore.norm, by = list(df$LABEL, df$locus), FUN = function(x) {do.call(sumfun, args = list(x, na.rm = TRUE))}, drop = FALSE)
  
  # best contig length
  dbestlength <- aggregate(df$bestlength, by = list(df$LABEL, df$locus), FUN = function(x) {do.call(sumfun, args = list(x, na.rm = TRUE))}, drop = FALSE)
  
  # alignment length in query (contig)
  dqaln <- aggregate(df$qaln, by = list(df$LABEL, df$locus), FUN = function(x) {do.call(sumfun, args = list(x, na.rm = TRUE))}, drop = FALSE)
  
  # alignment length in target (reference sequence)
  dtaln <- aggregate(df$taln, by = list(df$LABEL, df$locus), FUN = function(x) {do.call(sumfun, args = list(x, na.rm = TRUE))}, drop = FALSE)
  
  # alignment proportion in query (contig) (query alignment length / query length [bestlength])
  dqfrac <- aggregate(df$qfrac, by = list(df$LABEL, df$locus), FUN = function(x) {do.call(sumfun, args = list(x, na.rm = TRUE))}, drop = FALSE)
  
  # alignment proportion in target (reference sequence) (target alignment length / target length [reference sequence length])
  dtfrac <- aggregate(df$tfrac, by = list(df$LABEL, df$locus), FUN = function(x) {do.call(sumfun, args = list(x, na.rm = TRUE))}, drop = FALSE)
  
  # alignment length in query (contig) relative to target length (query alignment length / target length [reference sequence length])
  dtqfrac <- aggregate(df$tqfrac, by = list(df$LABEL, df$locus), FUN = function(x) {do.call(sumfun, args = list(x, na.rm = TRUE))}, drop = FALSE)
  
  # reference sequence GC content <tgc> (reference-sequence-based stat)
  dtgc <- aggregate(df$tgc, by = list(df$LABEL, df$locus), FUN = function(x) {do.call(sumfun, args = list(x, na.rm = TRUE))}, drop = FALSE)
  
  # reference sequence length <tlen> (reference-sequence-based stat)
  dtlen <- aggregate(df$tlen, by = list(df$LABEL, df$locus), FUN = function(x) {do.call(sumfun, args = list(x, na.rm = TRUE))}, drop = FALSE)
  
  stopifnot(all.equal(dptaxa$LABEL, dncontigs$Group.1),
            all.equal(dptaxa$locus, dncontigs$Group.2),
            all.equal(dptaxa$LABEL, dtlen$Group.1),
            all.equal(dptaxa$locus, dtlen$Group.2))
  
  ## Combine
  dagg <- data.frame(LABEL = dptaxa$LABEL,
                     ntaxa = dptaxa$ntaxa,
                     locus = dptaxa$locus,
                     ptaxa = dptaxa$ptaxa,
                     ncontigs = dncontigs$x,
                     ncontigscat = round(dncontigs$x),
                     bestscore = dbestscore$x,
                     bestscore.norm = dbestscore.norm$x,
                     bestlength = dbestlength$x,
                     qaln = dqaln$x,
                     taln = dtaln$x,
                     qfrac = dqfrac$x,
                     tfrac = dtfrac$x,
                     tqfrac = dtqfrac$x,
                     tgc = dtgc$x,
                     tlen = dtlen$x)
  dagg$LABEL <- factor(dagg$LABEL, levels = levels(df$LABEL))
  dagg$locus <- factor(dagg$locus, levels = levels(df$locus))
  dagg$ncontigscat[dagg$ncontigs >= XandMore ] <- paste0(">=", XandMore)
  dagg$ncontigscat <- factor(dagg$ncontigscat, levels = XandMorelevs)
  
  ## Return
  return(dagg)
}

# plot aggregated data
plot.agg <- function(df, x = "LABEL", y, fill = "LABEL", ylab = "", ybreaks = NULL, nloci.filt = NULL, sloci.filt = NULL, thr = NULL, alpha = 0.5) {
  if (!is.null(nloci.filt)) {
    stopifnot(!is.null(thr))
    comp <- ifelse(filterexts[filtervars == y] == "max", "above", "below")
    if (comp == "above") hjust <- "bottom" else hjust <- "top"
  }
  p <- ggplot(na.omit(df[,c(x, y, fill)]), aes_string(x = x, y = y, fill = fill)) +
    geom_boxplot(alpha = alpha) +
    geom_violin(alpha = alpha) +
    scale_fill_manual(guide = FALSE, values = grcolsf) +
    labs(x = "", y = paste(funlab, ylab)) +
    coord_flip(clip = "off") +
    theme_bw()
  if (!all(is.na(df[,y]))) {
    if (!is.null(ybreaks)) p <- p + scale_y_continuous(breaks = ybreaks)
    if (!is.null(nloci.filt)) {
      p <- p + geom_hline(aes(yintercept = thr), col = "tomato") +
        annotate("text", x = (1:ngr.passed) + 0.2, y = thr, 
                 label = paste0(nloci.filt, " (", round(100*nloci.filt/nloci,2), "%)"), 
                 hjust = hjust, vjust = "left") +
        ggtitle(paste0("Number of target loci: ", nloci, " ; ", comp, " ", thr, ": ", 
                       sloci.filt, " (", round(100*sloci.filt/nloci,2), "%)"))
    } else {
      p <- p + ggtitle(paste0("Number of target loci: ", nloci))
    }
  } else {
    p <- p + 
      ggtitle(paste0("Number of target loci: ", nloci, " ; all ", y, " are NA"))
  }
  invisible(p)
}

# plot pairwise data
plot.comp <- function(df, x, y, xlab = x, ylab = y, shape = "lpassed", linetype = "lpassed", color = "LABEL", alpha = 0.5, smooth = TRUE, boxes = FALSE) {
  p <- ggplot(na.omit(df[,c(x,y,shape,linetype,color)]), aes_string(x = x, y = y, shape = shape, linetype = linetype, color = color)) +
    ggtitle(paste0("Number of target loci: ", nloci)) +
    labs(x = xlab, y = ylab) +
    scale_color_manual(guide = FALSE, values = grcolsf) +
    scale_shape_manual(values = c(5,19)) +
    scale_linetype_manual(values = c(9,1)) +
    theme_bw()
  if (!all(is.na(df[,y]))) {
    p <- p + 
      facet_wrap(as.formula(paste("~", color))) +
      ggtitle(paste0("Number of target loci: ", nloci, " ; passed: ", lpassed))
    if (smooth) {
      p <- p + 
        geom_point(alpha = alpha) +
        geom_smooth(method = "lm", se = FALSE)
      # geom_smooth(method = "gam", se = FALSE, formula = y ~ s(x, bs = "cs"))
    }
    if (boxes) {
      p <- p +
        geom_boxplot(alpha = alpha) +
        geom_violin(alpha = alpha) 
    }
  } else {
    p <- p + 
      ggtitle(paste0("Number of target loci: ", nloci, " ; passed: ", lpassed, " ; all ", y, " are NA"))
  }
  invisible(p)
}

# create a matrix for heatmaps
do.mat <- function(df, var, NA.replacement = 0, fac = "locus", ids = "ID") {
  mat <- matrix(as.numeric(df[,var]), ncol = nlevels(df[,fac]))
  rownames(mat) <- unique(df[,ids])
  mat[is.na(mat)] <- NA.replacement
  return(mat)
}

# extrafun to heatmap.2
myfun <- function() {
  plot(0, 0, type = "n", xlab = "", ylab = "", axes = FALSE, bty = "n")
  legend("topleft", legend = levels(dmerge$GROUP), 
         pch = 15, 
         pt.cex = 3, 
         cex = 1, 
         col = grcols, 
         ncol = ceiling(nlevels(dmerge$GROUP)/2),
         bty = "n", xpd = TRUE)
}
myfunf <- function() {
  plot(0, 0, type = "n", xlab = "", ylab = "", axes = FALSE, bty = "n")
  legend("topleft", legend = levels(dmergesub.nonzero$GROUP), 
         pch = 15, 
         pt.cex = 3, 
         cex = 1, 
         col = grcolsf, 
         ncol = ceiling(nlevels(dmergesub.nonzero$GROUP)/2),
         bty = "n", xpd = TRUE)
}

# wrapper around heatmap.2
do.heatmap <- function(mat, max.nbreaks = 200, labRow, Colv = TRUE, reverse.colors = FALSE, highquant = 1, lmat, lwid, lhei, hclustmethod, RowSideColors, extrafun, key.title) {
  
  stopifnot(is.matrix(mat), highquant > 0)
  
  # set heatmap legend color breaks
  states <- unique(as.numeric(mat))
  nstates <- length(states)
  matmax <- max(states, na.rm = TRUE)
  matmin <- min(states, na.rm = TRUE)
  # if (matmax <= 1) start <- 0.001 else start <- 0.5
  if (matmax <= 1) start <- 0.001 else start <- matmin + 1/nstates
  end <- matmax #- start
  
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
  if (length(breaks) < length(col) + 1) breaks <- c(0, seq(start, matmax, length.out = length(col)))
  if (any(duplicated(breaks))) breaks <- c(0, seq(start, matmax, length.out = length(col)))
  
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
  if (!sum(states) == 0) {
    p <- try(heatmap.2(x = mat,
                       Rowv = FALSE, # no reordering of rows (individuals)
                       Colv = Colv, # would reorder columns (loci) according to decreasing fun value (but separately for each cluster)
                       labRow = labRow, labCol = NA,
                       add.expr = box(),
                       #hclustfun = function(x) hclust(x, method = hclustmethod), # here, x is going to be distfun(mat)
                       #distfun = dist,                                           # 
                       scale = "none",
                       margins = c(3,0), # bottom margin will be 3
                       trace = "none", 
                       symkey = FALSE, 
                       symbreaks = FALSE, 
                       dendrogram = "none",
                       density.info = "histogram",
                       denscol = "black", # line color of histogram
                       col = col,
                       breaks = breaks,
                       key.title = key.title,
                       key.xlab = "",
                       key.ylab = "",
                       key.xtickfun = key.xtickfun,
                       keysize = 1,
                       extrafun = extrafun,
                       RowSideColors = RowSideColors,
                       key.par = list(mar = c(2,0,2,0)), # margins around key (bottom, left, top, right)
                       lmat = lmat, lhei = lhei, lwid = lwid # to understand this, see show.layout(layout(lmat,lhei,lwid))
    ), silent = TRUE)
    if (inherits(p, "try-error")) {
      msg <- paste0("ERROR: heatmap (or part of it) of ", substitute(mat), " could not be plotted.\nCheck figure margins and/or conformity of key breaks (", length(breaks), ") and colors (", length(col), ").\n", p[1])
      plot.new()
      mtext(msg)
      cat(msg)
    }
  } else {
    msg <- paste0("heatmap of ", substitute(mat), " not plotted (all values are zero)\n")
    plot.new()
    mtext(msg)
  }
}

# hexplot
hexplot <- function(df, nbins = 50, x = "x", y = "y", z = "z", 
                    xlab = x, ylab = y, zlab = z,
                    legend.position = "top", plot.margin = unit(c(0.5, 1, 0.5, 0.5), "cm"),
                    FUN = length, print = TRUE) {
  p <- ggplot(df, aes_string(x = x, y = y, z = z)) +
    stat_summary_hex(bins = nbins, fun = function(x) do.call(FUN, list(x))) +
    xlab(xlab) +
    ylab(ylab) +
    scale_fill_continuous(name = zlab) +
    theme_bw() +
    theme(legend.position = legend.position, plot.margin = plot.margin)
  if (print) print(p)
  invisible(p)
}

# ggplot2 colors
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

###################################################################################

## Read contig stats
dc <- read.delim(file, check.names = FALSE)
names(dc)[1] <- "ID"

## Read reference fasta
ref <- read.FASTA(refseqs)
dgc <- numeric()
for (i in seq(length(ref))) {
  dgc <- c(dgc, GC.content(ref[i]))
}
dref <- data.frame("locus" = names(ref), "tlen" = lengths(ref), "tgc" = dgc)

## Read metadata
if (!is.na(meta)) {
  dmeta <- read.delim(meta)
  names(dmeta)[1] <- "ID"
  names(dmeta)[2] <- "GROUP"
}

## Merge data
cat("\nmerging data...\n")

# merge with taxon metadata
if (!is.na(meta)) {
  dmerge <- merge(dc, dmeta, by = "ID", all.x = TRUE, all.y = FALSE, sort = FALSE)
  dmerge$GROUP <- as.character(dmerge$GROUP)
} else {
  dmerge <- dc
  dmerge$GROUP <- "NA"
}
# merge with reference sequence data
if (!all(unique(as.character(dmerge$locus)) %in% as.character(dref$locus))) {
  stop("locus names in ", file, " do not match locus names in ", refseqs)
}
dmerge <- merge(dmerge, dref, by = "locus", all.x = TRUE, all.y = FALSE, sort = FALSE)

## Handle NA values as an own category (any missing taxon in meta or any taxon labelled as "NA")
dmerge$GROUP[is.na(dmerge$GROUP)] <- "NA"
dmerge$GROUP <- factor(dmerge$GROUP)

## Get number of individuals and groups, number of target loci
ntaxa <- nlevels(dmerge$ID) # all taxa
ngroups <- nlevels(dmerge$GROUP) # all mapped groups
nloci <- nlevels(dc$locus) # all target loci
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

## Add number of individuals per group
dmerge$LABEL <- as.character(dmerge$GROUP)
for (group in unique(dmerge$GROUP)) {
  dmerge$LABEL[dmerge$LABEL == group] <- paste0(names(grtab[group]), " (n = ", grtab[group], ")")
}
dmerge$LABEL <- as.factor(dmerge$LABEL)

## Add number of contigs category (ncontigscat)
XandMorelevs <- as.character(0:XandMore)
XandMorelevs[length(XandMorelevs)] <- paste0(">=", XandMore)
dmerge$ncontigscat <- as.character(dmerge$ncontigs)
dmerge$ncontigscat[dmerge$ncontigs >= XandMore ] <- paste0(">=", XandMore)
dmerge$ncontigscat <- factor(dmerge$ncontigscat, levels = XandMorelevs)

## Get alignment lengths and fractions (if available)
dmerge$qaln <- abs(dmerge$qend - dmerge$qstart) # query alignment length
dmerge$taln <- abs(dmerge$tend - dmerge$tstart) # target alignment length
dmerge$qfrac <- dmerge$qaln / dmerge$bestlength # aligned portion in query
dmerge$tfrac <- dmerge$taln / dmerge$tlen # aligned portion in target
dmerge$tqfrac <- dmerge$qaln / dmerge$tlen # query alignment length relative to target length

## Normalize bestscore for alignment length
# in the case of a Smith-Waterman alignment (affine:local, exhaustive yes), the raw scores are the sum of the substitution matrix scores and the gap penalties
# to normalize it, we divide the raw score by the target alignment length
# in the case of a Smith-Waterman alignment (affine:local, exhaustive yes), the upper limit of the normalized score is usually 5 (follows from the Smith-Waterman scoring algorithm)
#dmerge$bestscore.norm <- dmerge$bestscore / dmerge$length # if normalized for target length, normalized scores are low if the alignment length is small
dmerge$bestscore.norm <- dmerge$bestscore / dmerge$taln

## Correct ncontigs if alignment ranges are non-overlapping
## wait for contig-ranges.txt

## Order rows (prepares to get a heatmap matrix (wide format) of any variable (long format) and columns (first fnames, then rest)
dmerge <- dmerge[order(dmerge$locus, dmerge$GROUP, dmerge$ID), c(fnames, names(dmerge)[!names(dmerge) %in% fnames])]
rownames(dmerge) <- NULL

## Filter taxa based on min.ploci (minimum proportion of loci with data)
min.nloci <-  min.ploci * nloci # minimum required number of covered loci
dmerge.nonzero <- subset(dmerge, ncontigs > 0)

dploci <- aggregate(dmerge.nonzero$ncontigs, by = list(dmerge.nonzero$ID), FUN = length, drop = FALSE) # keeps all taxa in aggregation
names(dploci) <- c("ID", "covered_loci")
dploci$covered_loci[is.na(dploci$covered_loci)] <- 0 # all zero contigs in a taxon
dploci$required_loci <- min.nloci
dploci$total_loci <- nloci
dploci$ploci <- dploci$covered_loci / dploci$total_loci
dploci$tpassed <- factor(ifelse(dploci$ploci < min.ploci, 0, 1), levels = c("0", "1"))
dploci <- cbind(dploci, dmerge.nonzero[!duplicated(dmerge.nonzero$ID), c("GROUP","LABEL"), drop = F][order(dmerge.nonzero[!duplicated(dmerge.nonzero$ID),"ID"]),])
dploci <- dploci[order(dploci$LABEL, dploci$ID),c("ID", "GROUP", "LABEL", "covered_loci", "required_loci", "total_loci", "ploci", "tpassed")]
rownames(dploci) <- NULL

## Get number of passed (min.ploci) and considered (non-NA group) taxa and groups
# passed (min.ploci)
taxa.passed <- as.character(dploci$ID[dploci$tpassed == "1"]) # names of taxa passing min.ploci
tpassed <- length(taxa.passed) # number of taxa passing min.ploci
grtab.passed <-  table(dploci$GROUP, dploci$tpassed) # table of passed and filtered GROUPs
labtab.passed <- table(dploci$LABEL, dploci$tpassed) # table of passed and filtered LABELs
ngr.passed <- sum(grtab.passed[,"1", drop = FALSE] > 0) # number of GROUPs/LABELs passing min.ploci
gr.passed <- grtab.passed[grtab.passed[,"1", drop = FALSE] > 0, "1", drop = TRUE] # table of GROUPs passing min.ploci
lab.passed <- labtab.passed[labtab.passed[,"1", drop = FALSE] > 0, "1", drop = TRUE] # table of LABELs passing min.ploci
ntaxa.ploci <- labtab.passed[,"0", drop = TRUE] # table of LABELs not passing min.ploci

# considered (passed and non-NA group)
if (length(gr.passed) == 1) {
  # if there is 1 passed group, R will turn names(gr.passed) to rownames(gr.passed), which causes problems for the definition of (n)gr.cons and staxa.cons.
  # 1 passed group is either <NA> (no meta specified) or <NA> (only passed group out of >= 1 groups) or a single, non-NA group
  # in all these cases, that group is considered for locus filtering
  gr.passed <- grtab.passed[grtab.passed[,"1", drop = FALSE] > 0, "1", drop = FALSE] # table of GROUPs passing min.ploci
  lab.passed <- labtab.passed[labtab.passed[,"1", drop = FALSE] > 0, "1", drop = FALSE] # table of LABELs passing min.ploci
  ntaxa.ploci <- labtab.passed[,"0", drop = TRUE] # table of LABELs not passing min.ploci
  gr.cons <- gr.passed # table of GROUPs considered
  lab.cons <- lab.passed # table of LABELs considered
  names(gr.passed) <- rownames(gr.passed)
  names(lab.passed) <- rownames(lab.passed)
  names(ntaxa.ploci) <- rownames(ntaxa.ploci)
  names(gr.cons) <- rownames(gr.cons)
  names(lab.cons) <- rownames(lab.cons)
} else {
  gr.cons <- gr.passed[!names(gr.passed) == "NA"] # table of GROUPs considered
  lab.cons <- lab.passed[!names(gr.passed) == "NA"] # table of LABELs considered
}
ngr.cons <- length(gr.cons) # number of GROUPs/LABELs considered
staxa.cons <- sum(gr.cons) # number of taxa considered

cat(paste0("\nfound ", tpassed, "/", ntaxa, " taxa (", round(100*tpassed/ntaxa, 2), "%) from ", ngr.passed, "/", ngroups, " groups (", round(100*ngr.passed/ngroups, 2), "%) with data for ", min.nloci, "/", nloci, " (>=", 100*min.ploci, "%) target loci...\n"))

## Aggregate by group (passed taxa)
cat(paste0("\naggregating data from ", tpassed, " taxa from ", ngr.passed, " group(s) using ", tolower(funlab)," values ...\n"))
dmergesub.nonzero <- subset(dmerge.nonzero, ID %in% taxa.passed)
dmergesub.nonzero$ID <- droplevels(dmergesub.nonzero$ID)
dmergesub.nonzero$GROUP <- droplevels(dmergesub.nonzero$GROUP)
dmergesub.nonzero$LABEL <- droplevels(dmergesub.nonzero$LABEL)

# check that at least 1 taxon passed
if (tpassed == 0) {
  dmerge$tpassed <- 0
  write.table(dmerge, file = outtxt, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
  write.table(dploci, file = aggtax, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
  stop(paste0("No taxa passed the min. target locus coverage proportion (", min.ploci, "). Consider to decrease this filter.\n"))
}
dagg <- get.agg(df = dmergesub.nonzero, sumfun = sumfun)

## Update merged data with passed taxa
dmerge$tpassed <- factor(ifelse(dmerge$ID %in% taxa.passed, 1, 0))

## Get number of filtered loci per variable
for (filtervar in filtervars) {
  filterext <- filterexts[which(filtervars == filtervar)]
  filtervalue <- get(paste(filterext, filtervar, sep = "."))
  switch(filterext,
         min = {
           # number of filtered loci per group
           nloci.filtervar <- tapply(dagg[,filtervar], INDEX = dagg[,"LABEL"], FUN = function(x) {sum(x < filtervalue, na.rm = TRUE)})
           
           # ID of filtered loci
           loci.filtervar <- na.omit(unique(dagg[dagg[,filtervar] < filtervalue, "locus"]))
           
           # sum of filtered loci
           sloci.filtervar <- length(loci.filtervar)
         },
         max = {
           nloci.filtervar <- tapply(dagg[,filtervar], INDEX = dagg[,"LABEL"], FUN = function(x) {sum(x > filtervalue, na.rm = TRUE)})
           loci.filtervar <- na.omit(unique(dagg[dagg[,filtervar] > filtervalue, "locus"]))
           sloci.filtervar <- length(loci.filtervar)
         })
  assign(paste("nloci", filtervar, sep = "."), nloci.filtervar) # e.g. nloci.ncontigs
  assign(paste("loci", filtervar, sep = "."), loci.filtervar) # e.g. loci.ncontigs
  assign(paste("sloci", filtervar, sep = "."), sloci.filtervar) # e.g. sloci.ncontigs
  rm(nloci.filtervar, loci.filtervar, sloci.filtervar, filterext, filtervalue)
}

## Get variable matrices (all except ptaxa) for all taxa and loci
for (filtervar in filtervars[-1]) {
  mat.filtervar <- do.mat(df = dmerge, var = filtervar)
  assign(paste("mat", filtervar, sep = "."), mat.filtervar)
  rm(mat.filtervar)
}

## Order variable matrices according to bestscore (if sortmat = FALSE)
if (!sortmat) {
  neworder <- order(apply(mat.bestscore, 2, sumfun, na.rm = TRUE), decreasing = FALSE)
  for (filtervar in filtervars[-1]) {
    mat.filtervar <- get(paste("mat", filtervar, sep = "."))[,neworder]
    assign(paste("mat", filtervar, sep = "."), mat.filtervar)
  }
  rm(mat.filtervar, neworder)
}

## Select loci based on min.ptaxa, max.ncontigs (if available, also based on min.bestscore, min.bestlength, min.bestscore.norm, min.taln, min.tfrac)
cat(paste0("\nfinding target loci passing filters in ", staxa.cons, " taxa from ", ngr.cons, " considered group(s):\n\n"))
print(gr.cons)

if (all(is.na(dagg$taln))) {
  cat("\nWARNING: target (reference) alignment length not available, omitting filters for minimum best normalized alignment score and minimum target alignment length / proportion\n")
  min.bestscore.norm <- NA
  min.taln <- NA
  min.tfrac <- NA
}

if (all(is.na(dagg$bestlength))) {
  cat(paste0("\nWARNING: length of best matching contigs not available in <", file, ">, omitting filter for minimum median length of best contig (", min.bestlength, ")\n"))
}

list.passed <- list()
loci.passed <- character()
for (group in names(lab.cons)) {
  dfilt <- subset(dagg, LABEL == group & ptaxa >= min.ptaxa & ncontigs <= max.ncontigs & bestscore >= min.bestscore)
  if (!all(is.na(dfilt$bestlength))) {
    dfilt <- subset(dfilt, bestlength >= min.bestlength)
  }
  if (!all(is.na(dfilt$bestscore.norm))) {
    dfilt <- subset(dfilt, bestscore.norm >= min.bestscore.norm & taln >= min.taln & tfrac >= min.tfrac)
  }
  list.passed[[group]] <- dfilt$locus
  if (length(loci.passed) == 0) loci.passed <- dfilt$locus else loci.passed <- intersect(loci.passed, dfilt$locus)
}
nloci.passed <- sapply(list.passed, length) # number of passed loci per considered group
names(list.passed) <- paste0(names(list.passed), " (", "taxa passed: ", gr.cons, " [", round(100*gr.cons/grtab[names(gr.cons)],2), "%] ; loci passed: ", nloci.passed, " [", round(100*nloci.passed/nloci,2), "%])")
lpassed <- length(loci.passed)
cat(paste0("\nfound ", lpassed, "/", nloci, " (", round(100*lpassed/nloci,2), "%) target loci passing filters in ", ngr.cons, " considered group(s)\n"))

# check that at least 1 locus passed
if (lpassed == 0) {
  dmerge$lpassed <- factor(0, levels = c("0","1"))
  write.table(dmerge, file = outtxt, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
  write.table(dploci, file = aggtax, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
  write.table(dagg, file = aggtxt, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
  write.table(taxa.passed, file = taxtxt, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
  stop(paste0("No target locus passed the filters in all groups. Consider to decrease one or more of them.\n"))
}

## Update merged data with passed target loci
dmerge$lpassed <- factor(ifelse(dmerge$locus %in% loci.passed, 1, 0), levels = c("0","1"))
dagg$lpassed <- factor(ifelse(dagg$locus %in% loci.passed, 1, 0), levels = c("0","1"))

## Subset merged data
dmerge.passed <- subset(dmerge, tpassed == "1" & lpassed == "1")
dmerge.passed$locus <- droplevels(dmerge.passed$locus)

## Get variable matrices (all except ptaxa) for passed taxa and loci
for (filtervar in filtervars[-1]) {
  mat.filtervar <- do.mat(df = dmerge.passed, var = filtervar)
  assign(paste("matf", filtervar, sep = "."), mat.filtervar)
  rm(mat.filtervar)
}

## Order filtered variable matrices according to bestscore (if sortmat = FALSE)
if (!sortmat) {
  neworder <- order(apply(matf.bestscore, 2, sumfun, na.rm = TRUE), decreasing = FALSE)
  for (filtervar in filtervars[-1]) {
    mat.filtervar <- get(paste("matf", filtervar, sep = "."))[,neworder]
    assign(paste("matf", filtervar, sep = "."), mat.filtervar)
  }
  rm(mat.filtervar, neworder)
}

## Get statistics of filtered loci
# subset dagg and dmerge to get loci that are present in some but absent in other groups, etc.

###########
## Plots ##
###########

## Set heatmap labels and group colors for heatmaps and VENN diagram
labRow <- rownames(mat.ncontigs)
labRowf <- rownames(matf.ncontigs)

grcols <- gg_color_hue(ngroups)
grcolsf <- grcols[levels(dmerge$LABEL) %in% names(lab.passed)] # passed groups
grcolsv <- grcols[levels(dmerge$LABEL) %in% names(lab.cons)] # considered groups

RowSideColors <- grcols[dmerge[!duplicated(dmerge$ID),]$LABEL] # all taxa
RowSideColorsf <- grcols[dmerge.passed[!duplicated(dmerge.passed$ID),]$LABEL] # passed taxa

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

## Create plots
# Taxon filtering
p0 <- ggplot(dploci, aes(LABEL, covered_loci, fill = LABEL)) +
  geom_boxplot(alpha = 0.5) +
  geom_point(aes(color = LABEL)) +
  geom_hline(aes(yintercept = min.nloci), col = "tomato") +
  geom_hline(aes(yintercept = nloci)) +
  scale_fill_discrete(guide = FALSE) +
  scale_colour_discrete(guide = FALSE) +
  scale_y_continuous(sec.axis = sec_axis(~./nloci, breaks = seq(0, 1, by = 0.1), labels = seq(0, 1, by = 0.1)), limits = c(0, nloci)) +
  labs(x = "", y = "Number of target loci covered by a taxon") +
  annotate("text", x = (1:ngroups) + 0.2, y = min.nloci, 
           label = paste0(ntaxa.ploci, " (", round(100*ntaxa.ploci/ntaxa,2), "%)"), 
           hjust = "top", vjust = "left") +
  ggtitle(paste0("Number of target loci: ", nloci, " ; required: ", min.nloci, " (", round(100*min.ploci,2), "%)\n",
                 "Number of taxa: ", ntaxa, " ; passed: ", tpassed, " (", round(100*tpassed/ntaxa,2), "%)")) +
  coord_flip() +
  theme_bw()

# Aggregated data
p1 <- plot.agg(dagg, y = "ptaxa", ylab = "proportion of taxa covered by a target locus", nloci.filt = nloci.ptaxa, sloci.filt = sloci.ptaxa, thr = min.ptaxa)

p2 <- plot.agg(dagg, y = "ncontigs", ylab = "number of (non-zero) contigs", nloci.filt = nloci.ncontigs, sloci.filt = sloci.ncontigs, thr = max.ncontigs)

p3 <- plot.agg(dagg, y = "bestscore.norm", ylab = "best normalized alignment score", nloci.filt = nloci.bestscore.norm, sloci.filt = sloci.bestscore.norm, thr = min.bestscore.norm)

p4 <- plot.agg(dagg, y = "taln", ylab = "target alignment length", nloci.filt = nloci.taln, sloci.filt = sloci.taln, thr = min.taln)

p5 <- plot.agg(dagg, y = "tfrac", ylab = "target alignment fraction", nloci.filt = nloci.tfrac, sloci.filt = sloci.tfrac, thr = min.tfrac)

p6 <- plot.agg(dagg, y = "bestscore", ylab = "best alignment score", nloci.filt = nloci.bestscore, sloci.filt = sloci.bestscore, thr = min.bestscore)

p7 <- plot.agg(dagg, y = "bestlength", ylab = "best contig length", nloci.filt = nloci.bestlength, sloci.filt = sloci.bestlength, thr = min.bestlength)

# target GC vs. number of contigs
p8 <- plot.comp(dagg, "ncontigscat", "tgc", "Number of contigs", "Target GC content", smooth = F, box = T)

# target length vs. number of contigs
p9 <- plot.comp(dagg, "ncontigscat", "tlen", "Number of contigs", "Target length", smooth = F, box = T)

# target alignment length vs. number of contigs
p10 <- plot.comp(dagg, "ncontigscat", "taln", "Number of contigs", "Target alignment length", smooth = F, box = T)

# target alignment fraction vs. number of contigs
p11 <- plot.comp(dagg, "ncontigscat", "tfrac", "Number of contigs", "Target alignment fraction", smooth = F, box = T)

# best normalized alignment score vs. number of contigs
p12 <- plot.comp(dagg, "ncontigscat", "bestscore.norm", "Number of contigs", "Normalized best alignment score", smooth = F, box = T)

# target alignment length vs. target GC content
p13 <- plot.comp(dagg, "tgc", "taln", "Target GC content", "Target alignment length", smooth = T, box = F)

# target alignment fraction vs. target GC content
p14 <- plot.comp(dagg, "tgc", "tfrac", "Target GC content", "Target alignment fraction", smooth = T, box = F)

# best normalized alignment score vs. target GC content
p15 <- plot.comp(dagg, "tgc", "bestscore.norm", "Target GC content", "Normalized best alignment score", smooth = T, box = F)

## Write PDF
pdf(file = outpdf, width = plot.width, height = plot.height)

# taxon filtering
cat("\nwriting ggplots...\n")
print(p0) # dploci$covered_loci

# VENN diagram for loci that passed filters in all groups
if (length(list.passed) <= 5) {
  # load libraries
  suppressPackageStartupMessages(library(grid))
  suppressPackageStartupMessages(library(VennDiagram))
  
  # set plot margins and color function
  oldmar <- par()$mar
  par(mar = rep(0,4))
  
  # plot VENN diagram for passed loci
  plot.new()
  futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
  grid.draw(venn.diagram(list.passed, filename = NULL, 
                         height = 3000, width = 3000,
                         cat.dist = rep(0, length(grcolsv)),
                         cat.col = grcolsv,
                         col = grcolsv,
                         fill = grcolsv,
                         alpha = 0.6,
                         cex = 2, cat.cex = 0,
                         lty = rep(2, length(grcolsv)),
                         print.mode = c("raw")))
  legend("topleft", names(list.passed), text.col = grcolsv, bty = "n", cex = 1)
  par(mar = oldmar)
} else {
  plot.new()
  msg1 <- paste0("Number of passed target loci: ", lpassed, "/", nloci, " (", round(100*lpassed/nloci,2), "%)")
  legend("topleft", c(names(list.passed), "", msg1), text.col = c(grcolsv, 1, 1), bty = "n", cex = 1)
  msg <- paste("more than 5 non-NA groups present: no VENN diagram shown.")
  mtext(msg, side = 1)
}

# aggreated locus data
print(p1) # dagg$ptaxa
print(p2) # dagg$ncontigs
print(p3) # dagg$bestscore.norm
print(p4) # dagg$taln
print(p5) # dagg$tfrac
print(p6) # dagg$bestscore
print(p7) # dagg$bestlength

# relationships between target sequence length / GC content and filtering variables
print(p8) # dagg$tgc vs. ncontigscat
print(p9) # dagg$tlen vs. ncontigscat
print(p10) # dagg$taln vs. ncontigscat
print(p11) # dagg$tfrac vs. ncontigscat
print(p12) # dagg$bestscore.norm vs. ncontigscat
print(p13) # dagg$taln vs. tgc
print(p14) # dagg$tfrac vs. tgc
print(p15) # dagg$bestscore.norm vs. tgc

# heatmaps all loci
cat("\nwriting heatmaps...\n")
do.heatmap(mat.ncontigs, max.nbreaks = 200, labRow, Colv = sortmat, reverse.colors = TRUE, highquant = 0.999, lmat, lwid, lhei, hclustmethod, RowSideColors, myfun, paste0("Number of (non-zero) contigs (", ntaxa, " taxa, ", nloci, " loci)"))
do.heatmap(mat.bestscore.norm, max.nbreaks = 200, labRow, Colv = sortmat, reverse.colors = FALSE, highquant = 0.975, lmat, lwid, lhei, hclustmethod, RowSideColors, myfun, paste0("Best normalized alignment score (", ntaxa, " taxa, ", nloci, " loci)"))
do.heatmap(mat.taln, max.nbreaks = 200, labRow, Colv = sortmat, reverse.colors = FALSE, highquant = 0.975, lmat, lwid, lhei, hclustmethod, RowSideColors, myfun, paste0("Target alignment length (", ntaxa, " taxa, ", nloci, " loci)"))
do.heatmap(mat.tfrac, max.nbreaks = 200, labRow, Colv = sortmat, reverse.colors = FALSE, highquant = 0.75, lmat, lwid, lhei, hclustmethod, RowSideColors, myfun, paste0("Target alignment fraction (", ntaxa, " taxa, ", nloci, " loci)"))
do.heatmap(mat.bestscore, max.nbreaks = 200, labRow, Colv = sortmat, reverse.colors = FALSE, highquant = 0.975, lmat, lwid, lhei, hclustmethod, RowSideColors, myfun, paste0("Best alignment score (", ntaxa, " taxa, ", nloci, " loci)"))
do.heatmap(mat.bestlength, max.nbreaks = 200, labRow, Colv = sortmat, reverse.colors = FALSE, highquant = 0.975, lmat, lwid, lhei, hclustmethod, RowSideColors, myfun, paste0("Best contig length (", ntaxa, " taxa, ", nloci, " loci)"))

# heatmaps passed loci
do.heatmap(matf.ncontigs, max.nbreaks = 200, labRowf, Colv = sortmat, reverse.colors = TRUE, highquant = 0.999, lmat, lwid, lhei, hclustmethod, RowSideColorsf, myfunf, paste0("Number of (non-zero) contigs (", tpassed, " passed taxa, ", lpassed, " passed loci)"))
do.heatmap(matf.bestscore.norm, max.nbreaks = 20, labRowf, Colv = sortmat, reverse.colors = FALSE, highquant =0.975, lmat, lwid, lhei, hclustmethod, RowSideColorsf, myfunf, paste0("Best normalized alignment score (", tpassed, " passed taxa, ", lpassed, " passed loci)"))
do.heatmap(matf.taln, max.nbreaks = 200, labRowf, Colv = sortmat, reverse.colors = FALSE, highquant = 0.975, lmat, lwid, lhei, hclustmethod, RowSideColorsf, myfunf, paste0("Target alignment length (", tpassed, " passed taxa, ", lpassed, " passed loci)"))
do.heatmap(matf.tfrac, max.nbreaks = 200, labRowf, Colv = sortmat, reverse.colors = FALSE, highquant = 0.75, lmat, lwid, lhei, hclustmethod, RowSideColorsf, myfunf, paste0("Target alignment fraction (", tpassed, " passed taxa, ", lpassed, " passed loci)"))
do.heatmap(matf.bestscore, max.nbreaks = 200, labRowf, Colv = sortmat, reverse.colors = FALSE, highquant = 0.975, lmat, lwid, lhei, hclustmethod, RowSideColorsf, myfunf, paste0("Best alignment score (", tpassed, " passed taxa, ", lpassed, " passed loci)"))
do.heatmap(matf.bestlength, max.nbreaks = 200, labRowf, Colv = sortmat, reverse.colors = FALSE, highquant = 0.975, lmat, lwid, lhei, hclustmethod, RowSideColorsf, myfunf, paste0("Best contig length (", tpassed, " passed taxa, ", lpassed, " passed loci)"))

graphics.off()

#########
## LOG ##
#########
t2 <- Sys.time()
dlog <- c(
  paste0("=== Filter taxa and loci based on visualization of (aggregated) statistics ==="),
  paste0(""),
  paste0("Analysis started at: ", t1),
  paste0(""),
  paste0("Input data:"),
  paste0("-----------"),
  paste0("Locus stats file: ", file),
  paste0("Reference sequences: ", refseqs),
  paste0("Metadata file (mapping taxa to groups): ", meta),
  paste0(""),
  paste0("Thresholds for taxon filtering:"),
  paste0("-------------------------------"),
  paste0("Min. proportion of target loci covered by taxon:           ", min.ploci),
  paste0(""),
  paste0("Thresholds for locus filtering:"),
  paste0("-------------------------------"),
  paste0("Min. ", tolower(funlab), " proportion of taxa covered by target loci:     ", min.ptaxa),
  paste0("Max. ", tolower(funlab), " number of (non-zero) contigs per target locus: ", max.ncontigs),
  paste0("Min. ", tolower(funlab), " best normalized alignment score:               ", min.bestscore.norm),
  paste0("Min. ", tolower(funlab), " target locus alignment length:                 ", min.taln),
  paste0("Min. ", tolower(funlab), " target locus alignment fraction:               ", min.tfrac),
  paste0("Min. ", tolower(funlab), " best raw alignment score:                      ", min.bestscore),
  paste0("Min. ", tolower(funlab), " best contig length:                            ", min.bestlength),
  paste0(""),
  paste0("Other (internal) parameters:"),
  paste0("----------------------------"),
  paste0("Aggregation summary function: ", tolower(funlab)),
  paste0("Contig number summary categories: ", 0, " - >=", XandMore),
  paste0("Hierarchical clustering of loci in heatmaps: ", as.character(sortmat)),
  paste0("Hierarchical clustering method: ", hclustmethod),
  paste0("PDF height (inches): ", plot.height),
  paste0("PDF width (inches): ", plot.width),
  paste0(""),
  paste0("Taxon filtering results:"),
  paste0("------------------------"),
  paste0("Total number of taxa: ", ntaxa),
  paste0("Number of passed taxa: ", tpassed, " (", round(100*tpassed/ntaxa,2), "%)"),
  paste0("Number of considered taxa: ", staxa.cons, " (", round(100*staxa.cons/ntaxa,2), "%)"),
  paste(""),
  paste0("Total number of groups: ", ngroups),
  paste0("Number of passed groups: ", ngr.passed, " (", round(100*ngr.passed/ngroups,2), "%)"),
  paste0("Number of considered groups: ", ngr.cons, " (", round(100*ngr.cons/ngroups,2), "%)"),
  paste0(""),
  paste0("Locus filtering results:"),
  paste0("------------------------"),
  paste0("Total number of target loci: ", nloci),
  paste0("Number of passed target loci (passed in all non-NA groups): ", lpassed, " (", round(100*lpassed/nloci,2), "%)"),
  paste0("Min. number of target loci covered by taxon: ", ceiling(min.ploci*nloci), " (", round(100*(min.ploci*nloci)/nloci,2), "%)"),
  paste0(""),
  paste0("Number of target loci with : "),
  paste0("Min. ", tolower(funlab), " proportion of taxa covered by target loci below ", min.ptaxa, ": ", sloci.ptaxa, " (", round(100*sloci.ptaxa/nloci,2), "%)"),
  paste0("Max. ", tolower(funlab), " number of (non-zero) contigs per target locus above ", max.ncontigs, ": ", sloci.ncontigs, " (", round(100*sloci.ncontigs/nloci,2), "%)"),
  paste0("Min. ", tolower(funlab), " best normalized alignment score below ", min.bestscore.norm, ": ", sloci.bestscore.norm, " (", round(100*sloci.bestscore.norm/nloci,2), "%)"),
  paste0("Min. ", tolower(funlab), " target locus alignment length below ", min.taln, ": ", sloci.taln, " (", round(100*sloci.taln/nloci,2), "%)"),
  paste0("Min. ", tolower(funlab), " target locus alignment fraction below ", min.tfrac, ": ", sloci.tfrac, " (", round(100*sloci.tfrac/nloci,2), "%)"),
  paste0("Min. ", tolower(funlab), " best raw alignment score below ", min.bestscore, ": ", sloci.bestscore, " (", round(100*sloci.bestscore/nloci,2), "%)"),
  paste0("Min. ", tolower(funlab), " best contig length below ", min.bestlength, ": ", sloci.bestlength, " (", round(100*sloci.bestlength/nloci,2), "%)"),
  paste0(""),
  paste0("Output files:"),
  paste0("-------------"),
  paste0("Log file:                               ", logtxt),
  paste0("Visualization of (aggregated) stats:    ", outpdf),
  paste0("Merged target locus stats and metadata: ", outtxt),
  paste0("Aggregated locus statistics:            ", aggtxt),
  paste0("Passed target loci:                     ", loctxt),
  paste0("Taxon statistics:                       ", aggtax),
  paste0("Passed taxa:                            ", taxtxt),
  paste0(""),
  paste0("Analysis completed at: ", t2)
)

## Write TXT
cat("\nwriting results .txt and .log...\n")
write.table(dmerge, file = outtxt, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
write.table(dploci, file = aggtax, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
write.table(dagg, file = aggtxt, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
write.table(loci.passed, file = loctxt, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
write.table(taxa.passed, file = taxtxt, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
writeLines(text = dlog, con = logtxt)

## Print finish time
cat("\n")
print(t2)
print(t2 - t1)



