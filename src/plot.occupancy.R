#!/cluster/apps/r/3.6.0_openblas/x86_64/bin/Rscript

#########################################################
## Plot occupancy stats (per individual and per group) ##
#########################################################

## Usage
# run script without arguments to see the arguments it takes

## Author
# simon.crameri@env.ethz.ch, Mar 2019

## Load required library
suppressPackageStartupMessages(library(ggplot2))

## Get arguments
args = commandArgs(trailingOnly = TRUE)

## Check arguments
options(warning.length = 2000L)
options(width = 1000)
if (! length(args) %in% c(2, 3)) {
  stop("At least 2 arguments needed:
       REQUIRED
       1) <tab|CHR>: path to occupancy stats file (.occ) ;
       2) <thresh|NUM>: minimum occupancy per individual ;
        
       OPTIONAL
       3) <meta|CHR>: path to meta file that maps individuals to groups",
       call.=FALSE)
}

## Set arguments
tab <- as.character(args[1])
thresh <- as.numeric(as.character(args[2]))
meta <- as.character(args[3])

## Set arguments (for debugging)
# tab <- "test.occ"
# thresh <- 0
# meta <- "mapfile.txt"

## Additional arguments
thresh.lab <- 0.6     # inds with occupancy between 0 and thresh.lab will be displayed
occvar <- "occupancy" # name of occupancy field 
ycols <- c("#A6CEE3","#99CD91","#B89B74","#F06C45","#ED8F47","#825D99","#B15928") # group colors, will be recycled

## Check parameters
stopifnot(file.exists(tab), thresh <= 1, thresh >=0)
if (!is.na(meta)) stopifnot(file.exists(meta))

## Define helperfuncitons
# print keeping number and proportion
get.keeping.prop <- function(thresh, df, var = "occupancy", return.prop = TRUE) {
  n <- nrow(df)
  nkeep <- sum(as.numeric(as.character(df[,var])) >= thresh)
  pkeep <- 100*nkeep/n
  
  if (return.prop) {
    return(round(pkeep,1))
  } else {
    return(nkeep)
  }
}


## Read input
dd <- read.delim(file = tab)
names(dd)[1] <- "taxon"
names(dd)[2] <- occvar

## Read metadata
if (!is.na(meta)) {
  dmeta <- read.delim(meta)
  names(dmeta)[1] <- "taxon"
  names(dmeta)[2] <- "GROUP"
}

## Merge occupancy and meta data
if (!is.na(meta)) {
  if (!all(dd$taxon %in% dmeta$taxon)) {
    notfound <- dd$taxon[!dd$taxon %in% dmeta$taxon]
    stop(paste0("found ", length(notfound), " individuals in <", tab, "> (1st field) that could not be matched to <", meta, "> (1st field):\n", paste(notfound, collapse = ", "), "\n"))
  }
  dmerge <- merge(dd, dmeta, by = "taxon", all.x = TRUE, all.y = FALSE, sort = FALSE)
  dmerge$GROUP <- droplevels(dmerge$GROUP)
  
  dmerge$GROUP <- as.character(dmerge$GROUP)
  dmerge[is.na(dmerge$GROUP), "GROUP"] <- "NA"
  dmerge <- dmerge[order(dmerge[,occvar]),]
  dmerge$LABEL <- paste(dmerge$taxon, dmerge$GROUP)
  dmerge <- dmerge[,c("taxon","GROUP","LABEL", names(dmerge)[!names(dmerge) %in% c("taxon","GROUP","LABEL")])]
  
  cat("found", nrow(dmerge), "matches:")
  print(table(dmerge$GROUP))
} else {
  dmerge <- dd
  dmerge$GROUP <- "NA"
  dmerge$LABEL <- paste(dmerge$taxon, dmerge$GROUP)
  dmerge <- dmerge[,c("taxon","GROUP","LABEL", names(dmerge)[!names(dmerge) %in% c("taxon","GROUP","LABEL")])]
}

## Handle NA values as an own category
dmerge$GROUP[is.na(dmerge$GROUP)] <- "NA"
dmerge$GROUP <- factor(dmerge$GROUP)

## Order by group (highest to lowest occupancy)
dgr <- aggregate(dmerge[,occvar], by = list(dmerge$GROUP), FUN = median, na.rm = T)
levs <- dgr[order(dgr$x),1]
dmerge$GROUP <- factor(dmerge$GROUP, levels = levs)

## Be verbose
testseq <- sort(unique(c(seq(0, 1, by = 0.05), thresh)))
names(testseq) <- as.character(testseq)

cat("proportion kept:\n")
print(sapply(testseq, FUN = get.keeping.prop, df = dmerge, return.prop = TRUE))

cat("\nnumber of taxa kept:\n")
print(sapply(testseq, FUN = get.keeping.prop, df = dmerge, return.prop = F))

cat("\nnumber of taxa removed:\n")
print(nrow(dmerge)-sapply(testseq, FUN = get.keeping.prop, df = dmerge, return.prop = F))

cat("\n")

## Plot occupancy stats 
# per individual
n <- nrow(dmerge)
nkeep <- sum(dmerge[,occvar] >= thresh)
nrm <- n-nkeep
dlab <- dmerge[dmerge[,occvar] <= thresh.lab,]

p1 <- ggplot(dmerge, aes_string(rev(occvar))) +
  geom_histogram(aes(y = ..count..), bins = nrow(dmerge), color = "gray", alpha = 0.1) +
  geom_density(aes(y = ..density..), fill = "tomato", alpha = 0.5) +
  geom_vline(xintercept = thresh) +
  #lims(x = c(0,1)) +
  xlab("Occupancy (fraction of alignment sites with nucleotides)") +
  ggtitle(paste0("Threshold ", thresh, " keeps ", nkeep, " / ", n, " [", round(100*nkeep/n,2), "%], removes ", nrm, " / ", n, " [", round(100*nrm/n,2), "%] taxa")) +
  theme_bw()

if (nrow(dlab) > 0) p1 <- p1 + geom_text(aes_string(x = occvar, y = 1.1, label = "LABEL"), data = dlab, angle = 90, size = 2.5, vjust = 0, hjust = 0)

# per group
nygr <- nlevels(dmerge[,"GROUP"])
ycols <- rep(ycols, ceiling(nygr/length(ycols)))[1:nygr]

# add number of individuals per group (all samples, all groups)
grtab <- table(dmerge[!duplicated(dmerge[,"taxon"]),"GROUP"])
dmerge[,"LABEL2"] <- as.character(dmerge[,"GROUP"])
for (group in levels(dmerge[,"GROUP"])) {
  dmerge[,"LABEL2"][dmerge[,"LABEL2"] == group] <- paste0(names(grtab[group]), " (n = ", grtab[group], ")")
}
dgr <- aggregate(dmerge[,occvar], by = list(dmerge[,"LABEL2"]), FUN = median, na.rm = T)
levs <- dgr[order(dgr$x),1]
dmerge[,"LABEL2"] <- factor(dmerge[,"LABEL2"], levels = levs)

p2 <- ggplot(dmerge, aes_string("LABEL2", occvar, fill = "LABEL2")) +
  geom_point() +
  geom_boxplot(alpha = 0.5) +
  geom_violin(alpha = 0.5) +
  xlab("") +
  ylab("Occupancy (fraction of alignment sites with nucleotides)") +
  scale_fill_manual(guide = FALSE, values = ycols[sapply(as.character(levs), function(x) which(sort(as.character(levs)) == x))]) +
  lims(y = c(0,1)) +
  coord_flip() +
  theme_bw()

pdf(file = paste0(tab, ".pdf"), width = 7, height = 7)
suppressWarnings(print(p1))
print(p2)
graphics.off()

## Write merged data
dmerge <- dmerge[order(dmerge[,occvar]),]
write.table(dmerge, file = paste0(tab, ".txt"), quote = F, row.names = F, col.names = T, sep = "\t")

## Write taxa to be removed
taxrm <- dmerge[which(dmerge[,occvar] < thresh),]
taxrm <- taxrm[order(taxrm[,occvar]),]
if (nrow(taxrm) > 0) {
  write.table(taxrm[,1], file = paste0(tab, "_rm_", thresh), quote = F, row.names = F, col.names = F)
}

