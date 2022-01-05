#!/cluster/apps/r/3.6.0_openblas/x86_64/bin/Rscript

#####################################
## Plot contig numbers (per group) ##
#####################################

## Usage
# run script without arguments to see the arguments it takes

## Author
# simon.crameri@env.ethz.ch, Apr 2019

## Load required library
suppressPackageStartupMessages(library(ggplot2))

## Get arguments
args = commandArgs(trailingOnly = TRUE)

## Check arguments
if (! length(args) %in% c(1, 2)) {
  stop("At least 1 argument needed:
       REQUIRED
       1) <file|CHR>: path to collected contig numbers (loci_contignumbers.txt)

       OPTIONAL
       2) <meta|CHR>: path to metadata file mapping individuals (1st column) to groups (2nd column). Header and tab separation expected. More individuals in different order than in <file> are ok.
       ", call.=FALSE)
}

## Set arguments
file <- as.character(args[1])
meta <- as.character(args[2])

## Set arguments (for debugging)
# file = "loci_contignumbers.txt"
# meta <- "mapfile.orig.cpgroup.txt"

## Additional arguments
ycols = c("#A6CEE3","#99CD91","#B89B74","#F06C45","#ED8F47","#825D99","#B15928") # group colors, will be recycled

## Read contig numbers
dc <- read.delim(file, check.names = FALSE)
names(dc)[1] <- "ID"

## Read metadata
if (!is.na(meta)) {
  dmeta <- read.delim(meta)
  names(dmeta)[1] <- "ID"
  names(dmeta)[2] <- "GROUP"
}

## Merge it
if (!is.na(meta)) {
  dmerge <- merge(dc, dmeta, by = "ID")
  dmerge$GROUP <- as.character(dmerge$GROUP)
} else {
  dmerge <- dc
  dmerge$GROUP <- "NA"
}

## Handle NA values as an own category
dmerge$GROUP[is.na(dmerge$GROUP)] <- "NA"
dmerge$GROUP <- factor(dmerge$GROUP)

## Add number of individuals per group
grtab <- table(dmerge[!duplicated(dmerge$ID),"GROUP"])
cat("\nfound", sum(grtab), "matches:\n")
print(grtab)
dmerge$LABEL <- as.character(dmerge$GROUP)
for (group in unique(dmerge$GROUP)) {
  dmerge$LABEL[dmerge$LABEL == group] <- paste0(names(grtab[group]), " (n = ", grtab[group], ")")
}

## Recode column names
checknames <- names(dmerge)[seq(min(grep("ID", names(dmerge)))+1, min(grep("GROUP", names(dmerge)))-1, by = 1)]
sumnames <- gsub("exactlyatleast", "atleast", paste0("exactly", gsub("^>=", "atleast", checknames)))
names(dmerge)[names(dmerge) %in% checknames] <- sumnames
dmerge$nloci <- apply(dmerge[,sumnames], 1, sum)

## Reshape to long format for plotting
dpl <- reshape(dmerge[,c("ID","GROUP","LABEL", sumnames)], 
               idvar = "ID",
               varying = sumnames,
               times = sumnames,
               v.names = c("number"),
               direction = "long")
dpl$time <- factor(dpl$time, levels = sumnames)

## Plot it
nygr <- nlevels(dmerge$GROUP)
ycols <- rep(ycols, ceiling(nygr/length(ycols)))[1:nygr]

pdf(file = gsub(".txt$", ".pdf", file))
ggplot(dpl, aes(time, number, fill = LABEL)) +
  geom_boxplot(alpha = 0.75, width = 0.5) +
  geom_hline(yintercept = unique(dmerge$nloci)) +
  labs(y = "Number of regions",
       x = "Number of successful exonerate alignments",
       fill = "") +
  scale_x_discrete(labels = checknames) +
  scale_fill_manual(values = ycols) +
  ggtitle(paste("Total number of regions:", unique(dmerge$nloci))) +
  theme_bw()
graphics.off()


