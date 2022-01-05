#!/cluster/apps/r/3.6.0_openblas/x86_64/bin/Rscript

## Usage: plot.pcr.duplicates.R <sfile> <stats>

## Load libraries
suppressPackageStartupMessages(library(ggplot2)) # ggplot

## Author: simon.crameri@env.ethz.ch, Mar 2020

## Get arguments
args <- commandArgs(trailingOnly=TRUE)

## Check arguments
options(warning.length = 5000L)
if (! length(args) %in% c(2)) {
  stop("2 arguments needed:
       REQUIRED
       1) <sfile|CHR>:              path to samples file. Header and tab-separation expected.
                                    Sample IDs must be in the first column. Group IDs can be specified in the 
                                    second column (if not specified, all samples are assumed to constitute one group).
                                    The group ID is used to apply region filtering criteria 5-12 within all considered 
                                    groups, to determine regions passing the filtering criteria in all groups.
                                    Samples that do not belong to any specified group (second column empty or 'NA') 
                                    will be displayed in summary plots but will not be considerd during region filtering. 
                                    Additional columns are ignored.
       2) <stats|CHR>:              path to PCR duplication stats. Header and tab-separation expected.
                                    Sample IDs must be in the first column. Percent duplication must be in the second column.",
       call.=FALSE)
}

## Set arguments
sfile <- as.character(args[1])
stats <- as.character(args[2])

## Set arguments (for debugging)
#sfile = "mapfile.txt"
#stats = "percdup.txt"

## Additional arguments
idvar <- "ID"
grpvar1 <- "GROUP"
grpvar2 <- "LABEL"
navar <- "NA"

ycols <- c("#A6CEE3","#99CD91","#B89B74","#F06C45","#ED8F47","#825D99","#B15928") # group colors, will be recycled

statext <- paste0(".", rev(unlist(strsplit(stats, split = "[.]")))[1])
oname <- gsub(statext, ".pdf", stats)
height = 7
width = 7

## Read PCR duplication stats
dd <- read.delim(stats, header = TRUE)
names(dd) <- c(idvar,"PERCDUP")

## Read sample / group file
ds <- read.delim(sfile, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
if (ncol(ds) == 1) ds[,grpvar1] <- rep("Undefined", nrow(ds))
names(ds) <- c(idvar, grpvar1)
samples <- ds[,idvar]
ds[ds[,grpvar1] %in% c(NA, "NA", ""),grpvar1] <- navar # sets <NA> or <> group values to "NA"

## Merge data
dmerge <- merge(dd, ds, by = idvar, all.x = FALSE, all.y = TRUE, sort = FALSE)

## Add number of individuals per group (all samples, all groups)
grtab <- table(dmerge[!duplicated(dmerge[,idvar]),grpvar1])
dmerge[,grpvar2] <- as.character(dmerge[,grpvar1])
for (group in unique(dmerge[,grpvar1])) {
  dmerge[,grpvar2][dmerge[,grpvar2] == group] <- paste0(names(grtab[group]), " (n = ", grtab[group], ")")
}
dmerge[,grpvar2] <- as.factor(dmerge[,grpvar2])
dmerge[,grpvar1] <- as.factor(dmerge[,grpvar1])

## Plot data
nygr <- nlevels(dmerge[,grpvar2])
ycols <- (rep(ycols, ceiling(nygr/length(ycols)))[1:nygr])

p <- ggplot(dmerge, aes_string(x = grpvar2, y = "PERCDUP", fill = grpvar2)) +
  geom_boxplot(alpha = 0.5) +
  geom_violin(alpha = 0.5) +
  scale_fill_manual(guide = FALSE, values = ycols) +
  labs(x = "", y = "Fraction of duplicated reads") +
  lims(y = c(0,1)) +
  coord_flip() +
  theme_bw()

pdf(oname, height = height, width = width)
print(p)
graphics.off()

