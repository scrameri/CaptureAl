#!/cluster/apps/r/3.6.0_openblas/x86_64/bin/Rscript

## Usage: plot.readcounts.R <read.counts.txt>

## Value: barplot of read counts

## Author: simon.crameri@env.ethz.ch, Feb 2020

## Get arguments
args = commandArgs(trailingOnly = TRUE)

## Check arguments
options(warning.length = 2000L)
if (length(args) < 1) {
  stop("1 argument needed:
       REQUIRED
       1) <x|CHR>:                path to file with read counts (sample name in first column, read counts in second column, no header) 
       ",
       call.=FALSE)
}

## Set arguments
x <- as.character(args[1])

## Check arguments
stopifnot(file.exists(x))

## Get read counts
dd <- read.delim(x, header = FALSE)
stopifnot(is.numeric(dd$V2))
cat("found read counts for", nrow(dd), "samples\n")

## Get basic stats
# mior <- sum(dd$V2)/1E06
# summ <- summary(dd$V2)/1E06
scaleby=1E06
roundby=3
s <- paste0("Stats in Millions:\n",
            "Min: ", round(min(dd$V2)/scaleby, roundby), 
            "\n", "q25: ", round(quantile(dd$V2, 0.25)/scaleby, roundby), 
            "\n", "med: ", round(quantile(dd$V2, 0.5)/scaleby, roundby), 
            "\n", "q75: ", round(quantile(dd$V2, 0.75)/scaleby, roundby), 
            "\n", "max: ", round(max(dd$V2)/scaleby, roundby))

## Plot read counts
oname <- "read.counts.pdf"
pdf(oname)
omar <- par()$mar
par(mar = c(8, 6, 4, 2))
barplot(dd$V2, names.arg = dd$V1, las = 2, cex.names = 0.5, cex.axis = 0.5, ylab = "Read count")
abline(h = seq(0, 100E06, by = 1E06), lty = 2, col = "gray45")
text(x = 1, y = max(dd$V2), labels = s, xpd = TRUE, pos = 2, cex = 0.5)

## Plot density
if ("ggplot2" %in% rownames(installed.packages())) {
  library(ggplot2)
  ggplot(dd, aes(x = factor(0), V2)) +
    geom_boxplot(fill = "blue", alpha = 0.5) +
    geom_violin(fill = "blue", alpha = 0.5) +
    labs(x = "", y = "Read count") +
    theme_bw()
}
graphics.off()
par(mar = omar)
