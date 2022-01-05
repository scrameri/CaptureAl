#!/cluster/apps/r/3.6.0_openblas/x86_64/bin/Rscript

## Usage: plot.fasta.lengths.R <table of fasta lengths> <OPTIONAL second, tird, forth, etc. table of fasta lengths for comparison>

## Load libraries
suppressPackageStartupMessages(library(ggplot2))

## Author: simon.crameri@env.ethz.ch, May 2019

## Get arguments
args <- commandArgs(trailingOnly=TRUE)

## Check arguments
if (!length(args) >= 1) {
  stop("n arguments taken (1 needed):
       REQUIRED
       1) 1st file with fasta sequence lengths (1st column: fasta header, 2nd column: length, tab-delmited, no header) (will be labelled '1'); 
       
       OPTIONAL
       2) 2nd file with fasta sequence lengths (1st column: fasta header, 2nd column: length, tab-delmited, no header) (will be labelled '2');
       3) 3rd file with fasta sequence lengths (1st column: fasta header, 2nd column: length, tab-delmited, no header) (will be labelled '3');
       4), 5), ... n)",
       call.=FALSE)
}

## Get arguments
for (tab in args) {
  id <- which(args == tab)
  dd <- read.delim(tab, header = F, stringsAsFactors = F)
  names(dd) <- c("Locus", "Length")
  dd$FASTA <- id
  if (!exists("d.plot")) {
    d.plot <- dd
  } else {
    d.plot <- rbind(d.plot, dd)
  }
}

## Helperfunctions
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

## Get numer of sequences and median lengths
d.plot$FASTA <- factor(d.plot$FASTA)
meds <- tapply(X = d.plot$Length, INDEX = d.plot$FASTA, FUN = median)
lens <- tapply(X = d.plot$Length, INDEX = d.plot$FASTA, FUN = length)

## Get annotated fasta labels and median lengths
d.plot <- merge(d.plot, data.frame(FASTA = names(meds), med = meds, n = lens), by = "FASTA", all.x = T, sort = F)
d.plot$LABEL <- paste0(d.plot$FASTA, " (n=", d.plot$n, ", med=", d.plot$med, ")")

## Get total lengths
d.totlen <- aggregate(d.plot$Length, by = list(d.plot$LABEL), FUN = sum)
names(d.totlen) <- c("LABEL", "totlen")

## Get colors and caption
cols <- gg_color_hue(length(unique(d.plot$FASTA))) # color choice here
ptit <- paste0("sequence lengths min: ", min(d.plot$Length),
               " ; 1st Qu.: ", quantile(d.plot$Length, 0.25),
               " ; med: ", median(d.plot$Length), 
               " ; mean: ", round(mean(d.plot$Length)),
               " ; 3rd Qu.: ", quantile(d.plot$Length, 0.75),
               " ; max: ", max(d.plot$Length))

## Make density plots and boxplots
# total length
p1 <- ggplot(d.totlen, aes(x = factor(LABEL), y = totlen, fill = factor(LABEL))) +
  geom_bar(stat = "identity", alpha = 0.5) +
  theme_bw() +
  scale_fill_discrete(guide = F) +
  labs(x = "", y = "Total length in FASTA [bp]", caption = ptit)
  #coord_flip()
  
p2 <- ggplot(d.plot, aes(x = Length, fill = LABEL)) +
  #geom_histogram(aes(y = ..density..), alpha = 0.3, bins = 50, position = "identity") +
  geom_density(aes(colour = LABEL), alpha = 0.3) +
  xlab("Sequence length [bp]") +
  scale_fill_manual(name = "", values = cols) +
  scale_color_manual(name = "", values = cols) +
  labs(caption = ptit) +
  theme_bw() + 
  theme(legend.position = "top")

p3 <- p2 + scale_x_continuous(trans = "log10")

p4 <- ggplot(d.plot, aes(x = LABEL, y = Length, fill = LABEL)) +
  geom_boxplot(alpha = 0.5) +
  geom_violin(alpha = 0.5) +
  scale_fill_manual(name = "", values = cols) +
  labs(x = "", y = "Sequence length [bp]", caption = ptit) +
  theme_bw() + 
  theme(legend.position = "top")

p5 <- p4 + scale_y_continuous(trans = "log10")

## Plot histograms
oname <- paste0(args[1], ".pdf")
if (file.exists(oname)) system(paste("mv", oname, gsub(".pdf$", ".bak.pdf", oname)))
pdf(file = oname, width = 6+length(args), height = 6+length(args))
print(p1)
print(p2)
print(p3)
print(p4)
print(p5)
graphics.off()

