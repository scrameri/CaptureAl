#!/cluster/apps/r/3.6.0_openblas/x86_64/bin/Rscript

## Usage: plot.fasta.stats.R <fasta1> <OPTIONAL second, tird, forth, etc. fastas for comparison>

## Load libraries
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ape))

## Author: simon.crameri@env.ethz.ch, May 2019

## Get arguments
print(Sys.time())
args <- commandArgs(trailingOnly=TRUE)

## Check arguments
if (!length(args) >= 1) {
  stop("n arguments taken (1 needed):
       REQUIRED
       1) 1st fasta (will be labelled '1'); 
       
       OPTIONAL
       2) 2nd fasta (will be labelled '2');
       3) 3rd fasta (will be labelled '3');
       4), 5), ... n)",
       call.=FALSE)
}

## Helperfunctions
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
get.stats <- function(x)  {
  tab <- table(factor(x, levels = c("a","c","g","t","k","y","w","s","r","m","b","d","h","v","n","-","?")))
  nagct <- sum(tab[c("a","c","g","t")])
  gc <- sum(tab[c("g","c")])/nagct
  cons <- sum(tab[c("k","y","w","s","r","m","b","d","h","v")])/nagct
  ns <- sum(tab[c("n","?")])/nagct
  gaps <- sum(tab["-"])/nagct
  len <- length(x)
  df <- data.frame(length = len, gc = gc, cons = cons, gaps = gaps, ns = ns)
  return(df)
}

## Get arguments
for (tab in args) {
  id <- which(args == tab)
  cat("processing", id, "/", length(args), "fastas...\r")

  if (file.exists(tab)) {

    fas <- read.FASTA(tab)
  
    dd <- data.frame(as.data.frame(t(sapply(as.character(fas), get.stats))))
    for (i in seq(ncol(dd))) dd[,i] <- unlist(dd[,i])
    dd$FASTA <- id
    dd$file <- basename(tab)
  
    if (!exists("d.plot")) {
      d.plot <- dd
    } else {
      d.plot <- rbind(d.plot, dd)
    }
  } else {
  warning("file <", tab, "> not found, skipping...")
  }
}

cat("\ndone\n")

## Additional arguments
writetxt <- TRUE # if TRUE, will write d.plot to disk

## Get number of sequences
d.plot$FASTA <- factor(d.plot$FASTA, levels = as.character(seq(length(args))))
# gcmed <- tapply(X = d.plot$gc, INDEX = d.plot$FASTA, FUN = median)
lens <- tapply(X = d.plot$gc, INDEX = d.plot$FASTA, FUN = length)

## Get annotated fasta labels
d.plot <- merge(d.plot, data.frame(FASTA = names(lens), n = lens), by = "FASTA", all.x = T, sort = F)
d.plot$LABEL <- paste0(d.plot$FASTA, " (n=", d.plot$n, ")")

cols <- gg_color_hue(length(unique(d.plot$FASTA))) # color choice here
ptit <- paste0(1:length(args), ": ", args, collapse = ", ")
ptit <- gsub(", 11", ",\n11", gsub(", 6", ",\n6", ptit)) # line break after 5 or 10 fastas

## Get total length in FASTA
d.totlen <- aggregate(d.plot$length, by = list(d.plot$LABEL), FUN = sum)
names(d.totlen) <- c("LABEL","totlen")

## Make plots
cat("plotting...")

p1.1 <- ggplot(d.totlen, aes(x = factor(LABEL), y = totlen, fill = factor(LABEL))) +
  geom_bar(stat = "identity", alpha = 0.5) +
  theme_bw() +
  scale_fill_discrete(guide = F) +
  labs(x = "", y = "Total length in FASTA [bp]", caption = ptit)
  #coord_flip()

# lenghts
p1.2 <- ggplot(d.plot, aes(x = length, fill = LABEL)) +
  geom_density(aes(colour = LABEL), alpha = 0.3) +
  scale_fill_manual(name = "", values = cols) +
  scale_color_manual(name = "", values = cols) +
  labs(x = "Sequence length [bp]", caption = ptit) +
  theme_bw() +
  theme(legend.position = "top")

p1.3 <- ggplot(d.plot, aes(x = length, fill = LABEL)) +
  geom_density(aes(colour = LABEL), alpha = 0.3) +
  scale_fill_manual(name = "", values = cols) +
  scale_color_manual(name = "", values = cols) +
  scale_x_continuous(trans = "log10") +
  labs(x = "Sequence length [bp]", caption = ptit) +
  theme_bw() +
  theme(legend.position = "top")

p1.4 <- ggplot(d.plot, aes(x = LABEL, y = length, fill = LABEL)) +
  geom_boxplot(alpha = 0.5) +
  geom_violin(alpha = 0.5) +
  scale_fill_manual(name = "", values = cols) +
  scale_color_manual(name = "", values = cols) +
  labs(x = "", y = "Sequence length [bp]", caption = ptit) +
  theme_bw() +
  theme(legend.position = "top")

p1.5 <- ggplot(d.plot, aes(x = LABEL, y = length, fill = LABEL)) +
  geom_boxplot(alpha = 0.5) +
  geom_violin(alpha = 0.5) +
  scale_fill_manual(name = "", values = cols) +
  scale_color_manual(name = "", values = cols) +
  scale_y_continuous(trans = "log10") +
  labs(x = "", y = "Sequence length [bp]", caption = ptit) +
  theme_bw() +
  theme(legend.position = "top")

# gc content
p2.1 <- ggplot(d.plot, aes(x = gc, fill = LABEL)) +
  geom_density(aes(colour = LABEL), alpha = 0.3) +
  scale_fill_manual(name = "", values = cols) +
  scale_color_manual(name = "", values = cols) +
  labs(x = "GC content", caption = ptit) +
  theme_bw() + 
  theme(legend.position = "top")

p2.2 <- ggplot(d.plot, aes(x = LABEL, y = gc, fill = LABEL)) +
  geom_boxplot(alpha = 0.5) +
  geom_violin(alpha = 0.5) +
  scale_fill_manual(name = "", values = cols) +
  scale_color_manual(name = "", values = cols) +
  labs(x = "", y = "GC content", caption = ptit) +
  theme_bw() +
  theme(legend.position = "top")

# gc vs length
p3.1 <- ggplot(d.plot, aes(x = gc, y = length, fill = LABEL, color = LABEL)) +
  geom_point(alpha = 0.1) +
  geom_density_2d(aes(colour = LABEL), alpha = 1) +
  xlab("GC content") +
  ylab("Sequence length [bp]") +
  scale_fill_manual(name = "", values = cols) +
  scale_color_manual(name = "", values = cols) +
  labs(caption = ptit) +
  theme_bw() + 
  theme(legend.position = "top")

p3.2 <- ggplot(d.plot, aes(x = gc, y = length, fill = LABEL, color = LABEL)) +
  geom_point(alpha = 0.1) +
  geom_density_2d(aes(colour = LABEL), alpha = 1) +
  xlab("GC content") +
  ylab("Sequence length [bp]") +
  scale_fill_manual(name = "", values = cols) +
  scale_color_manual(name = "", values = cols) +
  scale_y_continuous(trans = "log10") +
  labs(caption = ptit) +
  theme_bw() +
  theme(legend.position = "top")

# consensus content
p4 <- ggplot(d.plot, aes(x = LABEL, y = cons, fill = LABEL)) +
  geom_boxplot(alpha = 0.5) +
  geom_violin(alpha = 0.5) +
  scale_fill_manual(name = "", values = cols) +
  labs(x = "", y = "Consensus base (k,y,w,s,r,m,b,d,h,v) content", caption = ptit) +
  theme_bw() + 
  theme(legend.position = "top")

# gap content
p5 <- ggplot(d.plot, aes(x = LABEL, y = gaps, fill = LABEL)) +
  geom_boxplot(alpha = 0.5) +
  geom_violin(alpha = 0.5) +
  scale_fill_manual(name = "", values = cols) +
  labs(x = "", y = "Gap (-) content", caption = ptit) +
  theme_bw() +
  theme(legend.position = "top")

# undet base content
p6 <- ggplot(d.plot, aes(x = LABEL, y= ns, fill = LABEL)) +
  geom_boxplot(alpha = 0.5) +
  geom_violin(alpha = 0.5) +
  scale_fill_manual(name = "", values = cols) +
  labs(x = "", y = "Undetermined base (-, ?) content", caption = ptit) +
  theme_bw() +
  theme(legend.position = "top")

cat("done\n")

## Write txt
ext <- paste0(".", rev(unlist(strsplit(args[1], split = "[.]")))[1], "$")
if (writetxt) {
  tfile <- paste0(gsub(ext, "", basename(args[1])), ".stats.txt")
  write.table(d.plot, file = tfile, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
}

## Plot histograms
cat("writing pdf...")

pfile <- paste0(gsub(ext, "", basename(args[1])), ".stats.pdf")
if (file.exists(pfile)) system(paste("mv", pfile, gsub(".pdf$", ".bak.pdf", pfile)))
pdf(file = pfile, width = 6+length(args), height = 6+length(args))
print(p1.1) # bar totlen
#print(p1.2) # density
#print(p1.3) # density log
print(p1.4) # violin
print(p1.5) # violin log
#print(p2.1) # density
print(p2.2) # violin
print(p3.1) # scatter, density2D
print(p3.2) # scatter, density2D
print(p4) # boxplot
print(p5) # boxplot
print(p6) # boxplot

graphics.off()
cat("done\n")
print(Sys.time())
