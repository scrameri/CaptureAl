############################################################
### PCA and NJ of Dalbergia population set .vcf variants ###
############################################################

# simon.crameri@usys.ethz.ch, Mar 2022

# libraries
library(vcfR) # v1.12.0
library(adegenet) # v2.1.5
library(ape) # v5.6.1
library(poppr) # v2.9.3
library(ggplot2) # v3.3.5, part of tidyverse package
source("HELPERFUNCTIONS.R")

## Create genind object from VCF file of filtered SNPs
# read vcf
input <- "Dalbergia_51_2396_116500.filtered.vcf" # 116,500 variants after variant filtering
dd <- read.vcfR(file = input)

# heterozygosity
dhet <- is_het(extract.gt(dd), na_is_false = FALSE)
phet <- apply(dhet, 2, function(x) sum(x, na.rm = T)/length(x[!is.na(x)]))
summary(phet)


## Convert vcfR to genind
# Rename LOCNAME (shorten and remove .)
# gi.all <- vcfR2genind(dd, return.alleles = T) # error due to "." in LOCNAME
dd@fix[,"CHROM"] <- gsub("[.]", "", sub(".*_ID_(.*$)", "\\1", dd@fix[,"CHROM"])) # sub("merged.merged", "merged", 
dd <- addID(dd)
gi.all <- vcfR2genind(dd, return.alleles = T)

# Rename samples (remove _L001 suffix and substitute - with _)
indNames(gi.all) <- gsub("-", "_", sub("_L001$", "", indNames(gi.all)))


## Filter SNPs
# subset biallelic SNPs
gi.bi <- gi.all[loc = which(gi@loc.n.all == 2), drop = TRUE]

# subset bi-allelic subset with no missing data
loc.nomiss <- unique(gsub("[.].*$", "", names(which(!apply(gi.bi@tab, 2, anyNA)))))
length(loc.nomiss) # 60204 loci with zero missingness
gi.pca <- gi.bi[loc = loc.nomiss, drop = TRUE]

## convert genind to genlight
# NOTE: vcfR2genlight(dd) is not equivalent because dd contains SNPs with missing genotypes
gl <- genind2genlight(gi.pca)


## PCA on genlight object (much faster than on genind object)
gl.pca <- glPca(gl, center = TRUE, scale = FALSE, parallel = TRUE, n.cores = 4, nf = 4)

# percent of explained variance
gl.pca$perc.var <- round(100*gl.pca$eig/sum(gl.pca$eig), 2)

# assign species to clusters
gl@other$cluster <- cutree(hclust(dist(gl.pca$scores[,1:2])), k = 3)
gl@pop <- factor(NA, levels = c("monticola","orientalis","spB"))
gl@pop[gl@other$cluster %in% gl@other$cluster[grep("SH611",indNames(gl))]] <- "monticola"
gl@pop[gl@other$cluster %in% gl@other$cluster[grep("EME018",indNames(gl))]] <- "orientalis"
gl@pop[gl@other$cluster %in% gl@other$cluster[grep("SH411",indNames(gl))]] <- "spB"
gi.pca@pop <- gl@pop
table(gl@pop)

# plot pca (note that PC2 was flipped to better match the NJ tree)
p.pca <- ggplot(data.frame(gl.pca$scores, taxon = gl@pop), aes(PC1, -PC2, color = taxon)) +
  geom_point() +
  labs(x = paste0("PC 1 (", gl.pca$perc.var[1], "%)"),
       y = paste0("PC 2 (", gl.pca$perc.var[2], "%)")) +
  scale_color_manual(values = c("#0099E6", "#FF6600", "#33004D")) +
  theme_bw()

pdf("PCA_populationset.pdf")
print(p.pca)
graphics.off()

## NJ tree
gi.dist <- nei.dist(gi.pca)
gi.nj <- nj(gi.dist)

# plot NJ tree
pdf("NJ_populationset.pdf")
plot(gi.nj, "unrooted", tip.col = c("#0099E6", "#FF6600", "#33004D")[gi.pca@pop])
graphics.off()

