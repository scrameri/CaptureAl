############################################################
### PCA and NJ of Dalbergia population set .vcf variants ###
############################################################

# simon.crameri@usys.ethz.ch, Mar 2022

# libraries
library(vcfR)
library(adegenet)
source("HELPERFUNCTIONS.R")

## Arguments
max.missing <- 0.05     # maximum fraction of missing genotypes for a SNP to be kept
min.mac <- 3            # minimum minor allele count for a SNP to be kept (use min.mac = 2 to remove singletons)
biallelic.only <- FALSE # if TRUE, only keep SNPs with exactly 2 allels
expected.only <- TRUE   # if TRUE, only keep SNPs with all(alleles %in% exp.char)
exp.char = c("A","T","G","C") # expected alleles (if expected.only = TRUE, any locus with other alleles is dropped)
n.snps <- 3             # number of SNPs randomly selected per target region
region.split <- "_"     # character used to split locNames(genind), which should contain the target region
region.index <- 1       # integer index of split string that represents the target region
na.int <- -9            # value used to code missing genotypes in STRUCTURE input file

## Create genind object from VCF file of filtered SNPs
# read vcf
input <- "data/Dalbergia_51_2396_116500.filtered.vcf" # 116,500 variants after variant filtering
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

# Reorder individuals and loci
gi.all <- gi.all[order(indNames(gi.all))]


## Filter SNPs
# subset SNPs for missingness
cat("subsetting SNPs...\n")
gistru <- filter.gi(gi = gi.all, biallelic.only = biallelic.only, expected.only = expected.only,
                    exp.char = exp.char, min.mac = min.mac, max.missing = max.missing)

# exclude monomorphic SNPs 
locvar <- apply(gistru@tab, 2, var, na.rm = TRUE)
gistru <- gistru[loc = unique(sapply(strsplit(names(locvar[locvar > 0]), split = "[.]"), "[", 1))] # removes loci with 2 alleles, all heterozygotes

# subset N SNPs 
set.seed(1503) # this seed may not result in the same random subset of SNPs as the one used in the article
gistrusub <- gi.filter.per.region(gi = gistru, n = n.snps, split = region.split, split.index = region.index)

# write in <onerowperind> format (2 successive columns represent the maternal and paternal haplotype, 1 = A, 2 = T, 3 = G, 4 = C)
df <- genind2df(gistrusub, oneColPerAll = TRUE, usepop = FALSE)
dfnum <- apply(df, 2, function(x) {as.numeric(factor(x, levels = exp.char))}) # STRUCTURE needs integers
rownames(dfnum) <- rownames(df)
d.stru <- apply(dfnum, 2, function(x) {x[is.na(x)] <- na.int ; return(x)})

ncol(d.stru)/2 # 7156 SNPs

write.table(d.stru, file = "data/Monticola.51.7156.Dalbergia_CH1.3_51_116500.str",
            row.names = TRUE, sep = " ", col.names = F, quote = F)



# use the .str file to run STRUCTURE analyses


## Plot STRUCTURE results
# get data
input <- "Chapter1.3_mapsnp-2396_51/results"
# load("data/d.uce_v1.rda") # hidden
# load("data/data.rda") # hidden coordinates

# set arguments
simple <- FALSE # TRUE does not work with Clumpak, but FALSE does (in this case Structure results need not to contain popoulation labels)
Kmax <- 10
R <- 10
write <- TRUE

##########################################################

# get helperfunctions
clean <- function(string) {
  repeat{
    string <- gsub("  ", " ", string)
    if (length(grep("  ", string)) == 0) {
      string <- gsub("^ ", "", string)
      string <- gsub(" $", "", string)
      break
    } 
  }
  return(string)
}

ReadAndClean <- function(file) {
  f <- suppressWarnings(readLines(file))
  start <- grep("^Inferred ancestry", f)[1]+2
  end <- grep("^Estimated Allele Frequencies", f)[1]-3
  fsub <- f[start:end]
  fdat <- sapply(fsub, clean, USE.NAMES = FALSE)
  return(fdat)
}

# create output dir
ofolder <- file.path(dirname(input),"Clumpakinput")
if (!dir.exists(ofolder)) dir.create(ofolder)

# get pop identifyer
files <- list.files(input, full.names = T)
first <- ReadAndClean(files[1])
ids <- sapply(strsplit(first, split = " "), "[", 2)

# sort (using non-public coordinates)
# d.uce$ID.Miseq.cleaned <- gsub("-","_",d.uce$ID.Miseq)
# rownames(d.uce) <- d.uce$ID.Miseq.cleaned
# rownames(data) <- data$ID_Lab
# 
# idlab <- d.uce[ids,c("ID_Lab")]
# 
# dlab <- data.frame(ID.Miseq.cleaned = d.uce[ids,c("ID.Miseq.cleaned")], data[idlab,c("ID_Lab","Species","LatitudeDecimal","LongitudeDecimal","MinimumElevation")])
# 
# # sort 
# xlat <- suppressWarnings(as.numeric(as.character(dlab$LatitudeDecimal)))
# xlat[is.na(xlat)] <- mean(xlat, na.rm = T)
# dlab <- dlab[order(as.character(dlab$Species), -xlat),]

# sort
labs <- c("SH0387","SH0412","SH0411","SH0377","SH0455","SH0479","SH0481", 
          "SH0299","SH0300","SH0290","SH0349","SH0394","SH0390","SH0424",
          "SH0401","SH0486","SH0482","H0060", "H0059", "RBE2282","H0082",
          "SH0565","SH0573","SH0591","SH0560","SH0564","SH0611","SH0609",
          "SH0612","RAF0014","SH0335","SH0039","SH0051","RBE2479-1","SH0226",   
          "H0074", "SH0015","CR7303","CR7296","RZK7699","CR7314","RBE2236",
          "SH0598","H0043","EME0024","EME0018","RZK8028","SFR0253","SFR0252",
          "RZK8037","SFR0258")

inds <- c("SH387","SH412","SH411","SH377_S3","SH455_S6","S3_S11","SH481",
          "SH299","S8_S16","SH290","SH349","SH394","S2_S10","SH424",
          "SH401","SH486","S7_S15","H060","H0059_R_S13","RBE2282","H082",
          "SH565_S21","SH573","SH591","S1_S9","SH564","SH611","609_65_S1",
          "SH612","RAF014","SH335_S9","SH039","SH051","RBE2479_1","C1_S6",
          "H074","SH015_S5","CR7303","CR7296","RZK7699_S3","CR7314","RBE2236_S7",
          "P8_22_S24","H043_S19","EME024","EME018","RZK8028","SFR253","SFR252",  
          "RZK8037","SFR258")

dlab <- data.frame(ID.Miseq.cleaned = inds, ID_Lab = labs,
           Species = factor(c(rep("sp.B",7),rep("monticola",22),rep("orientalis",22))))


# write mapfile
dmap <- data.frame(INDEX = as.numeric(dlab$Species), LABEL = dlab$Species)
if (write) write.table(dmap$INDEX, file = file.path(dirname(ofolder), "mapfile.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ")

dlabs <- dmap[!duplicated(dmap),]
dlabs <- dlabs[order(as.character(dlabs$LABEL)),]
# if (write) write.table(dlabs, file = file.path(dirname(ofolder), "labelfile.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ")

## Write ordered STRUCTURE results files to use as CLUMPAK input files
for (file in files) {
  dat <- data.frame(t(sapply(strsplit(ReadAndClean(file), split = " "), function(x) matrix(x))))
  K <- ncol(dat)-4
  names(dat) <- c("ROW","ID.Miseq.cleaned","XMIS","DOT",paste0("Q",1:K))
  
  # sort
  dat <- dat[dlab$ID.Miseq.cleaned,]
  rownames(dat) <- NULL
  
  if (simple) {
    df <- dat[,paste0("Q", 1:K)]
    write.table(df, file = file.path(ofolder, paste0(basename(file), ".txt")), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ")
  } else {
    dat$POP <- as.numeric(factor(dlab$Species))
    df <- dat[c("ID.Miseq.cleaned","XMIS","POP","DOT",paste0("Q", 1:K))]
    write.table(df, file = file.path(ofolder, paste0(basename(file), ".txt")), quote = FALSE, row.names = TRUE, col.names = FALSE, sep = " ")
  }
}
