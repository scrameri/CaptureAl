#!/cluster/apps/r/3.6.0_openblas/x86_64/bin/Rscript

## Usage: get.taxonstats.R <fasta alignment> <OPTIONAL: name of output file>

## Load libraries
suppressPackageStartupMessages(library(ape)) # read.fasta, GC.content
# suppressPackageStartupMessages(library(pegas)) # site.spectrum
# suppressPackageStartupMessages(library(ggplot2)) # ggplot

## Author: simon.crameri@env.ethz.ch, Mar 2020

## Get arguments
args <- commandArgs(trailingOnly=TRUE)

## Check arguments
options(warning.length = 5000L)
if (! length(args) %in% c(1,2)) {
  stop("1 arguments needed (2 taken):
       REQUIRED
       1) <ffile|CHR>:  path to FASTA alignment
       
       OPTIONAL
       2) <ofile|CHR>:  name of output file [DEFAULT: paste0(gsub(paste0(suffix, '$'), '', ffile), '.tstats')]",
       call.=FALSE)
}

## Set arguments
ffile <- as.character(args[1])

## Set arguments (for debugging)
# ffile <- "consDalbergia_NC_033812.1_LG_09_9957318_9957503_ID_5750.merged.all.aln.etr.itr.fasta"

## Additional arguments
suffix <- ".all.aln.etr.itr.fasta"
loc <- gsub(paste0(suffix, "$"), "", basename(ffile))
if (is.na(args[2])) ofile <- paste0(loc, ".tstats") else ofile <- as.character(args[2])

## Check arguments
stopifnot(file.exists(ffile))

## read FASTA
f <- as.matrix(ape::read.FASTA(ffile))

## get character matrix
df <- as.character(f)

# ## Site frequency spectrum
# sf <- site.spectrum(f, folded = F)
# plot(sf)

## get major allele
get.maj <- function(x, NA.char = c("-","N")) {
  names(which.max(table(x[!x %in% NA.char])))
}
maj <- apply(df, MARGIN = 2, FUN = get.maj)

## get SNP
get.SNP <- function(x, maj, NA.char = c("-","N")) {
  x[x %in% NA.char] <- maj[x %in% NA.char]
  ref <- maj[x != maj]
  alt <- x[x != maj]
  paste0(ref, ">", alt)
}
snp <- apply(df, MARGIN = 1, FUN = get.SNP, maj = maj)

## get transitions and transversions and missingness
get.trans.transv <- function(x) {
  # Rambaut et al. 2008 MBE: 
  # The dominant form of postmortem single-base modifications are C to T changes, 
  # which will be indistinguishable from G to A changes because damage can occur on either strand
  ag <- sum(x %in% c("a>g","t>c")) # trans
  ct <- sum(x %in% c("c>t","g>a")) # trans
  ac <- sum(x %in% c("a>c","t>g")) # transv
  at <- sum(x %in% c("a>t","t>a")) # transv
  gc <- sum(x %in% c("g>c","c>g")) # transv
  gt <- sum(x %in% c("g>t","c>a")) # transv
  c("trans" = sum(ag, ct), "transv" = sum(ac, at, gc, gt), ag = ag, ct = ct, ac = ac, at = at, gc = gc, gt = gt)
}

get.NA <- function(x, NA.char = c("-","N")) {
  sum(x %in% NA.char)
}

tr <- data.frame(ID = names(snp), t(sapply(snp, get.trans.transv)), 
                 mis = apply(df, MARGIN = 1, FUN = get.NA), 
                 len = apply(df, MARGIN = 1, FUN = length),
                 loc = loc)

## Write results
write.table(tr, file = ofile, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# ## Plot
# ggplot(tr) +
#   geom_bar(aes(x = ID, y = trans), stat="identity") +
#   coord_flip() +
#   theme_bw()
# 
# ggplot(tr) +
#   geom_bar(aes(x = ID, y = transv), stat="identity") +
#   coord_flip() +
#   theme_bw()
# 
# ggplot(tr) +
#   geom_bar(aes(x = ID, y = ct), stat="identity") +
#   coord_flip() +
#   theme_bw()
# 
# ggplot(tr) +
#   geom_bar(aes(x = ID, y = mis), stat="identity") +
#   coord_flip() +
#   theme_bw()

