#!/usr/bin/Rscript

## Load dependencies
library(vcfR)
library(ape)
library(adegenet)

## Get arguments
args = commandArgs(trailingOnly=TRUE)
t1 <- Sys.time()

## Check arguments
if (length(args) != 2) {
  stop("Two arguments needed: 
       REQUIRED
       1) <rdapaths|CHR>: file with paths to .rda files with genind objects ; 
       2) <fastafile|CHR>: path to reference .fasta", call.=FALSE)
}

## Set arguments
rdapaths = as.character(args[1])
ref = as.character(args[2])

## Check arguments
stopifnot(file.exists(rdapaths), file.exists(ref))

## Read paths file
f <- readLines(rdapaths)

## Read reference .fasta file
dna <- read.FASTA(ref)
dna <- dna[sort(names(dna))]

## Get regions
# start <- 1
# regions <- names(dna)[start:length(dna)]
regions <- tools::file_path_sans_ext(basename(f))
stopifnot(all(regions %in% names(dna)))
regions <- names(dna)[names(dna) %in% regions]

## Bind gi objects
count <- 0
for (region in regions) {
  
  cat(which(regions == region), "/", length(regions), "\r")
  
  # load gi object
  load(f[which(regions == region)]) # object name = region name
  dd <- get(region) ; rm(list=region)
  
  if (!exists("gi")) {
    gi <- dd$gi
  } else {
    stopifnot(all.equal(indNames(gi), indNames(dd$gi)))
    dp <- gi@other$DP
    qu <- gi@other$QUAL
    gi <- genind(tab = cbind(gi@tab, dd$gi@tab), pop = gi@pop, 
                 prevcall = gi@call, ploidy = gi@ploidy, type = gi@type, 
                 strata = gi@strata, hierarchy = gi@hierarchy)
    gi@other$DP <- c(dp, dd$gi@other$DP)
    gi@other$QUAL <- c(qu, dd$gi@other$QUAL)
  }
  
  # save gi (every 100th)
  count <- count + 1
  if (count %in% seq(100,length(dna),by=100)) {
    save(gi, file = paste0("gi_", nInd(gi), "_", nLoc(gi), ".rda"))
  }
}

## Save final binded genind object
save(gi, file = paste0("gi_", nInd(gi), "_", nLoc(gi), ".rda"))

## Finish
t2 <- Sys.time()
cat("\nFinish time:", as.character(t2), "\n")
print(t2-t1)

