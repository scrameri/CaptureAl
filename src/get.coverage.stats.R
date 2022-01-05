#!/cluster/apps/r/3.6.0_openblas/x86_64/bin/Rscript

## Usage: get.coverage.stats.R <bam> <refseqs> <ppaired=TRUE> <allatonce =FALSE>

## Load libraries
suppressPackageStartupMessages(library(ape)) # read.FASTA

## Author: simon.crameri@env.ethz.ch, Mar 2020

## Get arguments
args <- commandArgs(trailingOnly=TRUE)

## Check arguments
if (! length(args) %in% c(2:4)) {
  stop("2 arguments needed (4 taken):
       REQUIRED
       1) <bam|CHR>:           path to .bam file for which mapped length and average read coverage needs to be calculated.
       2) <refseqs|CHR>:       path to reference .fasta (used to map reads to). Will stop if some region IDs in <bam> do not match.
       
       OPTIONAL
       3) <ppaired|BOOLEAN>:   if TRUE, will use the 'samtools mpileup' option to determine regions (only regions with properly 
                               paired reads), if FALSE, will use the 'samtools mpileup -A' option to determine regions 
                               (all regions with any mapped reads) [DEFAULT: TRUE].
       4) <allatonce|BOOELAN>: if TRUE, will read coverage data for all regions in the .bam file to memory. This is quite
                               fast, but will use A LOT of memory, especially if this script is executed in parallel. 
                               If FALSE, coverage data is read to memory for each region consecutively, which saves memory at the
                               expense of computation time [DEFAULT: FALSE].",
       call.=FALSE)
}

## Set arguments
bam <- as.character(args[1])
refseqs <- as.character(args[2])

ppaired <- as.character(args[3])
if (is.na(ppaired)) ppaired <- TRUE
allatonce <- as.character(args[4])
if (is.na(allatonce)) allatonce <- FALSE

## Set arguments (for debugging)
# bam <- "depth.bam"
# refseqs <- "test.fasta"
# ppaired <- TRUE
# allatonce <- TRUE

# t1 <- Sys.time()
# paste0("Starting time: ", t1)

## Additional arguments
get.boolean <- function(x, truestrings = c("T","TRUE", "t","true", "True"), 
                        falsestrings = c("F","FALSE","f","false","False")) {
  if (x %in% truestrings) {
    return(TRUE)
  } else {
    if (x %in% falsestrings) {
      return(FALSE)
    } else stop("<", x, "> cannot be interpreted. Please specify one of ", 
                paste(c(truestrings, falsestrings), collapse = ", "), ".")
  }
}
ppaired <- get.boolean(ppaired) ; pchar <- ifelse(ppaired, "PROPERLY-PAIRED reads only", "ANY mapped read")
allatonce <- get.boolean(allatonce)
suffix <- ".coverage.txt"
verbose <- FALSE

###########################################3

## Check arguments
stopifnot(file.exists(bam), file.exists(refseqs), is.logical(ppaired), is.logical(allatonce))

## Read reference
ref <- ape::read.FASTA(refseqs) ; nreg.all <- length(ref)

## Get region names from .bam
# region names can be retrieved from a .bam using 'samtools mpileup':
# 'samtools mpileup - -A' returns all regions with any mapped reads (but no regions without any mapped read)
# 'samtools mpileup -'    returns only regions with properly paired mapped reads
if (verbose) cat(paste0("reading regions with ", pchar, "..."))
if (ppaired) {
  cmd <- paste("samtools view -b", bam, "2>/dev/null | samtools mpileup -  2>/dev/null | awk '{print $1}' | sort | uniq 2>/dev/null")
} else {
  cmd <- paste("samtools view -b", bam, "2>/dev/null | samtools mpileup - -A 2>/dev/null | awk '{print $1}' | sort | uniq 2>/dev/null")
}
regs <- sort(system(command = cmd, intern = TRUE)) # stores regions
nreg <- length(regs)
cat(paste0(nreg, " / ", nreg.all , " (", round(100*nreg/nreg.all,2), "%) region(s) found\n"))
ref <- ref[regs] # only considered regions, and in the same order

## Check conformity with reference .fasta
if (!all(regs %in% names(ref))) {
  mismatch <- as.character(regs[which(!regs %in% names(ref))])
  stop(length(mismatch), " regions in <", bam, "> not found in <", refseqs, ">:\n",
       paste(mismatch, collapse = ", "))
}

## Get Coverage stats
# get depth using 'genomeCoverageBed' (alias for 'bedtools genomecov')
# genomeCoverageBed documentation: http://gensoft.pasteur.fr/docs/bedtools/2.19.1/content/tools/genomecov.html
# genomeCoverageBed vs. samtools depth: https://www.biostars.org/p/67579/
# -> samtools depth skips secondary alignments and aberrand read pairs, and will therefore result in 
#    lower depths reported in some cases.
if (allatonce) {
  ## Reads coverage data to memory for all regios in the bam file at once. Fast but can only be done if there is sufficient memory
  # especially if this script is run in parallel
  if (verbose) cat("reading coverages for ALL regions at once...\n")
  cmd <- paste("samtools view -b", bam, "2>/dev/null |  genomeCoverageBed  -ibam - -dz")
  l <- strsplit(system(command = cmd, intern = TRUE), split = "\t") # stores all coverage data
  ls <- data.frame(region = sapply(l, "[", 1), POS = as.numeric(sapply(l, "[", 2)), DEPTH = as.numeric(sapply(l, "[", 3)))
  
  ## Only considered regions
  ls <- ls[ls$region %in% regs,] ; ls$region <- droplevels(ls$region)

  ## Get Coverage stats
  if (verbose) cat("calculating average coverages...\n")
  
  # length in bam
  linbam <- tapply(ls$POS, ls$region, length)
  
  # length in refseq
  stopifnot(all.equal(levels(ls$region), names(ref)))
  linref <- lengths(ref)
  
  # avg. coverage in bam
  cinbam <- tapply(ls$DEPTH, ls$region, mean)
  
  # avg. coverage in refseq
  bases <- tapply(ls$DEPTH, ls$region, sum)
  stopifnot(all.equal(names(bases), names(linref)))
  cinref <- bases/linref
  
  # bind
  dd <- data.frame(region = names(linbam), linbam = linbam, linref = linref, cinbam = cinbam, cinref = cinref)
} else {
  if (verbose) cat("reading and calculating coverages for region\n")
  dd <- array(data = NA, dim = c(0, 5), dimnames = list(NULL, c("region","linbam","linref","cinbam","cinref")))
  for (region in regs) { # names(ref)
    cat(which(regs == region), " / ", nreg, "\r")
    
    cmd <- paste("samtools view -b", bam, region, "2>/dev/null |  genomeCoverageBed  -ibam - -dz")
    l <- strsplit(system(command = cmd, intern = TRUE), split = "\t") # stores coverage data for one region
    ls <- data.frame(region = sapply(l, "[", 1), POS = sapply(l, "[", 2), DEPTH = sapply(l, "[", 3))
    
    # length in bam
    linbam <- nrow(ls)
    
    # length in refseq
    linref <- lengths(ref[region])
    
    # avg. coverage in bam
    dp <- as.numeric(as.character(ls$DEPTH))
    cinbam <- mean(dp)
    
    # avg. coverage in refseq
    cinref <- sum(dp)/linref
    
    # bind
    dd <- rbind(dd, data.frame(region = region, linbam = linbam, linref = linref, cinbam = cinbam, cinref = cinref))
  }
}

## Write output
if (verbose) cat("Writing output...\n")
ext <- paste0(".", tools::file_ext(bam))
write.table(dd, file = gsub(".bwa-mem.sorted", "", gsub(ext, suffix, bam)), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# t2 <- Sys.time()
# paste0("Finish time: ", t2)
