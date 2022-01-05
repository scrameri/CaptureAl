#!/cluster/apps/r/3.6.0_openblas/x86_64/bin/Rscript

## Usage: collect.exonerate.stats.R <sample file> <OPT: dir> <OPT: suffix>

## Value: Collects sample alignment statistics

## Author: simon.crameri@env.ethz.ch, Mar 2020

## Load required library
suppressPackageStartupMessages(library(data.table))

## Get arguments
args = commandArgs(trailingOnly = TRUE)

## Check arguments
if (!length(args) %in% c(1:3)) {
  stop("1 argument needed (3 taken): 
       REQUIRED:
       1) <sfile|CHR>:  path to sample file containing paths to mapping directories ;

       OPTIONAL (if one is given, all must be given in this order):
       2) <dir|CHR>:    folder with subdirectories for each sample [DEFAULT: current directory ]
       3) <suffix|CHR>: sample suffix present in <sfile> and absent in <dir> subdirectories  [DEFAULT: '']", call.=FALSE)
}

## Set parameters
sfile = as.character(args[1])
dir = as.character(args[2])
if (is.na(dir)) dir = getwd()
suffix = as.character(args[3])
if (is.na(suffix)) suffix = ""

## Set parameters (for debugging)
# sfile = "samples.txt"
# dir = "best.contigs.test"
# suffix = ""

## Check parameters
stopifnot(file.exists(sfile),
          dir.exists(dir))

## Additional parameters
XandMore = 3 # will summarize for 0, 1, 2, ..., >=XandMore contigs

tabname = "contig-table.txt"   # expects a header
ncontigname = "ncontigs"       # column name in <tabname> where number of contigs is stored
nlociname = "nloci"            # column name in <tabname> where number of loci is stored
outfile = "loci_contignumbers.txt" # name of output file

statname = "contig-stats.txt" # expects a header
statfields = c("ID","target","locid","npassed","nfailed","ncontigs","bestscore","bestscore.norm","bestlength","taln") # complete list of expected header
# numerics = c("npassed","nfailed","ncontigs","bestscore","bestscore.norm","bestlength","taln") # numeric columns in <statname>
outfile2 = "loci_stats.txt"   # name of second output file

rangename = "contig-ranges.txt" # expects a header
addnumerics = c("qstart","qend","tstart","tend","score") # numeric columns in <rangename>
outfile3 = "loci_ranges.txt"  # name of third output file

verbose = TRUE # prints loop progress

## Read sfile
samples <- as.character(read.delim(sfile, header = F)$V1)

## Remove old output files
move <- function(file) {
  if (file.exists(file)) system(paste("mv", file, paste0(file, ".bak")))
}
invisible(sapply(c(outfile,outfile2,outfile3), move)) # moves any existing output file(s) to *.bak

## Loop through files
dcounts <- data.frame(array(data = NA, dim = c(0, XandMore + 2), dimnames = list(NULL, c("ID", 0:XandMore))), check.names = FALSE)

# get contig counts and stats
for (sample in samples) {
  if (verbose) cat(sample, " ")
  tabfile <- file.path(dir, sample, tabname)
  statfile <- file.path(dir, sample, statname)
  rangefile <- file.path(dir, sample, rangename)
  
  stopifnot(file.exists(tabfile),
            file.exists(statfile))
  
  # tab
  dtab <- read.delim(tabfile)
  dcount <- data.frame(array(data = NA, dim = c(1, XandMore + 2), dimnames = list(NULL, c("ID", 0:XandMore))), check.names = FALSE)
  names(dcount)[ncol(dcount)] <- paste0(">=", XandMore)
  dcount[,"ID"] <- gsub(paste0(suffix, "$"), "", sample)
  
  for (j in 0:(XandMore-1)) {
    n = as.numeric(dtab[dtab[,ncontigname] == j, nlociname])
    if (length(n) == 1) dcount[,as.character(j)] <- n else dcount[,as.character(j)] <- 0
  }
  
  dcount[,ncol(dcount)] <- sum(dtab[dtab[,ncontigname] >= XandMore, nlociname])
  stopifnot(sum(dcount[1,-1]) == sum(dtab[,nlociname]))
  dcounts <- rbind(dcounts, dcount)
  
  # stat
  dstat <- fread(statfile)
  if (any(!statfields %in% names(dstat))) {
    addfields <- statfields[which(!statfields %in% names(dstat))]
    dadd <- array(NA, dim = c(nrow(dstat), length(addfields)), dimnames = list(NULL, addfields))
    dstat <- cbind(dstat, dadd)
  }
  dstat <- dstat[,..statfields]
  if (!file.exists(outfile2)) {
    fwrite(dstat, file = outfile2, row.names = F, col.names = T, sep = "\t", quote = F, na = "NA")
  } else {
    fwrite(dstat, file = outfile2, row.names = F, col.names = F, sep = "\t", quote = F, append = TRUE, na = "NA")
  }
  
  # alignment range
  if (file.exists(rangefile)) {
    cat("\n")
    drange <- fread(rangefile)
    if (! "ID" %in% names(drange)) drange <- data.table(ID = rep(sample, nrow(drange)), drange)
    if (!file.exists(outfile3)) {
      fwrite(drange, file = outfile3, row.names = F, col.names = T, sep = "\t", quote = F, na = "NA")
    } else {
      fwrite(drange, file = outfile3, row.names = F, col.names = F, sep = "\t", quote = F, append = TRUE, na = "NA")
    }
  } else {
    if (verbose) cat(" - ", rangefile, "does not exist, skipping\n")
  }
}

## Write binded results
write.table(dcounts, file = outfile, row.names = F, col.names = T, sep = "\t", quote = F)

