#!/cluster/apps/r/3.6.0_openblas/x86_64/bin/Rscript

## Usage: collect.coverage.stats.R <sample file> <Q> <dir>

## Needs: ggplot2, tidyr

## Value: Collects sample .bam coverage statistics (missing sample x region combinations are expanded)

## Author: simon.crameri@env.ethz.ch, Dec 2021

## Get arguments
args = commandArgs(trailingOnly = TRUE)

## Check arguments
if (!length(args) %in% c(1:4)) {
  stop("1 argument needed (4 taken): 
       REQUIRED:
       1) <sfile|CHR>:  path to sample file containing paths to sample directories with mapping results. No header expected. ;

       OPTIONAL (must be given in this order):
       2) <Q|NUM>:      mapping quality threshold. Will look for SAMPLE.Q${Q}.coverage.txt [DEFAULT: 10]
       3) <maxcov|NUM>: maximum coverage threshold. Will look for SAMPLE.Q${Q}.nodup.cov${maxcov}.coverage.txt [DEFAULT: 500]
       4) <dir|CHR>:    folder with subdirectories for each sample [DEFAULT: current directory ]",
       call=FALSE)
}

## Set parameters
sfile = as.character(args[1])

Q = as.character(args[2])
if (is.na(Q)) Q = 10
maxcov = as.character(args[3])
if (is.na(maxcov)) maxcov = 500
dir = as.character(args[4])
if (is.na(dir)) dir = getwd()


## Set parameters (for debugging)
# sfile = "samples.txt"
# Q = 10
# maxcov = 500
# dir = "."

## Check parameters
stopifnot(file.exists(sfile),
          Q >= 0,
          maxcov >= 0,
          dir.exists(dir))

## Additional parameters
# input paths (path to coverage stats file)
# will look for <covfileQ> in each sample subdirectory. 
# <SAMPLE>, <Q> and <maxcov> are part of the file path and will be handled using regular expressions.
covfileQ <- "__SAMPLE__.Q__MAPQ__.coverage.txt"
covfileNODUP <- "__SAMPLE__.Q__MAPQ__.nodup.cov__MAXCOV__.coverage.txt"

# variables expected in <covfileQ>
regvar <- "region"    # variable name in <covfileQ> denoting <region id>
lvar1 <- "linref"     # variable name in <covfileQ> denoting <length of target (reference) sequence>
pvar1 <- "linbam"     # variable name in <covfileQ> denoting <length of mapped target (reference) sequence>
cvar2 <- "cinref"     # variable name in <covfileQ> denoting <average coverage of target (reference) sequence>
pvar2 <- "cinbam"     # variable name in <covfileQ> denoting <average coverage of mapped part of target (reference) sequence>

# variables created in merged data.frame
idvar <- "ID"         # variable name denoting <sample ID>

# output file
cov_Q = paste0("coverage_stats.Q", Q, ".txt") # name of output file
cov_NODUP = paste0("coverage_stats.Q", Q, ".nodup.cov", maxcov, ".txt") # name of output file

##########################################################################################

## Helperfunctions
move <- function(file) {
  if (file.exists(file)) system(paste("mv", file, paste0(file, ".bak")))
}
read.cov <- function(covstat, verbose = FALSE) {
  covdata <- data.frame(array(NA, dim = c(0, 6), dimnames = list(NULL, c(idvar,regvar,pvar1,lvar1,pvar2,cvar2))))
  for (i in samples) {
    cat(which(samples == i), " / ", length(samples), "\r")
    p <- file.path(i, paste0(sub("__MAXCOV__", maxcov, sub("__MAPQ__", Q, sub("__SAMPLE__", i, covstat)))))
    if (!file.exists(p)) {
      message("file <", p, "> not found, did you specify the correct path (line 54-55 in script)?")
    } else {
      dd <- read.delim(p, header = TRUE, stringsAsFactors = FALSE)
      dd <- data.frame(ID = rep(i, nrow(dd)), dd, stringsAsFactors = FALSE)
      if (!regvar %in% names(dd)) {stop("variable <", regvar, "> not found in <", p, ">, please check line 58 in script")}
      if (!lvar1 %in% names(dd))  {stop("variable <", lvar1,  "> not found in <", p, ">, please check line 59 in script")}
      if (!pvar1 %in% names(dd))  {stop("variable <", pvar1,  "> not found in <", p, ">, please check line 60 in script")}
      if (!cvar2 %in% names(dd))  {stop("variable <", cvar2,  "> not found in <", p, ">, please check line 61 in script")}
      if (!pvar2 %in% names(dd))  {stop("variable <", pvar2,  "> not found in <", p, ">, please check line 62 in script")}
      covdata <- rbind(covdata, dd)
    }
  }
  names(covdata)[1] <- c(idvar)
  
  # fill up regions with no coverage (requires tidyr and ggplot2)
  if (nrow(covdata) != length(unique(covdata[,idvar]))*length(unique(covdata[,regvar]))) {
    if (verbose) cat("\n\nexpanding sample x region combinations...\n")

    covdata[,idvar] <- factor(covdata[,idvar])
    covdata[,regvar] <- factor(covdata[,regvar])
    covdata <- data.frame(tidyr::complete(covdata, !! ggplot2::sym(idvar), !! ggplot2::sym(regvar)))
  }

  # Sort coverage stats
  covdata <- covdata[order(covdata[,regvar], covdata[,idvar]),]
  return(covdata)
}

##########################################################################################

## Read sfile
samples <- as.character(read.delim(sfile, header = FALSE)$V1)

## Change to mapping directory
setwd(dir)

## Move existing oupt file
#invisible(sapply(c(cov_Q, cov_NODUP), move)) # moves any existing output file(s) to *.bak

## Read coverage data
cat(paste0("collecting / expanding coverage stats (Q", Q, ")...\n"))
covdataQ <- read.cov(covfileQ)
cat(paste0("collecting / expanding coverage stats (Q", Q, ".nodup.cov", maxcov, ")...\n"))
covdataNODUP <- read.cov(covfileNODUP)

## Write results
cat("\n\nwriting coverage stats...\n")
write.table(covdataQ, file = cov_Q, row.names = F, col.names = T, sep = "\t", quote = F)
write.table(covdataNODUP, file = cov_NODUP, row.names = F, col.names = T, sep = "\t", quote = F)
