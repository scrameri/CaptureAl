#!/gdc_home4/scrameri/bin/Rscript

## Usage: collect.coverage.stats.R <sample file> <mapQ> <dir>

## Value: Collects sample .bam coverage statistics (missing sample x region combinations are expanded)

## Author: simon.crameri@env.ethz.ch, Mar 2020

## Get arguments
args = commandArgs(trailingOnly = TRUE)

## Check arguments
if (!length(args) %in% c(1:3)) {
  stop("1 argument needed (3 taken): 
       REQUIRED:
       1) <sfile|CHR>:  path to sample file containing paths to sample directories with mapping results. No header expected. ;

       OPTIONAL (if one is given, all must be given in this order):
       2) <mapQ|CHR>:   mapping quality. Will look for ${SAMPLE}.bwa-mem.mapped.Q${MAPQ}.sorted.bam.coverage.txt [DEFAULT: 10]
       3) <dir|CHR>:    folder with subdirectories for each sample [DEFAULT: current directory ]",
       call=F)
}

## Set parameters
sfile = as.character(args[1])

mapQ = as.character(args[2])
if (is.na(mapQ)) mapQ = 10
dir = as.character(args[3])
if (is.na(dir)) dir = getwd()


## Set parameters (for debugging)
# sfile = "samples.txt"
# dir = "."
# suffix = ""

## Check parameters
stopifnot(file.exists(sfile),
          dir.exists(dir))

## Additional parameters
# input paths (path to coverage stats file)
# will look for <covstatfile> in each sample subdirectory. 
# <SAMPLE> and <MAPQ> can be part of the path and will be handled using regular expressions.
covstatfile <- "SAMPLE.bwa-mem.mapped.QMAPQ.sorted.coverage.txt" 

# variables expected in <covstatfile>
regvar <- "region"    # variable name in <covstatfile> denoting <region id>
lvar1 <- "linref"     # variable name in <covstatfile> denoting <length of target (reference) sequence>
pvar1 <- "linbam"     # variable name in <covstatfile> denoting <length of mapped target (reference) sequence>
cvar2 <- "cinref"     # variable name in <covstatfile> denoting <average coverage of target (reference) sequence>
pvar2 <- "cinbam"     # variable name in <covstatfile> denoting <average coverage of mapped part of target (reference) sequence>

# variables created in merged data.frame
idvar <- "ID"                   # variable name denoting <sample ID>

# output file
outfile = "coverage_stats.txt" # name of output file

## Helperfunctions
move <- function(file) {
  if (file.exists(file)) system(paste("mv", file, paste0(file, ".bak")))
}

##########################################################################################

## Move existing oupt file
invisible(sapply(c(outfile), move)) # moves any existing output file(s) to *.bak

## Read sfile
samples <- as.character(read.delim(sfile, header = FALSE)$V1)

## Read coverage data
cat("reading coverage data...\n")
covdata <- data.frame(array(NA, dim = c(0, 6), dimnames = list(NULL, c(idvar,regvar,pvar1,lvar1,pvar2,cvar2))))
for (i in samples) {
  cat(which(samples == i), " / ", length(samples), "\r")
  p <- file.path(i, paste0(gsub("MAPQ", mapQ, gsub("SAMPLE", i, covstatfile))))
  if (!file.exists(p)) warning("file <", p, "> not found,\ndid you specify the correct path (line 46 in script)?")
  dd <- read.delim(p, header = TRUE, stringsAsFactors = FALSE)
  dd <- data.frame(ID = rep(i, nrow(dd)), dd, stringsAsFactors = FALSE)
  if (!regvar %in% names(dd)) {stop("variable <", regvar, "> not found in <", p, ">, please check line 79 in script")}
  if (!pvar1 %in% names(dd))  {stop("variable <", pvar1,  "> not found in <", p, ">, please check line 80 in script")}
  if (!lvar1 %in% names(dd))  {stop("variable <", lvar1,  "> not found in <", p, ">, please check line 81 in script")}
  if (!pvar2 %in% names(dd))  {stop("variable <", pvar2,  "> not found in <", p, ">, please check line 82 in script")}
  if (!cvar2 %in% names(dd))  {stop("variable <", cvar2,  "> not found in <", p, ">, please check line 83 in script")}
  covdata <- rbind(covdata, dd)
}
names(covdata)[1] <- c(idvar)

## fill up regions with no coverage (requires tidyr and ggplot2)
if (nrow(covdata) != length(unique(covdata[,idvar]))*length(unique(covdata[,regvar]))) {
  cat("\n\nexpanding sample x region combinations...\n")
  
  # Load libraries
  suppressPackageStartupMessages(library(ggplot2)) # sym
  suppressPackageStartupMessages(library(tidyr))  # complete

  covdata[,idvar] <- factor(covdata[,idvar])
  covdata[,regvar] <- factor(covdata[,regvar])
  covdata <- data.frame(tidyr::complete(covdata, !! sym(idvar), !! sym(regvar)))
}

## Sort coverage stats
covdata <- covdata[order(covdata[,regvar], covdata[,idvar]),]

## Write binded results
cat("\nwriting coverage data...\n")
write.table(covdata, file = outfile, row.names = F, col.names = T, sep = "\t", quote = F)
