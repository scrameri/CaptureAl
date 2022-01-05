#!/cluster/apps/r/3.6.0_openblas/x86_64/bin/Rscript

#######################################################
## FILTER CONTIGS BY NORMALIZED BEST EXONERATE SCORE ##
#######################################################

## Usage
# run script without arguments to see the arguments it takes

## Author
# simon.crameri@env.ethz.ch, Apr 2019

## Load required library

## Get arguments
args = commandArgs(trailingOnly = TRUE)

## Check arguments
options(warning.length = 2000L)
if (! length(args) %in% c(2, 3)) {
  stop("3 arguments taken (2 required): 
       REQUIRED
       1) <indir|CHR>: path to folder with exonerate results ; 
       2) <locfile|CHR>: path to file with exonerate stats ;
       
       OPTIONAL
       3) <min.normbestscore|NUM>: minimum normalized best exonerate score; any contig below that score will be filtered (i.e., treated as missing) [DEFAULT: 2]", 
       call.=FALSE)
}

## Set arguments
indir <- as.character(args[1])
locfile <- as.character(args[2])

min.normbestscore <- as.numeric(as.character(args[3]))
if (is.na(min.normbestscore)) min.normbestscore <- 2

## Set arguments (for debugging)
# indir <- "best.contigs.9.6555/EME001_S1_L001"
# locfile <- "loci_stats-0.33-2-2.txt"
# min.normbestscore <- 2

## Additional arguments
id <- basename(indir)                 # taxon ID
col.normbestscore <- "bestscore.norm" # column name in <locfile> pointing to normalized best exonerate scores
col.passed <- "lpassed"                # column name in <locfile> pointing to passed loci
col.id <- "ID"                        # column name in <locfile> pointing to taxon ID
col.loc <- "locus"                    # column name in <locfile> pointing to locus ID
fsuffix <- ".bestScore.fasta"         # contig file suffix in <indir>
resultempty <- "--------------------" # any filtered contig will be replaced by this string
oprefix <- "loci_rm_min.normbestscore_" # prefix of output filtering directory in <indir>
odir <- file.path(indir, paste0(oprefix, min.normbestscore)) # filtered contigs will be move to this directory
verbose <- TRUE                       # if TRUE, prints number of filtered contigs to screen

## Read locfile
dloc <- read.delim(locfile)

## Get passed loci of taxon ID <id> 
dpassed <- dloc[dloc[,col.id] == id & dloc[,col.passed] == 1,]

## Get contigs that are below <min.normbestscore>
dfilt <- dpassed[dpassed[,col.normbestscore] < min.normbestscore & !is.na(dpassed[,col.normbestscore]),, drop = FALSE]
nfilt <- nrow(dfilt)
nfiltloc <- length(unique(dfilt[,col.loc]))
nloc <- length(unique(dpassed[,col.loc]))

## Be verbose
if (verbose) cat(paste0("contigs below min. normalized best exonerate score ", min.normbestscore, " in ", id, ": ", nfilt, " / ", nloc, " (", round(100*nfiltloc/nloc,2), "%) loci.\n"))

## Move any previously filtered contigs back and remove previous loci_rm_min.normbestscore.below_<min.ncormbestscore> folder
filtdir <- list.files(indir, pattern = oprefix, full.names = TRUE)
stopifnot(length(filtdir) <= 1)

if (length(filtdir) == 1) {
  if (length(list.files(path = filtdir)) > 0) {
    cmd <- paste("mv", file.path(filtdir, "*"), indir)
    system(cmd)
  }
  if (length(list.files(path = filtdir)) == 0) system(paste("rm -rf", filtdir))
}

## Move filtered contigs to loci_rm_min.normbestscore.below_<min.ncormbestscore> folder
if (nfilt > 0) {
  if (!dir.exists(odir)) dir.create(odir)
  for (locus in dfilt[,col.loc]) {
    fname <- paste0(id, ".", locus, fsuffix)
    fasta <- file.path(indir, fname)
    stopifnot(file.exists(fasta))
    cmd <- paste("mv", fasta, odir)
    newfas <- c(readLines(fasta, n = 1), resultempty) # header and <resultempty>
    
    if (!file.exists(file.path(odir, fname))) {
      # move filtered contig to <odir>
      system(command = cmd)
      
      # create empty contig in <indir>
      writeLines(newfas, con = fasta)
    } else {
      warning(fname, " already in ", odir, "\n")
    }
  }
}

