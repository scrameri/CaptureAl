#!/cluster/apps/r/3.6.0_openblas/x86_64/bin/Rscript

#########################
## Get Exonerate Stats ##
#########################

## Usage
# run script without arguments to see the arguments it takes

## Author
# simon.crameri@env.ethz.ch, Mar 2019

## Load required library
suppressPackageStartupMessages(library(data.table))

## Get arguments
args = commandArgs(trailingOnly = TRUE)

## Check arguments
# options(warning.length = 2000L)
# options(width = 1000)
if (! length(args) %in% c(1)) {
  stop("1 argument needed:
       REQUIRED
       1) <folder|CHR>: path to exonerate results folder", 
       call.=FALSE)
}

## Set arguments
folder <- as.character(args[1])

## Set arguments (for debugging)
#folder <- "best.contigs.test/sample/"

## Additional parameters
# prefixes and suffixes
suf.allscore <- ".allScore"     # suffix of the allScore file containing the contig names and exonerate scores of all successful exonerate alignments
suf.exonerate <- ".exonerate"   # suffix of the exonerate file containing the alignment sugar (query qrange target trange score)
# suf.bestscore <- ".bestScore"
# suf.fasta <- ".fasta"                            
locus.prefix <- paste0("^", basename(folder), ".") # string that comes before the locus name in an <*.suf.allscore> file

# grep contig length from contig name
get.length <- TRUE        # set to TRUE if the contig length can be grepped from the contig name
length.split <- "_"       # if get.lenght = TRUE, specifies a string by which to split the contigname to access the length info, e.g. <_> for contigname <3442_contig_937_length"> 
length.string <- "length" # if get.length = TRUE, specifies a unique string that comes before or after the length info, e.g. <length> for contigname <3442_contig_937_length"> 
length.where <- -1        # if get.length = TRUE, specifies an index relative to <length.string> that directs to the length info, e.g. <-1> for contigname <3442_contig_937_length">

# get alignment ranges from .exonerate files
get.range <- TRUE         # set to TRUE if there are .exonerate result files with alignment sugar
get.overlap <- TRUE       # if TRUE, will compute overlap length of alignments between a target and several queries (takes time)
max.contigs <- 20         # will summarize overlaps between up to <max.contigs> aligned contigs
format <- c("query", "qstart", "qend", "target", "tstart", "tend", "score") # sugar format (sequence of fields in .exonerate file)
tname <- "target"  # name of target sequence column in .exonerate file
tsname <- "tstart" # name of target alignment start in .exonerate file
tename <- "tend"   # name of target alignment end in .exonerate file 
q.split <- " - "   # splits query from target and / or target from score
t.split <- " [+] " # splits query from target and / or target from score

# name of output files
tabfile <- "contig-table.txt" # name of contig table file
ddfile <- "contig-stats.txt" # name of contig stats file
rangefile <- "contig-ranges.txt" # name of contig alignment ranges file

## Define Helperfunctions
# contig scores
get.scores <- function(scorefile) {
  dscores <- unlist(strsplit(readLines(scorefile), split = " "))
  if (length(dscores) > 1) {
    # exonerate might output > 1 alignment per query and target even if bestn is set to 1
    contigs <- dscores[seq(1, length(dscores), by = 2)]
    # this prevents duplicate alignments (same query and target) from being counted
    scores <- as.numeric(dscores[seq(2, length(dscores), by = 2)])[!duplicated(contigs)]
  } else {
    scores <- as.numeric(dscores)
  }
  return(scores)
}
get.lengths <- function(scorefile, split, string, where) {
  dscores <- unlist(strsplit(readLines(scorefile), split = " "))
  dlengths <- unlist(strsplit(dscores, split = split))
  if (length(dscores) > 1) {
    # exonerate might output > 1 alignment per query and target even if bestn is set to 1
    contigs <- dscores[seq(1, length(dscores), by = 2)]
    # this prevents duplicate alignments (same query and target) from being counted
    lengths <- as.numeric(dlengths[grep(string, dlengths) + where])[!duplicated(contigs)]
  } else {
    lengths <- as.numeric(dscores)
  }
  return(lengths)
}
get.ranges <- function(exoneratefile, format, split1, split2) {
  f <- readLines(exoneratefile)
  drange <- data.frame(array(NA, dim = c(max(c(length(f),1)), 7), 
                             dimnames = list(NULL, format)))
  if (length(f) > 0) {
    for (i in seq(length(f))) {
      sugar <- unlist(strsplit(unlist(strsplit(f[i], split = split1)), split = split2))
      
      qlist <- strsplit(sugar[1], split = " ")
      drange[i, 1] <- sapply(qlist, "[", 1)
      drange[i, 2] <- as.numeric(sapply(qlist, "[", 2))
      drange[i, 3] <- as.numeric(sapply(qlist, "[", 3))
      
      tlist <- strsplit(sugar[2], split = " ")
      drange[i, 4] <- sapply(tlist, "[", 1)
      drange[i, 5] <- as.numeric(sapply(tlist, "[", 2))
      drange[i, 6] <- as.numeric(sapply(tlist, "[", 3))
      
      drange[i, 7] <- as.numeric(sugar[3])
      }
  } else {
    drange[1, 4] <- gsub(paste0("^", basename(dirname(exoneratefile)), "."), "", 
                                    gsub(paste0(suf.exonerate, "$"), "", basename(exoneratefile)))
    drange[1, c(1:3,5:6)] <- NA
    drange[1, 7] <- 0
  }
  # this prevents duplicate alignments (same query and target) from being counted (only keeps the one with higher score)
  drange <- drange[order(drange[,4],drange[,7], method = "radix", decreasing = c(FALSE, TRUE)),]
  return(drange[!duplicated(drange[,1]),])
}
get.overlaps <- function(df, tsname, tename, max.contigs = 20) {
  for (j in 1:max.contigs) {
    k <- j + 1
    t1 <- sort(as.numeric(df[j,c(tsname, tename), with = FALSE]))
    repeat {
      t2 <- sort(as.numeric(df[k,c(tsname, tename), with = FALSE]))
      if (length(t2) == 2) {
        df[k,paste0("o", j)] <- length(intersect(t1[1]:t1[2], t2[1]:t2[2]))
        k <- k + 1
      } else {
        break
      }
    }
    if (j >= nrow(df)) {
      break
    }
  }
  return(df)
}

##########################################################################################

## Get file paths
p.allscore <- list.files(folder, pattern = suf.allscore, full.names = TRUE)
p.exonerate <- list.files(folder, pattern = suf.exonerate, full.names = TRUE)

## Get exonerate scores
# score
lscores <- list()
llengths <- list()
for (locus in p.allscore) {
  lscores[[locus]] <- get.scores(locus)
  if (get.length) {
    llengths[[locus]] <- get.lengths(locus, split = length.split, string = length.string, where = length.where)
  }
}
# lscores <- sapply(p.allscore, FUN = get.scores) # takes too long for >3000 loci
names(lscores) <- gsub(locus.prefix, "", gsub(paste0(suf.allscore, "$"), "", basename(names(lscores))))

# contig lengths
if (get.length) {
  # llengths <- sapply(p.allscore, FUN = get.lengths, split = length.split, string = length.string, where = length.where) # takes too long for >3000 loci
  names(llengths) <- names(lscores)
}

# alignment ranges
# t1 <- Sys.time()
if (length(p.exonerate) == 0) get.range <- FALSE
if (get.range) {
  if (get.overlap) {
    lranges <- data.table(array(NA, dim = c(0, 7+max.contigs),
                                dimnames = list(NULL, c(format, paste0("o", 1:max.contigs)))))
  } else {
    lranges <- data.table(array(NA, dim = c(0, 7),
                                dimnames = list(NULL, c(format))))
  }
  
  for (locus in p.exonerate) {
    # get ordered alignment range data table
    drange <- data.table(get.ranges(locus, format = format, split1 = q.split, split2 = t.split))
    setkey(drange, target, score)
    drange <- drange[order(target, -score)]
    
    if (get.overlap) {
      # get length of multi-contig alignment overlaps
      for (h in 1:max.contigs) {drange[,paste0("o", h) := as.numeric(NA)]}
      lranges <- rbindlist(list(lranges, get.overlaps(drange, tsname, tename, max.contigs)))
    } else {
      lranges <- rbindlist(list(lranges, drange))
    }
  }
  # order by target locus and decreasing score
  setkey(lranges, target, score)
  lranges <- lranges[order(target, -score)]

  # only keep best to merge with dd
  if (get.overlap) {
    bestranges <- lranges[!duplicated(lranges[, tname, with = FALSE])][,!paste0("o", 1:max.contigs), with = FALSE]
  } else {
    bestranges <- lranges[!duplicated(lranges[, tname, with = FALSE])]
  }
}
# t2 <- Sys.time()
# print(t2-t1) # for 6555 loci, took 0.01 secs if get.range=F, 16.2 secs if get.range=T and get.overlap=F, 1.17 mins if get.range=T and get.overlap=T

# combine
dd <- data.frame(ID = basename(folder), 
                 locus = names(lscores), # defined 'locus' as name of target
                 ncontigs = unlist(lapply(lscores, length)),
                 scores = unlist(lapply(lscores, paste, collapse = ", ")),
                 stringsAsFactors = FALSE,
                 row.names = NULL)
dd$whichbest <- as.numeric(sapply(strsplit(dd$scores, split = ", "), function(x) which.max(as.numeric(x))))
dd$bestscore <- as.numeric(sapply(seq(nrow(dd)), FUN = function(x) {sapply(strsplit(dd$scores[x], split = ", "), "[", dd$whichbest[x])}))

if (get.length) {
  dd$lengths <- unlist(lapply(llengths, paste, collapse = ", "))
  dd$bestlength <- as.numeric(sapply(seq(nrow(dd)), FUN = function(x) {sapply(strsplit(dd$lengths[x], split = ", "), "[", dd$whichbest[x])}))
} else {
  dd$lengths <- NA
  dd$bestlength <- NA
}

if (get.range) {
  names(bestranges)[names(bestranges) == tname] <- "locus"
  stopifnot(all(bestranges$locus %in% dd$locus))
  ddnames <- names(dd)
  dd <- merge(dd, bestranges, by = "locus", all.x = TRUE, all.y = FALSE, sort = FALSE)
  dd <- dd[,c(ddnames, names(dd)[!names(dd) %in% ddnames])]
} else {
  for (i in format[!format==tname]) dd[,i] <- NA
}

# Handle NAs
dd$ncontigs[dd$bestscore == 0] <- 0
dd[which(dd$ncontigs == 0), c("scores","bestscore","whichbest","lengths","bestlength")] <- NA

## Tabulate number of queries (contigs) where exonerate was able to align a query (contig) to a target (reference sequence)
contab <- table(dd$ncontigs)
tab <- data.frame("ncontigs" = names(contab), "nloci" = as.numeric(contab))

## Write results
write.table(tab, file = file.path(folder, tabfile), row.names = FALSE, quote = FALSE, sep = "\t")
write.table(dd, file = file.path(folder, ddfile), row.names = FALSE, quote = FALSE, sep = "\t")

if (get.range) {
  write.table(lranges, file = file.path(folder, rangefile), row.names = FALSE, quote = FALSE, sep = "\t")
}


