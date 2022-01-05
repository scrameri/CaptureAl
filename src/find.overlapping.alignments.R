#!/cluster/apps/r/3.6.0_openblas/x86_64/bin/Rscript

## Usage: find.overlapping.alignments.R <.blast.filtered> <...>

## Value: summary of overlapping fasta alignments in <.blast.filtered> results (a table in .overlap and list in .list)

## Author: simon.crameri@env.ethz.ch, May 2019

## Get arguments
args = commandArgs(trailingOnly = TRUE)

## Check arguments
options(warning.length = 2000L)
if (!length(args) %in% c(1:9)) {
  stop("1 argument needed (9 arguments taken):
       REQUIRED
       1) <bf|CHR>:                path to filtered blast vs. self results (*.blast.filtered) 
     
       OPTIONAL
       2) <check.lg|BOOLEAN>:      if TRUE, checks linkage group (LG) conformity based on query and subject identifiers, string.lg.before and string.lg.after [DEFAULT: TRUE]
       3) <string.lg.before|CHR>:  unique string in locus name that comes just BEFORE the linkage group ID [DEFAULT: 'LG_']
       4) <string.lg.after|CHR>:   string in locus name that comes just AFTER the linkage group ID [DEFAULT: '_']
       5) <check.overlap|BOOLEAN>: if TRUE, checks whether hits are at alignment ends based on q/s.start/end and query/subject.length and <tol>/<both.ends> [DEFAULT: TRUE]
       6) <tol|NUM>:               alignments ending at <tol> basepairs from a query/subject start/end end will still be considered as being at an alignment end [DEFAULT: 5]
       7) <both.ends|BOOLEAN>:      if TRUE, both query and subject alignments need to be at an alignment end ; if FALSE, one is enough [DEFAULT: FALSE]
       8) <min.alnlen|NUM>:        minimum alignment length for hit to be considered [DEFAULT: 0]
       9) <min.percid|NUM>:        minimum percent identity for hit to be considered [DEFAULT: 0]
       ",
       call.=FALSE)
}

## Set arguments
bf <- as.character(args[1])

check.lg <- as.character(args[2])
if (is.na(check.lg)) check.lg <- TRUE
string.lg.before <- as.character(args[3])
if (is.na(string.lg.before)) string.lg.before = "LG_"
string.lg.after <- as.character(args[4])
if (is.na(string.lg.after)) string.lg.after = "_"

check.overlap <- as.character(args[5])
if (is.na(check.overlap)) check.overlap <- TRUE
tol <- as.numeric(as.character(args[6]))
if (is.na(tol)) tol <- 5
both.ends <- as.character(args[7])
if (is.na(both.ends)) both.ends <- FALSE

min.alnlen <- as.numeric(as.character(args[8]))
if (is.na(min.alnlen)) min.alnlen <- 0
min.percid <- as.numeric(as.character(args[9]))
if (is.na(min.percid)) min.percid <- 0

## Get booleans
get.boolean <- function(x, truestrings = c("T","TRUE","t","true","True"), 
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
check.lg <- get.boolean(check.lg)
check.overlap <- get.boolean(check.overlap)
both.ends <- get.boolean(both.ends)

# ## Set arguments (for debugging)
# bf <- "test.vs.self.blast.filtered"
# check.lg = TRUE
# string.lg.before = "LG_"
# string.lg.after = "_"
# check.overlap = TRUE
# tol = 5
# both.ends = FALSE
# min.alnlen = 0
# min.percid = 0

## Additional arguments
# input file
expected <- c("query.id", "subject.id") # needed fields to remove duplicated hits
expected.over <- c("q.start", "q.end", "query.length", "s.start", "s.end", "subject.length") # needed fields if check.overlap = T

# output files
outtab <- gsub(".vs.self.blast.filtered", ".overlap", bf) # name of output table
outlist <- gsub(".vs.self.blast.filtered", ".list", bf)   # name of output list

# verbosity
verbose = TRUE

## Check arguments
stopifnot(file.exists(bf),
          is.logical(check.lg),
          is.logical(check.overlap),
          is.numeric(tol), tol >= 0,
          is.logical(both.ends),
          is.numeric(min.alnlen), min.alnlen >= 0,
          is.numeric(min.percid), min.percid >= 0, min.percid <= 100,
          is.logical(verbose))

## Be verbose
if (verbose) {
  cat("\n### Find overlapping alignments ###\n\n")
  cat("input file:  ", bf, "\n")
  cat("output table:", outtab, "\n")
  cat("output list: ", outlist, "\n\n")
  if (check.lg) {
    cat(paste0("will filter for hits within the same linkage group, assuming that all locus names have a <", string.lg.before, "xxx", string.lg.after, "> pattern\n"))
  }
  if (check.overlap) {
    cat(paste0("will filter for hits at alignment ends in ", ifelse(both.ends, "both", "at least one"), " involved alignment", ifelse(both.ends, "s", ""), " using a tolerance of ", tol, "bp\n"))
  }
}

## Define helperfunctions
# quiet stop
stopQuietly <- function(...) {
  blankMsg <- sprintf("\r%s\r", paste(rep(" ", getOption("width")-1L), collapse=" "))
  stop(simpleError(blankMsg))
}

# get linkage group
get.lg <- function(x, string.lg.before, string.lg.after) {
  sapply(strsplit(sapply(strsplit(x, split = string.lg.before), function(x) {rev(x)[1]}), split = string.lg.after), "[", 1)
}

# get end overlap
get.overlap <- function(x, tol = 0, both.ends = TRUE) {
  t1 <- as.numeric(x["q.start"]) %in% 1:tol
  t2 <- as.numeric(x["q.start"]) %in% as.numeric(x["query.length"]):(as.numeric(x["query.length"])-tol)
  t3 <- as.numeric(x["s.start"]) %in% 1:tol
  t4 <- as.numeric(x["s.start"]) %in% as.numeric(x["subject.length"]):(as.numeric(x["subject.length"])-tol)
  
  t5 <- as.numeric(x["q.end"]) %in% 1:tol
  t6 <- as.numeric(x["q.end"]) %in% as.numeric(x["query.length"]):(as.numeric(x["query.length"])-tol)
  t7 <- as.numeric(x["s.end"]) %in% 1:tol
  t8 <- as.numeric(x["s.end"]) %in% as.numeric(x["subject.length"]):(as.numeric(x["subject.length"])-tol)
  
  if (both.ends) {
    res <- any(c(t1 && t3, t1 && t4, t2 && t3, t2 && t4, 
                 t1 && t7, t1 && t8, t2 && t7, t2 && t8,
                 t5 && t3, t5 && t4, t6 && t3, t6 && t4, 
                 t5 && t7, t5 && t8, t6 && t7, t6 && t8))
  } else {
    res <- any(c(t1,t2,t3,t4,t5,t6,t7,t8))
  }
  return(res)
}

## Move any existing output files
if (file.exists(outtab)) system(paste("mv", outtab, paste0(outtab, ".bak")))
if (file.exists(outlist)) system(paste("mv", outlist, paste0(outlist, ".bak")))

## Read .blast.filtered
dbf <- read.delim(bf, stringsAsFactors = FALSE)
if (nrow(dbf) < 1) {
  cat(paste0("\nno blast hits found in <", bf, ">!\n"))
  stopQuietly()
}
if (!all(expected %in% names(dbf))) {
  cat("field names", paste(expected, collapse = ", "), "expected!\n")
  stopQuietly()
}

## Sort involved loci per hit
dbf$query.subject.sorted <- NA
for (i in seq(nrow(dbf))) {
  dbf$query.subject.sorted[i] <- paste(sort(c(dbf$query.id[i], dbf$subject.id[i])), collapse = " <-> ")
}

## Check whether locus pairs are from the same linkage group (putatively physically adjacent)
if (check.lg) {
  dbf$lg1 <- get.lg(dbf$query.id, string.lg.before, string.lg.after)
  dbf$lg2 <- get.lg(dbf$subject.id, string.lg.before, string.lg.after)
  dbf$samelg <- dbf$lg1 == dbf$lg2
} else {
  dbf$lg1 <- dbf$lg2 <- NA
  dbf$samelg <- TRUE
}

## Check whether locus pairs overlap at alignment ends (putatively physically adjacent)
if (check.overlap) {
  if (!all(expected.over %in% names(dbf))) {
    cat("field names", paste(expected.over, collapse = ", "), "expected if check.overlap = TRUE!\n")
    stopQuietly()
  }
  dbf$endoverlap <- apply(dbf, MARGIN = 1, FUN = get.overlap, tol = tol, both.ends = both.ends)
} else {
  dbf$endoverlap <- TRUE
}

## Check alignment length and percent identity
if (min.alnlen > 0) {
  if (!"alignment.length" %in% names(dbf)) {
    cat("field name <alignment.length> expected if min.alnlen > 0!\n")
    stopQuietly()
  }
  dbf$alnlen <- ifelse(dbf$alignment.length >= min.alnlen, TRUE, FALSE)
} else {
  dbf$alnlen <- TRUE
}
if (min.percid > 0) {
  if (!"perc.identity" %in% names(dbf)) {
    cat("field name <perc.identity> expected if min.percid > 0!\n")
    stopQuietly()
  }
  dbf$percid <- ifelse(dbf$perc.identity >= min.percid, TRUE, FALSE)
} else {
  dbf$percid <- TRUE
}

## Check for duplicated hits (A <> B & A <> B or A <> B & B <> A) and keep the best-scoring one
dbf$dup <- TRUE
for (i in unique(dbf$query.subject.sorted)) {
  ls <- dbf[dbf[,"query.subject.sorted"] == i,]
  dbf[names(which.max(apply(ls[,c("samelg","endoverlap","alnlen","percid")], MARGIN = 1, FUN = sum))), "dup"] <- FALSE
}

## Filter BLAST hits for non-duplicated and good-quality hits between alignment ends and within the same LG
# get non-duplicated hits between alignment ends and within the same LG passing min.alnlen and min.percid
dbf.uniq <-   subset(dbf, dup == FALSE & alnlen == TRUE & percid == TRUE)
dbf.passed <- subset(dbf, dup == FALSE & alnlen == TRUE & percid == TRUE & samelg == TRUE & endoverlap == TRUE)
rownames(dbf.passed) <- rownames(dbf.uniq) <- NULL

if (nrow(dbf.passed) > 0) {
  
  # get hits involving more than 1 other locus (twoormore == TRUE, e.g. A <-> X <-> B)
  o.loci <- c(dbf.passed$query.id, dbf.passed$subject.id)
  dup.o.loci <- o.loci[duplicated(o.loci)]
  dbf.passed$twoormore <- dbf.passed$query.id %in% dup.o.loci | dbf.passed$subject.id %in% dup.o.loci
  
  dbf.passed_2 <- subset(dbf.passed, twoormore == FALSE)
  dbf.passed_3 <- subset(dbf.passed, twoormore == TRUE)
  
  # get involved loci with hits to exactly 1 other locus
  if (nrow(dbf.passed_2) > 0) {
    # id of all overlapping loci (duplicates removed)
    l.2.uniq <- character()
    for (i in seq(nrow(dbf.passed_2))) {
      loci2 <- sort(c(dbf.passed_2[i,"query.id"], dbf.passed_2[i,"subject.id"]))
      l.2.uniq <- c(l.2.uniq, paste(loci2, collapse = " "))
    }
  } else {
    l.2.uniq <- NA
  }
  
  # get involved loci with hits to more than 1 other locus
  if (nrow(dbf.passed_3) > 0) {
    # id of all overlapping loci (duplicated for each locus)
    l.involved <- character()
    for (i in seq(nrow(dbf.passed_3))) {
      ids <- grep(dbf.passed_3[i,"subject.id"], dbf.passed_3$query.subject.sorted)
      involved <- sort(unique(c(dbf.passed_3[ids,c("query.id")], dbf.passed_3[ids,c("subject.id")])))
      l.involved[i] <- paste(involved, collapse = ", ")
    }
    
    # indices list of all overlapping loci (duplicated for each locus)
    keep.involved <- list()
    for (i in 1:length(l.involved)) {
      involved <- unlist(strsplit(l.involved[i], split = ", "))
      ids <- grep(paste(involved, collapse = "|"), l.involved)
      keep.involved[[i]] <- ids
    }
    
    # indices list of all overlapping loci (duplicates removed)
    keep.involved.uniq <- keep.involved[which(!duplicated(keep.involved))]
    keep.involved.uniq.order <- order(unlist(lapply(keep.involved.uniq, length)), decreasing = TRUE)
    
    # id of all overlapping loci (duplicates removed)
    l.3.uniq <- character()
    for (i in 1:length(keep.involved.uniq)) {
      loci3 <- sort(unique(unlist(strsplit(l.involved[keep.involved.uniq[[i]]], split = ", "))))
      l.3.uniq <- c(l.3.uniq, paste(loci3, collapse = " "))
      # print(length(keep.involved.uniq[[i]]))
      # print(length(loci3))
      # scan()
    }
    l.3.uniq <- l.3.uniq[keep.involved.uniq.order]
  } else {
    l.3.uniq <- NA
  }
} else {
  l.2.uniq <- l.3.uniq <- NA
}

## Get and check overlapping loci
l.res <- na.omit(c(l.3.uniq, l.2.uniq))
two <- sort(unlist(strsplit(l.2.uniq, split = " ")))
if (!is.na(l.3.uniq)) more <- sort(unlist(strsplit(l.3.uniq, split = " "))) else more <- NA
if (any(duplicated(two))) {
  cat("\nWARNING: duplicated loci detected and removed in output list (2)!\n")
  #stopQuietly()
}
if (any(duplicated(more))) {
  cat("\nWARNING: duplicated loci detected and removed in output list (>2)!\n")
  #stopQuietly()
}
if (sum(c(length(two), length(na.omit(more)))) != length(unlist(strsplit(l.res, split = " ")))) {
  cat("\nWARNING: duplicated loci detected and removed in output list!\n")
}

# be verbose
if (verbose) {
  cat("\nfound", nrow(dbf), "hit(s) in total\n")
  cat("removed", sum(dbf$dup), "duplicated hit(s) between the same query/subject pair\n")
  if (min.alnlen > 0) cat("removed", sum(!dbf$alnlen), "hit(s) with alignment length below", min.alnlen, "\n")
  if (min.percid > 0) cat("removed", sum(!dbf$percid), "hit(s) with percent identity below", min.percid, "\n")
  if (check.lg) cat("removed", sum(!dbf.uniq$samelg), "hit(s) between different linkage groups\n")
  if (check.overlap) cat("removed", sum(!dbf.uniq$endoverlap), "hit(s) that are not between alignment ends\n")
  cat("\n")
  print(table(samelg=dbf.uniq$samelg, endoverlap=dbf.uniq$endoverlap))
  cat("\nthis amounts to", nrow(dbf.passed), "overlap(s) within the same linkage groups and between alignment ends:\n")
  cat(" ", length(unlist(strsplit(l.res, split = " "))), "loci involved\n")
  cat(" ", nrow(dbf.passed_2), "overlap(s) involve exactly 2 loci\n")
  cat(" ", nrow(dbf.passed_3), "overlap(s) involve more than 2 loci\n")
  cat("\nthis amounts to", length(na.omit(c(l.3.uniq, l.2.uniq))), "region(s) with overlaps:\n")
  # cat(" ", length(na.omit(l.2.uniq)), "region(s) involve exactly 2 loci\n")
  cat(" ", length(na.omit(l.3.uniq)), "region(s) involve more than 2 loci\n")
  if (sum(dbf.passed$twoormore) > 0 & sum(dbf.passed$twoormore) < 10) {
    print(l.3.uniq)
  }
  cat("\n")
}

## Write results
if (nrow(dbf.passed) > 0) {
  write.table(dbf, file = outtab, row.names = F, col.names = T, quote = F, sep = "\t")
  write.table(l.res[!duplicated(l.res)], file = outlist, row.names = F, col.names = F, quote = F, sep = "\t")
} else {
  cat("no overlapping loci found\n")
}

