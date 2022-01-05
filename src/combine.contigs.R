#!/cluster/apps/r/3.6.0_openblas/x86_64/bin/Rscript

#####################################
## Combine non-overlapping Contigs ##
#####################################

## Usage
# combine.contigs.R <folder> <OPT: cpath> <OPT: fpath> <OPT: get.range> <OPT: combine> <OPT: min.aln> <OPT: min.norm.score>
# run script without arguments to see the arguments it takes

## Author
# simon.crameri@env.ethz.ch, Mar 2020

## Load required library
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ape))

## Get arguments
args = commandArgs(trailingOnly = TRUE)

## Check arguments
options(warning.length = 2000L)
setDTthreads(threads = 1) # prevent cluster to use more than 1 process per parallel execution
if (! length(args) %in% 1:7) {
  stop("1 argument needed (5 taken):
       REQUIRED
       1) <folder|CHR>:        path to exonerate results folder ; 
       
       OPTIONAL
       2) <cpath|CHR>:         path to contigs <e.g. path/to/consensus_contigs.fasta>, can include SAMPLE and LOCUS tags that are replaced using regex 
                               [DEFAULT: assemblies/SAMPLE.dipspades/extracted_reads_SAMPLE.fastq.LOCUS.ids.spades/dipspades/consensus_contigs.fasta] ;
       3) <fpath|CHR>:         path to best contig(s) <e.g. path/to/*bestScore.fasta>, can include SAMPLE and LOCUS tags that are replaced using regex 
                               [DEFAULT: file.path(folder, SAMPLE.LOCUS.bestScore.fasta)] ;
       4) <get.range|BOOLEAN>: if TRUE, will look for <*.exonerate> files in <folder> and identify overlapping (paralogous loci?) vs. non-overlapping 
                               (locus fractions?) contigs from the (required) alignment sugar. [DEFAULT: TRUE] ;
       5) <combine|BOOLEAN>:   if TRUE and get.range = TRUE, will combine non-overlapping contigs into a single supercontig with adequate number 
                               of Ns used as a spacer [DEFAULT: TRUE] ;
       6) <min.aln|NUMERIC>:   minimum target (reference) alignment length for kept / combined contigs. Only active if get.range = TRUE and combine = TRUE. 
                               [DEFAULT: 80] ;
       7) <min.norm.score|NUM>:minimum normalized alignment score for kept / combined contigs.
                               Only active if get.range = TRUE and combine = TRUE. [DEFAULT: 2] ",
       call.=FALSE)
}

## Set prefixes and suffixes
suf.allscore <- ".allScore"     # suffix of the allScore file containing the contig names and exonerate scores of all successful exonerate alignments
suf.exonerate <- ".exonerate"   # suffix of the exonerate file containing the alignment sugar (query qrange target trange score)
suf.bestscore <- ".bestScore"   # suffix of bestScore file containing the score of the best-matching contig (or zero if no contig aligned)
suf.fasta <- ".fasta"           # suffix for .fasta formatted files
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

## Set arguments
folder <- as.character(args[1])

cpath <- as.character(args[2])
if (is.na(cpath)) cpath <- "assemblies/SAMPLE.dipspades/extracted_reads_SAMPLE.fastq.LOCUS.ids.spades/dipspades/consensus_contigs.fasta"
fpath <- as.character(args[3])
if (is.na(fpath)) fpath <- file.path(folder, paste0("SAMPLE.LOCUS", suf.bestscore, suf.fasta))
get.range <- as.character(args[4])
if (is.na(get.range)) get.range <- TRUE else get.range <- get.boolean(get.range)
combine <- as.character(args[5])
if (is.na(combine)) combine <- TRUE else combine <- get.boolean(combine)
min.aln <- as.numeric(as.character(args[6]))
if (is.na(min.aln)) min.aln <- 80
min.norm.score <- as.numeric(as.character(args[7]))
if (is.na(min.norm.score)) min.norm.score <- 2

## Set arguments (for debugging)
# folder <- "best.contigs.63.2984/609-65_S1_L001"
# 
# # contig path (SAMPLE and LOCUS will be replace accordingly using regex)
# cpath <- "assemblies/SAMPLE.dipspades/extracted_reads_SAMPLE.fastq.LOCUS.ids.spades/dipspades/consensus_contigs.fasta"
# fpath <- file.path(folder, paste0("SAMPLE.LOCUS", suf.bestscore, suf.fasta))
# 
# # handle overlapping contigs
# get.range <- TRUE         # set to TRUE if there are .exonerate result files with alignment sugar.
# combine <- TRUE           # if TRUE and get.range = TRUE, will combine non-overlapping contigs
# 
# # filtering thresholds
# min.aln = 80              # minimum target alignment lenght for contig to be used or combined
# min.norm.score = 2        # minimum normalized alignment score for contig to be used or combined

## Check arguments
stopifnot(dir.exists(folder),
         is.logical(get.range), is.logical(combine),
         min.aln >= 0, min.norm.score >= 0)

## Additional parameters
# sample name (from folder basename)
sample <- basename(folder)  # string encoding a sample ID
locus.prefix <- paste0("^", sample, ".") # string that comes before the locus name in an <*.suf.allscore> file

# target locus covered by contigs (from .exonerate files)
max.contigs <- 20           # will summarize overlaps between up to <max.contigs> aligned contigs
rewrite <- TRUE             # if TRUE, will rewrite passed best contigs that are not merged to supercontigs (no changes made)

# contig merging thresholds
tol = 0                     # maximum alignment overlap allowed for two contigs to be considered as potentially non-overlapping
max.overlap = 10            # maximum overlap between contigs to be combined 
min.qlen = 80               # minimum query (contig) length for contig to be used or combined
mis = "--------------------" # reports this sequence if no contig passes thrsholds
NA.char = "-"               # this character will be used to fill gaps between merged contigs

# grep contig length from contig name
get.length <- FALSE         # set to TRUE if get.range = FALSE or get.range = TRUE and combine = FALSE, and the contig length can be grepped from the contig name. 
length.split <- "_"         # if get.lenght = TRUE and get.range = FALSE, specifies a string by which to split the contigname to access the length info, e.g. <_> for contigname <3442_contig_937_length"> 
length.string <- "length"   # if get.length = TRUE and get.range = FALSE, specifies a unique string that comes before or after the length info, e.g. <length> for contigname <3442_contig_937_length"> 
length.where <- -1          # if get.length = TRUE and get.range = FALSE, specifies an index relative to <length.string> that directs to the length info, e.g. <-1> for contigname <3442_contig_937_length">

# expected variable names
format <- c("query", "qstart", "qend", "target", "tstart", "tend", "score") # sugar format (sequence of fields in .exonerate file)
tname <- "target"           # name of target sequence column in .exonerate file
tsname <- "tstart"          # name of target alignment start in .exonerate file
tename <- "tend"            # name of target alignment end in .exonerate file 
q.split <- " - "            # splits query from target and / or target from score
t.split <- " [+] "          # splits query from target and / or target from score

# name of output files
tabfile <- "contig-table.txt"     # name of contig table file
dfile <- "contig-stats.txt"       # name of contig stats file
rangefile <- "contig-ranges.txt"  # name of contig alignment ranges file

##########################################################################################

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

# contig lengths
get.lengths <- function(scorefile, split, string, where) {
  dscores <- unlist(strsplit(readLines(scorefile), split = " "))
  dlengths <- unlist(strsplit(dscores, split = split))
  if (length(dscores) > 1) {
    # exonerate might output > 1 alignment per query and target even if bestn is set to 1
    contigs <- dscores[seq(1, length(dscores), by = 2)]
    # this prevents duplicate alignments (same query and target) from being counted
    lengths <- as.numeric(dlengths[grep(string, dlengths) + where])[!duplicated(contigs)]
  } else {
    lengths <- numeric() #as.numeric(dscores)
  }
  return(lengths)
}

# contig names
get.names <- function(scorefile) {
  dscores <- unlist(strsplit(readLines(scorefile), split = " "))
  if (length(dscores) > 1) {
    # exonerate might output > 1 alignment per query and target even if bestn is set to 1
    contigs <- dscores[seq(1, length(dscores), by = 2)]
    names <- contigs[!duplicated(contigs)]
  } else {
    names <- character()
  }
  return(names)
}

# contig ranges
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

# contig overlaps
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

# multicontig stats
get.nonoverlaps <- function(dl, locid = "locid", score = "score", o1 = "o1", ID = "ID",
                            tol = 0, all = FALSE, stats = FALSE, plot = FALSE) {
  
  # nb of additional (discarded) contigs
  s <- tapply(X = dl[[o1]], INDEX = dl[[locid]], FUN = function(x) {length(x[!is.na(x)])})
  
  # nb of additional (discarded) contigs that are not overlapping with the selected contig (if ovar = o1)
  sno <- tapply(X = dl[[o1]], INDEX = dl[[locid]], FUN = function(x) {sum(x <= tol, na.rm = T)})
  
  # score of best (selected) contig
  dlsub <- dl[is.na(o1)]
  bs <- tapply(X = dlsub[[score]], INDEX = dlsub[[locid]], FUN = function(x) {identity(x)})
  
  # maximum score of additional (discarded) contigs
  dlsub <- dl[o1 <= tol]
  msc <- tapply(X = dlsub[[score]], INDEX = dlsub[[locid]], FUN = function(x) {max(x, na.rm = T)})
  
  # create results data.frame
  d.pno <- data.table(names(s), s, sno, sno/s, bs, msc, msc/bs) # nsamples x nloci rows
  names(d.pno) <- c(locid, "nb.contig.surp", "nb.contig.surp.no", "no.contig.ratio", "best.score", "best.score.no", "score.ratio")
  d.pno[,ID] <- dl[match(d.pno[,locid], dl[,locid]),ID]
  d.pno <- d.pno[order(d.pno[,locid]),]
  d.pno.add <- d.pno[nb.contig.surp.no > 0]
  # rownames(d.pno) <- rownames(d.pno.add) <- NULL
  
  # aggregated statistics per individual
  if (stats) {
    # using data.table
    A <- d.pno.add[,.(count=.N),by=ID]
    a <- d.pno.add[,.(mean=mean(nb.contig.surp)),by=ID]
    n <- d.pno.add[,.(mean=mean(nb.contig.surp.no)),by=ID]
    s <- d.pno.add[,.(mean=mean(best.score.no)),by=ID]
    r <- d.pno.add[,.(mean=mean(score.ratio)),by=ID]
    d.added <- data.table(A, a$mean, n$mean, s$mean, r$mean)
    names(d.added) <- c(ID, "mean.nloc.surp.no", "mean.nb.surp", "mean.nb.surp.no", "mean.score", "mean.score.ratio")
    
  }
  
  # plot
  if (plot) {
    # number of loci
    nums <- d.pno.add[,.(count=.N),by=ID]
    dnums <- data.frame(ID = nums[[ID]], N = nums$count, 
                        x = min(d.pno.add[,"best.score"]), y = max(d.pno.add[,"best.score"]),
                        label = paste0("n = ", nums$count))
    d.pno.add[[ID]] <- factor(d.pno.add[[ID]], levels = dnums[order(-dnums$N),ID])
    
    # ggplots
    library(ggplot2)
    p1 <- ggplot(d.added, aes_string(x = "mean.nloc.surp.no", y = "mean.score", size = "mean.score.ratio", colour = ID)) +
      geom_point() +
      labs(x = "Mean number of loci with non-overlapping additional contigs",
           y = "Mean alignment score of best non-overlapping additional contig",
           size = "Mean ratio: best additional score / best score",
           colour = "") +
      theme_bw()
    
    p2 <- ggplot(d.pno.add, aes_string(x = "best.score", y = "best.score.no", colour = ID)) + # size = "score.ratio",
      geom_point() +
      geom_density_2d(colour = "black", alpha = 0.5) +
      labs(x = "Alignment score of best contig",
           y = "Alignment score of best non-overlapping additional contig",
           size = "Ratio: best additional score / best score",
           colour = "") +
      geom_abline(slope = 1, intercept = 0, colour = "gray45", linetype = 2) +
      geom_text(data = dnums, inherit.aes = FALSE, 
                aes(x = x, y = y, label = label),
                vjust = 0.5, hjust = -0.5) +
      scale_colour_discrete(guide = FALSE) +
      xlim(min(d.pno.add$best.score.no), max(d.pno.add$best.score)) +
      ylim(min(d.pno.add$best.score.no), max(d.pno.add$best.score)) +
      facet_wrap(formula(paste("~",ID))) +
      theme_bw()
  }
  
  # return
  res <- list()
  res[["added"]] <- d.pno.add
  if (all) res[["all"]] <- d.pno
  if (stats) res[["stats"]] <- d.added
  if (plot) res[["p1"]] <- p1
  if (plot) res[["p2"]] <- p2
  return(res)
}

#  get non-overlapping contigs
get.added.contigs <- function(dlocid, ID = "ID", o1 = "o1", query = "query", 
                              max.contigs = 20, tol = 0, verbose = FALSE) {
  if (verbose) cat(paste0(unique(dlocid[,ID]), "\r"))
  dlocid[,query] <- as.character(dlocid[,query])
  res <- character()
  idx <- 1
  for (i in 1:min(nrow(dlocid),max.contigs)) {
    t <- dlocid[i,o1]
    q <- dlocid[i,get(paste0("o",rev(idx)[1]))]
    if (is.na(t)) {
      res <- sort(unique(c(res, dlocid[i,query])))
    }
    if (!is.na(t) & t <= tol & q <= tol) { # !is.na(q) &
      res <- sort(unique(c(res, dlocid[c(i,idx),query])))
      idx <- c(idx, i)
    }
  }
  r <- dlocid[query %in% res]
  return(r)
}

# get reverse complement
get.rc <- function (seq) {
  fwd <- unlist(strsplit("ACGTMRWSYKVHDBNI", split = ""))
  rev <- unlist(strsplit("TGCAKYWSRMBDHVNI", split = ""))
  isupper <- any(seq %in% fwd)
  rv <- function(x) {gsub(x, rev[which(x==fwd)], x)}
  r <- unname(rev(sapply(toupper(seq), FUN = rv)))
  if (isupper) res <- r else res <- tolower(r)
  return(res)
}

# read contigs
read.contigs <- function(folder = "best.contigs.63.2984/H121/", prefix = paste0(basename(folder), "[.]"), suffix = ".bestScore.fasta") {
  con <- list.files(folder, pattern = suffix, full.names = TRUE)
  z <- as.list(ape::as.DNAbin("ATG"))
  for (i in con) {
    z <- c(z, as.list(ape::read.FASTA(i)))
    names(z)[length(z)] <- gsub(paste0("^", prefix), "", gsub(paste0(suffix, "$"), "", basename(i)))
  }
  return(z[-1])
}

# write combined contigs
write.combined.contigs <- function(dl, NA.char = "-", score = "score", o1 = "o1", ID = "ID", locid = "locid",
                                   query = "query", target = "target", tstart = "tstart", tend = "tend",
                                   cpath, fpath, sample, max.contigs = 20, tol = 0, min.aln = 80, min.qlen = 80, max.overlap = 10,
                                   min.norm.score = 2, verbose = FALSE, overwrite = TRUE) {
  
  ## Debugging
  # score = "score"; o1 = "o1"; ID = "ID"; locid = "locid";query = "query"; target = "target"; tstart = "tstart"; tend = "tend" ; verbose = T ; overwrite = T
  
  ## Helperfunction
  get.spacer <- function(dp, i, b) {
    if (i > b) {
      d0 <- dp[b,"dir"]
      if (d0 == 1 & d1 == 1) {
        # last contig in + direction and this contig in + direction
        e0 <-  dp[b,"qlen"] - dp[b,"qend"]     # e0 is 3' overhang of former contig (ACTG---)
        e1 <- dp[i,"qstart"]                   # e1 is 5' overhang of this contig   (---ACTG)
      }
      if (d0 == 1 & d1 == -1) {
        # last contig in + direction and this contig in - direction
        e0 <-  dp[b,"qlen"] - dp[b,"qend"]     # e0 is 3' overhang of former contig (ACTG---)
        e1 <- dp[i,"qlen"] - dp[i,"qstart"]    # e1 is 5' overhang of this contig   (---ACTG)
      }
      if (d0 == -1 & d1 == 1) {
        # last contig in - direction and this contig in + direction
        e0 <- dp[b,"qend"]                     # e0 is 3' overhang of former contig (ACTG---)
        e1 <- dp[i,"qstart"]                   # e1 is 5' overhang of this contig   (---ACTG)
      }
      if (d0 == -1 & d1 == -1) {
        # last contig in - direction and this contig in - direction
        e0 <- dp[b,"qend"]                     # e0 is 3' overhang of former contig (ACTG---)
        e1 <- dp[i,"qlen"] - dp[i,"qstart"]    # e1 is 5' overhang of this contig   (---ACTG)
      }
      if (length(c) > 0) {
        spacing <- as.numeric(unlist(dp[i,"tstart"] - e1) - (dp[b,"tend"] + e0))
      } else {
        spacing <- 0
      }
    } else {
      spacing <- 0
    }
    return(spacing)
  }
  
  ## Order by target locus and decreasing score
  setkey(dl, target, score)
  setorder(dl, target, -score)

  ## Add sample x locus combinations
  # dl <- data.table(rep(sample, nrow(dl)), dl) ; names(dl)[1] <- ID
  # dl[,locid] <- interaction(unlist(dl[,"target"]), unlist(dl[,"ID"]))
  # 
  # ## Add query and target (reference) alignment length
  # dl[,"qaln"] <- abs(dl[,"qend"] - dl[,"qstart"])
  # dl[,"taln"] <- abs(dl[,"tend"] - dl[,"tstart"])
  # 
  # ## Add normalized alignment score
  # dl[,"score.norm"] <- dl[,"score"] / dl[,"taln"]
  
  ## Get loci with non-overlapping contigs
  d.add <- get.nonoverlaps(dl, locid = locid, score = score, o1 = o1, ID = ID)$added
  da <- rbindlist(sapply(unique(d.add[,locid]), FUN = function(x) {
    get.added.contigs(dlocid = dl[locid == x], ID = ID, o1 = o1, query = query, max.contigs = max.contigs, tol = tol, verbose = verbose)},
    simplify = FALSE))
  da[,"dir"] <- ifelse(da[,qend]-da[,qstart] > 0, 1, -1)
  
  ## Loop through loci
  cnonoverlap <- sort(unique(da[[target]]))
  crest <- sort(unique(dl[! target %in% cnonoverlap,target]))
  dl[,"ok"] <- as.integer(rep(NA_integer_, nrow(dl)))
  dl[,"order"] <- as.integer(rep(NA_integer_, nrow(dl)))
  dl[,"start"] <- as.character(rep(NA, nrow(dl))) 
  dl[,"end"] <- as.character(rep(NA, nrow(dl)))   
  dl[,"qlen"] <- as.integer(rep(NA_integer_, nrow(dl)))
  dl[,"start2"] <- as.character(rep(NA, nrow(dl)))
  dl[,"end2"] <- as.character(rep(NA, nrow(dl)))  
  dl[,"qlen2"] <- as.integer(rep(NA_integer_, nrow(dl)))   
  dl[,"spacer"] <- as.integer(rep(NA_integer_, nrow(dl)))  
  dl[,"dir"] <- as.integer(ifelse(dl[,qend]-dl[,qstart] > 0, 1, -1))
  
  ## Loop through loci without non-overlapping contigs (just one best contig)
  for (loc in crest) {
    if (verbose) cat(loc, "\n")
    
    # get contig paths
    contigs <- gsub("LOCUS", loc, gsub("SAMPLE", sample, cpath))
    best <- gsub("LOCUS", loc, gsub("SAMPLE", sample, fpath))
    
    if (file.exists(contigs)) {
      dc <- as.character(ape::read.FASTA(contigs))
      
      names <- unlist(dl[target == loc & ID == sample,"query"])
      
      if (!all(is.na(names))) {
        ## Case I: contig exits, alignment exists
        dl[target == loc & ID == sample,"ok"] <- 0
        dl[target == loc & ID == sample,"order"] <- NA
        dl[target == loc & ID == sample,"start"] <- dl[target == loc & ID == sample,"start2"] <- unlist(lapply(dc, FUN = function(x) {paste(x[1:8], collapse = "")}))[names]
        dl[target == loc & ID == sample,"end"] <- dl[target == loc & ID == sample,"end2"] <- unlist(lapply(dc, FUN = function(x) {paste(rev(rev(x)[1:8]), collapse = "")}))[names]
        dl[target == loc & ID == sample,"qlen"] <- dl[target == loc & ID == sample,"qlen2"] <- unlist(lapply(dc, FUN = function(x) {length(x)}))[names]
        
        dp <- dl[target == loc & ID == sample]
        
        # do not use contigs if criteria are not met
        if (dp[1,qlen] >= min.qlen & dp[1,"taln"] >= min.aln & dp[1,"score.norm"] >= min.norm.score) {
          dp[1,"ok"] <- 1
          dl[target == loc & ID == sample & query == names[1], "ok"] <- 1
          dl[target == loc & ID == sample & query == names[1], "order"] <- 1
        }
      } else {
        ## Case II: contig exists, no alignment
        dp <- dl[target == loc & ID == sample]
        dp[1,"ok"] <- 0
      }

      # write new contig file (only if not ok)
      if (dp[1,"ok"] == 0) {
        # move old contig file if overwrite == FALSE
        if (!overwrite) {
          system(command = paste("mv", best, gsub("bestScore.fasta$", "bestScore.fasta.bak", best)))
        }
        writeLines(c(paste0(">", sample), mis), con = best)
      }
      if (dp[1,"ok"] == 1 & rewrite) {
        c <- as.matrix(t(dc[[unlist(dp[1,"query"])]]))
        rownames(c) <- sample
        ape::write.FASTA(ape::as.DNAbin(c), file = best)
      }
    }
  }
  
  ## Loop through loci with non-overlapping contigs (best contig plus additional, non-overlapping contigs)
  for (loc in cnonoverlap) {
    if (verbose) cat(loc, "\n")
    
    # get contig paths
    contigs <- gsub("LOCUS", loc, gsub("SAMPLE", sample, cpath))
    best <- gsub("LOCUS", loc, gsub("SAMPLE", sample, fpath))
    
    # get all contigs
    if (file.exists(contigs)) {
      dc <- as.character(ape::read.FASTA(contigs))
      
      dl[target == loc & ID == sample,"ok"] <- 0
      dl[target == loc & ID == sample,"start"] <- unlist(lapply(dc, FUN = function(x) {paste(x[1:8], collapse = "")}))[unlist(dl[target == loc & ID == sample,"query"])]
      dl[target == loc & ID == sample,"end"] <- unlist(lapply(dc, FUN = function(x) {paste(rev(rev(x)[1:8]), collapse = "")}))[unlist(dl[target == loc & ID == sample,"query"])]
      dl[target == loc & ID == sample,"qlen"] <- unlist(lapply(dc, FUN = function(x) {length(x)}))[unlist(dl[target == loc & ID == sample,"query"])]
      
      # paste non-overlapping contigs together
      dp <- da[ID == sample & target == loc,]
      
      if (nrow(dp) > 0) {
        # sort according to tstart
        dp <- dp[order(dp[,tstart]),]
        
        # add contig lengths
        dp[,"qlen"] <- lengths(dc)[na.omit(match(names(dc), dp[,query]))]
        dp[,"ok"] <- as.numeric(rep(NA, nrow(dp)))
        
        # concatenate contigs
        spacing <- NA
        bestscore <- max(as.numeric(unlist(dp[,"score"])))
        bestname <- as.character(unlist(dp[score == bestscore,"query"]))
        
        for (i in seq(nrow(dp))) {
          d1 <- unlist(dp[i,"dir"])
          name <- unlist(dp[i,"query"])
          
          # determine contig spacer: first (no spacer defined)
          if (i == 1) {
            if (length(dc[[name]]) >= min.qlen & dp[i,"taln"] >= min.aln & dp[i,"score.norm"] >= min.norm.score) {
              if (d1 == 1) {
                # first contig in + direction
                c <- dc[[name]]
              } else {
                # first contig in - direction
                c <- get.rc(dc[[name]])
              }
              dp[i,"ok"] <- 1
            } else {
              c <- character()
              dp[i,"ok"] <- 0
            }
            assign(paste0("c", i), c)
          }
          
          # determine contig spacer: second to last (spacer needs to be defined)
          if (i > 1)  {
            
            bb <- which(unlist(dp[1:(i-1),"ok"]) == 1)
            
            # CASE: no previous ok contig (spacer set to zero)
            if (length(bb) == 0) {
              if (length(dc[[name]]) >= min.qlen & dp[i,"taln"] >= min.aln & dp[i,"score.norm"] >= min.norm.score) {
                if (d1 == 1) {
                  # first contig in + direction
                  c <- dc[[name]]
                } else {
                  # first contig in - direction
                  c <- get.rc(dc[[name]])
                }
                dp[i,"ok"] <- 1
              } else {
                c <- c
                dp[i,"ok"] <- 0
              }
              assign(paste0("c", i), c)
            }
            
            # CASE: at least one previous ok contig (spacer not zero)
            if (length(bb) > 0) {
              b <- max(bb)
              spacing <- get.spacer(dp, i, b)
              
              # case: spacer positive (no trimming needed)
              if (spacing >= 0) {
                if (length(c) > 0) {r <- rep(NA.char, spacing)} else {r <- NULL}
                if (d1 == 1) {
                  a <- c(r, dc[[name]])
                }
                if (d1 == -1) {
                  a <- c(r, get.rc(dc[[name]]))
                }
                
                # do not combine contigs if criteria are not met
                if (length(dc[[name]]) >= min.qlen & dp[i,"taln"] >= min.aln & dp[i,"score.norm"] >= min.norm.score) {
                  c <- c(c, a)
                  dp[i,"ok"] <- 1
                } else {
                  c <- c
                  dp[i,"ok"] <- 0
                }
                assign(paste0("c", i), c)
              }
              
              # case: spacer negative (trimming needed, or discarding of contig)
              if (spacing < 0) {
                if (length(c) > 0) {r <- rep(NA.char, 1)} else {r <- NULL}
                if (d1 == 1) {
                  a <- c(r, dc[[name]][-c(1:(abs(spacing)+1))])
                }
                if (d1 == -1) {
                  a <- c(r, get.rc(dc[[name]])[-c(1:abs(spacing))])
                }
                
                # do not combine contigs if criteria are not met
                if (length(dc[[name]]) >= min.qlen & dp[i,"taln"] >= min.aln & abs(spacing) <= max.overlap & dp[i,"score.norm"] >= min.norm.score) {
                  c <- c(c, a)
                  dp[i,"ok"] <- 1
                } else {
                  # criteria not met, but best score needs to be kept and preceding ok contig removed
                  if (dp[i,"score"] == bestscore) {
                    dp[b,"ok"] <- 0
                    dl[target == loc & ID == sample & query == unlist(dp[b,"query"]),"ok"] <- 0
                    bb <- which(unlist(dp[1:(i-1),"ok"]) == 1)
                    if (d1 == 1) {
                      a <- dc[[name]]
                    }
                    if (d1 == -1) {
                      a <- get.rc(dc[[name]])
                    }
                    
                    # no previous ok contig (no spacer needed)
                    if (length(bb) == 0) {
                      if (length(dc[[name]]) >= min.qlen & dp[i,"taln"] >= min.aln & dp[i,"score.norm"] >= min.norm.score) {
                        c <- a
                        dp[i,"ok"] <- 1
                      } else {
                        c <- character()
                        dp[i,"ok"] <- 0
                      }
                    }
                   
                    # at least one previous ok contig (new spacer needs to be defined)
                    if (length(bb) > 0) {
                      b <- max(bb)
                      spacing <- get.spacer(dp, i, b) 
                      stopifnot(spacing >= 0) # will always be positive
                      
                      c <- get(paste0("c", b))
                      if (length(c) > 0) {r <- rep(NA.char, spacing)} else {r <- NULL}
                      c <- c(c, r, a)
                      dp[i,"ok"] <- 1
                    }
                  } else {
                    # case overlapping contig is not best contig
                    c <- c
                    dp[i,"ok"] <- 0
                  }
                }
                assign(paste0("c", i), c)
              }
            }
          }
          
          # fill dl
          dl[target == loc & ID == sample & query == name,"ok"] <- dp[i,"ok"]
          dl[target == loc & ID == sample & query == name,"spacer"] <- spacing
          
        }
        
        # fill dl
        dl[target == loc & ID == sample & ok == 1, "order"][order(dl[target == loc & ID == sample & ok == 1, tstart])] <- 1:nrow(dl[target == loc & ID == sample & ok == 1])
        dl[target == loc & ID == sample & ok == 1, "start2"] <- paste(c[1:8], collapse = "")
        dl[target == loc & ID == sample & ok == 1, "end2"] <- paste(rev(rev(c)[1:8]), collapse = "")
        dl[target == loc & ID == sample & ok == 1, "qlen2"] <- length(c)
    
        # move old contig file if overwrite == FALSE
        if (!overwrite) {
          system(command = paste("mv", best, gsub("bestScore.fasta$", "bestScore.fasta.bak", best)))
        }
        
        # write new contig file
        if (any(dp[,"ok"] == 1)) {
          # do not combine contigs if criteria (best-scored contig kept) are not met
          if (any(dp[ok == 1,score] == bestscore)) {
            c <- as.matrix(t(c))
            rownames(c) <- sample
            ape::write.FASTA(ape::as.DNAbin(c), file = best)
          } else {
            cat(paste0(sample, ".", loc, ": contig with best raw score (", bestname, ") discared, replaced with contig(s) with better normalized scores (", paste(unlist(dp[ok == 1,query]), collapse = ", "), ")\n"))
            c <- as.matrix(t(c))
            rownames(c) <- sample
            ape::write.FASTA(ape::as.DNAbin(c), file = best)
          }
        } else {
          writeLines(c(paste0(">", sample), mis), con = best)
        }
      }
    }
  }
  
  # return results View(dl[,-c(10:29)])
  return(dl)
}

# quick sequence alignment
aln <- function(char1, char2, gapOpening = -2, gapExtension = -8, type = "global") {
  require(Biostrings)
  sigma <- nucleotideSubstitutionMatrix(match = 2, mismatch = -1, baseOnly = FALSE)
  aln <- pairwiseAlignment(paste(gsub("", "",toupper(char1)), collapse = ""), 
                           paste(gsub("", "",toupper(char2)), collapse = ""),
                           substitutionMatrix = sigma, gapOpening = gapOpening,
                           gapExtension = gapExtension, scoreOnly = FALSE,
                           type = type)
  aln
}

# print quick sequence alignment
printPairwiseAlignment <- function(alignment, chunksize=60, returnlist=FALSE) {
  require(Biostrings)  
  seq1aln <- pattern(alignment)
  seq2aln <- subject(alignment)
  alnlen  <- nchar(seq1aln)
  starts  <- seq(1, alnlen, by=chunksize)
  n       <- length(starts)
  seq1alnresidues <- 0
  seq2alnresidues <- 0
  for (i in 1:(n)) {
    chunkseq1aln <- substring(seq1aln, starts[i], starts[i]+chunksize-1)
    chunkseq2aln <- substring(seq2aln, starts[i], starts[i]+chunksize-1)
    gaps1 <- countPattern("-",chunkseq1aln)
    gaps2 <- countPattern("-",chunkseq2aln)
    seq1alnresidues <- seq1alnresidues + chunksize - gaps1
    seq2alnresidues <- seq2alnresidues + chunksize - gaps2
    if (returnlist == 'FALSE')
    {
      print(paste(chunkseq1aln,seq1alnresidues))
      print(paste(chunkseq2aln,seq2alnresidues))
      print(paste(' '))
    }
  }
  if (returnlist == 'TRUE')
  {
    vector1 <- s2c(substring(seq1aln, 1, nchar(seq1aln)))
    vector2 <- s2c(substring(seq2aln, 1, nchar(seq2aln)))
    mylist <- list(vector1, vector2)
    return(mylist)
  }
}

# alignment wrapper
align <- function(ref = "consDalbergia_4c_2984.fasta", folder = "best.contigs.63.2984/H121/", loc = "ID_658.merged", type = "global-local", winsize = 200) {
  fas <- ape::read.FASTA(ref)
  con <- list.files(folder, pattern = paste0(loc, "[.]"), full.names = TRUE)
  con <- con[grep(".bestScore.fasta$", con)]
  x <- as.character(ape::read.FASTA(con))[[1]]
  y <- as.character(fas)[[grep(paste0(loc, "$"),names(fas))]]
  a <- aln(x,y,type = type) # type %in% "global", "local", "overlap", "global-local", and "local-global" 
  printPairwiseAlignment(a, winsize)
  res <- list()
  res[["ref"]] <- toupper(paste(y, collapse = ""))
  res[["tar"]] <- toupper(paste(x, collapse = ""))
  invisible(res)
}


##########################################################################################

## Case I: no suf.exonerate files with alignment sugar (get.range = FALSE)
if (!get.range) {
  
  ## Get exonerate scores
  p.allscore <- list.files(folder, pattern = suf.allscore, full.names = TRUE)
  
  # score
  lscores <- list()
  lnames <- list()
  llengths <- list()
  for (locus in p.allscore) {
    lscores[[locus]] <- get.scores(locus)
    lnames[[locus]] <- get.names(locus)
    if (get.length) {
      llengths[[locus]] <- get.lengths(locus, split = length.split, string = length.string, where = length.where) # ordered by decreasing score
    }
  }
  # lscores <- sapply(p.allscore, FUN = get.scores) # takes too long for >3000 loci
  names(lscores) <- names(lnames) <- gsub(locus.prefix, "", gsub(paste0(suf.allscore, "$"), "", basename(names(lscores))))
  
  ## Get contig lengths
  if (get.length) {
    # llengths <- sapply(p.allscore, FUN = get.lengths, split = length.split, string = length.string, where = length.where) # takes too long for >3000 loci
    names(llengths) <- names(lscores)
  } else {
    # read lengths from contigs
    for (locus in names(lscores)) {
      c <- gsub("LOCUS", locus, gsub("SAMPLE", sample, gsub(locus.prefix, "", gsub(paste0(suf.allscore, "$"), "", cpath))))
      if (file.exists(c)) {
        llengths[[locus]] <- unname(lengths(ape::read.FASTA(c))[lnames[[locus]]]) # ordered by decreasing length, reordered by decreasing score
        if (all(is.na(llengths[[locus]]))) llengths[[locus]] <- numeric()
      } else {
        llengths[[locus]] <- numeric()
      }
    }
  }
  
  ## Combine
  d <- data.frame(ID = basename(folder),
                  locus = names(lscores), # defined 'locus' as name of target
                  ncontigs = unlist(lapply(lscores, length)),
                  scores = unlist(lapply(lscores, paste, collapse = ", ")),
                  stringsAsFactors = FALSE,
                  row.names = NULL)
  d$whichbest <- as.numeric(sapply(strsplit(d$scores, split = ", "), function(x) which.max(as.numeric(x))))
  d$bestscore <- as.numeric(sapply(seq(nrow(d)), FUN = function(x) {sapply(strsplit(d$scores[x], split = ", "), "[", d$whichbest[x])}))
  
  d$lengths <- unlist(lapply(llengths, paste, collapse = ", "))
  d$lengths[d$lengths == "NA"] <- NA
  d$bestlength <- as.numeric(sapply(seq(nrow(d)), FUN = function(x) {sapply(strsplit(d$lengths[x], split = ", "), "[", d$whichbest[x])}))
  for (i in format[!format==tname]) d[,i] <- NA
  
  # Handle NAs
  d$ncontigs[d$bestscore == 0] <- 0
  d[which(d$ncontigs == 0), c("scores","bestscore","whichbest","lengths","bestlength")] <- NA
  
}


## Case II: alignment sugar in suf.exonerate files
if (get.range) {
  p.exonerate <- list.files(folder, pattern = suf.exonerate, full.names = TRUE)
  if (length(p.exonerate) == 0) get.range <- FALSE
  
  lranges <- data.table(array(NA, dim = c(0, 7+max.contigs),
                              dimnames = list(NULL, c(format, paste0("o", 1:max.contigs)))))
  for (locus in p.exonerate) {
    # get ordered alignment range data table
    drange <- data.table(get.ranges(locus, format = format, split1 = q.split, split2 = t.split))
    setkey(drange, target, score)
    drange <- drange[order(target, -score)]
    
    # get length of multi-contig alignment overlaps
    for (h in 1:max.contigs) {drange[,paste0("o", h) := as.numeric(NA)]}
    lranges <- rbindlist(list(lranges, get.overlaps(drange, tsname, tename, max.contigs)))
  }
  
  # add sample x locus combinations
  lranges <- data.table(rep(sample, nrow(lranges)), lranges) ; names(lranges)[1] <- "ID"
  lranges[,"locid"] <- interaction(unlist(lranges[,"target"]), unlist(lranges[,"ID"]))
  
  # add query and target (reference) alignment length
  lranges[,"qaln"] <- abs(lranges[,"qend"] - lranges[,"qstart"])
  lranges[,"taln"] <- abs(lranges[,"tend"] - lranges[,"tstart"])
  
  # add normalized alignment score
  lranges[,"score.norm"] <- lranges[,"score"] / lranges[,"taln"]
  
  # order by target locus and decreasing score
  setkey(lranges, target, score)
  lranges <- lranges[order(target, -score)]

  ## Case II.1: Combine non-overlappig contigs passing filters and replace the best contig with the supercontig.
  if (combine) {
    res <- write.combined.contigs(dl = lranges, NA.char = NA.char, cpath = cpath, fpath = fpath, sample = sample, max.contigs = max.contigs, 
                                  tol = tol, min.aln = min.aln, min.qlen = min.qlen, max.overlap = max.overlap, min.norm.score = min.norm.score)
    
    # number of used (single or combined) contigs
    x1 <- res[,.(npassed = sum(ok == 1, na.rm = T)), by = target]
    # table(x1$npassed, useNA= "a")
    
    # number of discarded (potentially paralogous) contigs
    x2 <- res[,.(nfailed = sum(ok == 0, na.rm = T)), by = target]
    # table(x2$nfailed, useNA= "a")
    
    # bestscore
    x3 <- res[ok == 1,.(bestscore = mean(score, na.rm = T)), by = target]
    # range(x3$bestscore) 
    
    # bestscore.norm
    x4 <- res[ok == 1,.(bestscore.norm = mean(score.norm, na.rm = T)), by = target]
    # range(x4$bestscore.norm)
    
    # bestlength
    x5 <- res[ok == 1,.(bestlength = mean(qlen2)), by = target]
    # range(x5$bestlength)
    
    # best taln
    x6 <- res[ok == 1,.(taln = sum(taln)), by = target]
    # range(x6$taln)
    
    d <- merge(data.table(ID = rep(sample, nrow(x1)), x1, x2[,2]), data.table(x3, x4[,2], x5[,2], x6[,2]), all = T, sort = F)
    d$ncontigs <- d$nfailed + ifelse(d$npassed > 0, 1, 0) # passed contigs are counted as one combined contig
    d$locid <- interaction(d$ID, d$target)
    setkey(d, target)
    setcolorder(d, c("ID","target","locid","npassed","nfailed","ncontigs","bestscore","bestscore.norm","bestlength","taln"))
    
  } else {
    ## Case II.2 Do not combine non-overlappig contigs, only summarize alignment sugar and contig overlaps
    res <- lranges
    
    if (get.length) {
      dlengths <- unlist(strsplit(unlist(res[,"query"]), split = length.split))
      res[,"qlen"] <- as.numeric(rep(NA, nrow(res)))
      res[!is.na(query),"qlen"] <- as.numeric(dlengths[grep(length.string, dlengths) + length.where])
    } else {
      # read lengths from contigs
      res[,"qlen"] <- as.numeric(rep(NA, nrow(res)))
      for (locus in unique(unlist(res[,"target"]))) {
        c <- gsub("LOCUS", locus, gsub("SAMPLE", sample, gsub(locus.prefix, "", gsub(paste0(suf.allscore, "$"), "", cpath))))
        if (file.exists(c)) {
          res[target == locus,"qlen"] <- unname(lengths(ape::read.FASTA(c))[unlist(res[target == locus,"query"])]) # ordered by decreasing length, reordered by decreasing score
        } else {
          res[target == locus,"qlen"] <- NA
        }
      }
    }
    
    res[,"passed"] <- ifelse(res[,"taln"] >= min.aln & res[,"qlen"] >= min.qlen & res[,"score.norm"] >= min.norm.score, 1, 0)
    res[,"ok"] <- res[,.(ifelse(score == max(score, na.rm = TRUE) & !all(score == 0), 1, ifelse(all(score==0), as.numeric(NA), 0))), by = target]$V1
    
    # if there are just non-ok contigs for a locus, ok is set to false but the contig is still present in .bestScore.fasta!
    
    x1 <- res[,.(npassed = sum(ok == 1, na.rm = T)), by = target]
    x2 <- res[,.(nfailed = sum(ok == 0, na.rm = T)), by = target]
    x3 <- res[ok == 1,.(bestscore = mean(score, na.rm = T)), by = target]
    x4 <- res[ok == 1,.(bestscore.norm = mean(score.norm, na.rm = T)), by = target]
    x5 <- res[ok == 1,.(bestlength = mean(qlen)), by = target]
    x6 <- res[ok == 1,.(taln = sum(taln)), by = target]
    d <- merge(data.table(ID = rep(sample, nrow(x1)), x1, x2[,2]), data.table(x3, x4[,2], x5[,2], x6[,2]), all = T, sort = F)
    d$ncontigs <- d$nfailed + ifelse(d$npassed > 0, 1, 0) # passed contigs are counted as one combined contig
    d$locid <- interaction(d$ID, d$target)
    setkey(d, target)
    setcolorder(d, c("ID","target","locid","npassed","nfailed","ncontigs","bestscore","bestscore.norm","bestlength","taln"))
  }
}

## Tabulate number of queries (contigs) where exonerate was able to align a query (contig) to a target (reference sequence)
contab <- table(d$ncontigs)
tab <- data.frame("ncontigs" = names(contab), "nloci" = as.numeric(contab))

## Check results
# x <- unique(res[res$ok != res3$ok]$target)
# x <- unname(unlist(d[npassed > 1,"target"]))
# 
# for (i in seq(length(x))) {print(res[target == x[i]]) ; cat("\n") ; align(loc = x[i], type = "global-local") ; scan(quiet = T)}

## Write results
write.table(tab, file = file.path(folder, tabfile), row.names = FALSE, quote = FALSE, sep = "\t")
write.table(d, file = file.path(folder, dfile), row.names = FALSE, quote = FALSE, sep = "\t")

if (get.range) {
  write.table(res, file = file.path(folder, rangefile), row.names = FALSE, quote = FALSE, sep = "\t")
}

