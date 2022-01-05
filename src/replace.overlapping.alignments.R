#!/cluster/apps/r/3.6.0_openblas/x86_64/bin/Rscript

## Usage: replace.overlapping.alignments.R

## Value: replaces old files in <alndir> with new files ind <mrgdir> according to <mrglist>

## Author: simon.crameri@env.ethz.ch, Apr 2019

## Get arguments
args = commandArgs(trailingOnly = TRUE)

## Check arguments
if (!length(args) %in% c(3,4)) {
  stop("3 arguments required (4 taken): 
       REQUIRED
       1) <alndir|CHR>: path to directory with overlapping fasta alignments to be replaced (mafft.NIND.NLOC.Y) ; 
       2) <mrgdir|CHR>: path to directory with corresponding merged (and filtered) fasta alignments ('mafft.overlapping.X.Y') ;
       3) <mrglist|CHR>: path to overlapping loci list ('mafft.NIND.NLOC.Y.consZ.vs.self.blast.filtered.list') ;
       
       OPTIONAL
       4) <revert|BOOLEAN>: if TRUE, will revert any previous filtering [DEFAULT: TRUE]", 
       call.=FALSE)
}

## Set arguments
alndir <- as.character(args[1])
mrgdir <- as.character(args[2])
mrglist <- as.character(args[3])

revert <- as.logical(as.character(args[4]))
if (is.na(revert)) revert <- TRUE

## Set arguments (for debugging)
# alndir = "testfolder"
# mrgdir = "mafft.overlapping.3"
# mrglist = "testfolder.cons-0.5-0.vs.self.blast.filtered.list"
# revert <- TRUE

## Check arguments
stopifnot(dir.exists(alndir),
          dir.exists(mrgdir),
          file.exists(mrglist))

## Additional arguments
# paths
odir <- file.path(alndir, "loci_overlapping") # name of output directory with replaced FASTA files
pdfdir <- paste0(alndir,".viz")
logdir <- paste0(alndir,".logs")
olog <- "replaced.txt" # file that keeps track of the alignment replacement

# infixes and suffixes
mrginfix <- ".merged" # file infix of merging (corresponds to <infix> argument in align.contigs.of.overlapping.loci.sh)
mrgsuffix <- unique(sapply(strsplit(list.files(mrgdir), split = mrginfix), function(x) {rev(x)[1]})) # file extension of files in mrgdir, something like ".all.aln.etr.itr.fasta" 
stopifnot(length(mrgsuffix) == 1)
mrginsuffix <- paste0(mrginfix, mrgsuffix) # complete suffix of merged alignments in mrgdir

# other
noverlap.max <- 5 # maximum number of overlapping loci (corresponds to implementation in align.contigs.of.overlapping.loci.sh)
verbose <- TRUE

## Define helperfunctions
do.revert <- function(dir, odir, olog) {
  if (dir.exists(dir)) {
    odirs <- list.dirs(dir)
    rmdir <- odirs[grep(odir, odirs)] # remove this directory after files have been moved
    if (length(rmdir) == 1) {
      # remove log file <olog>
      rmfiles <- readLines(file.path(rmdir, olog))
      
      # revert previous filtering (move filterd files up one directory
      mvfiles <- list.files(rmdir)[grep(olog, list.files(rmdir), invert = TRUE)]
      for (file in rmfiles) {
        file.remove(file.path(dir, file))
      }
      for (file in mvfiles) {
        system(paste("mv", file.path(rmdir, file), dir))
      }
      if (list.files(rmdir) == olog) {system(paste("rm -r", rmdir))}
    }
  }
}

do.move <- function(dir, odir, mrgdir, oldfiles, newfiles, olog) {
  if (dir.exists(dir)) {
    if (!dir.exists(odir)) dir.create(odir)
    oldfiles <- unique(oldfiles)
    newfiles <- unique(newfiles)
    if (length(oldfiles) > 0) {
      for (oldfile in oldfiles) {
        cmd <- paste("mv", file.path(dir, oldfile), odir)
        system(cmd)
      } 
    }
    if (length(newfiles) > 0) {
      for (newfile in newfiles) {
        cmd <- paste("cp", file.path(mrgdir, newfile), dir)
        system(cmd)
      }
    }
    write.table(newfiles, file.path(odir, olog), quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
}

## Check for duplicated alignments A <-> B and B <-> A
dmrg <- readLines(con = mrglist)
lmrg <- strsplit(dmrg, split = " ")
if (exists("df.mrg")) rm(df.mrg)
for (i in 1:noverlap.max) {
  if (!exists("df.mrg")) {
    df.mrg <- data.frame(sapply(lmrg, "[", i))
  } else {
    df.mrg <- cbind(df.mrg, data.frame(sapply(lmrg, "[", i)))
  }
  names(df.mrg)[i] <- paste0("locus", i)
  df.mrg[,i] <- as.character(df.mrg[,i])
}
df.mrg$mrgfile <- paste0(df.mrg$locus1, mrginsuffix)
df.mrg$original <- dmrg
df.mrg$ordered <- NA
for (i in seq(nrow(df.mrg))) {
  df.mrg$ordered[i] <- paste(sort(na.omit(as.character(df.mrg[i,1:(ncol(df.mrg)-3)]))), collapse = " ")
}
dups <- which(duplicated(df.mrg$ordered))
if (length(dups) > 0) df.mrg.uniq <- df.mrg[-dups,] else df.mrg.uniq <- df.mrg
if (verbose) {
  cat("\nfound", length(dups), "duplicated alignment(s):", 
      paste(df.mrg$ordered[dups], collapse = ", "), "\n")
  cat("\n")
}

## Get (filtered) merged alignments (w/o duplicates) to be used as replacement
newfiles <- df.mrg.uniq$mrgfile
if (!all(newfiles %in% list.files(mrgdir))) {
  notfound <- df.mrg.uniq$mrgfile[!df.mrg.uniq$mrgfile %in% list.files(mrgdir)]
  if (verbose) {
   cat(paste0("\nfound ", length(notfound), " missing merged alignments in ", mrgdir, ":\n"))
   print(notfound)
   cat("these should correspond to the filtered merged alignments\n\n")
  }
  newfiles <- newfiles[!newfiles %in% notfound]
}

## Get overlapping alignments to be replaced
alnfas <- list.files(alndir, pattern = mrgsuffix)
if (!length(alnfas) > 1) stop(paste0("No files with extension <", mrgsuffix, "> (suffix in ", mrgdir, ") found in ", alndir, "!"))
df.mrg.uniq.filt <- subset(df.mrg.uniq, mrgfile %in% newfiles)
oldfiles <- character()
for (i in seq(nrow(df.mrg.uniq.filt))) {
  oldfiles <- c(oldfiles, paste0(na.omit(as.character(df.mrg.uniq.filt[i,1:(ncol(df.mrg.uniq.filt)-3)])), mrgsuffix))
}

## Reverse any previous filtering
if (revert) {
  # LOGS
  do.revert(dir = logdir, odir = file.path(logdir, "loci_overlapping"), olog = olog)

  # PDF
  do.revert(dir = pdfdir, odir = file.path(pdfdir, "loci_overlapping"), olog = olog)
  
  # FASTA
  do.revert(dir = alndir, odir = file.path(alndir, "loci_overlapping"), olog = olog)
}

## Move oldfiles to <dir>/odir and newfiles to <dir>
if (verbose) {
  cat(paste0("replacing ", length(unique(oldfiles)), " overlapping alignments with ", length(unique(newfiles)), " merged alignments\n\n"))
}

# LOGS
do.move(dir = logdir, odir = file.path(logdir, "loci_overlapping"), mrgdir = paste0(mrgdir, ".logs"), olog = olog,
        oldfiles = c(gsub(".fasta$", ".err", oldfiles),gsub(".fasta$", ".log", oldfiles)),
        newfiles = c(gsub(".fasta$", ".err", newfiles),gsub(".fasta$", ".log", newfiles)))

# PDF
do.move(dir = pdfdir, odir = file.path(pdfdir, "loci_overlapping"), mrgdir = paste0(mrgdir, ".viz"), olog = olog,
        oldfiles = gsub(".fasta$", ".pdf", oldfiles),
        newfiles = gsub(".fasta$", ".pdf", newfiles))
        
# FASTA
do.move(dir = alndir, odir = file.path(alndir, "loci_overlapping"), mrgdir = mrgdir, olog = olog,
        oldfiles = oldfiles,
        newfiles = newfiles)
