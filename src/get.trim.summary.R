#!/cluster/apps/r/3.6.0_openblas/x86_64/bin/Rscript

## Usage: get.trim.summary.R <trimming directory>

## Author: simon.crameri@env.ethz.ch, May 2019

## Get arguments
args = commandArgs(trailingOnly = TRUE)

## Check arguments
if (length(args)!=1) {
  stop("1 argument required:
       REQUIRED
       1) <logsdir|CHR>: path to trimmed alignments log directory", 
       call.=FALSE)
}

## Set arguments
logsdir <- as.character(args[1])

## Set arguments (for debugging)
#logsdir <- "testfolder.logs"

## Check arguments
stopifnot(dir.exists(logsdir))

## Read log files
logs <- list.files(logsdir, pattern = ".log")
logs <- logs[logs!="params.log"]
if (length(logs) == 0) stop("No alignment trimming *.log files found. You need to specify a trimmed alignments log directory.")
params <- readLines(file.path(logsdir, "params.log"))

loglist <- list()
for (log in logs) {
  loglist[[gsub(".log$", "", log)]] <- readLines(file.path(logsdir, log))
}

## Get stats
# determine trimming type
keptfracs <- lapply(lapply(lapply(loglist, FUN = function(x) {x[grep("^kept", x)]}), FUN = function(x) {sapply(strsplit(sapply(strsplit(x, split = "[(]"), "[", 2), split = "%"), "[", 1)}), FUN = function(x) {as.numeric(x)/100})
if (all(lengths(keptfracs) == 2)) itr <- TRUE else itr <- FALSE

# fraction of kept bases
pkepts <- lapply(keptfracs, prod)

# fracton of masked bases (pbasestrim), fraction of sites with masking (psitestrim)
trimlist <- lapply(loglist, FUN = function(x) {xtr <- x[grep("^trimmed", x)] ; xtr2 <- xtr[!xtr %in% xtr[grep("first|last", xtr)]] ; return(xtr2)})
trimsplit <- lapply(trimlist, FUN = strsplit, split = "[(]")
pbasestrim <- lapply(trimsplit, FUN = function(x) {0.01*as.numeric(sapply(strsplit(sapply(x, "[", 2), split = "[%)]"), "[", 1))})
psitestrim <- lapply(trimsplit, FUN = function(x) {0.01*as.numeric(sapply(strsplit(sapply(x, "[", 3), split = "[%)]"), "[", 1))})

# fraction of samples with masking (pseqstrim) 
if (itr) {
  pseqstrim <- lapply(trimsplit, FUN = function(x) {0.01*as.numeric(sapply(strsplit(sapply(x, "[", 4), split = "[%)]"), "[", 1))})
} else {
  pseqstrim <- lapply(keptfracs, FUN = function(x) {x[0]})
}

## Summary
cat("\n")
cat("Fraction of kept alignment sites:\n")
print(summary(unlist(pkepts)))
cat("\n")
cat("Fraction of trimmed bases:\n")
print(summary(unlist(pbasestrim)))
cat("\n")
cat("Fraction of sites with at least one trimmed base:\n")
print(summary(unlist(psitestrim)))
cat("\n")
if (itr) {
  cat("Fraction of sequences with at least one trimmed base:\n")
  print(summary(unlist(pseqstrim)))
  cat("\n")
}
