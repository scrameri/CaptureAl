#!/cluster/apps/r/3.6.0_openblas/x86_64/bin/Rscript

## Usage: Rscript filter.blast.vs.self.R <.blast output>
## Author: simon.crameri@env.ethz.ch, Apr 2019

## Define and sanity-check arguments
args = commandArgs(trailingOnly=T)

## Check arguments
if (! length(args)  == 1) {
  stop("1 argument needed:
       REQUIRED
       1) <file|CHR>: .vs.self.blast file",
       call. = FALSE)
}

##
# warn if file extension != ".blast"
filesplit <- strsplit(args[1], split = "[.]")
fileext <- sapply(filesplit, "[", length(unlist(filesplit)))
if (fileext != "blast") warning("File input does not have the .blast extension, might not be in correct format!")

## Set arguments
file <- args[1]

## Set arguments (for debugging)
# file <- "test.vs.self.blast"

## Additional arguments
verbose <- TRUE

## Read data
tit <- read.delim(file, nrows = 1, header = F)
tit <- unlist(strsplit(as.character(tit$V1), split = ", "))
tit <- gsub("[.][.]", ".", gsub(" ", ".", tit))
tit[1] <- "query.id"
dat <- read.delim(file, skip = 1, header = F, stringsAsFactors = F)
names(dat) <- tit
if (verbose) cat("found", nrow(dat), "hits\n")

## Filter for non-self hits
nonself <- numeric()
for (i in 1:nrow(dat)) {
  ls <- dat[i,]
  if (ls[,"query.id"] != ls[,"subject.id"]) nonself <- c(nonself, i) 
}
dat2 <- dat[nonself,]
if (verbose) cat("filtered", nrow(dat)-nrow(dat2), "self hits\n")

## Filter for non-reciprocal hits (remove subject-query <-> query-subject redundancy) (keep hit with longest alignment length)
dat3 <- dat2
nonreciprocal <- numeric()
for (i in 1:nrow(dat3)) {
  # cat(i,"/",nrow(dat3),"\n")
  ls <- dat3[i,]
  if (all(is.na(dat3[i,]))) break else {
    if (paste(ls[,"query.id"], ls[,"subject.id"]) %in% paste(dat3[,"subject.id"], dat3[,"query.id"])) {
      dtest <- dat3[c(i, which(paste(dat3[,"subject.id"], dat3[,"query.id"]) == paste(ls[,"query.id"], ls[,"subject.id"]))),]
      dat3 <- dat3[-which(rownames(dat3) %in% rownames(dtest)[-which.max(dtest[,"alignment.length"])]),]
    }
  }
}
if (verbose) cat("filtered", nrow(dat2)-nrow(dat3), "reciprocal hits\n")
if (verbose) cat("\n", nrow(dat3), "hits remain!\n\n")

## order filtered data according to query length
dat3 <- dat3[order(dat3$query.id),]

## Output filtered .blast table
write.table(dat3, file = paste0(file, ".filtered"), quote = F, row.names = F, sep = "\t")

