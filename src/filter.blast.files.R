#!/cluster/apps/r/3.6.0_openblas/x86_64/bin/Rscript

#########################
### Filter BLAST hits ###
#########################

## Usage
# run script without arguments to see the arguments it takes

## Details
# RULE to select best out of multiple hits per query
# 1) alignment length must be at least 0.01*min.perc.alignment.length*query.length AND at least min.alignment.length
# 2) subject length must be at least min.subject.length AND at most max.subject.length
# 3) perc.identity of alingment must be at least min.perc.identity
# 4) if there are still >= 2 hits remaining, take the hit with the lower evalue

## Author
# simon.crameri@env.ethz.ch, Apr 2019

## Load required library
suppressPackageStartupMessages(library("data.table")) # enables fast handling of large data frames

## Get arguments
args = commandArgs(trailingOnly=TRUE)

## Check arguments
options(warning.length = 2000L)
if (!length(args) %in% c(1,6)) {
  stop("At least one argument needed: 
       REQUIRED
       1) f                           character         path to .blast file
       
       OPTIONAL
       2) min.perc.alignment.length   numeric [0, 100]  minimum percentage of query.length in alignment.length [DEFAULT: 0]
       3) min.alignment.length        numeric [1, ]     minimum length of aligned portion of subject and query [DEFAULT: 80]
       4) min.subject.length          numeric [1, ]     minimum subject length [DEFAULT: 100]
       5) max.subject.length          numeric [1, ]     maximum subject length [DEFAULT: 10000]
       6) min.perc.identity           numeric [0, 100]  minimum percentage of identity between aligned portion of subject and query [DEFAULT: 80]",
       call.=FALSE)
}

## Set arguments
f <- as.character(args[1])
min.perc.alignment.length <- as.numeric(as.character(args[2])) ; if (is.na(min.perc.alignment.length)) min.perc.alignment.length <- 0
min.alignment.length <- as.numeric(as.character(args[3])) ; if (is.na(min.alignment.length)) min.alignment.length <- 80
min.subject.length <- as.numeric(as.character(args[4])) ; if (is.na(min.subject.length)) min.subject.length <- 100
max.subject.length <- as.numeric(as.character(args[5])) ; if (is.na(max.subject.length)) max.subject.length <- 10000
min.perc.identity <- as.numeric(as.character(args[6])) ; if (is.na(min.perc.identity)) min.perc.identity <- 80

## Echo arguments
cat("\n### FILTER BLAST HITS ###\n")
cat("\nlooking for .blast results in", f, "\n")
cat("minimum percentage of query length in alignment length set to", min.perc.alignment.length, "\n")
cat("minimum alignment length set to", min.alignment.length, "\n")
cat("minimum subject length set to", min.subject.length, "\n")
cat("maximum subject length set to", max.subject.length, "\n")
cat("minimum percentage of identity between aligned portion of subject and query set to", min.perc.identity, "\n")

## Check arguments
stopifnot(file.exists(f), 
          is.numeric(min.perc.alignment.length), min.perc.alignment.length >= 0, min.perc.alignment.length <= 100,
          is.numeric(min.alignment.length), min.alignment.length > 0,
          is.numeric(min.subject.length), min.subject.length > 0,
          is.numeric(max.subject.length), max.subject.length > 0,
          is.numeric(min.perc.identity), min.perc.identity >= 0, min.perc.identity <= 100)

## Process blast data (expects output generated from <blast.fasta.seqs.sh>)
blastnames <- c("query.id", "query.length", "subject.id", "subject.length", 
                "alignment.length", "perc.identity", "q.start", "q.end", 
                "s.start", "s.end", "evalue", "bit.score", "score", "subject.title")

## Define helperfunction
# select passed and best contigs
select.contig <- function(x, choice = "evalue",
                          min.perc.align.len = min.perc.alignment.length, min.align.len = min.alignment.length, 
                          min.subject.len = min.subject.length, max.subject.len = max.subject.length, 
                          min.perc.ident = min.perc.identity) {
                            
    contigs <- t[i=which(t$id == x), ]
    ls <- which(contigs$alignment.length >= 0.01*min.perc.align.len*unique(contigs$query.length) &
                 contigs$alignment.length >= min.align.len &
                  contigs$subject.length >= min.subject.len & 
                   contigs$subject.length <= max.subject.len & 
                    contigs$perc.identity >= min.perc.ident)
    
    # set passed
    passed <- numeric(nrow(contigs))
    passed[ls] <- 1
    set(t, i=which(t$id==x), j = "passed", value = passed)
    
    # set best
    # if (length(ls) > 1) ls <- ls[which.min(contigs[,get(choice)][ls])] # old way
    if (length(ls) > 1) {
      choose <- which(contigs[ls,get(choice)] == min(contigs[ls,get(choice)]))
      if (length(choose) > 1) {
        ls <- ls[which.max(contigs[choose,subject.length])]
      } else {
        ls <- ls[choose]
      }
    }
    best <- numeric(nrow(contigs))
    best[ls] <- 1
    set(t, i=which(t$id==x), j = "best", value = best)    
}

## Process file
dat <- data.table(array(NA, dim = c(0, length(blastnames)+4), dimnames = list(NULL, c("ID", "id", blastnames, "best", "passed")))) # add "Freq" here if computed
 
## Read BLAST output
t <- data.table(read.table(f, stringsAsFactors = F, header = F))
colnames(t) <- blastnames

## Attach metadata from filename (expects metadata in <f>, evalue after <contigbase> and sample (s) after '.on.')
# f looks like: "${blastdir}/${sample}/${refbase}.on.${sample}.${contigbase}.${evalue}.blast"
flist <- strsplit(f, split = "[.]")
evalue <- sapply(flist, function(x) {rev(x)[2]})
s <- sapply(flist, function(x) {x[which(x == "on")+1]})
#evalue <- gsub(".blast", "", sapply(strsplit(f, split = paste0(".", contigbase, "."), fixed = TRUE), "[", 2))
#s <- gsub(paste0(".", contigbase, ".", evalue, ".blast"), "", sapply(strsplit(f, split = ".on.", fixed = TRUE), "[", 2)))

t$ID <- s
t$id <- interaction(t$ID, t$query.id)
colnames(t) <- c(blastnames, "ID", "id")

## Add number of hits per query
# d.hits <- data.table(tapply(X = t$query.id, INDEX = t$ID, FUN = table)[[1]])
# d.hits$ID <- s
# d.hits$id <- interaction(d.hits$ID, d.hits$V1)
# colnames(d.hits) <- c("query.id", "Freq", "ID", "id")
# t <- merge(t, d.hits[,c("id","Freq")], by = c("id"), all.x = T, sort = F)
 
## Select best out of multiple hits per query
# set order by id and evalue
setcolorder(t, c("ID",  "id", blastnames)) # add "Freq" here if computed above
setorderv(t, cols = c("id", "alignment.length", "subject.length", "perc.identity", "evalue"), order = c(1, -1, 1, -1, 1))
t$best <- numeric(nrow(t))
t$passed <- numeric(nrow(t))

# select contigs according to the rule defined in select.contig function (see above)
s <- sapply(unique(as.character(t$id)), FUN = select.contig)

## Bind data
dat <- rbindlist(list(dat, t))

## Clean up
rm(list=c("t", "s")) # add d.hits here if computed

## Write results
obase1 <- paste(gsub(paste0(".", evalue, ".blast"), "", f), evalue, min.perc.alignment.length, min.alignment.length, min.subject.length, max.subject.length, min.perc.identity, sep = ".")
write.table(subset(dat, best == 1), file = paste0(obase1, ".besthits.txt"), quote = FALSE, row.names = FALSE, sep = "\t")
write.table(subset(dat, passed == 1), file = paste0(obase1, ".passedhits.txt"), quote = FALSE, row.names = FALSE, sep = "\t")
write.table(dat, file = paste0(obase1, ".allhits.txt"), quote = FALSE, row.names = FALSE, sep = "\t")

cat("Done!\n")
