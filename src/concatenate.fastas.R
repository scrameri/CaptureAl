#!/cluster/apps/r/3.6.0_openblas/x86_64/bin/Rscript

############################################################
## Concatenate alignments in .fasta format to .phy format ##
############################################################

## Usage
# run script without arguments to see the arguments it takes

## Author
# simon.crameri@env.ethz.ch, Apr 2019

## Load required library
suppressPackageStartupMessages(library(ape))

## Get arguments
args = commandArgs(trailingOnly = TRUE)

## Check arguments
options(warning.length = 2000L)
if (! length(args) %in% c(1)) {
  stop("1 argument needed:
       REQUIRED
       1) <folder|CHR>: path to folder with .fasta alignments"
       , call.=FALSE)
}

## Set arguments
folder <- as.character(args[1])

## Set arguments (for debugging)
# folder <- "testfolder"

## Additional arguments
suffix <- ".fasta"     # will look for files with this suffix in <folder>
NA.char <- c("-")      # will count the proportion of <NA.char> relative to all sites to get the individual-based supermatrix occupancy (<folder.occ> output)
package.size = 1E5     # will procude matrices up to <package.size> alignment size and paste all matrices together in a final step
model.prefix <- "DNA"  # will attach this string before basename(fasta) in <folder.model> output
model.suffix <- "gene" # will attach this string after basename(fasta) in <folder.model> output
oname <- paste0(basename(folder), ".phy") # name of main output file (.phylip format)
coordfile <- gsub(".phy$", ".coords", oname)
modelfile <- gsub(".phy$", ".model", oname)
occfile <- gsub(".phy$", ".occ", oname)
logfile <- gsub(".phy$", ".log", oname)

## Check output file
if (file.exists(oname)) stop("output file ", oname, " already exists!")

## Define helperfunctions
create.labs <- function(x, maxlablen) {
  paste0(x, paste0(rep(" ", maxlablen - nchar(x) + 1), collapse = ""), collapse = "")
}
get.occ <- function(x, NA.char = "-") {
  length(x[!x %in% NA.char])/length(x)
}

## Loop through files
fastas <- list.files(folder, pattern = suffix, full.names = TRUE)

## Get first alignment
# read
cat("\nreading alignments from", folder, "...\n\n")
# cat(paste0(basename(fastas[1]), " [", which(fastas == fastas[1]), " / ", length(fastas), "]\r"))
fas <- read.FASTA(fastas[1])

# check that it is alignmed data
stopifnot(length(unique(lengths(fas))) == 1)
fas <- as.character(as.matrix(fas))

# create label and space
maxlablen <- max(nchar(rownames(fas)))
labs <- as.character(sapply(rownames(fas), FUN = create.labs, maxlablen))
origlabs <- rownames(fas)
nseqs <- nrow(fas)
rownames(fas) <- labs

# initiate alignment length and package number
dlen <- data.frame(locus = basename(fastas[1]), length = ncol(fas))
pck.nb <- 1

for (fasta in fastas[-1]) {
  
  cat(paste0(" [", which(fastas == fasta), " / ", length(fastas), "] [", pck.nb, "]\r"))
  
  # read
  fas2add <- read.FASTA(fasta)
  
  # check that it is alignmed data
  stopifnot(length(unique(lengths(fas2add))) == 1)
  fas2add <- as.character(as.matrix(fas2add))
  
  # check that labels match
  stopifnot(all.equal(rownames(fas2add), origlabs))
  
  # add alignment length
  dlen <- rbind(dlen, data.frame(locus = basename(fasta), length = ncol(fas2add)))
  
  # bind
  if ((ncol(fas) + ncol(fas2add)) > package.size) {
    assign(paste0("fas", pck.nb), fas)
    fas <- fas2add
    pck.nb <- pck.nb + 1
  } else {
    fas <- cbind(fas, fas2add)
  }
}

assign(paste0("fas", pck.nb), fas)

## Write .model file
dlen$to <- cumsum(dlen$length)
dlen$from <- dlen$to - dlen$length + 1
dlen <- dlen[,c("locus", "length", "from", "to")]
dlen$model <- paste0(model.prefix, ",", gsub(suffix, "", dlen$locus), "_", model.suffix, seq(nrow(dlen)), "=", dlen$from, ",", dlen$to)
write.table(dlen, file = coordfile, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(dlen$model, file = modelfile, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

## Paste everything together and write output .phy
# header
ncols <- 0
fas.all <- array(NA, dim = c(nseqs, 0), dimnames = list(labs, NULL))
cat("\n\n")
for (i in 1:(pck.nb)) {
  
  cat(paste0("binding package ", i, " / ", pck.nb, "\r"))
  
  ncols <- ncols + ncol(get(paste0("fas", i)))
  fas.all <- cbind(fas.all, get(paste0("fas", i)))
}

write.table(paste(nseqs, ncols), file = oname, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "")

# labels and all sequences in a single row
cat("\n\nwriting data to", oname, "... ")
write.table(fas.all, file = oname, quote = FALSE, row.names = TRUE, col.names = FALSE, sep = "", append = TRUE)
cat("done!\n")

## Write .occ file
docc <- apply(X = fas.all, MARGIN = 1, FUN = get.occ, NA.char = NA.char)
write.table(data.frame(taxon = origlabs, occupancy = docc), file = occfile, 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

## Be verbose and write .log file
l1 <- paste("Supermatrix dimensions:", nseqs, "taxa,", length(fastas), "loci and", ncols, "aligned columns\n")
l2 <- paste("Overall matrix occupancy:", paste0(round(100*mean(docc),2), "%\n\n"))
cat("\n\n")
cat(l1)
cat("\n")
cat(l2)
writeLines(text = c(l1, l2), con = logfile)
            
## Old way: does not work because there seems to be no way to paste columns without any delimiter
# newfile <- paste0(oname, ".add")
# if (file.exists(oname)) {
#   write.table(fas, file = newfile, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "")
#   system(paste("paste", oname, newdat, ">", paste0(oname, ".tmp")))
#   system(paste("mv", paste0(oname, ".tmp"), oname))
# } else {
#   write.table(fas, file = oname, quote = FALSE, row.names = TRUE, col.names = FALSE, sep = "")
# }


