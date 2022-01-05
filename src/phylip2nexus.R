#!/cluster/apps/r/3.6.0_openblas/x86_64/bin/Rscript

## Usage: phylip2nexus.R <.phy concatenated alignment file> <.model partition file contaiting limits between loci>

## Author: simon.crameri@env.ethz.ch, Nov 2018

## Get arguments
args = commandArgs(trailingOnly = TRUE)

## Check arguments
if (length(args)!=2) {
  stop("Two arguments required:
        REQUIRED
        1) path to concatenated alignment file <.phy>
        2) path to partition file <.model>", call.=FALSE)
}

## Set arguments
f.phy <- as.character(args[1])
f.mod <- as.character(args[2])

## Set arguments (for debugging)
# f.phy = "mafft.108.1483.c0.5.d0.25.trim-4-20-8-1.phy"
# f.mod = "mafft.108.1483.c0.5.d0.25.trim-4-20-8-1.model"

## Additional arguments
gapchar <- "?"
mischar <- "-"
cat("using", gapchar, "to encode gaps and", mischar, "to encode missing data\n")
ext <- paste0(".", sapply(strsplit(basename(f.phy), split = "[.]"), function(x) {rev(x)[1]}))
phybase <- gsub(paste0(ext, "$"), "", basename(f.phy))
ofile <- paste0(phybase, ".nex") # name of output file

## Check arguments
stopifnot(file.exists(f.phy), file.exists(f.mod))

## Read <.model> file containing character partitions (e.g. limits between loci)
mod <- readLines(f.mod)

## Read <.phy> file containing the concatenated alignments
phy <- readLines(f.phy)

## Check and adjust <.model> file
mod.adj <- paste0(gsub(",", "-", gsub("=", ":", gsub("DNA,", "", mod))), ",")

## Check and adjust <.phy> file
phyinfo <- unlist(strsplit(phy[1], split = " "))
split <- unlist(strsplit(phy[2], split = " "))
stopifnot(all.equal(length(phy)-1, as.numeric(phyinfo[1])),
          all.equal(nchar(split[length(split)]), as.numeric(phyinfo[2])))
physeq <- phy[-1]

phylist <- vector("list", length(phy)-1)
for (i in seq(length(physeq))) {
  split <- unlist(strsplit(physeq[i], split = " "))
  phylist[[i]] <- split[split != ""]
}
stopifnot(unique(unlist(lapply(phylist, length))) == 2)
phylabels <- sapply(phylist, "[", 1)
phybegin <- max(nchar(phylabels)) + 2
physpace <- phybegin - nchar(phylabels)

for (i in seq(length(physeq))) {
  phylist[[i]][3] <- paste0(gsub("-", "_", phylist[[i]][1]), paste(rep(" ", physpace[i]), collapse = ""), phylist[[i]][2])
}
physeqs <- sapply(phylist, "[", 3)

## Create NEXUS output
ntax <- as.numeric(phyinfo[1])
nchar <- as.numeric(phyinfo[2]) 
datatype <- "DNA"
missing <- mischar
gap <- gapchar

# DATA chunk
d.start <- "BEGIN DATA;"
d.dim <- paste0("\tDIMENSIONS NTAX=", ntax, " NCHAR=", nchar, ";")
d.format <- "\tFORMAT"
d.dtype <- paste0("\t\tDATATYPE=", datatype)
d.missing <- paste0("\t\tMISSING=", missing)
d.gap <- paste0("\t\tGAP=", gap, ";")
d.matrix <- "\tMATRIX"
# physeqs here
d.end <- ";\nEND;\n"

# SETS chunk
s.start <- "BEGIN SETS;"
s.cpart <- paste0("\tCHARPARTITION LOCI=")
# mod.adj here
s.end <- c(";", "END;") # needs to have this ending

## Write lines
string <- c("#NEXUS", d.start, d.dim, d.format, d.dtype, d.missing, d.gap, d.matrix, physeqs, d.end, s.start, s.cpart, mod.adj, s.end)
writeLines(string, con = ofile)
