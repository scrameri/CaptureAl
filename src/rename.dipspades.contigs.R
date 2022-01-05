#!/cluster/apps/r/3.6.0_openblas/x86_64/bin/Rscript

### Rename dipspades contigs ###

## Usage: Rscript rename.dipspades.contigs.R <individual folder> <suffix to be omitted in fasta headers>
## Author: simon.crameri@env.ethz.ch, 18.2.2020

## Get arguments
args = commandArgs(trailingOnly=TRUE)

## Check argument
if (length(args)<1) {
  stop("At least one argument needed: 1) name of folder (one per individual) with subfolders (one for each locus) with contigs ; 2) OPTIONAL: individual folder suffix", call.=FALSE)
}

## Set argument
id <- args[1]
sample.suffix <- args[2] # suffix present in <id> but not desired in contig header

## Additional arguments
contig.path <- "/dipspades/consensus_contigs.fasta" # <id>/*/PATHTO/contig.fasta

## Define renaming function
rename.contigs <- function(dir, ind, sample.suffix, contig.path) {
  
  ## Read .fasta
  id <- basename(gsub(paste0(sample.suffix, "$"), "", ind))
  path <- paste0(dir, contig.path)

  if (file.exists(path)) {
	
    if (file.size(path)>0) {

	## Read .fasta
	Sys.chmod(path, mode = "0777", use_umask = TRUE)
	fas <- read.delim(path, header = F, stringsAsFactors = F)
  
	## Rename .fasta
	header.i <- grep("^>", fas$V1)
	fas$V1[header.i] <- paste0(">", id, "__", 1:length(header.i))

	## Write .fasta
	write.table(x = fas, file = path, row.names = F, col.names = F, quote = F)
    } 
  }
}

## Execute renaming function
dirs <- list.dirs(path = id, recursive = F)
sapply(dirs, rename.contigs, ind = id, sample.suffix = sample.suffix, contig.path = contig.path)
