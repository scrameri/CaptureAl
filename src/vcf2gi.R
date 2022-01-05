#!/usr/bin/Rscript

## Load dependencies
suppressPackageStartupMessages(library(vcfR))
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(adegenet))

## Get arguments
args = commandArgs(trailingOnly=TRUE)

## Get arguments
if (length(args) != 2) {
  stop("Two arguments needed: 
       REQUIRED
       1) <vcf|CHR>: path to .vcf file ; 
       2) <ref|CHR>: path to reference .fasta", call.=FALSE)
}

## Set arguments
vcf = as.character(args[1])
ref = as.character(args[2])

## Check arguments
stopifnot(file.exists(ref))

## Read reference .fasta file
dna <- read.FASTA(ref)
dna <- dna[sort(names(dna))]
# regions <- names(dna)[start:length(dna)]

## Additional Arguments
ofolder <- tools::file_path_sans_ext(vcf)
region <- basename(ofolder)

exp.chars = c("a", "c", "t", "g")
thr.het = 1
filter1 = FALSE
filter2 = FALSE
verbose = TRUE

abbreviate = TRUE
abbreviate.fun <- function(x = vcf@fix[,"CHROM"]) {
  gsub("[.]", "_", sapply(strsplit(x, split = "_ID_"), "[", 2))
}

################################################################################

## Define function
vcf2gi <- function(vcffile, dna, ofolder, abbreviate = FALSE, abbreviate.fun = NULL,
                   exp.chars = c("a","c","t","g"),
                   thr.het = 1, filter1 = FALSE, filter2 = FALSE,
                   verbose = TRUE) {
  
  ## USAGE
  # vcffile         path to .vcf file
  # dna             DNAbin object (full reference sequences)
  # ofolder         path to output folder (where .fasta and .rda files will be written, and .vcf file will be moved)
  # abbreviate      replace region names [vcf@fix[,"CHROM"]] by ouput of abbreviate.fun(vcf@fix[,"CHROM"]) 
  # abbreviate.fun  function used to abbreviate vcf@fix[,"CHROM"]
  # exp.chars       expected characters in DNAbin alignment
  # thr.het         if  < 1, will use reference allele instead of ambiguity codes at sites with > thr.het heterozygous genotypes
  # filter1         if TRUE, will use reference allele instead of ambiguity codes at sites without homozygous genotypes for minor allele (= pseudo-monomorphic sites)
  # filter2         if TRUE, will use reference allele instead of alternate alleles at sites with > 2 alleles (= non-biallelic sites)
  # verbose         if TRUE, prints more verbose output
  
  ## Read .vcf
  stopifnot(file.exists(vcffile))
  vcf <- try(read.vcfR(file = vcffile, verbose = verbose), silent = T)
  if (inherits(vcf, "try-error")) {
    vcf <- new(Class = "vcfR") ; vcf@gt <- array(dim=c(0,1))
  }
  nsnp <- nrow(vcf@gt)
  
  if (nsnp > 0) {
    cat("found", nsnp, "variants\n")
    
    # check reference regions
    regs <- gsub("contig=ID=", "",
                 lapply(queryMETA(vcf, element = "contig"), "[", 1))
    
    stopifnot(all(regs %in% names(dna)))
    
    region <- gsub(".vcf$", "", basename(vcffile))
    
    ## Create chromR
    chrom <- suppressWarnings(create.chromR(name = region, 
                                            vcf = vcf, 
                                            seq = dna[region], 
                                            verbose = verbose))
    
    ## Generate DNAbin of SNPs
    dnabin <- vcfR2DNAbin(chrom, extract.indels = T, unphased_as_NA = T,
                          extract.haps = F, consensus = T,
                          ref.seq = dna[region], start.pos = 1,
                          verbose = verbose)
    dnabin.mat <- as.character(dnabin)
    
    # all SNP posistions
    snppos <- which(apply(dnabin.mat, 2, function(x) {!all(is.na(x))}))
    if(length(snppos) != dim(vcf)["variants"]) cat(paste0(dim(vcf)["variants"], " variants expected but only ", length(snppos), " SNPs found (indels in variants?)\n"))
    stopifnot(ncol(dnabin) == lengths(dna[region]))
    stopifnot(nrow(dnabin) %in% c(dim(vcf)["gt_cols"]-1, 2*(dim(vcf)["gt_cols"]-1)))
    
    # replace 'n' with NA
    dnabin.mat <- apply(dnabin.mat, 2, function(x) {x[x=="n"] <- NA ; return(x)})
    
    # exclude sites where > thr.het of individuals are heterozygous
    if (thr.het < 1) {
      get.het.freq <- function(x, het.chars) {
        if (!all(is.na(x))) {
          freq <- sum(x %in% het.chars) / sum(!is.na(x))
        } else {
          freq <- 0
        }
        return(freq)
      }
      all.chars <- sort(names(table(dnabin.mat)))
      het.chars <- all.chars[!all.chars %in% exp.chars]
      het.freqs <- apply(dnabin.mat, 2, get.het.freq, het.chars = het.chars)
      pos.excl1 <- which(het.freqs > thr.het)
      snppos.filtered1 <- snppos[!snppos %in% pos.excl1]
      cat(paste0("excluded ", length(pos.excl1), " heterozygous sites, ", length(snppos.filtered1), 
                 " (", round(100*length(snppos.filtered1)/length(snppos),1), "%) remain\n"))
    } else {
      pos.excl1 <- numeric()
      snppos.filtered1 <- snppos
    }
    if (length(pos.excl1)>0) dnabin.mat[,pos.excl1] <- NA
    
    # exclude sites where only an expected character and an ambiguitiy code remain
    if (filter1) {
      get.polymorphic.sites <- function(x, mat, exp.chars) {
        if (!all(is.na(mat[,x]))) {
          chars <- names(table(mat[,x]))
          nexp <- length(chars[chars %in% exp.chars])
          incl <- ifelse(nexp < 2, F, T)
        } else {
          incl <- T
        }
        return(incl)
      }
      pos.excl2 <- which(!sapply(1:ncol(dnabin.mat), FUN = get.polymorphic.sites, 
                                 mat = dnabin.mat, exp.chars = exp.chars))
      snppos.filtered2 <- snppos.filtered1[!snppos.filtered1 %in% pos.excl2]
      cat(paste0("excluded ", length(pos.excl2), " pseudo-monomorphic sites, ", length(snppos.filtered2), 
                 " (", round(100*length(snppos.filtered2)/length(snppos),1), "%) remain\n"))
    } else {
      pos.excl2 <- numeric()
      snppos.filtered2 <- snppos.filtered1
    }
    if (length(pos.excl2)>0) dnabin.mat[,pos.excl2] <- NA
    
    # exclude triallelic SNPs
    if (filter2) {
      get.biallelic <- function(x, mat, exp.chars) {
        if (!all(is.na(mat[,x]))) {
          chars <- names(table(mat[,x]))
          nexp <- length(chars[chars %in% exp.chars])
          incl <- ifelse(nexp <= 2, T, F)
        } else {
          incl <- T
        }
        return(incl)
      }
      pos.excl3 <- which(!sapply(1:ncol(dnabin.mat), FUN = get.biallelic, 
                                 mat = dnabin.mat, exp.chars = exp.chars))
      snppos.filtered3 <- snppos.filtered2[!snppos.filtered2 %in% pos.excl3]
      cat(paste0("excluded ", length(pos.excl3), " non-biallelic sites, ", length(snppos.filtered3), 
                 " (", round(100*length(snppos.filtered3)/length(snppos),1), "%) remain\n"))
    } else {
      pos.excl3 <- numeric()
      snppos.filtered3 <- snppos.filtered2
    }
    if (length(pos.excl3)>0) dnabin.mat[,pos.excl3] <- NA
    
    # insert reference allele at invariant / no-SNP posistions
    stopifnot(all.equal(ncol(dnabin.mat), length(unlist(as.character(dna[region])))))
    
    insert.refbase <- function(mat, ref, NA.char = "N") {
      stopifnot(is.matrix(mat), class(ref) == "DNAbin")
      refseq <- unlist(as.character(ref))
      stopifnot(ncol(mat) == length(refseq))
      
      ins.ref <- function(mat, ref, index) {
        if (all(is.na(mat[,index]))) {
          mat[,index] <- refseq[index]
        }
        return(mat[,index])
      }
      
      mat.ins <- sapply(1:ncol(mat), FUN = ins.ref, mat = mat, ref = ref)
      mat.ins[is.na(mat.ins)] <- NA.char
      return(mat.ins)
    }
    dnabin.mat <- insert.refbase(mat = dnabin.mat, ref = dna[region])
    dnabin2 <- as.DNAbin(dnabin.mat)
    # image.DNAbin(dnabin2
    
    ## Create genind
    # with ambiguity codes for heterozygous individuals (not meaningful)
    # only retains sites where there are at least 2 expected bases (a,c,t,g)
    # gi <- DNAbin2genind(dnabin2, polyThres = 0, exp.char = all.chars)
    
    # with allele counts for heteozygous genotypes
    # rename loci (no . in locus names allowed)
    if (abbreviate) {
      vcf@fix[,"CHROM"] <- abbreviate.fun(vcf@fix[,"CHROM"])
    }
    
    # create genind (code heterozygous individuals)
    gi <- suppressWarnings(vcfR2genind(vcf))
    
    # re-include entirely non-type individual(s) 
    ids <- colnames(vcf@gt)[-1]
    im <- ids[!ids%in%indNames(gi)]
    im.tab <- array(data = integer(), dim = c(length(im), ncol(gi@tab)), dimnames = list(im, colnames(gi@tab)))
    gi <-  genind(tab = rbind(gi@tab, im.tab), pop = gi@pop, 
                  prevcall = gi@call, ploidy = rep(2, length(ids)), type = gi@type, 
                  strata = gi@strata, hierarchy = gi@hierarchy)
    gi <- gi[ids]
    
    # get allele infos
    fix <- data.frame(vcf@fix[,c("CHROM","POS","ID","REF","ALT","QUAL"),drop=F], stringsAsFactors = FALSE)
    dp <- extract.info(vcf, element = "DP")
    fixdp <- function(x) {as.numeric(unique(unlist(strsplit(x, split = ","))))}
    fix$DP <- sapply(dp, FUN = fixdp)
    
    # match all.names to actual alleles
    gi.an <- sapply(1:nLoc(gi), FUN = function(x) {
      unlist(strsplit(paste(fix[x,c("REF","ALT")], collapse=","), split = ","))
    }, simplify = FALSE)
    names(gi.an) <- names(gi@all.names)
    gi@all.names <- gi.an
    colnames(gi@tab) <- paste(as.character(gi@loc.fac), 
                              unlist(sapply(seq(nLoc(gi)), FUN = function(x) {
                                gi@all.names[[x]][1:gi@loc.n.all[[x]]]
                              })),
                              sep = ".")
    
    # add DP and QUAL info in @other slot
    gi@other$DP <- list(fix$DP)
    gi@other$QUAL <- list(as.numeric(fix$QUAL))
    names(gi@other$QUAL) <- names(gi@other$DP) <- unique(fix$CHROM)
    
    #table(gi@loc.n.all)
    # gi@all.names
  } else {
    dnabin2 <- gi <- NULL # if no variants detected
  }
  
  ## Generate results
  res <- list()
  # res[["chrom"]] <- chrom
  res[["dnabin"]] <- dnabin2
  res[["gi"]] <- gi
  # res[["vcf"]] <- vcf
  return(res)
  
}

genind2genlight <- function(gi, verbose = TRUE) {
  
  # check gi
  stopifnot(class(gi)[1] == "genind")
  mat <- gi@tab
  
  # get factors for locus and position
  loc.all <- strsplit(colnames(mat), split = "[.]")
  loc.names <- sapply(loc.all, "[", 1)
  loc.list <- strsplit(loc.names, split = "_")
  loc.id <- sapply(loc.list, "[", 1)
  fac.id <- factor(loc.id, levels = unique(loc.id))
  fac.pos <- factor(loc.names, levels = unique(loc.names))
  
  if (verbose) cat("filtering for biallelic loci...\n")
  d.loc <- data.frame(loc = sapply(loc.all, "[", 1),
                      all = sapply(loc.all, "[", 2),
                      stringsAsFactors = F)
  d.loc$loc <- factor(d.loc$loc, levels = unique(d.loc$loc))
  loctab <- table(d.loc$loc, d.loc$all)
  biloc <- rownames(loctab)[which(rowSums(loctab) == 2)]
  biloc.perc <- 100*(length(biloc) / nlevels(d.loc$loc))
  d.loc.bi <- subset(d.loc, loc %in% biloc)
  locsel <- paste(d.loc.bi$loc, d.loc.bi$all, sep = ".")
  
  if (verbose) {
    cat(paste0(length(biloc), " / ", nlevels(d.loc$loc),
               " (", round(biloc.perc, 0), "%) are biallelic\n"))
  }
  
  # take every second allele
  matsub <- mat[,locsel[seq(2,length(locsel),by=2)],drop=F]
  matsubnameslist <- strsplit(colnames(matsub), split = "[.]")
  
  fac.pos <- sapply(matsubnameslist, "[", 1)
  fac.pos <- factor(fac.pos, levels = unique(fac.pos))
  
  fac.id <- sapply(strsplit(sapply(matsubnameslist, "[", 1), split = "_"), "[", 1)
  fac.id <- factor(fac.id, levels = unique(fac.id))
  
  locselloc <- d.loc.bi$loc[seq(1,nrow(d.loc.bi),by=2)]
  locselref <- locsel[seq(1,length(locsel),by=2)]
  locselalt <- locsel[seq(2,length(locsel),by=2)]
  
  # store reference and alternate allele
  ref <- toupper(sapply(strsplit(locselref, split = "[.]"), "[", 2))
  alt <- toupper(sapply(strsplit(locselalt, split = "[.]"), "[", 2))
  
  # create new locus name with reference and alternate allele
  locselrefalt <- paste0(locselloc, ".", ref, "/", alt)
  colnames(matsub) <- locselrefalt
  
  # generate gl (genlight) object
  if (verbose) cat("generating genlight object...\n")
  l.snp <- apply(matsub, 1, function(x) {l <- list(x) ; return(l)})
  l.snp <- lapply(l.snp, unlist)
  gl <- new(Class = "genlight", l.snp)
  gl@ploidy <- gi@ploidy
  locNames(gl) <- locselrefalt
  
  # attach pop and Clade
  gl@pop <- gi@pop
  if (!is.null(gi@other)) gl@other <- gi@other
  return(gl)
  cat("Done!\n")
}

################################################################################

## Create output log
t1 <- Sys.time()
cat("============================\n")
cat("======= vcf2gi.R LOG =======\n")
cat("============================\n")
cat("\nStarting time:", as.character(t1), "\n")
cat(R.version.string, "\n")
cat(paste0("ape version = ", packageVersion("ape")), "\n")
cat(paste0("adegenet version = ", packageVersion("adegenet")), "\n")
cat(paste0("vcfR version = ", packageVersion("vcfR")), "\n")
cat("\nArguments:\n")
cat(paste0("vcf = ", vcf, "\n"))
cat(paste0("reference = ", ref, "\n"))
cat(paste0("ofolder = ", ofolder, "\n"))
cat(paste0("exp.chars = ", paste(exp.chars, collapse = ", "), "\n"))
cat(paste0("thr.het = ", thr.het, "\n"))
cat(paste0("filter1 = ", filter1, "\n"))
cat(paste0("filter2 = ", filter2, "\n"))
cat(paste0("verbose = ", verbose, "\n"))
cat(paste0("abbreviate = ", abbreviate, "\n"))
cat(paste0("abbreviate.fun = "))
print(abbreviate.fun)
cat("\n")

## Create output folder
if (!dir.exists(ofolder)) dir.create(ofolder)

## Debugging
# load("Dalbergia_582_249059.filtered/gi_582_144484.rda")
# region=names(dna)[1]
# vcffile=paste0(file.path(vcffolder, region), ".vcf")
# exp.chars = c("a","c","t","g")
# thr.het=1
# verbose=T

################################################################################

## Process one region
# get input and output file
locname <- abbreviate.fun(region)
cat(locname, "[", which(names(dna) == region), "/", length(dna), "]\n")

odir <- ofolder
ifile <- vcf

# move ofile back to dirname(vcf)
ofile <- file.path(odir, basename(ifile))
if (file.exists(ofile) & !file.exists(ifile)) {
  system(command = paste("mv", ofile, dirname(vcf)))
}
stopifnot(file.exists(vcf))

# wait for input file to be written
repeat(
  if (!file.exists(ifile)) {
    # cat(".")
    Sys.sleep(5)
  } else {
    break
  }
)

if (!dir.exists(odir)) dir.create(odir)

# get genind object
dd <- vcf2gi(vcffile = ifile, dna = dna, ofolder = ofolder,
             abbreviate = abbreviate, abbreviate.fun = abbreviate.fun,
             exp.chars = exp.chars,
             thr.het = thr.het, filter1 = filter1, filter2 = filter2,
             verbose = FALSE)

if (!is.null(dd$gi)) {
  # order individuals
  dd$gi <- dd$gi[order(indNames(dd$gi))]
  
  # rename loci
  # locNames(dd$gi) <- paste0(locname, "_", locNames(dd$gi))
  
  # get genlight object  
  # dd$gl <- genind2genlight(dd$gi, verbose = FALSE)
  # dd$gl <- suppressWarnings(vcfR2genlight(dd$vcf))
  
  # save as .rda
  assign(region, dd)
  save(list = region, file = file.path(odir, paste0(region, ".rda")))
  
  # move .vcf to output directory  
  system(command = paste("mv", ifile, odir))
  
  # write FASTA
  write.FASTA(dd$dnabin, file = file.path(odir, paste0(region, ".fasta")))
} else {
  cat("ERROR: no variants read / passed filters\n")
}

# clean up
rm(list=c("dd","region","odir","ifile",region))

## Finish
t2 <- Sys.time()
cat("\nFinish time:", as.character(t2), "\n")
print(t2-t1)

