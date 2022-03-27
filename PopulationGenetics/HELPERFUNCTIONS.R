##########################################################
### HELPERFUNCTIONS for PCA, NJ and STRUCTURE analyses ###
##########################################################

# helperfunctions
genind2genlight <- function(gi, verbose = TRUE) {
  
  ## Usage
  # gi       genind object
  # verbose  if TRUE, prints information on progress and biallelic loci
  
  ## Author
  # simon.crameri@usys.ethz.ch, Mar 2022
  
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
  biloc.perc <- 100*(length(biloc) / length(unique(d.loc$loc)))
  d.loc.bi <- subset(d.loc, loc %in% biloc)
  locsel <- paste(d.loc.bi$loc, d.loc.bi$all, sep = ".")
  
  if (verbose) {
    cat(paste0(length(biloc), " / ", length(unique((d.loc$loc))),
               " (", round(biloc.perc, 0), "%) are biallelic\n"))
  }
  
  # take every second allele
  matsub <- mat[,locsel[seq(2,length(locsel),by=2)]]
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
}
gi.filter.per.region <- function(gi, n = 2, split = "_", split.index = 1) {
  
  ## Usage
  # gi          genind object
  # n           integer number of kept loci per region
  # split       character string that is used to split locNames(gi) into region and position
  #             e.g., the locus name "2325merged_570" will be interpreted as belonging to 
  #             region "2325merged" and position "570" if split = "_" and split.index = 1
  # split.index integer number denoting the index of the split locus name that represents the region, 
  #             e.g. sapply(strsplit(locNames(gi), split = split), "[", split.index) returns the region
  
  ## Author
  # simon.crameri@usys.ethz.ch, Mar 2022
  
  gi@other$locus <- factor(sapply(strsplit(locNames(gi), split = split), "[", split.index))
  idx <- unlist(sapply(levels(gi@other$locus), FUN = function(x) {
    snps <- which(gi@other$locus == x)
    if (length(snps) < n) {
      snps
    } else {
      sample(snps, n, replace = FALSE)
    }
  }
  ))
  gisub <- gi[loc = idx, drop = TRUE]
  gisub@other$locus <- factor(sapply(strsplit(locNames(gisub), split = split), "[", split.index))
  return(gisub)
}

filter.gi <- function(gi, biallelic.only = TRUE, expected.only = TRUE, exp.char = c("A","T","C","G"), 
                      min.mac = 2, max.missing = 0.05) {
  
  ## Usage
  # gi                genind object
  # biallelic.only    if TRUE, will only keep loci with exactly 2 alleles
  # expected.only     if TRUE, will only keep loci with <exp.char> alleles 
  # exp.char          character string with the expected alleles
  # min.mac           minimum minor allele count (mac): will only keep loci with mac >= <min.mac>
  # max.missing       will only keep loci with missingness <= <max.missing>
  #                   max.missing >= 1 is interpreted as number of missing genotypes
  #                   max.missing <1 is interpreted as the fraction of missing genotypes
  
  ## Author
  # simon.crameri@usys.ethz.ch, Mar 2022
  
  stopifnot(class(gi)[1] == "genind")
  
  ## Reduce to biallelic SNPs only
  if (biallelic.only) {
    locbi <- which(gi@loc.n.all == 2)
    gibi <- gi[loc = locbi, drop = TRUE]
    cat(paste0("excluded ", nLoc(gi)-length(locbi), " (", round(100*(nLoc(gi)-length(locbi))/nLoc(gi),2), "%) non-biallelic SNPs, ", nLoc(gibi), " remain\n"))
  } else {
    gibi <- gi
  }
  
  ## Exclude SNPs with other than exp.char
  if (expected.only) {
    weird <- sapply(gibi@all.names, function(x) any(!x %in% exp.char))
    gigood <- gibi[loc = which(!weird), drop = TRUE]
    cat(paste0("excluded ", sum(weird), " (", round(100*sum(weird)/nLoc(gibi),2), "%) SNPs with characters deviating from ", paste(exp.char, collapse = ", "), ", ", nLoc(gigood), " remain\n"))
  } else {
    gigood <- gibi
  }
  
  ## Remove SNPs with minor allele count (MAC) less than mac (mac = 2 removes singletons)
  ac <- apply(gigood@tab, 2, sum, na.rm = TRUE)
  locmac <- unique(sapply(strsplit(names(which(ac < min.mac)), split = "[.]"), "[", 1))
  gimin <- gigood[loc = locNames(gigood)[!locNames(gigood) %in% locmac], drop = TRUE]
  cat(paste0("excluded ", length(locmac), " (", round(100*length(locmac)/nLoc(gigood),2), "%) SNPs with minor allele count smaller than ", min.mac, ", ", nLoc(gimin), " remain\n"))
  
  ## Exclude SNPs with > max.missing missingness
  # missingness >= 1 is interpreted as number of missing genotypes
  # missingness <1 is interpreted as percentage of missing genotypes
  if (max.missing >= 1) max.missing <- max.missing / nInd(gimin)
  dmis <- apply(gimin@tab, 2, function(x) {sum(is.na(x))/length(x)})
  locsel <- unique(sapply(strsplit(names(which(dmis <= max.missing)), split = "[.]"), "[", 1))
  if (length(locsel) == 0) {
    print(summary(dmis))
    stop("max.missing is set too low, no variants remain. See missingness summary above.")
  }
  gi.cleaned <- gimin[loc = locsel, drop = TRUE]
  cat(paste0("excluded ", nLoc(gimin)-length(locsel), " (", round(100*(nLoc(gimin)-length(locsel))/nLoc(gimin),2), "%) SNPs with missingness above ", round(100*max.missing, 2), "%, ", nLoc(gi.cleaned), " remain\n"))
  
  ## Return filtered gi
  return(gi.cleaned)
}
