#######################################################################################
### Generate 6555 reference sequences for analysis of legume target enrichment data ###
#######################################################################################

# simon.crameri@usys.ethz.ch ; sfcrameri@gmail.com ; Mar 2022

library(ape)

## Define arguments
# Cajanus cajan genome and annotation downloaded from https://www.ncbi.nlm.nih.gov/assembly/GCF_000340665.1
genomefile <- "ncbi-genomes-2019-01-28/GCF_000340665.1_C.cajan_V1.0_genomic.fna"
gfffile <- "ncbi-genomes-2019-01-28/GCF_000340665.1_C.cajan_V1.0_genomic.gff"

# original 7201 reference sequences and position stats (provided by S. Zoller)
reffile <- "fasta/orig.fragments.of.all.probes.dalb.pset.two.sp.v1.fasta"
stats <- "fasta/orig.fragments.of.all.probes.dalb.pset.two.sp.v1.csv" # this was produced based on the files above, in lines 162-290 (takes long!)


## Define helperfunctions
# get physical distance between probes
get.physical.dist <- function(df, LG) {
  d.lg <- df[df$LG == LG, ]
  d.lg <- d.lg[order(d.lg$start),]
  if (nrow(d.lg)>1) {
    d.dist <- sapply(2:nrow(d.lg), function(x, df.lg = d.lg) {df.lg[x,"start"]-df.lg[x-1,"end"]-1})
    d.dist <- c(d.dist, NA)
  } else {
    d.dist <- NA
  }
  return(d.dist)
}

# merge physically close reference sequences
merge.close.seqs <- function(df, refseqs, ref, thr = 1, nchar = "n", verbose = TRUE) {
  
  # author: simon.crameri@env.ethz.ch, Mar 2019
  # arguments:
  #df = annotations for refseqs: "ID","LG","start","end","name","seqid","g.start","g.end","type","gene","product","length","disttonext"
  #refseqs = DNAbin list containing reference sequences (subset of ref)
  #ref = DNAbin list containing reference genome sequences (used to create refseqs)
  #thr = positive integer specifying which refseqs should be merged. Any refseq that is physically nearer than thr basepairs from the previous refseq will be merged with the previous refseq
  
  # identify probes above threshold
  find.probes <- function(dd, thr) {
    d.near <- which(dd$disttonext < thr)
    dd[d.near,"ID"]
  }
  
  # remove leading ns
  remove.leading.n <- function(x, nchar = "n") {
    repeat {
      if (x[1] %in% nchar) x <- x[-1] else break
    }
    return(x)
  }
  
  remove.trailing.n <- function(x, nchar = "n") {
    repeat {
      if (x[length(x)] %in% nchar) x <- x[-length(x)] else break
    }
    return(x)
  }
  
  # initialize
  df.merged <- df
  df.merged$mergedwith <- NA
  refseqs.merged <- refseqs
  continue <- FALSE
  
  # for debugging
  # df.merged <- d.probes.merged
  # refseqs.merged <- refseqs
  # continue <- FALSE
  
  # check input
  stopifnot(all.equal(nrow(df.merged), length(refseqs.merged)),
            all.equal(df.merged$ID, names(refseqs.merged)),
            all(tapply(df.merged$start, INDEX = df.merged$LG, FUN = function(x) {all.equal(x, sort(x))})))
  
  # repeat merging until criterion <thr> is met
  repeat{
    if (!continue) {
      ids <- find.probes(dd = df.merged, thr = thr)
      if (length(ids)==0) break
      cat(ids[1], "[", length(ids), "]\n")
    } else {
      if (length(ids)==0) break
      cat("repeat", ids[1], "[", length(ids), "]\n")
    }
    
    df1 <- df.merged[df.merged$ID == ids[1],]
    df2 <- df.merged[which(df.merged$ID == ids[1])+1,]
    
    seq1 <- as.character(as.matrix(refseqs.merged[df1$ID]))[1,]
    seq2 <- as.character(as.matrix(refseqs.merged[df2$ID]))[1,]
    
    # check whether seq2 is fully contained in seq1
    if (df1$disttonext < 0 & abs(df1$disttonext)+1 >= length(seq2)) {
      fully <- TRUE
    } else {
      fully <- FALSE
    }
    
    # get reference scaffold (read only if it changes)
    if (!exists("seqid")) {
      seqid <- names(ref)[grep(df1$seqid, names(ref))]
      refg <- ref[seqid]
      stopifnot(length(refg)==1)
      refseq <- as.character(as.matrix(refg))[1,]
    }
    if (is.na(df1$seqid)) {
      newseqid <- unique(na.omit(subset(df, LG == df1$LG)$seqid))
      if (length(newseqid) == 0) newseqid <- sapply(strsplit(names(ref)[grep(paste(unique(df1$LG, df2$LG), collapse = "|"),names(ref))], split = " "), "[", 1)
      df1$seqid <- df2$seqid <- newseqid
    }
    if (seqid != names(ref)[grep(df1$seqid, names(ref))]) {
      seqid <- names(ref)[grep(df1$seqid, names(ref))]
      refg <- ref[seqid]
      stopifnot(length(refg)==1)
      refseq <- as.character(as.matrix(refg))[1,]
    }
    
    # merge physically close sequences
    if (!fully) {
      seq.merged <- remove.trailing.n(remove.leading.n(refseq[df1$start:df2$end], nchar = nchar), nchar = nchar)
      stopifnot(all.equal(seq2, seq.merged[(length(seq.merged)-df2$length+1):length(seq.merged)]),
                all.equal(seq1, seq.merged[1:df1$length]))
      dd.merged <- data.frame(ID = df1$ID, LG = df1$LG, start = df1$start, end = df2$end,
                              name = gsub(df1$end, df2$end, df1$name),seqid = df1$seqid,
                              g.start = paste(unique(c(unlist(strsplit(df1$g.start, split = ", ")), unlist(strsplit(df2$g.start, split = ", ")))), collapse = ", "), 
                              g.end = paste(unique(c(unlist(strsplit(df1$g.end, split = ", ")), unlist(strsplit(df2$g.end, split = ", ")))), collapse = ", "), 
                              type = paste(unique(c(unlist(strsplit(df1$type, split = ", ")), unlist(strsplit(df2$type, split = ", ")))), collapse = ", "), 
                              gene = paste(unique(c(unlist(strsplit(df1$gene, split = ", ")), unlist(strsplit(df2$gene, split = ", ")))), collapse = ", "),
                              product = paste(unique(c(unlist(strsplit(df1$product, split = ", ")), unlist(strsplit(df2$product, split = ", ")))), collapse = ", "),
                              length = length(seq.merged), disttonext = df2$disttonext,
                              # mergedwith = df1$mergedwith, gsub("(.*)bp_([0-9]+)_.*", "\\2", df2$name),
                              mergedwith = gsub("^NA,", "", paste(df1$mergedwith, gsub("(.*)bp_([0-9]+)_.*", "\\2", df2$name), sep = ",")),
                              stringsAsFactors = FALSE) ##
    } else {
      seq.merged <- remove.leading.n(seq1, nchar = nchar)
      dd.merged <- df1
      if (is.na(df2$disttonext)) dd.merged$disttonext <- NA else dd.merged$disttonext <- df.merged[which(df.merged$ID == df2$ID)+1,"start"]-dd.merged$end-1
      dd.merged$mergedwith <- NA ##
    }
    
    # update df.merged and refseqs.merged
    df.merged <- df.merged[-which(df.merged$ID == df2$ID),]
    df.merged[which(df.merged$ID == df1$ID),] <- dd.merged
    stopifnot(nrow(df.merged) > 0)
    
    seq.dnabin <- as.list(as.DNAbin(seq.merged))
    names(seq.dnabin) <- dd.merged$ID
    refseqs.merged <- refseqs.merged[-which(names(refseqs.merged) == df2$ID)]
    refseqs.merged[df1$ID] <- seq.dnabin
    
    # prepare for next iteration
    if (!is.na(dd.merged$disttonext)) continue <- ifelse(dd.merged$disttonext < thr, TRUE, FALSE) else continue <- FALSE
  }
  return(list(df.merged = df.merged, refseqs.merged = refseqs.merged))
}

################################################################################

## Read genome
ref <- read.FASTA(genomefile)

# ## Read gff annotation file and split attributes (takes long)
# load("d.probes.rda") # d.probes, d.probes.lg, d.probes.near
# gff <- read.gff(gfffile, na.strings = c(".", "?"))
# gffattr <- strsplit(gff$attributes, split = ";")
# gffattr2 <- lapply(gffattr, function(x) {y <- strsplit(x, split = "=") ; n <- sapply(y, "[", 1) ; v = sapply(y, "[", 2) ; names(v) <- n ; return(v)})
# gffattrnames <- unique(unlist(lapply(gffattr2, function(x) {names(x)})))
# for (name in gffattrnames) {
#   cat(name, "\n")
#   ls <- unlist(lapply(gffattr2, "[", name))
#   stopifnot(length(ls) == nrow(gff))
#   gff[,name] <- ls
#   rm(list=c("name","ls"))
# }
# save(gff, file = "gff.rda")
# 
# ## Add LG and unplaced genomic scaffolds
# genome <- read.FASTA(genomefile)
# regions <- sapply(strsplit(names(genome), split = " "), "[", 1)
# 
# gff$LG <- NA
# lgs <- na.omit(unique(gff$chromosome))
# for (i in lgs) {
#   cat(i,"\n")
#   id <- which(lgs == i)
#   start <- which(gff$chromosome == i)
#   end <- (which(gff$chromosome == lgs[id+1])-1)[1]
#   if (length(start) == 1) {
#     range <- start:end
#     gff$LG[range] <- ifelse(nchar(i) == 1, paste0("0", i), i)
#   } else {
#     for (j in start) {
#       id <- which(start == j)
#       cat(id, "/", length(start), "\n")
#       begin <- j
#       end <- start[id+1]-1
#       if (is.na(end)) end <- nrow(gff)
#       range <- begin:end
#       scaffold <- sapply(strsplit(sapply(strsplit(names(genome)[which(regions == gff[begin,"seqid"])], split = "Scaffold"), "[", 2), split = ","), "[", 1)
#       gff$LG[range] <- paste0("Scaffold", scaffold)
#     }
#   }
# }
# save(gff, file = "gff.rda")
# 
# ## Match regions in d.probes to annotations in gff (any gene annotation that overlaps with any probe's start:end will be saved to a list <res>)
# load("gff.rda")
# 
# res <- list()
# for (i in seq(nrow(d.probes))) {
#   cat(i, "/", nrow(d.probes), "\n")
#   df <- d.probes[i,]
#   dg <- subset(gff, LG == df$LG & ! type %in% c("region","match") & start <= df$end & end >= df$start)
#   if (nrow(dg) > 0) {
#     dg$IDprobe <- df$ID
#     res[[df$ID]] <- dg
#   } else {
#     res[[df$ID]] <- NA
#   }
# }
# save(res, file = "gffprobes.rda")
# 
# ## Number of hits per probe
# load("gffprobes.rda")
# nhits <- lapply(res, function(x) {if (length(x)==1) 0 else nrow(x)})
# hist(unlist(nhits), breaks = 100, xlim = c(0,max(unlist(nhits))), main = "Number of gff entries overlapping with probes", xlab = "n")
# sum(unlist(nhits)==0) # this many regions have no hit
# sum(unlist(nhits)==1) # this many regions have 1 hit (pseudogene)
# sum(unlist(nhits)>=1) # this many regions have 1 or more hits (genes + pseudogenes with exons)
# 
# ## Gene, Product names, etc.
# repr <- data.frame(t(repr))
# repr$ID <- as.character(repr$ID)
# repr[is.na(repr[,"ID"]),"ID"] <- names(repr[is.na(repr[,"ID"]),"seqid"])
# for (i in seq(ncol(repr))) repr[,i] <- as.character(repr[,i])
# sum(is.na(repr$ID))
# 
# ## Merge and write to disk
# # merge
# d.probes.merged <- merge(d.probes, repr, by = "ID", all = TRUE, sort = F)
# 
# # map LG to seqids (updates seqid where no annotation overlapped with a probe)
# maplg2seqid <- sapply(unique(d.probes.merged$LG), FUN = function(x) {if (length(grep("Scaffold", x) == 1)) sapply(strsplit(names(ref)[grep(paste0("\\", x, "\\b"), names(ref))], split = " "), "[", 1) else NA})
# maplg2seqid <- lapply(maplg2seqid, function(x) {if (length(x) == 0) NA else x})
# stopifnot(all(lengths(maplg2seqid)==1))
# tomap <- unlist(maplg2seqid)[is.na(maplg2seqid)]
# tomap <- names(tomap)[grep("Scaffold", names(tomap), invert = TRUE)]
# maplg2seqid[tomap] <- as.list(sapply(tomap, FUN = function(x) unique(na.omit(d.probes.merged[d.probes.merged$LG == x, "seqid"]))))
# stopifnot(all(lengths(maplg2seqid)==1))
# library(plyr)
# seqids <- data.frame(old = d.probes.merged$seqid, 
#                      new = mapvalues(d.probes.merged$LG, from = unique(d.probes.merged$LG), to = unlist(maplg2seqid)), 
#                      stringsAsFactors = FALSE)
# stopifnot(length(which(seqids$new != seqids$old)) == 0)
# d.probes.merged$seqid <- seqids$new
# 
# # some stats
# sum(is.na(d.probes.merged$seqid))
# sum(!d.probes.merged$product %in% c("", NA)) # 6616/7201 with gene product
# length(unique(d.probes.merged$LG)) # 719 unique LG
# length(unique(d.probes.merged$seqid)) # 715/719 mapped
# sum(d.probes.merged$product %in% "") # 35/7201 without product (pseudogene)
# sum(is.na(d.probes.merged$product)) # 550/7201 without annotation
# d.probes.merged <- d.probes.merged[order(d.probes.merged$LG, d.probes.merged$start),]
# 
# ## Add average physical distance between probes
# l.dist <- lapply(unique(d.probes.merged$LG), function(x) {get.physical.dist(df = d.probes.merged, LG = x)})
# length(unlist(l.dist)) # 7201 loci
# d.probes.merged$disttonext <- unlist(l.dist)
# d.probes.merged$length <- d.probes.merged$end - d.probes.merged$start + 1
# 
# ## Check that actual probe length confirms <end>-<start>+1
# d.length <- lengths(refseqs)
# stopifnot(all.equal(names(refseqs), d.probes.merged$ID))
# 
# # 17 probes had Ns and are shorter than <end>-<start>+1
# d.probes.n <- data.frame(d.probes.merged[d.length != d.probes.merged$length,c("ID","start","end","length")], "actuallength" = d.length[d.length != d.probes.merged$length])
# nrow(d.probes.n)
# save(d.probes.n, file = "d.probes.n.rda")
# nrow(d.probes.n)
# 
# # fix lengths
# d.probes.merged$length <- d.length

# ## save results
# save(d.probes.merged, file = "d.probes.merged.rda")
# write.csv2(d.probes.merged, file = "d.probes.merged.csv", row.names = F, quote = F)

## Read original 7201 reference sequences
refseqs <- read.FASTA(reffile)
names(refseqs) <- sapply(strsplit(sapply(strsplit(names(refseqs), split = "__"), "[", 1), split = "_"),"[",2)

## Read stats of original 7201 reference sequences
d.probes.merged <- read.csv2(stats) # data.frame of dim [7201, 13]
d.probes.merged$ID <- as.character(d.probes.merged$ID)

## Average physical distance between probes
mean(tapply(d.probes.merged$disttonext, INDEX = d.probes.merged$LG, FUN = mean, na.rm = TRUE), na.rm = TRUE) # mean of per-LG average physical distance
mean(d.probes.merged$disttonext, na.rm = T) # mean physical distance

summary(d.probes.merged$length)
summary(d.probes.merged$disttonext) # there are overlapping target loci
sum(d.probes.merged$disttonext < 0, na.rm = T) # 378 overlapping target loci
sum(d.probes.merged$disttonext <= 100, na.rm = T) # 648 <= 100 bp
sum(d.probes.merged$disttonext <= 1000, na.rm = T) # 1073 <= 1000 bp
sum(d.probes.merged$disttonext <= 5000, na.rm = T) # 2927 <= 5000 bp

## Average coverage of genome
length(ref) # number of linkage groups and scaffolds: 36536
length(unique(d.probes.merged$LG)) # number of linkage groups and scaffolds covered by probes: 719
length(unique(d.probes.merged$LG[grep("^Scaffold", d.probes.merged$LG, invert = T)])) # number of linkage groups covered by probes: 1
length(unique(d.probes.merged$LG[grep("^Scaffold", d.probes.merged$LG)])) # number of scaffolds covered by probes: 708

## Merge probes that are < 100 bp physical distance apart 
#physically overlapping probes are merged, non-overlapping probes are complemented (up to thr-1 bases) with the reference genome sequence

# set threshold: any reference sequence that lies closer than <thr> to a neighboring reference sequence will be merged
# using the C. cajan genome as a spacer
thr = 100
res <- merge.close.seqs(df = d.probes.merged, refseqs = refseqs, ref = ref, thr = thr)

stopifnot(!any(duplicated(res$df.merged$name)),
          all.equal(sapply(strsplit(sapply(strsplit(res$df.merged$name, split = "__"), "[", 1), split = "_"), "[", 2), names(res$refseqs.merged)))
names(res$refseqs.merged) <- res$df.merged$id <- paste(res$df.merged$seqid, "LG", res$df.merged$LG, res$df.merged$start, res$df.merged$end, "ID", res$df.merged$ID, sep = "_")

## Average physical distance between probes
mean(res$df.merged$disttonext, na.rm = T) # 52'802.42 bp mean physical distance
mean(tapply(res$df.merged$disttonext, INDEX = res$df.merged$LG, FUN = mean, na.rm = TRUE), na.rm = TRUE) # mean of per-LG average physical distance: 35'415.12 bp

summary(res$df.merged$length) ; nrow(res$df.merged) # 100: 87-1424 bp (6555) ; 250: 87-2204 (6262) ; 500: 87-3037 (5956) ; 1000: 87-6673 (5498)
summary(res$df.merged$disttonext) # minimum threshold should be >= thr: 100: 100 - 1'306'862 ; 250: 250 - 1'306'862 ; 500: 502-1'306'862 ; 1000: 1003 - 1'306'862

## Physical distance between probes distribution
sum(res$df.merged$disttonext < 0, na.rm = T)     # 100: 0      ;  250: 0      300: 0    500: 0      1000: 0     overlapping target loci
sum(res$df.merged$disttonext <= 100, na.rm = T)  # 100: 2      ;  250: 0      300: 0    500: 0      1000: 0     <= 100 bp
sum(res$df.merged$disttonext <= 500, na.rm = T)  # 100: 599    ;  250: 306    300: 247  500: 0      1000: 0     <= 500 bp
sum(res$df.merged$disttonext <= 600, na.rm = T)  # 100: 699    ;  250: 406    300: 347  500: 100    1000: 0     <= 1000 bp
sum(res$df.merged$disttonext <= 1000, na.rm = T) # 100: 1057   ;  250: 764    300: 705  500: 458    1000: 0     <= 1000 bp
sum(res$df.merged$disttonext <= 2000, na.rm = T) # 100: 1622   ;  250: 1329   300: 1270 500: 1023   1000: 565   <= 2000 bp
sum(res$df.merged$disttonext <= 2500, na.rm = T) # 100: 1791   ;  250: 1498   300: 1439 500: 1192   1000: 734   <= 2500 bp
sum(res$df.merged$disttonext <= 5000, na.rm = T) # 100: 2281   ;  250: 1988   300: 1929 500: 1682   1000: 1224  <= 5000 bp

## Write results
save(res, file = paste0("d.probes.merged", thr, ".rda"))
write.FASTA(res$refseqs.merged, file = "fasta/Cajanus_cajan_6555reg.fasta")
write.csv2(res$df.merged, file = "fasta/Cajanus_cajan_6555reg.csv", row.names = F, quote = F)
