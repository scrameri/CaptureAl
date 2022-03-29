###################################################################
### Extract probes of Fabaceae1005 and Dalbergia2396 probe sets ###
###################################################################

# simon.crameri@usys.ethz.ch ; scrameri@gmail.com ; Mar 2022

library(ape)

## Define helperfunction
attach.probeseqs <- function(df, df.merged, df.probes, f, gff, verbose = TRUE) {
  
  # check input
  stopifnot("ID" %in% names(df.probes),
            "start" %in% names(df.probes),
            "end" %in% names(df.probes),
            "ID" %in% names(df.merged),
            "mergedwith" %in% names(df.merged),
            "pname" %in% names(df.probes),
            "seq" %in% names(df.probes),
            inherits(f, "DNAbin"),
            "gene" %in% names(gff),
            "type" %in% names(gff),
            "attributes" %in% names(gff),
            all.equal(nrow(df), length(f)))

  # gather genomic regions
  d.f <- data.frame(refseq_name = names(f), refseq_ID = as.numeric(gsub(".merged", "", sapply(strsplit(names(f), split = "_ID_"), "[", 2))))
  d.f <- d.f[match(df[,"ID"], d.f$refseq_ID),]
  stopifnot(all.equal(d.f$refseq_ID, df[,"ID"]))
  
  # attach ID of regions merged in pipeline step 6 (merged_step6)
  idx.merged <- grep("merged", d.f$refseq_name)
  df.merged[,"merged_step6"] <- as.character(NA)
  for (i in idx.merged) {
    
    # get ID mergedwith (7201 -> 6555)
    id <- df[i,"ID"] # df ID
    di <- df.merged[df.merged[,"ID"] %in% id, ] # df.merged line
    idm <- c(id, as.numeric(na.omit(unlist(strsplit(di[,"mergedwith"], split = ",|, "))))) # initial ids: id + mergedwith
    
    # find ID associated with gene(s) in merged region 
    genes <- unique(na.omit(unlist(strsplit(di[,"gene"], split = ",|, "))))
    ids7 <- unique(df.merged[df.merged[,"LG"] %in% di[,"LG"] & df.merged[,"gene"] %in% genes,"ID"])
    
    # find ID and mergedwith ID associated with a merged region
    idm7 <- unique(c(idm, as.numeric(na.omit(unlist(strsplit(df.merged[df.merged[,"ID"] %in% c(idm, ids7),c("mergedwith")], split = ",|, "))))))
    ids <- unique(c(ids7, idm7))
    
    # update df.merged
    df.merged[df.merged[,"ID"] %in% id, "merged_step6"] <- paste(ids, collapse = ", ") # may give redundant merged_step6
    
  }
    
  # loop over (merged) region ids and attach corresponding probe sequences and annotation
  df$probe_number <- df$probe_sequences <- df$probe_names <- df$ID_number <- df$probe_IDs <- NA
  df$ID_merged <- df$merged_step6 <- NA
  df$attribute_exon <- df$attribute_cds <- NA
  for (id in df$ID) {
    if (verbose) cat(which(df$ID == id), "/", nrow(df), "\r")
    
    # row with all ID
    df.ids <- df.merged[df.merged[,"ID"] %in% id, c("ID","mergedwith","merged_step6")]
    
    # all ID (ID, mergedwith [7201 -> 6555], merged_step6)
    ids <- suppressWarnings(as.numeric(na.omit(as.numeric(unique(unname(unlist(strsplit(apply(df.ids, 1, paste, collapse = ", "), split = ",|, "))))))))
    
    # subsets
    di <- df.probes[df.probes[,"ID"] %in% ids,] # subset of 12049 probes (length(ids) rows with all probe sequences)
    # dj <- df.merged[df.merged[,"ID"] %in% ids,] # subset of 6555 target regions (1 row for every involved ID)
    stopifnot(length(unique(di[,"LG"])) == 1)
    
    # attach probes (may be redundant)
    # df[df$ID %in% id, "probe_IDs"] <- paste(unique(ids), collapse = ", ")
    df[df$ID %in% id, "probe_IDs"] <- paste(unique(di$ID), collapse = ", ")
    df[df$ID %in% id, "ID_number"] <- length(unique(di$ID))
    df[df$ID %in% id, "probe_names"] <- paste(di[,"pname"], collapse = ", ")
    df[df$ID %in% id, "probe_sequences"] <- paste(di[,"seq"], collapse = ", ")
    df[df$ID %in% id, "probe_number"] <- nrow(di)
    
    # attach genomic region (reference sequence)
    df[df$ID %in% id, "genomic_region"] <- d.f[d.f$refseq_ID %in% id,"refseq_name"]
    df[df$ID %in% id,"length_region"] <- lengths(f[df[df$ID %in% id,"genomic_region"]])

    # attach merged_initial and merged_step6
    df[df$ID %in% id,"ID_merged"] <- sub("^$", NA, paste(ids[!ids %in% df.merged$ID], collapse = ", "))
    df[df$ID %in% id,"merged_step6"] <- as.numeric(grepl("merged$", df[df$ID %in% id,"genomic_region"]))

    # attach LG, start, end, length_conserved, gene    
    df[df$ID %in% id,"LG"] <- unique(di$LG)
    df[df$ID %in% id,"start"] <- min(di$start)
    df[df$ID %in% id,"end"] <- max(di$end)
    df[df$ID %in% id,"length_conserved"] <- sum(unlist(lapply(tapply(di$seq, INDEX = di$ID, FUN = nchar), FUN = function(x) {sum(x) - 50*(length(x)-1)})))
  
    # attach annotation
    df[df$ID %in% id,"gene"] <- paste(genes <- unique(unlist(strsplit(df.merged[df.merged[,"ID"] %in% ids,"gene"], split = ",|, "))), collapse = ", ")
    for (gene in genes) {
      ls <- gff[gff$gene %in% gene,]
      
      # attach annotation
      if (!gene %in% c("NA",NA)) {
        # attribute exon
        if (any(ls$type == "exon")) {
          le <- paste(unique(gsub("ID=[A-Za-z0-9]+;Parent=[A-Za-z0-9]+;(.*)", "\\1", unique(ls[ls$type == "exon","attributes"]))), collapse = " | ")
          if (is.na(df[df$ID %in% id, "attribute_exon"])) {
            df[df$ID %in% id, "attribute_exon"] <- le
          } else {
            df[df$ID %in% id, "attribute_exon"] <- paste(df[df$ID %in% id, "attribute_exon"], le, sep = " | ")
          }
        }
        
        # attribute CDS
        if (any(ls$type == "CDS")) {
          lc <- paste(unique(gsub("ID=[A-Za-z0-9]+;Parent=[A-Za-z0-9]+;(.*)", "\\1", unique(ls[ls$type == "CDS","attributes"]))), collapse = " | ")
          if (is.na(df[df$ID %in% id, "attribute_cds"])) {
            df[df$ID %in% id, "attribute_cds"] <- lc
          } else {
            df[df$ID %in% id, "attribute_cds"] <- paste(df[df$ID %in% id, "attribute_cds"], lc, sep = " | ")
          }
        }
      }
    }

  }
  
  # return new df (updates ID_number and probe variables) and new df.merged (mergedwith)
  return(res = list(df = df, df.merged = df.merged))
}


## Read 12049 probe sequences (used for target capture in the lab)
p <- read.FASTA("fasta/Cajanus_cajan_12049probes_6555reg.fasta")
df.merged <- df.merged.orig <- read.csv2("fasta/Cajanus_cajan_6555reg.csv")

## Read 1005 / 2396 reference sequences (used to map against)
f.subfamily <- read.FASTA("fasta/consFabaceae_4c_1005.fasta")
f.species <- read.FASTA("fasta/consDalbergia_4c_2396.fasta")

## Read Supplementary Tables (will attach probe sequences to these tables)
d.subfamily <- d.subfamily.orig <- data.frame(readxl::read_excel("../Tables/SupplementaryTables_S6-S9.xlsx", sheet = "TableS6", skip = 3))
d.species <- d.species.orig <- data.frame(readxl::read_excel("../Tables/SupplementaryTables_S6-S9.xlsx", sheet = "TableS7", skip = 3))

## Read genome annotation (takes time)
# gff <- read.gff("ncbi-genomes-2019-01-28/GCF_000340665.1_C.cajan_V1.0_genomic.gff", na.strings = c(".", "?"))
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
# save(gff, file = "ncbi-genomes-2019-01-28/gff.rda") # much smaller file
load("ncbi-genomes-2019-01-28/gff.rda") # gff


## Harmonize names
# df.probes
df.probes <- data.frame(ID = as.numeric(sub("dalb_([0-9]+)__.*", "\\1", sub("_p_[0-9]+$", "", names(p)))),
                        pname = names(p),
                        n = sub("_p_[0-9]+$", "", names(p)),
                        seq = sapply(lapply(as.character(p), FUN = toupper), function(x) paste(x, collapse = "")))

se <- sapply(strsplit(df.probes$n, split = "__"), "[", 3)
df.probes$start <- as.numeric(sapply(strsplit(se, split = "_"), "[", 1))
df.probes$end <- as.numeric(sapply(strsplit(se, split = "_"), "[", 2))
df.probes$LG <- sub("^CcLG", "", sapply(strsplit(df.probes$pname, split = "__"), "[", 2))
df.probes <- df.probes[order(df.probes$LG, df.probes$start, df.probes$pname),]

# harmonize
df.merged$n <- gsub("base.ccajan.extend.ccajan.1h.100i.100c.100bp_([0-9]+)__", "dalb_\\1__", df.merged$name)

## Attach probe sequences and gene annotations
old <- c("merged","probe_number","probe_names","probe_sequences","attribute_cds","attribute_exon")
res.subfamily <- attach.probeseqs(df = d.subfamily[,!names(d.subfamily) %in% old], df.merged = df.merged, df.probes = df.probes, f = f.subfamily, gff = gff)
res.species <- attach.probeseqs(df = d.species[,!names(d.species) %in% old], df.merged = df.merged, df.probes = df.probes, f = f.species, gff = gff)

scan()

## Write Tables S4-S5
col.order <- c("ID","genomic_region","length_region","merged_step6","ID_merged",
               "LG","start","end","length_conserved",
               "probe_IDs","ID_number","probe_number","probe_names","probe_sequences",
               "gene","attribute_cds","attribute_exon")

readr::write_tsv(res.subfamily$df[,col.order], file = "TableS6_Fabaceae1005reg.txt")
readr::write_tsv(res.species$df[,col.order], file = "TableS7_Dalbergia2396reg.txt")

## Write probe sets
probe_names.subfamily <- unique(unlist(strsplit(res.subfamily$df$probe_names, split = ", ")))
probe_names.species <- unique(unlist(strsplit(res.species$df$probe_names, split = ", ")))

stopifnot(all(probe_names.subfamily %in% names(p)),
          all(probe_names.species %in% names(p)))

length(p.subfamily <- p[probe_names.subfamily]) # 3273 (3139 not considering merging in step 0 ; 2558 not considering merging in step 7)
length(p.species <- p[probe_names.species]) # 6190 (6036 not considering merging in step 0 ; 5097 not considering merging in step 7)

write.FASTA(p.subfamily, paste0("fasta/Fabaceae1005_", length(p.subfamily), "probes_", nrow(res.subfamily$df), "reg.fasta"))
write.FASTA(p.species, paste0("fasta/Dalbergia2396_", length(p.species), "probes_", nrow(res.species$df), "reg.fasta"))

## Overlap
(m.over <- length(intersect(res.species$df$ID, res.subfamily$df$ID))) # 726 overlap in merged regions

(m.over / length(unique(res.subfamily$df$ID))) # 726/1005 (72.24%)
(m.over / length(unique(res.species$df$ID))) # 726/2396 (30.30%)

length(id.subfamily <- unique(as.numeric(unlist(strsplit(res.subfamily$df$probe_IDs, split = ", "))))) # 1612
length(id.species <-   unique(as.numeric(unlist(strsplit(  res.species$df$probe_IDs, split = ", "))))) # 3331
(p.over <- length(intersect(id.subfamily, id.species))) # 1267 overlap in probe set regions (877 not considering merging in step 7)

(p.over / length(id.subfamily)) # 1267/1612 (78.60%) # 877/1217 (72.06%) not considering merging in step 7
(p.over / length(id.species))   # 1267/3331 (38.04%) # 877/2699 (32.49%) not considering merging in step 7


## Check
# probe_IDs in same LG
all(sapply(strsplit(res.subfamily$df$probe_IDs, split = ", "), function(x) length(unique(df.merged[df.merged$ID %in% as.numeric(x),"LG"]))) == 1)
all(sapply(strsplit(  res.species$df$probe_IDs, split = ", "), function(x) length(unique(df.merged[df.merged$ID %in% as.numeric(x),"LG"]))) == 1)

