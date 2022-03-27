#########################################################
### BLASTN search against Angiosperms353 supercontigs ###
#########################################################

# simon.crameri@usys.ethz.ch, sfcrameri@gmail.com, Mar 2022

## helperfunctions
library(ape)
read.blast <- function(file, field.names, lev.ID, lev.ORTHOGROUP) {
  dd <- read.delim(file, header = F, skip = 1)
  names(dd) <- field.names
  
  # ID (Crameri et al. 2022)
  id <- grep("_ID_", as.character(dd[1,]))
  IDm <- sub("(.*)_ID_(.*)", "\\2", dd[,id])
  dd$ID <- factor(gsub("[m.]", "", IDm), levels = lev.ID)
  # dd$ID_merged <- sub("^[.]", "", gsub("[0-9]+", "", IDm))
  dd$ID_merged <- sub("[.]$", "", paste(dd$ID, sub("^[.]", "", gsub("[0-9]+", "", IDm)), sep = "."))
  
  # ORTHOGROUP (Johnson et al. 2018)
  id2 <- grep("^P0|_combined", as.character(dd[1,]))
  dd$ORTHOGROUP <- factor(gsub("|([A-Za-z]+)_([A-Za-z]+)-([0-9]+)", "\\3", gsub("(.*)-(.*)-(.*)", "\\3", dd[,id2])), levels = lev.ORTHOGROUP)
  
  return(dd)
}
filter.blast <- function(df, min.aln.len = 100, max.evalue = 1E-04) {
  df[df$alignment_length >= min.aln.len & df$evalue <= max.evalue,]
}
reciprocal.hits <- function(blast1, blast2, query_id1 = "ID", query_id2 = "ORTHOGROUP") {
  
  # column order
  query_ids <- paste(query_id1, query_id2, sep = "_")
  cols <- unique(c(c(query_id1,"ID_merged",query_id2,query_ids,"FILTER","AMBIGUOUS1","AMBIGUOUS2","AMBIGUOUS"), names(blast1)))
  
  # find reciprocal BLAST hits
  df.rec <- unique(data.frame(query_id1 = c(blast1[,query_id1], blast2[,query_id1]),
                              query_id2 = c(blast1[,query_id2], blast2[,query_id2])))
  names(df.rec) <- c(query_id1, query_id2)
  
  # bind reciprocal BLAST hits
  d.rec <- rbind(blast1[blast1[,query_id1] %in% df.rec[,query_id1] & blast1[,query_id2] %in% df.rec[,query_id2],],
                 blast2[blast2[,query_id1] %in% df.rec[,query_id1] & blast2[,query_id2] %in% df.rec[,query_id2],])
  if (identical(nrow(blast1)+nrow(blast2), nrow(d.rec))) {
    message("all hits are reciprocal!")
  } else {
    message("non-reciprocal hits detected!")
  }
  
  # remove redundancy (for every pair of query/subject, keep the hit with lowest evalue only)
  d.rec[,query_ids] <- paste(d.rec[,query_id1], d.rec[,query_id2], sep = "_")
  d.rec$FILTER <- 1
  d.rec$AMBIGUOUS1 <- d.rec$AMBIGUOUS2 <- d.rec$AMBIGUOUS <- NA
  d.rec <- d.rec[,cols] # reorder columns
  for (i in unique(d.rec[,query_ids])) {
    idx.match <- which(d.rec[,query_ids] %in% i)
    idx.beste <- which.min(d.rec[idx.match,"evalue"])
    
    idx.best.reciprocal <- which(
      d.rec[idx.match,"bit_score"] == d.rec[idx.match[idx.beste],"bit_score"] &
        d.rec[idx.match,"score"] == d.rec[idx.match[idx.beste],"score"])
    idx.best <- idx.best.reciprocal[grep("_ID_", d.rec[idx.match[idx.best.reciprocal],"query_id"])]
    
    if (length(idx.best) > 1) idx.best <- idx.best[1]
    if (length(idx.best) == 0) idx.best <- idx.beste
    d.rec[idx.match[idx.best], "FILTER"] <- 0
    amb1 <- d.rec[d.rec$FILTER == 0 & d.rec[,query_id1] %in% d.rec[idx.match[idx.best],query_id1], "AMBIGUOUS1"]
    amb2 <- d.rec[d.rec$FILTER == 0 & d.rec[,query_id2] %in% d.rec[idx.match[idx.best],query_id2], "AMBIGUOUS2"]
    if (length(amb1) > 1) {
      d.rec[d.rec$FILTER == 0 & d.rec[,query_id1] %in% d.rec[idx.match[idx.best],query_id1], "AMBIGUOUS1"] <- query_id2
    }
    if (length(amb2) > 1) {
      d.rec[d.rec$FILTER == 0 & d.rec[,query_id2] %in% d.rec[idx.match[idx.best],query_id2], "AMBIGUOUS2"] <- query_id1
    }
  }
  
  # collapse AMBIGUOUS1 and AMBIGUOUS2
  d.rec$AMBIGUOUS <- gsub("^NA, |, NA$", "", apply(d.rec[,c("AMBIGUOUS1","AMBIGUOUS2")],
                                                   MARGIN = 1,
                                                   FUN = function(x) paste(x, collapse = ", ")))
  d.rec$AMBIGUOUS[d.rec$AMBIGUOUS %in% "NA"] <- NA
  d.rec <- d.rec[,!names(d.rec) %in% c("AMBIGUOUS1","AMBIGUOUS2")]
  
  # check
  if (!identical(length(unique(d.rec[,query_ids])), sum(d.rec$FILTER == 0))) {
    message("number of unique ", query_ids, " is NOT equal to number of unfiltered hits")
  }
  
  # return results
  return(d.rec)
}


# ## shorten fasta names
# d.1005 <- read.FASTA("consFabaceae_4c_1005.fasta")
# d.2396 <- read.FASTA("consDalbergia_4c_2396.fasta")
# names(d.1005) <- gsub("merged","m", gsub("_LG_Scaffold", "_S", names(d.1005)))
# names(d.2396) <- gsub("merged","m", gsub("_LG_Scaffold", "_S", names(d.2396)))
# write.FASTA(d.1005, file = "consFabaceae_4c_1005_shortnames.fasta")
# write.FASTA(d.2396, file = "consDalbergia_4c_2396_shortnames.fasta")

# ## subset supercontigs353 (Fabaceae - 322 supercontigs)
# https://datadryad.org/stash/dataset/doi:10.5061/dryad.s3h9r6j
# # read supercontings
# p.353 <- "doi_10.5061_dryad.s3h9r6j__v1/sequences/all_supercontigs353.fasta"
# d.353 <- read.FASTA(p.353)
# 
# # subset supercontigs (keep Fabaceae supercontigs only)
# d.322 <- d.353[grep("^P03-10G", names(d.353))] # P03-10G = Gilbertiodendron ecoukense (Pellegr.) Burgt / Fabaceae, Detarioideae
# write.FASTA(d.322, file = "doi_10.5061_dryad.s3h9r6j__v1/sequences/P03-10G_supercontigs322.fasta")

## arguments
bnames <- c("query_id","query_length","subject_id","subject_length","alignment_length","perc_identity",
            "q_start","q_end","s_start","s_end","evalue","bit_score","score","subject_title")

## read annotations
# Angiosperm353 ORTHOGROUPS (Johnson et al. 2018)
dA <- data.frame(readxl::read_excel("Johnson_etal_2018_Extensive_Annotation_FixedFields_SC.xlsx")) # may have multiple lines

# collapse multiple lines (concatenate GO_slims)
d.A <- data.frame(array(NA, dim = c(0, 6), dimnames = list(NULL, c("ORTHOGROUP","Locus","Gene_Model","GO_term","annotations","details"))))
for (i in unique(dA$GO_term)) {
  tocollapse <- c("cat","code","GO_slim","GO_slim2","GO_slim3","GO_slim4","Reference","Made_by","date_last_modified")
  d <- dA[dA$GO_term %in% i,]
  d$GO_slims <- apply(d[,tocollapse], 1, FUN = function(x) {gsub("  ", " ", sub("^; ", "", paste(paste0(paste0("; ", tocollapse, ": "), x), collapse = "")))})
  d$annotations <- paste(unique(sub("^cat: ", "", sapply(strsplit(d$GO_slims, split = "; "), "[", 1))), collapse = " | ")
  d$details <- paste(d$GO_slims, collapse = " | ")
  # View(d[,!names(d) %in% tocollapse])
  d.A <- rbind(d.A, unique(d[,!names(d) %in% c(tocollapse, "GO_slims")]))
}

ddA <- data.frame(array(NA, dim = c(0, 5), dimnames = list(NULL, c("ORTHOGROUP","Locus","Gene_Model","GO_term","annotations"))))
for (i in unique(dA$ORTHOGROUP)) {
  tocollapse <- c("Gene_Model","GO_term","annotations")
  d <- d.A[d.A$ORTHOGROUP %in% i,]
  d$Gene_Models <- paste(unique(d$Gene_Model), collapse = " || ")
  d$GO_terms <- paste(unique(d$GO_term), collapse = " || ")
  d$Annotations <- paste(unique(d$annotations), collapse = "||")
  ddA <- rbind(ddA, unique(d[,!names(d) %in% c(tocollapse, "details")]))
}
lev.ORTHOGROUP <- unique(ddA$ORTHOGROUP)

# 6555 target regions (Crameri et al. 2022)
ddF <- read.csv2("orig.fragments.merged100_dropped0_6555.csv")
ID.2396 <- as.integer(gsub(".m", "", sub("(.*)_ID_(.*)", "\\2", names(read.FASTA("consDalbergia_4c_2396_shortnames.fasta")))))
ID.1005 <- as.integer(gsub(".m", "", sub("(.*)_ID_(.*)", "\\2", names(read.FASTA("consFabaceae_4c_1005_shortnames.fasta")))))

## read BLAST output
files <- list.files(path = "blast_results", pattern = "P03-10G|all_supercontigs", full.names = TRUE)
fnames <- c("d.AAD","d.AAF","d.DAA","d.DA","d.FAA","d.FA","d.AD","d.AF")
for (file in files) {
  
  # choose query levels
  lev.ID <- if (grepl("Dalbergia", file)) ID.2396 else ID.1005
  val <- read.blast(file = file, field.names = bnames,
                    lev.ID = lev.ID, lev.ORTHOGROUP = lev.ORTHOGROUP)
  
  assign(fnames[which(files == file)], value = val)
}


## filter for alignment length and score (aln.len >= 100, evalue <= 1E-04)
min.aln.len <- 100
max.evalue <- 1E-04

df.DA <- filter.blast(d.DA, min.aln.len = min.aln.len, max.evalue = max.evalue)
df.FA <- filter.blast(d.FA, min.aln.len = min.aln.len, max.evalue = max.evalue)
df.AD <- filter.blast(d.AD, min.aln.len = min.aln.len, max.evalue = max.evalue)
df.AF <- filter.blast(d.AF, min.aln.len = min.aln.len, max.evalue = max.evalue)

df.DAA <- filter.blast(d.DAA, min.aln.len = min.aln.len, max.evalue = max.evalue)
df.FAA <- filter.blast(d.FAA, min.aln.len = min.aln.len, max.evalue = max.evalue)
df.AAD <- filter.blast(d.AAD, min.aln.len = min.aln.len, max.evalue = max.evalue)
df.AAF <- filter.blast(d.AAF, min.aln.len = min.aln.len, max.evalue = max.evalue)

## Filter for reciprocal BLASTN hits
# these contain all reciprocal hits (all Angiosperm343 taxa with hits)
dd.subfamily <- reciprocal.hits(blast1 = df.FAA, blast2 = df.AAF)
dd.species <- reciprocal.hits(blast1 = df.DAA, blast2 = df.AAD)

# these are reduced to a single line per unique ID_ORTHOGROUP hit (only Angiosperm343 taxon with best hit)
d.subfamily <- dd.subfamily[dd.subfamily$FILTER == 0,]
d.species <- dd.species[dd.species$FILTER == 0,]

## Merge
# with gene annotations
dm.subfamily <- merge(merge(d.subfamily, ddF[,c("ID","LG","start","end","gene")], by = "ID", all.x = TRUE, all.y = FALSE, sort = FALSE), ddA, by = "ORTHOGROUP", all.x = TRUE, all.y = FALSE, sort = FALSE)
dm.species <- merge(merge(d.species, ddF[,c("ID","LG","start","end","gene")], by = "ID", all.x = TRUE, all.y = FALSE, sort = FALSE), ddA, by = "ORTHOGROUP", all.x = TRUE, all.y = FALSE, sort = FALSE)

# with taxon names
spec353 <- data.frame(readxl::read_excel("Voucher_table.xlsx"))[,1:5]
names(spec353)[1] <- "Accession"

dm.subfamily$Accession <- gsub("-[0-9]+$", "", dm.subfamily$subject_id)
dm.subfamily <- merge(dm.subfamily, spec353, by = "Accession", all.x = TRUE, all.y = FALSE, sort = FALSE)

dm.species$Accession <- gsub("-[0-9]+$", "", dm.species$subject_id)
dm.species <- merge(dm.species, spec353, by = "Accession", all.x = TRUE, all.y = FALSE, sort = FALSE)

## Sort columns and rows
col.order <- c("ID","ORTHOGROUP","AMBIGUOUS",
  "query_length","subject_length","alignment_length","perc_identity","q_start","q_end","s_start","s_end","evalue","bit_score","score",
  "query_id","LG","start","end","gene",
  "subject_id","Accession","APG.IV.Order","APG.IV.Family","Species","Author",
  "Locus","Gene_Models","GO_terms","Annotations")
dm.subfamily <- dm.subfamily[order(dm.subfamily$LG, dm.subfamily$start, as.numeric(as.character(dm.subfamily$ORTHOGROUP))),col.order]
dm.species <- dm.species[order(dm.species$LG, dm.species$start, as.numeric(as.character(dm.species$ORTHOGROUP))),col.order]

## Write output
readr::write_tsv(dm.species, file = "blast_results/BLASTN_Dalbergia2396_vs_Angiosperm353.txt")
readr::write_tsv(dm.subfamily, file = "blast_results/BLASTN_Fabaceae1005_vs_Angiosperm353.txt")

