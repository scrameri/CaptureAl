#!/cluster/apps/r/3.6.0_openblas/x86_64/bin/Rscript

##############################
## Get BLAST stats per grop ##
##############################

## Usage
# run script without arguments to see the arguments it takes

## Author
# simon.crameri@env.ethz.ch, Mar 2019

## Load required library
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(ggplot2))

## Get arguments
args = commandArgs(trailingOnly = TRUE)

## Check arguments
options(warning.length = 2000L)
if (! length(args) %in% c(6)) {
  stop("6 arguments needed:
       REQUIRED
       1) <tab|CHR>:                path to blast results file (output of <filter-blast-files-parallel.sh>) ;
       2) <meta|CHR>:               path to meta file that maps individuals (1st column) to groups (2nd column). File is expected to have a header. Individuals mapping to <NA> will not be considered for overlap calculation but will be displayed in plots. ;
       3) <ref>|CHR>:               path to reference .fasta file (used by <blast.contigs.parallel.sh>) ;
       4) <min.prop.hit|NUM>:       minimum proportion of indivivuals in GROUP having one or more BLAST hit(s) at a retained locus ;
       5) <max.mean.mult.hits|NUM>: maximum mean number of non-zero BLAST hits in GROUP at a retained locus  ;
       6) <min.prop.passed|NUM>:    minimum proportion of individuals in GROUP having one or more BLAST hits that meet the requirements as stated in blastres ;
       ", call.=FALSE)
}

## Set parameters
cat("\n=== BLAST Stats ===\n")
tab <- as.character(args[1])
meta <- as.character(args[2])
ref <- as.character(args[3])
min.prop.hit = as.numeric(args[4])
max.mean.mult.hits = as.numeric(args[5])
min.prop.passed = as.numeric(args[6])

## Set paramters (debugging)
# tab <- "Chapter1.1_assembleblast_6555_9/blastres.1e-04.0.80.100.5000.80.txt"
# meta <- "mapfile.orig.subfamily.txt"
# ref <- "orig.fragments.merged100_6555.fasta" #fabaceae_cons9_796.fasta"
# min.prop.hit = 0.5
# max.mean.mult.hits = 1.5
# min.prop.passed = 0.5

## Additional parameters
plot.all = FALSE # if TRUE, plots 26 plots, if FALSE, plots the 11 most important plots

parms <- unlist(strsplit(basename(tab), split = "[.]"))[-c(1,8)]
e=as.numeric(parms[1]) # evalue threshold (selects the corresponding .blast files)
p=as.numeric(parms[2]) # minimum percentage of query.length in BLAST alignment length
a=as.numeric(parms[3]) # minimum BLAST alignment length
n=as.numeric(parms[4]) # minimum contig length
x=as.numeric(parms[5]) # maximum contig length
i=as.numeric(parms[6]) # minimum percent identity of BLAST alignment

## Print parameters
cat("\nblast results used:", tab, "\n")
cat("   interpreted filters for passed hits: e =", e, "; p =", p, "; a =", a, "; n =", n, "; x =", x, "; i =", i, "\n")
cat("\nreference sequences used:", ref, "\n")
cat("\nmetadata file used:", meta, "\n")
cat("\noverlap parameters: min.prop.hit =", min.prop.hit, "; max.mean.mult.hits =", max.mean.mult.hits, "; min.prop.passed =", min.prop.passed, "\n")

thr.xormore <- 8
# cat("\nmultiple hits will be summarized up to", thr.xormore, "hits\n")

## Check parameters
stopifnot(file.exists(tab), file.exists(meta), file.exists(ref), 
          min.prop.hit >= 0, min.prop.hit <= 1,
          max.mean.mult.hits >= 1,
          min.prop.passed >= 0, min.prop.passed <= 1,
          is.numeric(thr.xormore), thr.xormore > 0)

## Define helperfunctions
# colors
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# aggregate variables
do.aggregate <- function(df, byvars, n) {
  
  ## Check df
  stopifnot(all(c("ID", "LABEL", "query.id", "nhits", "prop.passed", "prop.best", "passed") %in% names(df)))
  
  ## Get by list
  by = list()
  for (i in seq(length(byvars))) {
    by[[i]] <- df[,byvars[i]]
  }
  
  ## Aggregate
  # number of non-zero hits
  dhit     <- aggregate(df$nhits,       by = by, FUN = function(x) {length(x[x>0])}, drop = FALSE)
  
  # number of single-hits
  dsingle  <- aggregate(df$nhits,       by = by, FUN = function(x) {length(x[x==1])}, drop = FALSE)
  
  # number of non-zero multiple hits
  dmult    <- aggregate(df$nhits,       by = by, FUN = mean, drop = FALSE)
  
  # proportion of hits passed
  dpass    <- aggregate(df$prop.passed, by = by, FUN = mean, drop = FALSE)
  
  # proportion of best hits
  dbest    <- aggregate(df$prop.best,   by = by, FUN = mean, drop = FALSE)
  
  ## Get identifyers
  d.hits.ids <- data.frame(array(data = NA, dim = c(nrow(dhit), 0)))
  for (i in seq(length(byvars))) {
    d.hits.ids[,byvars[i]] <- factor(dhit[,paste0("Group.", i)])
  }
  
  ## Get numbers for proportions
  if (is.null(n)) {
    if (! "Group.2" %in% names(dhit)) {
      stop("n is NULL, attempt to set n from second element in by, but only 1 by element found!")
    }
    n <- as.numeric(gsub(")$", "", sapply(strsplit(as.character(dhit$Group.2), split = " \\(n ="), "[", 2)))
    if (anyNA(n)) stop("n is NULL but could not be grepped from second element in by!")
  }
  
  ## Return results
  d.hits <- data.frame(d.hits.ids,
                       n = n,
                       nb.hit =    dhit$x,    p.hits =  dhit$x / n,
                       nb.1hit =   dsingle$x, p.1hits = dsingle$x / n,
                       nb.hits =   dmult$x,
                       p.passed =  dpass$x,
                       p.best =    dbest$x)
  return(d.hits)
}

# aggregate variables (all and passed), and merge
get.aggr.data <- function(df, byvars, n) {
  
  ## Check df
  stopifnot(all(c("ID", "LABEL", "query.id", "nhits", "prop.passed", "prop.best", "passed") %in% names(df)))
  
  ## Get passed hits
  df.passed <- subset(df, passed > 0)
  
  ## Aggregate over byvars    
  d.hits.a <- do.aggregate(df = df,        byvars = byvars, n = n)
  d.hits.p <- do.aggregate(df = df.passed, byvars = byvars, n = n)
  
  d.hits <- merge(d.hits.a, d.hits.p, by = c(byvars,"n"), all = T, sort = F)
  names(d.hits) <- gsub(".y$", ".p", gsub(".x$", "", names(d.hits)))
  d.hits$p.hit.p <- d.hits$nb.hit.p / d.hits$nb.hit
  d.hits$p.1hit.p <- d.hits$nb.1hit.p / d.hits$nb.1hit
  #d.hits <- d.hits[,c(byvars,"n", sort(names(d.hits)[-(1:length(byvars)+1)]))]
  
  # replace NA values with zero
  d.hits[is.na(d.hits)] <- 0
  
  ## merge everything
  d.merge <- merge(d.hits, 
                   df[-which(duplicated(df[,byvars[length(byvars)]])), sort(unique(c(byvars[length(byvars)], "GROUP", "LABEL")))], 
                   by = byvars[length(byvars)])
  d.merge <- d.merge[,c("GROUP", "LABEL", names(d.hits)[which(!names(d.hits) %in% c("GROUP", "LABEL"))])]
  d.merge <- d.merge[order(d.merge$GROUP, d.merge[,byvars[1]]),]
  rownames(d.merge) <- seq(nrow(d.merge))
  return(d.merge)
}


#####################################################################################################################

## Read reference
dref <- read.FASTA(ref)
refnames <- names(dref)

## Read blast data
dd <- read.delim(tab)
dbcols <- c("ID","id","query.id","query.length","subject.id","subject.length","alignment.length",
            "perc.identity","q.start","q.end","s.start","s.end","evalue","bit.score","score",
            "subject.title","best","passed")
stopifnot(all.equal(names(dd), dbcols))

## Remove duplicated hits (keep better hit (evalue) if there are two contigs hitting the same refseq)
dd.orig <- dd[order(dd$id, dd$subject.id, -dd$best),]
rownames(dd.orig) <- seq(nrow(dd.orig))
dd <- dd.orig[!duplicated(dd.orig[,c("id","subject.id")]),]
stopifnot(all.equal(sum(dd$best), sum(dd.orig$best)))

## Add missing refnames levels
dd$query.id <- factor(dd$query.id, levels = refnames)
stopifnot(!anyNA(dd$query.id))

## Read metadata
dmeta <- read.delim(meta, header = T)
names(dmeta)[1] <- "ID"
names(dmeta)[2] <- "GROUP"
dmeta$GROUP <- as.character(dmeta$GROUP)
dmeta[is.na(dmeta$GROUP),"GROUP"] <- "NA"
dmeta$GROUP <- factor(dmeta$GROUP)

## Merge blast and meta data
stopifnot(all(dd$ID %in% dmeta$ID))
dd <- merge(dd, dmeta, by = "ID", all.x = T, all.y = F, sort = F)
dd$GROUP <- as.character(dd$GROUP)
dd[is.na(dd$GROUP), "GROUP"] <- "NA"

cat(paste0("\nfound ", length(unique(dd$GROUP)), " groups in ", meta, " matching individuals in ", tab, ":\n"))
grtab <- table(dd[-which(duplicated(dd$ID)),"GROUP"])
print(grtab)

## Add number of individuals per group
dd$LABEL <- dd$GROUP
for (group in unique(dd$GROUP)) {
  dd$LABEL[dd$LABEL == group] <- paste0(names(grtab[group]), " (n = ", grtab[group], ")")
}

## Get number of BLAST hits per locus
# get locus tables (best = contigs that passed and showed low evalue / long length ; passed = contigs that met epanxi thresholds)
dtab <- table(dd$id, dd$best)
dtabpassed <- table(dd$id, dd$passed)
stopifnot(all.equal(names(dtab[,1]), names(dtabpassed[,1])))

# compile number of hits
dp <- data.frame(id = names(dtab[,1]), 
                 worse = as.numeric(dtab[,"0"]), best = as.numeric(dtab[,"1"]), 
                 filtered = as.numeric(dtabpassed[,"0"]), passed = as.numeric(dtabpassed[,"1"]))
dp <- merge(dp, dd[-which(duplicated(dd$id)),c("id", "ID", "LABEL", "GROUP", "query.id", "query.length")], by = "id", all.x = T, all.y = F, sort = F)
stopifnot(all.equal(dp$worse + dp$best, dp$filtered + dp$passed))
dp$nhits <- dp$worse + dp$best

# group by number-of-hits-categories
dp$nhits2 <- factor(as.character(dp$nhits), levels = c(1:max(dp$nhits), paste0(">=", thr.xormore)))
dp$nhits2[dp$nhits >= thr.xormore] <- paste0(">=", thr.xormore)
dp$nhits2 <- droplevels(dp$nhits2)

# proportion of hits passed / used for alignment
dp$prop.passed <- dp$passed / dp$nhits
dp$prop.best <- dp$best / dp$nhits

## Stats per group
# Aggregate by individual (witin GROUP)
d.hits.IND <- get.aggr.data(df = dp, byvars = c("ID"), n = nlevels(dp[,"query.id"]))

# Aggregate by target locus (witin GROUP)
d.hits.LOC <- get.aggr.data(df = dp, byvars = c("query.id", "LABEL"), n = NULL)

# update thr.xormore (only used if categories are used)
if (thr.xormore >= max(d.hits.LOC$nb.hits, na.rm = T)-1) {
  thr.xormore <- ceiling(max(d.hits.LOC$nb.hits, na.rm = T)) #ceiling(max(d.hits.LOC$nb.hits))-1
}
cat("\nmultiple hits will be summarized up to", thr.xormore, "hits\n")

# number of BLAST hits by category
nb.hit.cat.max <- c(0:thr.xormore)
nb.hit.cat <- gsub("-1<x", "x", paste0(nb.hit.cat.max-1, "<x<=", nb.hit.cat.max))
d.hits.LOC$nb.hits.cat <- d.hits.LOC$nb.hits.p.cat <- factor(NA, levels = nb.hit.cat)
for (nb in seq(length(nb.hit.cat))) {
  if (nb == 1) lower.limit <- -0.001 else lower.limit <- nb.hit.cat.max[nb] - 1
  if (nb == length(nb.hit.cat)) upper.limit <- nb.hit.cat.max[nb] + 0.001 else upper.limit <- nb.hit.cat.max[nb]
  d.hits.LOC$nb.hits.cat[d.hits.LOC$nb.hits > lower.limit & d.hits.LOC$nb.hits <= upper.limit] <- levels(d.hits.LOC$nb.hits.cat)[nb]
  d.hits.LOC$nb.hits.p.cat[d.hits.LOC$nb.hits.p > lower.limit & d.hits.LOC$nb.hits.p <= upper.limit] <- levels(d.hits.LOC$nb.hits.p.cat)[nb]
}

## Filter loci per group (loci meeting min.prop.hit, max.mean.mult.hits, min.prop.passed)
d.list <- sapply(unique(d.hits.LOC$GROUP), FUN = function(x) {
  subset(d.hits.LOC, GROUP == x & p.hits >= min.prop.hit & nb.hits <= max.mean.mult.hits & p.hit.p >= min.prop.passed) # nb.hits or nb.hits.p
}, simplify = F)

## Get locus overlap
cat("\ncomputing overlap in all non-NA groups:\n")
if (any(names(grtab) == "NA")) nonnagroups <- names(grtab)[-which(names(grtab) == "NA")] else nonnagroups <- names(grtab)
grtab.nonna <- grtab[nonnagroups]
print(grtab.nonna)

d.list.nonna <- d.list[which(names(d.list) %in% nonnagroups)]
overlap <- d.list.nonna[[nonnagroups[1]]]$query.id
for (group in nonnagroups[-1]) {
  overlap <- intersect(overlap, d.list.nonna[[group]]$query.id)
}  
cat("\nfound", length(overlap), "overlapping target loci in all non-NA groups\n")

##############################################################################################

## Plot data
cat("\nplotting stats...\n")

pdf(file = gsub(".txt$", paste0("-", round(min.prop.hit,2),"-",round(max.mean.mult.hits,2),"-",round(min.prop.passed,2),".pdf"), tab))

## Plot VENN diagram
# rename with number of loci passing filters in each groups
d.list.venn <- lapply(d.list.nonna, FUN = function(x) unique(x$query.id))
names(d.list.venn) <- paste0(names(d.list.venn), " (n = ", lapply(d.list.venn, length), ")")

if (length(d.list.venn) <= 5) {
  
  # set plot margins and color function
  oldmar <- par()$mar
  par(mar = rep(0,4))
  
  # set colors
  venn.cols <- gg_color_hue(length(grtab))[which(names(grtab) %in% nonnagroups)]
  
  # plot VENN diagram
  suppressPackageStartupMessages(library(grid))
  suppressPackageStartupMessages(library(VennDiagram))
  plot.new()
  futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
  grid.draw(venn.diagram(d.list.venn, filename = NULL, 
                         height = 3000, width = 3000,
                         cat.dist = rep(0, length(venn.cols)),
                         # cat.pos = rep(0, length(venn.cols)),
                         cat.col = venn.cols,
                         col = venn.cols,
                         fill = venn.cols,
                         alpha = 0.6,
                         cex = 2, cat.cex = 0,
                         # direct.area = T,
                         # area.vector = rep(1, 10),
                         # euler.d = F,
                         lty = rep(2, length(venn.cols)),
                         # print.mode = c("raw","percent"), sigdigs = 3))
                         print.mode = c("raw")))
  legend("topleft", names(d.list.venn), text.col = venn.cols, bty = "n", cex = 1)
  par(mar = oldmar)
} else {
  cat("\nmore than 5 non-NA groups present: no VENN diagram shown.\n")
}

## Plot contig and alignment stats
# subject length (contig length)
p1 <- ggplot(dd, aes(LABEL, subject.length, fill = LABEL)) +
  geom_boxplot(alpha = 0.5) +
  geom_violin(alpha = 0.5) +
  geom_hline(aes(yintercept = n), col = "tomato") +
  geom_hline(aes(yintercept = x), col = "tomato") +
  coord_flip() +
  xlab("") +
  ylab("Contig length (bp)") +
  ggtitle("Contigs") +
  scale_fill_discrete(guide = F) +
  scale_y_continuous(limits = c(1, max(dd$subject.length)), trans = "log10") +
  theme_bw()

p2 <- ggplot(subset(dd, passed == 1), aes(LABEL, subject.length, fill = LABEL)) +
  geom_boxplot(alpha = 0.5) +
  geom_violin(alpha = 0.5) +
  geom_hline(aes(yintercept = n), col = "tomato") +
  geom_hline(aes(yintercept = x), col = "tomato") +
  coord_flip() +
  xlab("") +
  ylab("Contig length (bp)") +
  ggtitle("Contigs passed") +
  scale_fill_discrete(guide = F) +
  scale_y_continuous(limits = c(1, x), trans = "identity") +
  theme_bw()

p3 <- ggplot(subset(dd, query.id %in% overlap & passed == 1), aes(LABEL, subject.length, fill = LABEL)) +
  geom_boxplot(alpha = 0.5) +
  geom_violin(alpha = 0.5) +
  geom_hline(aes(yintercept = n), col = "tomato") +
  geom_hline(aes(yintercept = x), col = "tomato") +
  coord_flip() +
  xlab("") +
  ylab("Contig length (bp)") +
  ggtitle("Contigs passed from loci passed") +
  scale_fill_discrete(guide = F) +
  scale_y_continuous(limits = c(1, x), trans = "identity") +
  theme_bw()

p4 <- ggplot(subset(dd, query.id %in% overlap & best == 1), aes(LABEL, subject.length, fill = LABEL)) +
  geom_boxplot(alpha = 0.5) +
  geom_violin(alpha = 0.5) +
  geom_hline(aes(yintercept = n), col = "tomato") +
  geom_hline(aes(yintercept = x), col = "tomato") +
  coord_flip() +
  xlab("") +
  ylab("Contig length (bp)") +
  ggtitle("Best contigs") +
  scale_fill_discrete(guide = F) +
  scale_y_continuous(limits = c(1, x), trans = "identity") +
  theme_bw()

# alignment length 
p5 <- ggplot(dd, aes(LABEL, alignment.length, fill = LABEL)) +
  geom_boxplot(alpha = 0.5) +
  geom_violin(alpha = 0.5) +
  geom_hline(aes(yintercept = a), col = "tomato") +
  coord_flip() +
  xlab("") +
  ylab("Alignment length (bp)") +
  ggtitle("Contigs") +
  scale_fill_discrete(guide = F) +
  scale_y_continuous(limits = c(1, max(dd$alignment.length)), trans = "identity") +
  theme_bw()

p6 <- ggplot(subset(dd, passed == 1), aes(LABEL, alignment.length, fill = LABEL)) +
  geom_boxplot(alpha = 0.5) +
  geom_violin(alpha = 0.5) +
  geom_hline(aes(yintercept = a), col = "tomato") +
  coord_flip() +
  xlab("") +
  ylab("Alignment length (bp)") +
  ggtitle("Contigs passed") +
  scale_fill_discrete(guide = F) +
  scale_y_continuous(limits = c(1, max(dd$alignment.length)), trans = "identity") +
  theme_bw()

p7 <- ggplot(subset(dd, query.id %in% overlap & passed == 1), aes(LABEL, alignment.length, fill = LABEL)) +
  geom_boxplot(alpha = 0.5) +
  geom_violin(alpha = 0.5) +
  geom_hline(aes(yintercept = a), col = "tomato") +
  coord_flip() +
  xlab("") +
  ylab("Alignment length (bp)") +
  ggtitle("Contigs passed from loci passed") +
  scale_fill_discrete(guide = F) +
  scale_y_continuous(limits = c(1, max(dd$alignment.length)), trans = "identity") +
  theme_bw()

p8 <- ggplot(subset(dd, query.id %in% overlap & best == 1), aes(LABEL, alignment.length, fill = LABEL)) +
  geom_boxplot(alpha = 0.5) +
  geom_violin(alpha = 0.5) +
  geom_hline(aes(yintercept = a), col = "tomato") +
  coord_flip() +
  xlab("") +
  ylab("Alignment length (bp)") +
  ggtitle("Best contigs") +
  scale_fill_discrete(guide = F) +
  scale_y_continuous(limits = c(1, max(dd$alignment.length)), trans = "identity") +
  theme_bw()

# percent identity
p9 <- ggplot(dd, aes(LABEL, perc.identity, fill = LABEL)) +
  geom_boxplot(alpha = 0.5) +
  geom_violin(alpha = 0.5) +
  geom_hline(aes(yintercept = i), col = "tomato") +
  coord_flip() +
  xlab("") +
  ylab("Percent identity") +
  ggtitle("Contigs") +
  scale_fill_discrete(guide = F) +
  scale_y_continuous(limits = c(min(dd$perc.identity), 100), trans = "identity") +
  theme_bw()

p10 <- ggplot(subset(dd, passed == 1), aes(LABEL, perc.identity, fill = LABEL)) +
  geom_boxplot(alpha = 0.5) +
  geom_violin(alpha = 0.5) +
  geom_hline(aes(yintercept = i), col = "tomato") +
  coord_flip() +
  xlab("") +
  ylab("Percent identity") +
  ggtitle("Contigs passed") +
  scale_fill_discrete(guide = F) +
  scale_y_continuous(limits = c(min(dd$perc.identity), 100), trans = "identity") +
  theme_bw()

p11 <- ggplot(subset(dd, query.id %in% overlap & passed == 1), aes(LABEL, perc.identity, fill = LABEL)) +
  geom_boxplot(alpha = 0.5) +
  geom_violin(alpha = 0.5) +
  geom_hline(aes(yintercept = i), col = "tomato") +
  coord_flip() +
  xlab("") +
  ylab("Percent identity") +
  ggtitle("Contigs passed from loci passed") +
  scale_fill_discrete(guide = F) +
  scale_y_continuous(limits = c(min(dd$perc.identity), 100), trans = "identity") +
  theme_bw()

p12 <- ggplot(subset(dd, query.id %in% overlap & best == 1), aes(LABEL, perc.identity, fill = LABEL)) +
  geom_boxplot(alpha = 0.5) +
  geom_violin(alpha = 0.5) +
  geom_hline(aes(yintercept = i), col = "tomato") +
  coord_flip() +
  xlab("") +
  ylab("Percent identity") +
  ggtitle("Best contigs") +
  scale_fill_discrete(guide = F) +
  scale_y_continuous(limits = c(min(dd$perc.identity), 100), trans = "identity") +
  theme_bw()

# proportion of target loci with >=1 BLAST hit(s)
p13 <- ggplot(d.hits.IND, aes(LABEL, p.hits, fill = LABEL)) + 
  geom_boxplot(alpha = 0.5) + 
  geom_violin(alpha = 0.5) +
  geom_hline(aes(yintercept = min.prop.hit), col = "black") +
  coord_flip() +
  xlab("") +
  ylab("Proportion of target loci with >= 1 BLAST hit(s)") +
  scale_y_continuous(limits = c(0, 1), trans = "identity") +
  scale_fill_discrete(guide = F) +
  ggtitle(paste0("Per individual within group ; Number of target loci: ", length(refnames))) +
  theme_bw()

p14 <- ggplot(d.hits.IND, aes(x = LABEL, y = p.hit.p, fill = LABEL)) +
  geom_boxplot(alpha = 0.5) + 
  geom_violin(alpha = 0.5) +
  geom_hline(aes(yintercept = min.prop.passed), col = "black") +
  coord_flip() +
  xlab("") +
  ylab("Proportion of target loci with >= 1 BLAST hit(s) passed") +
  scale_y_continuous(limits = c(0, 1)) +
  scale_fill_discrete(guide = F) +
  ggtitle(paste0("Per individual within group ; Number of target loci: ", length(refnames))) +
  theme_bw()

p15 <- ggplot(d.hits.LOC, aes(LABEL, p.hits, fill = LABEL)) + 
  geom_boxplot(alpha = 0.5) + 
  geom_violin(alpha = 0.5) +
  geom_hline(aes(yintercept = min.prop.hit), col = "tomato") +
  coord_flip() +
  xlab("") +
  ylab("Proportion of individuals with >= 1 BLAST hit(s)") +
  scale_y_continuous(limits = c(0, 1), trans = "identity") +
  scale_fill_discrete(guide = F) +
  ggtitle(paste0("Per target locus within group ; Number of target loci: ", length(refnames))) +
  theme_bw()

p16 <- ggplot(d.hits.LOC, aes(x = LABEL, y = p.hit.p, fill = LABEL)) +
  geom_boxplot(alpha = 0.5) + 
  geom_violin(alpha = 0.5) +
  geom_hline(aes(yintercept = min.prop.passed), col = "black") +
  coord_flip() +
  xlab("") +
  ylab("Proportion of individuals with >= 1 BLAST hit(s) passed") +
  scale_y_continuous(limits = c(0, 1)) +
  scale_fill_discrete(guide = F) +
  ggtitle(paste0("Per target locus within group ; Number of target loci: ", length(refnames))) +
  theme_bw()

# number of (non-zero) BLAST hits
p17 <- ggplot(d.hits.IND, aes(x = LABEL, y = nb.hits, fill = LABEL)) +
  geom_boxplot(alpha = 0.5) + 
  geom_violin(alpha = 0.5) +
  geom_hline(aes(yintercept = max.mean.mult.hits), col = "black") +
  coord_flip() +
  xlab("") +
  ylab("Mean number of BLAST hits") +
  scale_y_continuous(limits = c(0, NA)) +
  scale_fill_discrete(guide = F) +
  ggtitle(paste0("Per individual within group ; Number of target loci: ", length(refnames))) +
  theme_bw()

p18 <- ggplot(d.hits.IND, aes(x = LABEL, y = nb.hits.p, fill = LABEL)) +
  geom_boxplot(alpha = 0.5) + 
  geom_violin(alpha = 0.5) +
  geom_hline(aes(yintercept = max.mean.mult.hits), col = "black") +
  coord_flip() +
  xlab("") +
  ylab("Mean number of passed BLAST hits") +
  scale_y_continuous(limits = c(0, NA)) +
  scale_fill_discrete(guide = F) +
  ggtitle(paste0("Per individual within group ; Number of target loci: ", length(refnames))) +
  theme_bw()

p19 <- ggplot(d.hits.LOC, aes(x = LABEL, y = nb.hits, fill = LABEL)) +
  geom_boxplot(alpha = 0.5) + 
  geom_violin(alpha = 0.5) +
  geom_hline(aes(yintercept = max.mean.mult.hits), col = "tomato") +
  coord_flip() +
  xlab("") +
  ylab("Mean number of BLAST hits") +
  scale_y_continuous(limits = c(0, NA)) +
  scale_fill_discrete(guide = F) +
  ggtitle(paste0("Per target locus within group ; Number of target loci: ", length(refnames))) +
  theme_bw()

p20 <- ggplot(d.hits.LOC, aes(x = LABEL, y = nb.hits.p, fill = LABEL)) +
  geom_boxplot(alpha = 0.5) + 
  geom_violin(alpha = 0.5) +
  geom_hline(aes(yintercept = max.mean.mult.hits), col = "black") +
  coord_flip() +
  xlab("") +
  ylab("Mean number of passed BLAST hits") +
  scale_y_continuous(limits = c(0, NA)) +
  scale_fill_discrete(guide = F) +
  ggtitle(paste0("Per target locus within group ; Number of target loci: ", length(refnames))) +
  theme_bw()

# proportion of single BLAST hit
p21 <- ggplot(d.hits.IND, aes(x = LABEL, y = p.1hits, fill = LABEL)) +
  geom_boxplot(alpha = 0.5) + 
  geom_violin(alpha = 0.5) +
  # geom_hline(aes(yintercept = X), col = "black") +
  coord_flip() +
  xlab("") +
  ylab("Proportion of target loci with single BLAST hit") +
  scale_y_continuous(limits = c(0, 1)) +
  scale_fill_discrete(guide = F) +
  ggtitle(paste0("Per individual within group ; Number of target loci: ", length(refnames))) +
  theme_bw()

p22 <- ggplot(d.hits.IND, aes(x = LABEL, y = p.1hit.p, fill = LABEL)) +
  geom_boxplot(alpha = 0.5) + 
  geom_violin(alpha = 0.5) +
  # geom_hline(aes(yintercept = X), col = "black") +
  coord_flip() +
  xlab("") +
  ylab("Proportion of target loci with single BLAST hit passed") +
  scale_y_continuous(limits = c(0, 1)) +
  scale_fill_discrete(guide = F) +
  ggtitle(paste0("Per individual within group ; Number of target loci: ", length(refnames))) +
  theme_bw()

p23 <- ggplot(d.hits.LOC, aes(x = LABEL, y = p.1hits, fill = LABEL)) +
  geom_boxplot(alpha = 0.5) + 
  geom_violin(alpha = 0.5) +
  # geom_hline(aes(yintercept = X), col = "black") +
  coord_flip() +
  xlab("") +
  ylab("Proportion of individuals with single BLAST hit") +
  scale_y_continuous(limits = c(0, 1)) +
  scale_fill_discrete(guide = F) +
  ggtitle(paste0("Per target locus within group ; Number of target loci: ", length(refnames))) +
  theme_bw()

p24 <- ggplot(d.hits.LOC, aes(x = LABEL, y = p.1hit.p, fill = LABEL)) +
  geom_boxplot(alpha = 0.5) + 
  geom_violin(alpha = 0.5) +
  # geom_hline(aes(yintercept = X), col = "black") +
  coord_flip() +
  xlab("") +
  ylab("Proportion of individuals with single BLAST hit passed") +
  scale_y_continuous(limits = c(0, 1)) +
  scale_fill_discrete(guide = F) +
  ggtitle(paste0("Per target locus within group ; Number of target loci: ", length(refnames))) +
  theme_bw()

# per number-of-hit category
p25 <- ggplot(d.hits.LOC, aes(x = LABEL, fill = nb.hits.cat)) +
  geom_bar(aes(y = ..count..), position = "dodge") +
  labs(x = "", y = "Number of target loci", fill = "BLAST hits") +
  geom_hline(aes(yintercept = length(refnames))) +
  coord_flip() +
  scale_fill_discrete(drop = FALSE) +
  ggtitle(paste0("All hits ; Number of target loci: ", length(refnames))) +
  theme_bw()
  
p26 <- ggplot(d.hits.LOC, aes(x = LABEL, fill = nb.hits.p.cat)) +
  geom_bar(aes(y = ..count..), position = "dodge") +
  labs(x = "", y = "Number of target loci", fill = "BLAST hits") +
  geom_hline(aes(yintercept = length(refnames))) +
  coord_flip() +
  scale_fill_discrete(drop = FALSE) +
  ggtitle(paste0("Passed hits ; Number of target loci: ", length(refnames))) +
  theme_bw()

## Write output
cat("\nwriting output...\n\n")

# overlapping loci
write.table(sort(overlap), file = gsub(".txt$", paste0("-", round(min.prop.hit,2),"-",round(max.mean.mult.hits,2),"-",round(min.prop.passed,2),"_overlap.txt"), tab),
            quote = F, row.names = F, col.names = F)

# plots
print(p1) #
print(p2) #
if (plot.all) print(p3)
if (plot.all) print(p4)

print(p5) #
print(p6) #
if (plot.all) print(p7)
if (plot.all) print(p8)

print(p9)  #
print(p10) #
if (plot.all) print(p11)
if (plot.all) print(p12)

if (plot.all) print(p13)
if (plot.all) print(p14)
print(p15) #
if (plot.all) print(p16)

if (plot.all) print(p17)
if (plot.all) print(p18)
print(p19) #
if (plot.all) print(p20)

print(p21) #
if (plot.all) print(p22)
print(p23) #
if (plot.all) print(p24)

print(p25) #
if (plot.all) print(p26)

dev.off()

# save output
save(dd, dp, d.hits.IND, d.hits.LOC, d.list, d.list.venn, 
     file = gsub(".txt$", paste0("-", round(min.prop.hit,2),"-",round(max.mean.mult.hits,2),"-",round(min.prop.passed,2),".rda"), tab))

cat("done!\n")


