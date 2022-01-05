#!/cluster/apps/r/3.6.0_openblas/x86_64/bin/Rscript

## Usage: plot.fastqc.R <file with paths> <pdf.name="fastqc.pdf"> <pdf.height=12> <pdf.width=12>

## Needs: fastqcr and ggplot2 packages

## Value: PDF with facetted FastQC plots for many samples

## Author: simon.crameri@usys.ethz.ch, Dec 2021

## Get arguments
args = commandArgs(trailingOnly = TRUE)

## Check arguments
options(warning.length = 2000L)
if (length(args) < 1 | length(args) > 4) {
  stop("1 argument needed (3 more are optional):
       REQUIRED
       1) <file|CHR>:       path to file with paths to fastqc.zip results 
       
       OPTIONAL
       2) <pdf.name|CHR>    output PDF file name (default: fastqc.pdf)
       3) <pdf.height|NUM>  height of output PDF file
       4) <pdf.width|NUM>   width of output PDF file
       ",
       call.=FALSE)
}

## Set arguments
# required
p <- as.character(args[1])

# optional
if (is.na(args[2])) pdf.name <- "fastqc.pdf" else pdf.name <- args[2]
if (is.na(args[3])) pdf.height <- 12 else pdf.height <- as.numeric(args[3])
if (is.na(args[4])) pdf.width <- 12 else pdf.width <- as.numeric(args[4])

# additional
col.sample <- "tomato"
col.test <- c("green4","orange","tomato")
col.atgc <- c("tomato","darkblue","orange","green4")
lev.mod <- c("Basic Statistics", "Per base sequence quality",
             "Per tile sequence quality","Per sequence quality scores",
             "Per base sequence content","Per sequence GC content",
             "Per base N content","Sequence Length Distribution",
             "Sequence Duplication Levels","Overrepresented sequences",
             # Overrepresented Sequences,
             "Adapter Content","Kmer Content")

## Set file paths
f <- read.table(p, header = FALSE, sep = "\t", stringsAsFactors = FALSE)[,1]

# debugging
# qc.dir = "."
# qc.pattern = ".fastqc.zip"
# f <- list.files(path = qc.dir, pattern = pattern, full.names = TRUE)
# pdf.name <- "fastqc.pdf"
# pdf.height <- 12
# pdf.width <- 12

################################################################################

## patch for qc_read_collection
# Error: Can't combine `Sample1$Max Obs/Exp Position` <double> and 
# `Sample2$Max Obs/Exp Position` <character>.
# -> happens if there are samples with zero rows in kmer_content statistics
qc_read_collection <- function(files, sample_names, modules = "all", verbose=T) 
{
  module_data <- lapply(files, qc_read, modules = modules, 
                        verbose = verbose)
  if (missing(sample_names) || length(sample_names) != length(files)) {
    sample_names <- lapply(module_data, function(x) unique(x$summary))
    sample_names <- unlist(sample_names)
  }
  names(module_data) <- sample_names
  module_names <- unique(unlist(lapply(module_data, names)))
  res <- list()
  for (i in seq_along(module_names)) {
    res[[i]] <- lapply(module_data, function(x) as.data.frame(x[[module_names[i]]]))
  }
  names(res) <- module_names
  
  ##<##<## begin patch
  res <- lapply(res, function(x) {
    dcl <- dplyr::bind_rows(lapply(x, function(y) {
      if (nrow(y) > 0) sapply(y, class)
      }))
    cl <- apply(dcl, 2, function(z) ifelse(any(z=="character"),"character",z[1]))
    lapply(x, function(w) {
      if (nrow(w) > 0) {for (i in names(w)) {class(w[,i]) <- cl[i]} ; w}
    })
  })
  ##<##<## end patch
  
  res <- lapply(res, dplyr::bind_rows, .id = "sample")
  res <- structure(res, class = c("list", "qc_read_collection"))
  res
}

################################################################################

## Load libraries
suppressPackageStartupMessages(require(fastqcr))
suppressPackageStartupMessages(require(ggplot2))

## Check arguments
stopifnot(any(file.exists(f)),
          is.numeric(pdf.height),
          is.numeric(pdf.width),
          pdf.height > 0,
          pdf.width > 0)

if (any(!file.exists(f))) {
  n <- f[which(!file.exists(f))]
  warning("These paths do not exist:\n", paste(n, collapse = "\n"))
  f <- f[!f%in%n]
}

## Read sample fastqc results
sample_names <- sub("_fastqc.zip$", "", basename(f))
qcrs <- qc_read_collection(files = f, verbose = FALSE, sample_names = sample_names)

## Read aggregated fastqc results
# status
qc.all <- qc_aggregate(qc.dir = unique(dirname(f))[1], progressbar = FALSE)
qc.all$tot.seq <- as.numeric(qc.all$tot.seq)
qc <- subset(qc.all, sample %in% sample_names)
qc$sample <- factor(qc$sample, levels = rev(sample_names))
qc$module <- factor(qc$module, levels = rev(lev.mod))

# status table
d.qc <- as.data.frame(table(qc$module, qc$status))
colnames(d.qc) <- c("module","test","Freq")
d.qc$test <- factor(d.qc$test, levels = c("PASS","WARN","FAIL"))
d.qc <- subset(d.qc, module != "Basic Statistics")

# nb_fail, nb_pass, nb_warn, failed, warned
suppressWarnings(s <- summary(qc))

# pct.dup, pct.gc, tot.seq, seq.length (also see qc_fails()/warns/problems)
stats <- qc_stats(qc)
ll <- strsplit(stats$seq.length, split = "-")
if (all(lengths(ll) == 1)) stats$seq.length <- as.numeric(stats$seq.length)
stats$seq.minlength <- as.numeric(sapply(ll, "[", 1))
stats$seq.maxlength <- as.numeric(sapply(ll, "[", 2))
if (all(is.na(stats$seq.maxlength))) stats$seq.maxlength <- stats$seq.minlength

if (nrow(qcrs[["basic_statistics"]]) > 0) {
  dseq.poor <- subset(qcrs[["basic_statistics"]],
                      Measure == "Sequences flagged as poor quality")
  dseq.poor$seq.poor <- as.numeric(dseq.poor$Value)
  stats <- merge(stats, dseq.poor[,c("sample","seq.poor")], by = "sample")
} else {
  stats$seq.poor <- 0
}

## Plot warnings and errors BARPLOT
pe <- ggplot(d.qc, aes(x = module, y = Freq, fill = test)) +
  geom_col(position = "dodge", alpha = 0.8) +
  labs(x = "", y = "Number of Samples", fill = "") +
  scale_fill_manual(values = col.test) +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "bottom") +
  ggtitle("FastQC WARNINGS and ERRORS")

## Plot warnings and errors HEATMAP
dstatus <- qcrs[[1]]
if (nrow(dstatus) > 0) {
  dstatus$status <- (factor(dstatus$status, levels = c("PASS", "WARN", "FAIL")))
  dstatus$module <- factor(dstatus$module, levels = (lev.mod))
  dstatus$sample <- factor(dstatus$sample, levels = rev(sample_names))
  dstatus <- subset(dstatus, module != "Basic Statistics")
  ps <- ggplot(dstatus, aes(x = module, y = sample, fill = status)) +
    geom_tile(alpha = 0.8) +
    scale_fill_manual(values = col.test) +
    labs(x = "", y = "", fill = "") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          legend.position = "bottom") +
    ggtitle("FastQC WARNINGS and ERRORS")
} else {
  ps <- NULL
}


## Plot simple stats
p1 <- ggplot(stats, aes(x = sample, y = tot.seq)) +
  geom_col(fill = "#825D99") +
  labs(x = "", y = "total number of reads") +
  scale_y_log10(breaks = c(0,1E01,1E02,1E03,1E04,1E05,1E06,1E07,1E08,1E09)) +
  coord_flip() +
  theme_bw() +
  ggtitle("Number of Reads")
p2 <- ggplot(stats, aes(x = sample, y = pct.dup)) +
  geom_col(fill = "#F06C45") +
  labs(x = "", y = "% of duplicate reads") +
  scale_y_continuous(limits = c(0, 100)) +
  coord_flip() +
  theme_bw() +
  ggtitle("Duplicate Reads")
p3 <- ggplot(stats, aes(x = sample, y = seq.poor)) +
  geom_col(fill = "#6F9E4C") +
  labs(x = "", y = "number of poor quality sequences") +
  coord_flip() +
  theme_bw() +
  ggtitle("Sequences Flagged as Poor Quality")
p4 <- ggplot(stats, aes(x = sample, y = pct.gc)) +
  geom_col(fill = "#B89B74") +
  labs(x = "", y = "% of GC content") +
  scale_y_continuous(limits = c(0, 100)) +
  coord_flip() +
  theme_bw() +
  ggtitle("GC Content")
if (all(lengths(ll) == 1)) {
  p5 <- ggplot(stats, aes(x = sample)) +
    geom_col(aes(y = seq.length), fill = "#A6CEE3")
} else {
  p5 <- ggplot(stats, aes(x = sample)) +
    geom_errorbar(aes(ymin = seq.minlength, ymax = seq.maxlength), width = .2,
                  colour = "#A6CEE3", size = 2, position = position_dodge(0.05))
}
p5 <- p5 +
  labs(x = "", y = "sequence length (bp)") +
  scale_y_log10(breaks = c(seq(0,100,by=20),150,300,1E03,1E04,1E05,1E06,1E07,1E08,1E09,1E10)) +
  coord_flip() +
  theme_bw() +
  ggtitle("Sequence length")

## Complex stats
# Per base sequence quality
d1 <- qcrs$per_base_sequence_quality
if (nrow(d1) > 0) {
  d1$sample <- factor(d1$sample, levels = sample_names)
  d1$Base <- factor(d1$Base, levels = d1$Base[!duplicated(d1$Base)])
  names(d1) <- gsub("90th", "X90", gsub("10th", "X10", gsub(" ", "", names(d1))))
  
  pl1 <- ggplot(d1, aes(x = Base, y = Mean, group = sample)) +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 0, ymax = 20),
              fill = col.test[3], alpha = 0.01) +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 20, ymax = 28),
              fill = col.test[2], alpha = 0.01) +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin =28, ymax = Inf),
              fill = col.test[1], alpha = 0.01) +
    geom_path() +
    geom_errorbar(aes_string(ymin="X10Percentile", ymax="X90Percentile"), width=.2,
                  position=position_dodge(0.05)) +
    scale_y_continuous(breaks = c(seq(0,20,by=5),seq(22,36,by=2),37),
                       limits = c(0,37)) +
    labs(x = "Position in read (bp)", y = "") +
    facet_wrap(~sample) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle("Per base sequence quality")
} else {
  pl1 <- NULL
}

# Per sequence quality scores
d2 <- qcrs$per_sequence_quality_scores
if (nrow(d2) > 0) {
  d2$sample <- factor(d2$sample, levels = sample_names)
  
  pl2.1 <- ggplot(d2, aes(x = Quality, y = Count)) +
    geom_path(colour = col.sample, alpha = 0.5) +
    scale_x_continuous(breaks = c(seq(0,36,by=2),37)) +
    labs(x = "Mean Sequence Quality (Phred Score)",
         y = "Average Quality per read") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle("Per sequence quality scores - All")
  
  pl2.2 <- ggplot(d2, aes(x = Quality, y = Count)) +
    geom_path(colour = col.sample) +
    scale_x_continuous(breaks = c(seq(0,36,by=2),37)) +
    labs(x = "Mean Sequence Quality (Phred Score)",
         y = "Average Quality per read") +
    facet_wrap(~sample) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle("Per sequence quality scores")
} else {
  pl2 <- NULL
}

# Per base sequence content
d3 <- qcrs$per_base_sequence_content
if (nrow(d3) > 0) {
  d3$sample <- factor(d3$sample, levels = sample_names)
  d3$Base <- factor(d3$Base, levels = d3$Base[!duplicated(d3$Base)])
  d3 <- reshape::melt(d3, id = c("sample","Base"))
  d3$variable <- factor(d3$variable, levels = c("A","T","G","C"))
  
  pl3 <- ggplot(d3, aes(x = Base, y = value, group = variable, colour = variable)) +
    geom_path() +
    labs(x = "Position in read (bp)", y = "Sequence content across all bases (%)",
         colour = "") +
    facet_wrap(~sample) +
    scale_colour_manual(values = col.atgc) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          legend.position = "bottom") +
    ggtitle("Per base sequence content")
} else {
  pl3 <- NULL
}

# Per sequence GC content
d4 <- qcrs$per_sequence_gc_content
if (nrow(d4) > 0) {
  d4$sample <- factor(d4$sample, levels = sample_names)
  names(d4) <- gsub(" ", "", names(d4))
  
  pl4.1 <- ggplot(d4, aes(x = GCContent, y = Count)) +
    geom_path(colour = col.sample, alpha = 0.5) +
    labs(x = "Mean GC content (%)",
         y = "GC count per read") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle("Per sequence GC content - All")
  
  pl4.2 <- ggplot(d4, aes(x = GCContent, y = Count)) +
    geom_path(colour = col.sample) +
    labs(x = "Mean GC content (%)",
         y = "GC count per read") +
    facet_wrap(~sample) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle("Per sequence GC content")
} else {
  pl4 <- NULL
}

# Per base N content
d5 <- qcrs$per_base_n_content
if (nrow(d5) > 0) {
  d5$sample <- factor(d5$sample, levels = sample_names)
  d5$Base <- factor(d5$Base, levels = d5$Base[!duplicated(d5$Base)])
  names(d5) <- gsub("-", "", names(d5))
  
  pl5 <- ggplot(d5, aes(x = Base, y = NCount, group = sample)) +
    geom_path() +
    scale_y_sqrt(breaks = c(0, 0.1, 0.33, 0.67, 1:4, seq(5, 25, by = 5),
                            seq(50, 100, by = 25)), limits = c(0, 100)) +
    labs(x = "Position in read (bp)", y = "Percent N (%)") +
    facet_wrap(~sample) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          legend.position = "bottom") +
    ggtitle("Per base N content")
} else {
  pl5 <- NULL
}

# Sequence Length Distribution
d6 <- qcrs$sequence_length_distribution
if (nrow(d6) > 0) {
  d6$sample <- factor(d6$sample, levels = sample_names)
  ll <- strsplit(as.character(d6$Length), split = "-")
  d6$seq.low <- sapply(ll, "[", 1)
  d6$seq.high <- sapply(ll, "[", 2)
  if (all(is.na(d6$seq.high))) d6$seq.high <- d6$seq.low
  d6$seq.high[is.na(d6$seq.high)] <- d6$seq.low[is.na(d6$seq.high)]
  d6$seq.range <- paste(d6$seq.low, d6$seq.high, sep = "-")
  d6 <- subset(d6, Count > 0)
  o <- order(as.numeric(d6$seq.low))
  d6$seq.range <- factor(d6$seq.range, levels = unique(d6$seq.range[o]))
  
  pl6 <- ggplot(d6, aes(x = seq.range, y = Count, group = sample)) +
    geom_point(colour = col.sample) +
    geom_path(colour = col.sample) +
    scale_y_log10(breaks = c(0,1E01,1E02,1E03,1E04,1E05,1E06,1E07,1E08,1E09)) +
    labs(x = "Sequence length (bp)",
         y = "Count") +
    facet_wrap(~sample) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle("Sequence length distribution")
} else {
  pl6 <- NULL
}

# Sequence Duplication Levels
d7 <- qcrs$sequence_duplication_levels
if (nrow(d7) > 0) {
  d7$sample <- factor(d7$sample, levels = sample_names)
  names(d7)[names(d7) == "Duplication Level"] <- "DL"
  d7$DL <- factor(d7$DL, levels = d7$DL[!duplicated(d7$DL)])
  d7 <- reshape::melt(d7, id = c("sample","DL"))
  
  pl7 <- ggplot(d7, aes(x = DL, y = value, group = variable, colour = variable)) +
    geom_path() +
    lims(y = c(0, 100)) +
    labs(colour = "", x = "Sequence Duplication Level", 
         y = "Percent of sequences remaining if deduplicated (%)") +
    facet_wrap(~sample) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          legend.position = "bottom") +
    ggtitle("Sequence duplication levels")
} else {
  pl7 <- NULL
}

# Adapter Content
d8 <- qcrs$adapter_content
if (nrow(d8) > 0) {
  d8$sample <- factor(d8$sample, levels = sample_names)
  d8$Position <- factor(d8$Position, levels = d8$Position[!duplicated(d8$Position)])
  d8 <- reshape::melt(d8, id = c("sample","Position"))
  
  pl8 <- ggplot(d8, aes(x = Position, y = value, group = variable, colour = variable)) +
    geom_path() +
    scale_y_sqrt(breaks = c(0, 0.1, 0.33, 0.67, 1:4, seq(5, 25, by = 5),
                            seq(50, 100, by = 25)), limits = c(0, 100)) +
    labs(colour = "", x = "Position in read (bp)", 
         y = "Percent Adapters (%)") +
    facet_wrap(~sample) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          legend.position = "bottom") +
    ggtitle("Adapter content")
} else {
  pl8 <- NULL
}

# Kmer Content
d9 <- qcrs$kmer_content
if (nrow(d9) > 0) {
  d9$sample <- factor(d9$sample, levels = sample_names)
  d9$Count <- factor(d9$Count)
  
  pl9 <- ggplot(d9, aes(x = Count)) +
    geom_bar() +
    labs(colour = "", x = "Kmer Count", y = "Count") +
    facet_wrap(~sample, scales = "free_x") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          legend.position = "bottom") +
    ggtitle("Kmer content")
} else {
  pl9 <- NULL
}

### Write simple stats
write.table(stats, file = sub(".pdf$", ".txt", pdf.name), sep = "\t", row.names = F, quote = F)

### Create PDF
pdf(pdf.name, width = pdf.width, height = pdf.height)

## Errors and warnings
print(pe)
print(ps)

## Simple stats
print(p1)
print(p2)
print(p3)
print(p4)
print(p5)

## Complex stats
# qc_plot_collection(qc = qcrs, modules = lev.mod[!lev.mod %in% "Basic Statistics"])
pl1
pl2.1
pl2.2
pl3
pl4.1
pl4.2
pl5
suppressMessages(print(pl6))
pl7
pl8
pl9

graphics.off()
