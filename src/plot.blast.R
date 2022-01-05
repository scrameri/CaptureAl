#!/cluster/apps/r/3.6.0_openblas/x86_64/bin/Rscript

########################
## PLOT BLAST RESULTS ##
########################

## Usage
# plot.blast.R <blast results> <OPT: number of subjects to draw>

## Author
# simon.crameri@env.ethz.ch, Apr 2020

## Load required library
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(data.table))

## Get arguments
args = commandArgs(trailingOnly = TRUE)

## Check arguments
options(warning.length = 2000L)
if (! length(args) %in% c(1:2)) {
  stop("1 arguments needed (2 taken): 
       REQUIRED
       1) <input|CHR>: path to blast results ; 
       OPTIONAL
       2) <n|NUM>:     number of subjects to draw [DEFAULT: 30]",
       call.=FALSE)
}

## Set arguments
input <- as.character(args[1])
n <- as.numeric(as.character(args[2]))
if (is.na(n)) n <- 30

## Additional arguments
ext <- paste0(".", rev(unlist(strsplit(input, split = "[.]")))[1])
oname <- gsub(paste0(ext, "$"), ".pdf", input)

width = 30
height = 30

## Check arguments
stopifnot(file.exists(input), n > 0, n < 1000)

## Read BLAST file
b <- fread(input)

## Plot
# select subjects
todraw <- sort(unique(b$subject.length), decreasing = TRUE)[1:n]
b2 <- b[b$subject.length %in% todraw,]

# order for subject
b2 <- b2[order(b2$subject.id, b2$s.start),]
b2$order <- NA_integer_

for (i in unique(b2$subject.id)) {
  b2[subject.id == i,"order"] <- 1:sum(b2$subject.id == i)
}

b2$subject.id <- factor(b2$subject.id, levels = as.character(unique(b2$subject.id[order(b2$subject.length, decreasing = TRUE)])))
len.que <- sum(unlist(b[!duplicated(b$query.id),"query.length"]))
len.sub <- sum(unlist(b[!duplicated(b$subject.id),"subject.length"]))
p1 <- ggplot(b2) +
  geom_segment(aes(x = 1, xend = subject.length, y = 0, yend = 0)) +
  geom_point(aes(x = s.start, y = order), col = "tomato", size = 3, shape = 3) +
  geom_point(aes(x = s.start + query.length, y = order), col = "tomato", size = 3, shape = 3) +
  geom_segment(aes(x = s.start, xend = s.start + query.length, y = order, yend = order), 
                            col = "tomato", size = 3) + 
  labs(x = "position in subject", y = "Number of queries aligned to subject") +
  facet_wrap(~subject.id, scales = "free") + 
  theme_bw() +
  ggtitle(paste0("Total subjects length: ", len.sub, " (n = ", length(unique(b$subject.id)), ") ; ",
                 "Total queries length: ", len.que, " (n = ", length(unique(b$query.id)), ", ", 
                 round(100*len.que/len.sub, 4), "%)"))

b2$subject.id <- factor(b2$subject.id, levels = rev(unique(b2$subject.id)))
p2 <- ggplot(b2, aes(subject.id, alignment.length/query.length)) +
  geom_boxplot(fill = "lightblue", alpha = 0.5) +
  geom_violin(fill = "lightblue", alpha = 0.5) +
  geom_point() +
  labs(x = "") +
  coord_flip() +
  ylim(0,NA) +
  theme_bw()

p3 <- ggplot(b, aes(alignment.length, perc.identity, size = query.length)) +
  geom_point() +
  geom_smooth(method = 'loess', formula = 'y ~ x') +
  theme_bw()

p4 <- ggplot(b) +
  geom_density(aes(query.length), alpha = 0.5, fill = "tomato") +
  geom_density(aes(alignment.length), alpha = 0.5, fill = "lightblue") +
  geom_point(data = data.frame(x = rep(NA_integer_,2), y = rep(NA_integer_,2), fac = c("1","2")),
             aes(x, y, colour = fac), size = 5, na.rm = TRUE) +
  labs(x = "length") +
  scale_colour_manual(name = "", values = c("tomato","lightblue"), 
                      labels = c("query.length","alignment.length")) +
  theme_bw() +
  theme(legend.position = "top")

pdf(oname, width, height)
print(p1)
print(p2)
print(p3)
print(p4)
graphics.off()


