#!/cluster/apps/r/3.6.0_openblas/x86_64/bin/Rscript

## Usage: add.to.terminal.branches.R <.tre> <OPTIONAL: length added to terminal branches> <OPTIONAL: whether to extend existing terminal branches>

## Author: simon.crameri@env.ethz.ch, Mar 2020

## Motivation for this script
# by Siavash Mirarab: https://github.com/smirarab/ASTRAL/blob/master/astral-tutorial.md#viewing-results-of-astral
# ASTRAL does not generate terminal branch lengths.
#
# by Liam J. Revell:  http://www.phytools.org/eqg/Exercise_3.2/
# The matrix edge contains the beginning and ending node number for all the nodes and tips in the tree. 
# By convention, the tips of the tree are numbered 1 through n for n tips; 
# and the nodes are numbered n + 1 through n + m for m nodes

## Load libraries
suppressPackageStartupMessages(library(ape)) # read.tree

## Get arguments
args <- commandArgs(trailingOnly=TRUE)

## Check arguments
if (!length(args) %in% 1:3) {
  stop("3 arguments taken (1 needed):
       REQUIRED
       1) <tre>     CHR     path to newick-formatted tree file; 

       OPTIONAL
       2) <bl>      NUM     branch length added to terminal branches (branches leading to tips) [DEFAULT: 1]
       3) <addtobl> BOOLEAN if TRUE, will add <bl> to terminal branches with existing lengths; if FALSE, will only add <bl> to terminal branches without existing lengths [DEFAULT: FALSE]",
       call.=FALSE)
}

## Set arguments
get.boolean <- function(x, truestrings = c("T","TRUE", "t","true", "True"), 
                        falsestrings = c("F","FALSE","f","false","False")) {
  if (x %in% truestrings) {
    return(TRUE)
  } else {
    if (x %in% falsestrings) {
      return(FALSE)
    } else stop("<", x, "> cannot be interpreted. Please specify one of ", 
                paste(c(truestrings, falsestrings), collapse = ", "), ".")
  }
}
tre <- as.character(args[1])
bl <- as.numeric(as.character(args[2]))
if (is.na(bl)) bl <- 1
addtobl <- as.character(args[3])
if (is.na(addtobl)) {addtobl <- FALSE} else {addtobl <- get.boolean(addtobl)}

## Set arguments (for debugging)
# tre <- "test.spectree"
# bl <- 1
# addtobl = FALSE

## Additional arguments
suffix = ".blen" # output file: suffix will be added to input file name

## Check arguments
stopifnot(file.exists(tre), bl >= 0)

## Read tree
# use ape to read the tip.labels
myphy <- ape::read.tree(file = tre)

# use readLines to get plain text (because ape::read.tree and ape::write.tree appear not to handle full ASTRAL node annotation correctly)
mytree <- readLines(tre)

## Add <bl> to non-internal edge.lengths
newtrees <- list()
ntrees <- length(mytree)
for (tree in 1:ntrees) {
 
  if (ntrees > 1) {
    newtree <- mytree[[tree]]
    oldtree <- myphy[[tree]]
  } else {
    newtree <- mytree
    oldtree <- myphy
  }
 
  # index of tips (leaves) in edge matrix of phylo object
  leaves <- which(oldtree$edge[,2] %in% 1:length(oldtree$tip.label))
 
  # tips (leaves) with terminal branch length
  havebl <- oldtree$tip.label[which(!is.na(oldtree$edge.length[leaves]))]
 
  # terminal branch lengths
  thebl <- na.omit(oldtree$edge.length[leaves])
 
  for (i in oldtree$tip.label) {
 
   # check if branch lengths for tips are already included
    hasbl <- i %in% havebl

    # check if quoting is necessary (- in tip label)
    addquote <- length(grep("-", i)) == 1

    if (hasbl) {
      if (addtobl) {
        newbl <- thebl[havebl == i] + bl
      } else {
        newbl <- thebl[havebl == i]
      }
    } else {
      newbl <- bl
    }

    if (addquote) {
      newid <- paste0("'", i, "'")
    } else {
      newid <- i
    }

    new <- paste0(newid, ":", newbl)

    # tip labels are adjacent to either ) or , or : (if they have branch lengths)
    newtree <- gsub(paste0(i, ")"), paste0(new, ")"), newtree)
    newtree <- gsub(paste0(i, ","), paste0(new, ","), newtree)
    if (addtobl) newtree <- gsub(paste0(i,":[0-9]+.[0-9]+"), paste0(new, ""), newtree)

    newtrees[[tree]] <- newtree
  }
}

## Save tree
writeLines(text = unlist(newtrees), con = paste0(tre, suffix))

