###########################
## Plot ASTRAL node pies ##
###########################

## sfcrameri@gmail.com

## packages
library("ggtree") # BiocManager::install("ggtree", version = "3.14", force = TRUE) for R 4.1
library("treeio") # tested on v1.18.1


## load plot.astral function and example tree
# load("plot.astral.rda") # plot.astral
load("plot.treedata.rda") # generic function plot.treedata() -> plot(treedata-object) works

# plot.treedata is the function I wrote to plot ASTRAL annotations as node-pie-charts
# treedata is the only required argument, the others may need more documentation
args(plot.treedata)


## try plot.astral on rooted tree
load("treedata/mytree.rda")

# inspect
# this example tree is rooted, and 0.1 coalescent units have been added to the zero-length tips
# it has a @phylo component with the tree and a @data component with ASTRAL annotations
mytree

# dtip is an optional data.frame, but helpful for tip label customization
# the data.frame can have a <color>, <alpha>, <size>, and/or <label> column with tip-specific settings
# either the row.names correspond to mytree@phylo$tip.label or rows are in tip order
head(dtip)

# by default, ASTRAL support values and quartet proportions are plotted
plot(treedata = mytree)

# simple tree
plot(treedata = mytree, add.pies = FALSE, add.support = FALSE)

# use size.support and alpha.support arguments
plot(treedata = mytree, add.pies = FALSE, size.support = 3, alpha.support = 0.65)

# use width.pie and alpha.pie arguments
plot(treedata = mytree, add.support = FALSE, width.pie = 0.05, alpha.pie = 0.9)

# if you specify the dt argument, you can adjust tip labels using a data.frame
plot(treedata = mytree, dt = dtip, col.edge = "blue", size.edge = 1, alpha.edge = 0.5,
     width.pie = 0.06, alpha.pie = 0.85, legend.pos.support = c(0.1,0.8))

# a low xmax uses more x axis space, a higher xmax gives space for long labels
# set width.pie and size.support to deal with different tree sizes and x.margin
plot(treedata = mytree, xmax = 1.2, width.pie = 0.1, size.support = 1.5,
     dt = data.frame(label = rep(paste(rep(LETTERS,2),collapse=""),Ntip(mytree@phylo))))



################################################################################

## FIGURE 1
# read using treeio:read.astral (read.tree also works but gives too long node labels)
mytree1 <- treeio::read.astral("treedata/mafft.110.1020.c0.5.d0.25.c0.4.n0.5.BS10.single.spectree") #raw ASTRAL output tree

# root the tree using treeio::root
mytree1@phylo <- treeio::root(phy = mytree1@phylo, outgroup = c("SC0124_R_S94","Polyg_S24_L001","SC0125_S95","SC0128_R_S96"),
                              edgelabel = TRUE, resolve.root = TRUE)

# this association of @phylo and @data slots gets messed up after re-rooting @phylo
plot(mytree1) # wrong association of pies and nodes, the backbone is highly supported in reality

# the @data slot needs to be fixed similarly to this (but it depends on the tree and rooting!)
mytree1@data[mytree1@data$node %in% c(113,115,120),-ncol(mytree1@data)] <- mytree1@data[mytree1@data$node %in% c(115,120,132),-ncol(mytree1@data)]
mytree1@data[mytree1@data$node %in% 132:218,-ncol(mytree1@data)] <- mytree1@data[mytree1@data$node %in% 133:219,-ncol(mytree1@data)]
mytree1@data[mytree1@data$node %in% 219,-ncol(mytree1@data)] <- mytree1@data[mytree1@data$node %in% 133,-ncol(mytree1@data)]

# plot
plot(mytree1, root.add = 0.5, tip.add = 0.5)


## FIGURE 2
# read using treeio:read.astral (read.tree also works but gives too long node labels)
mytree2 <- treeio::read.astral("treedata/mafft.63.2407.c0.5.d0.25.c0.4.n0.5.BS10.single.spectree") #raw ASTRAL output tree

# root the tree using treeio::root
mytree2@phylo <- treeio::root(phy = mytree2@phylo, outgroup = c("MBG030_S4_L001","MBG031_S10_L001"),
                              edgelabel = TRUE, resolve.root = TRUE)

# this association of @phylo and @data slots gets messed up after re-rooting @phylo
plot(mytree2) # wrong association of pies and nodes, the backbone is highly supported in reality

# the @data slot needs to be fixed similarly to this (but it depends on the tree and rooting!)
mytree2@data[mytree2@data$node %in% 65:124,-ncol(mytree2@data)] <- mytree2@data[mytree2@data$node %in% 66:125,-ncol(mytree2@data)]
mytree2@data[mytree2@data$node %in% 125,-ncol(mytree2@data)] <- mytree2@data[mytree2@data$node %in% 70,-ncol(mytree2@data)]

# plot
plot(mytree2, root.add = 0.5, tip.add = 0.5) # 0.1 coalescent units are added to zero-length terminal or root branches



## FIGURE S1
mytree1.2 <- treeio::read.astral("treedata/RAxML_bipartitions.mafft.110.1020.c0.5.d0.25.c0.4.n0.5.catBS100.nex") #raw RAxML output tree, no node labels read
mytree1.2@phylo <- treeio::read.tree(mytree1.2@file) # node labels read

# If resolve.root = TRUE, root adds a zero-length branch below the MRCA of the ingroup
mytree1.2@phylo <- treeio::root(phy = mytree1.2@phylo, outgroup = c("SC0124_R_S94","Polyg_S24_L001","SC0125_S95","SC0128_R_S96"),
                                edgelabel = TRUE, resolve.root = TRUE)

# plot
plot(treedata = mytree1.2, root.add = 0.1) # 0.1 are added to zero-length terminal or root branches


## FIGURE S2
mytree2.2 <- treeio::read.astral("treedata/RAxML_bipartitions.mafft.63.2407.c0.5.d0.25.c0.4.n0.5.catBS100.nex") #raw RAxML output tree, no node labels read
mytree2.2@phylo <- treeio::read.tree(mytree2.2@file) # node labels read

# If resolve.root = TRUE, root adds a zero-length branch below the MRCA of the ingroup
mytree2.2@phylo <- treeio::root(phy = mytree2.2@phylo, outgroup = c("MBG030_S4_L001","MBG031_S10_L001"),
                                edgelabel = TRUE, resolve.root = TRUE)

# plot
plot(treedata = mytree2.2, root.add = 0.1) # 0.1 are added to zero-length terminal or root branches
