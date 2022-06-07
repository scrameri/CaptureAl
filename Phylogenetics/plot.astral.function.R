## Define root.treedata helper functions by Guangchuang Yu
# https://rdrr.io/github/GuangchuangYu/treeio/src/R/utilities.R
# https://rdrr.io/github/YuLab-SMU/treeio/src/R/method-drop-tip.R
root.treedata <- function(phy, outgroup, node = NULL, edgelabel = TRUE, ...){
  # helpers
  old_new_node_mapping <- function(oldtree, newtree){
    treelab1 <- oldtree %>% 
      as_tibble() %>%
      dplyr::select(c("node", "label"))
    treelab2 <- newtree %>% 
      as_tibble() %>%
      dplyr::select(c("node", "label"))
    node_map <- dplyr::inner_join(treelab1, treelab2, by="label") %>%
      dplyr::select(c("node.x", "node.y")) %>%
      dplyr::rename(c(old="node.x", new="node.y"))
    return(node_map)
  }
  
  build_new_labels <- function(tree){
    node2label_old <- tree %>% as_tibble() %>% dplyr::select(c("node", "label")) 
    if (inherits(tree, "treedata")){
      tree <- tree@phylo
    }
    tree$tip.label <- paste0("t", seq_len(Ntip(tree)))
    tree$node.label <- paste0("n", seq_len(Nnode(tree)))
    node2label_new <- tree %>% as_tibble() %>% dplyr::select(c("node", "label")) 
    old_and_new <- node2label_old %>% 
      dplyr::inner_join(node2label_new, by="node") %>%
      dplyr::rename(old="label.x", new="label.y") 
    return (list(tree=tree, node2old_new_lab=old_and_new))
  }
  
  build_new_tree <- function(tree, node2old_new_lab){
    # replace new label with old label
    treeda <- tree %>% as_tibble()
    treeda1 <- treeda %>%
      dplyr::filter(.data$label %in% node2old_new_lab$new)
    treeda2 <- treeda %>%
      dplyr::filter(!(.data$label %in% node2old_new_lab$new))
    # original label
    treeda1$label <- node2old_new_lab[match(treeda1$label, node2old_new_lab$new), "old"] %>%
      unlist(use.names=FALSE)
    treeda <- rbind(treeda1, treeda2)
    tree <- treeda[order(treeda$node),] %>% as.phylo() 
    return (tree)
  }
  
  # root
  if (!missing(outgroup) && is.character(outgroup)){
    outgroup <- match(outgroup, phy@phylo$tip.label)
  }
  if (!edgelabel){
    ## warning message
    message("The use of this method may cause some node data to become incorrect (e.g. bootstrap values) if 'edgelabel' is FALSE.")
  }
  #object <- phy
  # generate node old label and new label map table.
  res <- build_new_labels(tree=phy)
  tree <- res$tree
  node2oldnewlab <- res$node2old_new_lab
  # reroot tree
  re_tree <- root(tree, outgroup = outgroup, node = node,
                  edgelabel = edgelabel, ...)
  
  node_map <- old_new_node_mapping(tree, re_tree)
  n.tips <- Ntip(re_tree)
  
  # replace new label with old label
  phy@phylo <- build_new_tree(tree=re_tree, node2old_new_lab=node2oldnewlab)
  
  # update data or extraInfo function
  update_data <- function(data, node_map) {
    cn <- colnames(data)
    cn <- cn[cn != "node"]
    data <- dplyr::inner_join(data, node_map, by=c("node"="old")) %>%
      dplyr::select(c("new", cn)) %>% 
      dplyr::rename(node=.data$new)
      
    # clear root data
    root <- data$node == (n.tips + 1)
    data[root,] <- NA
    data[root,'node'] <- n.tips + 1
    return(data)
  }
  if (nrow(phy@data) > 0) {
    phy@data <- update_data(phy@data, node_map)
  }    
  if (nrow(phy@extraInfo) > 0){
    phy@extraInfo <- update_data(phy@extraInfo, node_map)
  }
  
  return(phy)
}
  
## plot.treedata
# ## debug
# dt = NULL; astral = TRUE; ladderize = TRUE; print = TRUE; verbose = TRUE;
# root = FALSE; outgroup = NULL;
# root.add = 0.1; tip.add = 0.1; col.edge = "black"; size.edge = 0.25; alpha.edge = 1;
# col.tip = "black"; size.tip = 2.5; adj.tip = -0.02; xmax = 0.3;
# add.support = TRUE; size.support = 2; alpha.support = 1; legend.pos.support = "left";
# supportvar = if (all(pievars %in% colnames(treedata@data))) "pp1" else "support"; supportbreaks = c(0, 50, 80, 95, 99, 101);
# supportlabs = c("< 50%","[50-80)%","[80-95)%", "[95-99)%","[99-100]%");
# supportcols = c("red","orange","yellow","#0096FF","black");
# supportname = if (all(pievars %in% colnames(treedata@data))) "ASTRAL\nPosterior Support" else "Bootstrap Support";
# add.pies = TRUE; width.pie = 0.075; height.pie = width.pie; alpha.pie = 1;
# hjust.pie = 0.05; vjust.pie = 0.05;
# pievars = c("q1","q2","q3"); piecols = c("#4567AD","#EB6533","#CCC0B1")

## function to plot RAxML and ASTAL trees (with node pies)
plot.treedata <- function(
  treedata, dt = NULL, ladderize = TRUE, print = TRUE, verbose = TRUE,
  root = FALSE, outgroup = NULL,
  root.add = 0.1, tip.add = 0.1, col.edge = "black", size.edge = 0.25, alpha.edge = 1,
  col.tip = "black", size.tip = 2.5, adj.tip = -0.02, xmax = 0.3,
  add.support = TRUE, size.support = 2, alpha.support = 1, legend.pos.support = "left",
  supportvar = if (all(pievars %in% colnames(treedata@data))) "pp1" else "support", supportbreaks = c(0, 50, 80, 95, 99, 101),
  supportlabs = c("< 50%","[50-80)%","[80-95)%", "[95-99)%","[99-100]%"),
  supportcols = c("red","orange","yellow","#0096FF","black"),
  supportname = if (all(pievars %in% colnames(treedata@data))) "ASTRAL\nPosterior Support" else "Bootstrap Support",
  add.pies = TRUE, width.pie = 0.075, height.pie = width.pie, alpha.pie = 1,
  hjust.pie = 0.05, vjust.pie = 0.05,
  pievars = c("q1","q2","q3"), piecols = c("#4567AD","#EB6533","#CCC0B1")) {
  
  ## Usage
  # treedata   input treedata object (as obtained from treeio::read.astral(treefile))
  # dt         data.frame with information for tips (color field)
  # additional plotting parameters (colors, sizes, labels, adjustments)
  
  ## Value
  # ggtree plot (class "ggree", "gg", "ggplot")
  
  ## Author
  # sfcrameri@gmail.com, Feb 2022
  
  ## Load libraries
  require("ggtree") # tested on v3.2.1
  require("treeio") # tested on v1.18.1
  
  ## Define root.treedata helper functions by Guangchuang Yu
  # https://rdrr.io/github/GuangchuangYu/treeio/src/R/utilities.R
  # https://rdrr.io/github/YuLab-SMU/treeio/src/R/method-drop-tip.R
  root.treedata <- function(phy, outgroup, node = NULL, edgelabel = TRUE, ...){
    if (!missing(outgroup) && is.character(outgroup)){
      outgroup <- match(outgroup, phy@phylo$tip.label)
    }
    if (!edgelabel){
      ## warning message
      message("The use of this method may cause some node data to become incorrect (e.g. bootstrap values) if 'edgelabel' is FALSE.")
    }
    #object <- phy
    # generate node old label and new label map table.
    res <- build_new_labels(tree=phy)
    tree <- res$tree
    node2oldnewlab <- res$node2old_new_lab
    # reroot tree
    re_tree <- root(tree, outgroup = outgroup, node = node,
                    edgelabel = edgelabel, ...)
    
    node_map <- old_new_node_mapping(tree, re_tree)
    n.tips <- Ntip(re_tree)
    
    # replace new label with old label
    phy@phylo <- build_new_tree(tree=re_tree, node2old_new_lab=node2oldnewlab)
    
    # update data or extraInfo function
    update_data <- function(data, node_map) {
      cn <- colnames(data)
      cn <- cn[cn != "node"]
      data <- dplyr::inner_join(data, node_map, by=c("node"="old")) %>%
        dplyr::select(c("new", cn)) %>% 
        dplyr::rename(node=.data$new)
      
      # clear root data
      root <- data$node == (n.tips + 1)
      data[root,] <- NA
      data[root,'node'] <- n.tips + 1
      return(data)
    }
    if (nrow(phy@data) > 0) {
      phy@data <- update_data(phy@data, node_map)
    }    
    if (nrow(phy@extraInfo) > 0){
      phy@extraInfo <- update_data(phy@extraInfo, node_map)
    }
    
    return(phy)
  }
  
  old_new_node_mapping <- function(oldtree, newtree){
    treelab1 <- oldtree %>% 
      as_tibble() %>%
      dplyr::select(c("node", "label"))
    treelab2 <- newtree %>% 
      as_tibble() %>%
      dplyr::select(c("node", "label"))
    node_map <- dplyr::inner_join(treelab1, treelab2, by="label") %>%
      dplyr::select(c("node.x", "node.y")) %>%
      dplyr::rename(c(old="node.x", new="node.y"))
    return(node_map)
  }
  
  build_new_labels <- function(tree){
    node2label_old <- tree %>% as_tibble() %>% dplyr::select(c("node", "label")) 
    if (inherits(tree, "treedata")){
      tree <- tree@phylo
    }
    tree$tip.label <- paste0("t", seq_len(Ntip(tree)))
    tree$node.label <- paste0("n", seq_len(Nnode(tree)))
    node2label_new <- tree %>% as_tibble() %>% dplyr::select(c("node", "label")) 
    old_and_new <- node2label_old %>% 
      dplyr::inner_join(node2label_new, by="node") %>%
      dplyr::rename(old="label.x", new="label.y") 
    return (list(tree=tree, node2old_new_lab=old_and_new))
  }
  
  build_new_tree <- function(tree, node2old_new_lab){
    # replace new label with old label
    treeda <- tree %>% as_tibble()
    treeda1 <- treeda %>%
      dplyr::filter(.data$label %in% node2old_new_lab$new)
    treeda2 <- treeda %>%
      dplyr::filter(!(.data$label %in% node2old_new_lab$new))
    # original label
    treeda1$label <- node2old_new_lab[match(treeda1$label, node2old_new_lab$new), "old"] %>%
      unlist(use.names=FALSE)
    treeda <- rbind(treeda1, treeda2)
    tree <- treeda[order(treeda$node),] %>% as.phylo() 
    return (tree)
  }
  
  
  ## Handle node support values for non-ASTRAL trees
  stopifnot(inherits(treedata, "treedata"))
  astral <- all(pievars %in% colnames(treedata@data)) # guess if it's an ASTRAL tree
  
  # make support values (node label) numeric 
  if (!astral & is.character(treedata@phylo$node.label)) {
    treedata@phylo$node.label <- suppressWarnings(as.numeric(treedata@phylo$node.label))
  }
  
  # transfer support values (node label) to treedata@data
  if (!astral & !supportvar %in% names(treedata@data)) {
    support <- treedata@phylo$node.label
    node <- ggtree(treedata)$data[!ggtree(treedata)$data$label %in% treedata@phylo$tip.label,c("node")]
    treedata@data <- tidytree::tibble(support = support, node)
  }
  
  
  ## Create plot-related data.frames
  # data.frame of node labels (support values / ASTRAL pies) -> dn
  dn <- if (astral) {
    data.frame(treedata@data[,c(supportvar, pievars, "node")])
  } else {
    if (identical(Nnode(treedata), nrow(treedata@data))) {
      data.frame(treedata@data[,c(supportvar, "node")])
    } else {
      # experimental, please check
      data.frame(rbind(tidytree::tibble(s = NA, node = NA) %>% 
                         tidytree::rename(!!supportvar := s),
                       treedata@data[,c(supportvar, "node"),drop=F]))
    }
  }
  
  # data.frame of tips -> dt
  tl <- treedata@phylo$tip.label
  if (is.null(dt)) {
    dt <- data.frame(color = rep(col.tip, Ntip(treedata@phylo)),
                     label = tl)
  } else {
    if (!"color" %in% colnames(dt)) {
      dt$color <- col.tip
    }
    if (!"label" %in% colnames(dt)) {
      dt$label <- tl
    }
    if (all(tl %in% rownames(dt))) {
      dt <- dt[tl,, drop = FALSE]
    }
  }
  
  ## Re-root tree if indicated
  if (root & !is.null(outgroup)) {
    # message("rooting tree")
    treedata <- root.treedata(phy = treedata, outgroup = outgroup, edgelabel = TRUE)
  }
  
  ## Check (gives FALSE if ATRAL trees are rooted using root.treedata())
  # stopifnot(all.equal(Nnode(treedata), nrow(dn)),
  #           all.equal(Ntip(treedata), nrow(dt)))
  
  ## Add branch lengths
  # treedata@phylo$edge.length[is.nan(treedata@phylo$edge.length)] <- 0
  
  # if the root is resolved, the rooted ASTRAL tree has zero branch length for the outgroup
  p <- treedata@phylo
  
  # rn.idx <- which(p$edge[,1] == rootnode(treedata) & p$edge[,2] == Ntip(p) + Nnode(p)) # wrong?
  rn.idx <- which(p$edge[,1] == rootnode(treedata)) # check: gives more than 1?
  rn <- rn.idx[is.nan(treedata@phylo$edge.length[rn.idx])]
  
  if (length(rn) > 0) {
    if (verbose) cat("adding root.add =", root.add, 
                     "branch units to the root, which has zero branch length\n")
    treedata@phylo$edge.length[rn] <- root.add
  }
  
  # raw ASTRAL output has no edge.length on tips
  tn <- which(is.nan(treedata@phylo$edge.length))
  
  if (length(tn) > 0) {
    if (verbose) cat("adding tip.add =", tip.add, 
                     "branch units to tips, which have zero branch length\n")
    treedata@phylo$edge.length[tn] <- tip.add
  }
  
  ## Assing tip labels
  treedata@phylo$tip.label <- dt$label
  
  ## Basic tree
  gg <- ggtree(treedata, color = col.edge, root.position = 0,
               ladderize = ladderize, alpha = alpha.edge,
               size = size.edge) +
    geom_tiplab(size = size.tip, color = dt$color, hjust = adj.tip) +
    geom_treescale(0, treedata@phylo$Nnode, offset = 2, color = col.edge) +
    ggplot2::scale_color_discrete(guide = "none") +
    theme_tree(bgcolor = "white")
  
  # increase margin on the right (total.width)
  xlims <- max(ggplot2::ggplot_build(gg)$data[[1]]$x, na.rm = T)
  width <- xlims
  hjust <- -2*width
  total.width <- xlims+xmax*xlims
  
  gg <- gg + ggplot2::xlim(0, total.width)

  ## Create ASTRAL pies
  if (astral & add.pies) {
    pies <- nodepie(dn, cols = pievars, alpha = alpha.pie)
    pies <- lapply(pies, function(g) g + scale_fill_manual(values = piecols))
    gg <- suppressWarnings(
      gg + geom_inset(pies, width = width.pie, height = height.pie, 
                      hjust = hjust.pie, vjust = vjust.pie))
  }
  
  ## Map support values onto phylogeny
  if (add.support) {
    support.multiplyer <- if (astral) 100 else 1
    
    # classify support values using supportbreaks and supportlabs
    gg$data <- gg$data %>%
      tidytree::mutate(labelclass =
        factor(
          cut(support.multiplyer*gg$data[[supportvar]],
              breaks = supportbreaks, labels = supportlabs, right = FALSE),
          levels = supportlabs))
    
    # plot support values
    gg <- gg +
      geom_nodepoint(aes(fill = labelclass), shape = 21, size = size.support,
                     alpha = alpha.support) +
      scale_fill_manual(name = supportname, values = supportcols, drop = F) +
      theme(legend.position = legend.pos.support)
  }
  
  ## Print
  if (print) print(gg)

  ## Return
  invisible(gg)
}

save(plot.treedata, file = "plot.treedata.rda")

