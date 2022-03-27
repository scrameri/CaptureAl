# ## debug
# dt = NULL; astral = TRUE; ladderize = TRUE; print = TRUE; verbose = TRUE;
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
  
  ## Set defaults
  stopifnot(inherits(treedata, "treedata"))
  astral <- all(pievars %in% colnames(treedata@data))
  
  if (!astral & is.character(treedata@phylo$node.label)) {
    treedata@phylo$node.label <- suppressWarnings(as.numeric(treedata@phylo$node.label))
  }
  if (!astral & !supportvar %in% names(treedata@data)) {
    support <- treedata@phylo$node.label
    node <- ggtree(treedata)$data[!ggtree(treedata)$data$label %in% treedata@phylo$tip.label,c("node")]
    treedata@data <- tidytree::tibble(support = support, node)
    # if (identical(Nnode(treedata), nrow(treedata@data))) {
    #   treedata@data <- tidytree::tibble(support = treedata@phylo$node.label, treedata@data[,"node"])
    # } else {
    #   treedata@data <- tidytree::tibble(support = treedata@phylo$node.label, rbind(tidytree::tibble(node = NA), treedata@data[,"node"]))
    # }
  }
  
  tl <- treedata@phylo$tip.label
  dn <- if (astral) {
    data.frame(treedata@data[,c(supportvar, pievars, "node")])
  } else {
    if (identical(Nnode(treedata), nrow(treedata@data))) {
      data.frame(treedata@data[,c(supportvar, "node")])
    } else {
      data.frame(rbind(tidytree::tibble(s = NA, node = NA) %>% 
                         tidytree::rename(!!supportvar := s),
                       treedata@data[,c(supportvar, "node"),drop=F]))
    }
  }
  
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
  
  ## Check input
  stopifnot(all.equal(Nnode(treedata), nrow(dn)),
            all.equal(Ntip(treedata), nrow(dt)))
  
  ## Add branch lengths
  # if the root is resolved, the rooted ASTRAL tree has zero branch length for the outgroup
  p <- treedata@phylo
  rn <- which(p$edge[,1] == rootnode(treedata) & p$edge[,2] == Ntip(p) + Nnode(p))
  if (length(rn) != 0) {
    if (treedata@phylo$edge.length[rn] == 0) {
      if (verbose) cat("adding root.add =", root.add, 
                       "branch units to the root, which has zero branch length\n")
      treedata@phylo$edge.length[rn] <- root.add
    }
  }
  
  # raw ASTRAL output has no edge.length on tips
  if (sum(is.nan(treedata@phylo$edge.length)) == Ntip(treedata)) {
    if (verbose) cat("adding tip.add =", tip.add, 
                     "branch units to tips, which have zero branch length\n")
    treedata@phylo$edge.length[is.nan(treedata@phylo$edge.length)] <- tip.add
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
                      hjust = hjust.pie, vjust = vjust.pie) 
    )
  }
  
  ## Map support values onto phylogeny
  if (add.support) {
    support.multiplyer <- if (astral) 100 else 1
    
    gg$data <- gg$data %>%
      tidytree::mutate(labelclass = c(
        factor(rep(NA, nrow(gg$data)-nrow(dn)), levels = supportlabs),
        cut(support.multiplyer*dn[,supportvar],
            breaks = supportbreaks, labels = supportlabs, right = FALSE)))
    
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
