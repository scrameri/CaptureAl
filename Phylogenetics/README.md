# Post CaptureAl Workflows for Phylogenetics

## Phylogenetic trees including ASTRAL-III pie charts using R
- [plot_astral.R](plot_astral.R)

  - Read tree files outputted from RAxML (Stamatakis, 2014) or ASTRAL-III (Mirarab et al., 2014; Zhang et al., 2018) and convert it to a *treedata* object
  - Annotate the tree
  - Visualize the tree as shown below (ASTRAL-III quartet proportion as pie charts)

![Phylogenetics.png](https://raw.githubusercontent.com/scrameri/CaptureAl/master/Phylogenetics/Phylogenetics.png)

## Further Reading
- https://guangchuangyu.github.io/ggtree-book/chapter-ggtree.html

## References
- Mirarab, S., Reaz, R., Bayzid, M. S., Zimmermann, T., Swenson, M. S., & Warnow, T. (2014). ASTRAL: genome-scale coalescent-based species tree estimation. *Bioinformatics*, 30(17), i541–i548. [DOI](https://doi.org/10.1093/bioinformatics/btu462).
- Stamatakis, A. (2014). RAxML version 8: a tool for phylogenetic analysis and post-analysis of large phylogenies. *Bioinformatics*, 30(9), 1312–1313. [DOI](https://doi.org/10.1093/bioinformatics/btu033).
- Zhang, C., Rabiee, M., Sayyari, E., & Mirarab, S. (2018). ASTRAL-III: polynomial time species tree reconstruction from partially resolved gene trees. BMC *Bioinformatics*, 19, 153. [DOI](https://doi.org/10.1186/s12859-018-2129-y).
