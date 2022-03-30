# Post CaptureAl Workflows for Population Genetics

## Analysis of SNP data using R
- [vcf_PCA_NJ.R](vcf_PCA_NJ.R)

  - Create a *genind* and a *genlight* object from a `*.vcf` file of filtered SNPs
  - Perform Principal Component Analysis (PCA) 
  - Perform Neighbor-Joining (NJ) tree Analysis

## Subsetting SNPs for Analysis using STRUCTURE
- [vcf_STRUCTURE.R](vcf_STRUCTURE.R)

  - Randomly select a number of SNPs from a `*.vcf` file of filtered SNPs, with low missingness per genomic region
  - Create a corresponding input file for [STRUCTURE](https://web.stanford.edu/group/pritchardlab/structure.html) (Pritchard et al., 2000)
  - Sort STRUCTURE output based on population structure / geo-coordinates for visualization using [CLUMPAK](http://clumpak.tau.ac.il/)

## Further Reading
- https://grunwaldlab.github.io/Population_Genetics_in_R/Introduction.html

## References
- Jombart, T. (2008). adegenet: a R package for the multivariate analysis of genetic markers. *Bioinformatics*, 24(11), 1403–1405. [DOI](https://doi.org/10.1093/bioinformatics/btn129).
- Jombart, T., & Ahmed, I. (2011). adegenet 1.3-1: new tools for the analysis of genome-wide SNP data. *Bioinformatics*, 27(21), 3070–3071. [DOI](https://doi.org/10.1093/bioinformatics/btr521).
- Kamvar, Z. N., Tabima, J. F., & Grünwald, N. J. (2014). Poppr: an R package for genetic analysis of populations with clonal, partially clonal, and/or sexual reproduction. *PeerJ*, 2, e281. [DOI](https://doi.org/10.7717/peerj.281).
- Knaus, B. J., & Grünwald, N. J. (2017). VCFR: a package to manipulate and visualize variant call format data in R. *Molecular Ecology Resources*, 17(1), 44–53. [DOI](https://doi.org/10.1111/1755-0998.12549).
- Paradis, E., & Schliep, K. (2018). ape 5.0: an environment for modern phylogenetics and evolutionary analyses in R. *Bioinformatics*, 35(3), 526–528. [DOI](https://doi.org/10.1093/bioinformatics/bty633).
- Pritchard, J. K., Stephens, M., & Donnelly, P. (2000). Inference of population structure using multilocus genotype data. *Genetics*, 155(2), 945–959.
- Wickham, H., Averick, M., Bryan, J., Chang, W., McGowan, L., François, R., Grolemund, G., Hayes, A., Henry, L., Hester, J., Kuhn, M., Pedersen, T., Miller, E., Bache, S., Müller, K., Ooms, J., Robinson, D., Seidel, D., Spinu, V., Takahashi, K., Vaughan, D., Wilke, C., Woo, K., & Yutani, H. (2019). Welcome to the Tidyverse. *The Journal of Open Source Software*, 4(43), 1686. [DOI](https://doi.org/10.21105/joss.01686).
