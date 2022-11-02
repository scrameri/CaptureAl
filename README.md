# CaptureAl
Bioinformatic pipeline streamlining analysis of **target enrichment sequencing** data.

## How to Cite
Please cite the following paper when using CaptureAl in published work:
- [Crameri S., Fior S., Zoller S. & Widmer A. (2022) A target capture approach for phylogenomic analyses at multiple evolutionary timescales in rosewoods (*Dalbergia* spp.) and the legume family (Fabaceae). *Molecular Ecology Resources* 22(8):3087–3105.](https://doi.org/10.1111/1755-0998.13666)

## Before you start
1) Download/Clone the [CaptureAl repository](https://github.com/scrameri/CaptureAl) 
2) Follow the [installation instructions](https://github.com/scrameri/CaptureAl/blob/master/Install.md).
3) Determine how to implement [parallel computing](https://github.com/scrameri/CaptureAl/blob/master/Parallelize.md) in your computing environment (*GNU parallel* or *bsub*)

## Get started
Once the installation is complete, follow the [tutorial](https://github.com/scrameri/CaptureAl/blob/master/tutorial/) step by step. You can also skip preprocessing of reads and jump-start to CaptureAl **STEP 1** ([parallel](https://github.com/scrameri/CaptureAl/blob/master/tutorial/parallel/Step1_Read_Mapping.md) or [bsub](https://github.com/scrameri/CaptureAl/blob/master/tutorial/bsub/Step1_Read_Mapping.md)) if your reads are already quality-trimmed.

## Post-process Workflows for
- [Phylogenetics](Phylogenetics)
- [Population Genetics](PopulationGenetics)

## Pipeline overview
The pipeline is divided into **seven steps**, shown below. For each step, there is a limited number of scripts that have to be executed, and these scripts often execute other scripts located in the repository. Parts of these graph are also in the tutorial pages.

<br />

![CaptureAl.png](https://raw.githubusercontent.com/scrameri/CaptureAl/master/CaptureAl.png)

## Pipeline authors
Simon Crameri & Stefan Zoller
