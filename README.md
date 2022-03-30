# CaptureAl
Bioinformatic pipeline streamlining analysis of **target enrichment sequencing** data.

## Context
This repository accompanies a manuscript in review (Crameri et al. submitted), available [here](https://www.biorxiv.org/content/10.1101/2021.12.07.471551v1).

## Before you start
1) Download/Clone the [CaptureAl repository](https://github.com/scrameri/CaptureAl) 
2) Follow the [installation instructions](https://github.com/scrameri/CaptureAl/blob/master/Install.md).
3) Determine how to implement [parallel computing](https://github.com/scrameri/CaptureAl/blob/master/Parallelize.md) in your computing environment (*GNU parallel* or *bsub*)

## Get started
Once the installation is complete, follow the [tutorial](https://github.com/scrameri/CaptureAl/blob/master/tutorial/) step by step. You can also jump-start to **STEP 1** ([parallel](https://github.com/scrameri/CaptureAl/blob/master/tutorial/parallel/Step1_Read_Mapping.md) or [bsub](https://github.com/scrameri/CaptureAl/blob/master/tutorial/bsub/Step1_Read_Mapping.md)) if your reads are already quality-trimmed.

## Post-process Workflows for
- [Phylogenetics](Phylogenetics)
- [Population Genetics](PopulationGenetics)

## Pipeline overview
The pipeline is divided into **seven steps**, shown below. For each step, there is a limited number of scripts that have to be executed, and these scripts often execute other scripts located in the repository. Parts of these graph are also in the tutorial sites.

<br />

![CaptureAl.png](https://raw.githubusercontent.com/scrameri/CaptureAl/master/CaptureAl.png)

## Pipeline authors
Simon Crameri & Stefan Zoller
