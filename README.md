# CaptureAl
Bioinformatic pipeline streamlining analysis of **target enrichment sequencing** data. This repository accompanies a manuscript in review (Crameri et al. submitted), available [here](https://www.biorxiv.org/content/10.1101/2021.12.07.471551v1).

### Download the CaptureAl repository
To use this pipeline, clone the [CaptureAl repository](https://github.com/scrameri/CaptureAl) and follow the tutorial provided below. For each step, there is a limited number of scripts that have to be executed, and these scripts often execute other scripts located in the repository.

### Before you start
1) Follow the [installation instructions](https://github.com/scrameri/CaptureAl/blob/master/Install.md).
2) Determine how to implement [parallel computing](https://github.com/scrameri/CaptureAl/blob/master/Parallelize.md) in your computing environment (*GNU parallel* or *bsub*)

### Get started
Once the installation is complete, follow the tutorial step by step by viewing the numbered `.md` files, starting with **Sequence Quality Control** ([parallel](https://github.com/scrameri/CaptureAl/blob/master/tutorial/parallel/Step0.1_Sequence_Quality_Control.md) or [bsub](https://github.com/scrameri/CaptureAl/blob/master/tutorial/bsub/Step0.1_Sequence_Quality_Control.md)). Each step contains a link to the previous and next step.

### Graphical overview
The pipeline is divided into **seven steps**, shown below. Parts of these graph are also in the tutorial sites.

![CaptureAl.png](https://raw.githubusercontent.com/scrameri/CaptureAl/master/CaptureAl.png)
