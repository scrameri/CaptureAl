# Dependencies

## Software
This repository provides scripts (see [src](https://github.com/scrameri/CaptureAl/tree/master/src) directory) that use third-party software. These sofware tools must be executable from your computing environment:
- R and corresponding Rscript (tested on version 3.6 and 4.0)
- python (tested on version 2.7.11)
- java (tested on version 1.8.0_73)
- perl (tested on version 5.16.3)
- fastqc (tested on version 0.11.4)
- trimmomatic (tested on version 0.35)
- bwa (tested on version 0.7.17)
- sambamba (tested on version 0.8.0)
- samtools (tested on version 1.2)
- bedtools (tested on version 2.28.0)
- picard-tools (tested on version 2.23.8)
- spades or dipspades (tested on version 3.14.1)
- exonerate (tested on version 2.4.0)
- freebayes (tested on version 1.3.4)

## R packages
Make sure that the following R packages are installed in your computing environment and available through 'Rscript':

- ape       # used to manage FASTA files
- ggplot2   # used for graphics
- hexbin    # used for alignment assessment
