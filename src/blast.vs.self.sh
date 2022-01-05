#!/bin/bash

## Usage: blast.vs.self.sh <file.fasta> 

## Needs: blast.fasta.seqs.sh, filter.blast.vs.self.R

## Get arguments
fasta=$1

## Additional arguments
fbase=$(basename $fasta .fasta)

## Move any existing output files
if [ -f ${fbase}.vs.self.blast ] ; then mv ${fbase}.vs.self.blast ${fbase}.vs.self.bak.blast ; fi
if [ -f ${fbase}.vs.self.blast.filtered ] ; then mv ${fbase}.vs.self.blast.filtered ${fbase}.vs.self.bak.blast.filtered ; fi

## BLAST
blast.fasta.seqs.sh -q $fasta -d $fasta -e 1e-04
mv ${fbase}.on.${fbase}.blast ${fbase}.vs.self.blast

## Filter
echo "### Filtering BLAST results ###"
echo
filter.blast.vs.self.R ${fbase}.vs.self.blast
