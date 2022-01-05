#!/bin/bash

## Usage: rename.overlapping.consensusfiles.sh <dirbase> <consensus.fasta>

## Get arguments
dirbase=$(basename $1)
consbase=$(basename $2 .fasta)

## Additional arguments
suffix=".overlap"

## Move files
# dirbase
if [ -f ${dirbase}.fasta.check ] ; then mv ${dirbase}.fasta.check ${dirbase}${suffix}.fasta.check ; fi
if [ -f ${dirbase}.fasta.headers ] ; then mv ${dirbase}.fasta.headers ${dirbase}${suffix}.fasta.headers ; fi
if [ -f ${dirbase}.lengths ] ; then mv ${dirbase}.lengths ${dirbase}${suffix}.lengths ; fi
if [ -f ${dirbase}.lengths.pdf ] ; then mv ${dirbase}.lengths.pdf ${dirbase}${suffix}.lengths.pdf ; fi

# consbase
if [ -f ${consbase}.lengths ] ; then mv ${consbase}.lengths ${consbase}${suffix}.lengths ; fi
if [ -f ${consbase}.lengths.pdf ] ; then mv ${consbase}.lengths.pdf ${consbase}${suffix}.pdf ; fi

if [ -d ${consbase} ] ; then mv ${consbase} ${consbase}${suffix} ; fi
if [ -d ${consbase}.viz ] ; then mv ${consbase}.viz ${consbase}${suffix}.viz ; fi
if [ -d ${consbase}.logs ] ; then mv ${consbase}.logs ${consbase}${suffix}.logs ; fi
if [ -f ${consbase}.fasta ] ; then mv ${consbase}.fasta ${consbase}${suffix}.fasta ; fi
if [ -f ${consbase}.vs.self.blast ] ; then mv ${consbase}.vs.self.blast ${consbase}.vs.self${suffix}.blast ; fi
if [ -f ${consbase}.vs.self.blast.filtered ] ; then mv ${consbase}.vs.self.blast.filtered ${consbase}.vs.self${suffix}.blast.filtered ; fi
if [ -f ${consbase}.vs.self.blast.filtered.overlap ] ; then mv ${consbase}.vs.self.blast.filtered.overlap ${consbase}.vs.self${suffix}.blast.filtered.overlap ; fi
if [ -f ${consbase}.vs.self.blast.filtered.list ] ; then mv ${consbase}.vs.self.blast.filtered.list ${consbase}.vs.self${suffix}.blast.filtered.list ; fi
