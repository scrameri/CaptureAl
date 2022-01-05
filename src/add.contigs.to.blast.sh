#!/bin/bash

## Usage: add.contigs.to.blast.sh -s <name of added sample> -f <.fasta with sequences to be added (headers must be unique)> -d <folder with contigs to be blasted> -c <name of added contigs file> 

## Define arguments 
while getopts s:f:d:c: opts
do
        case "${opts}"
		in
 			s) addname=${OPTARG};;
 			f) fasta=${OPTARG};;
 			d) blastdir=${OPTARG};;
  			c) contigspath=${OPTARG};;
		esac
done

## Check arguments
if [ ! $addname ] ; then echo "added taxon name not specified (-s option)" ; exit 0 ; fi
if [ ! $fasta ] ; then echo "added sequences (.fasta format) not specified (-f option)" ; exit 0 ; fi
if [ ! $blastdir ] ; then echo "path to blast directory not specified (-d option), setting to 'blast.contigs'" ; blastdir='blast.contigs' ; fi
if [ ! $contigspath ] ; then echo "name of added contigs file not specified (-c option)" ; exit 0 ; fi

if [ ! -f $fasta ] ; then echo "added sequences <$fasta> not found" ; exit 0 ; fi
if [ ! -d $blastdir ] ; then echo "blast directory <$blastdir> not found" ; exit 0 ; fi

## Create new sample dir in $blastdir
if [ -d ${blastdir}/${addname} ] ; then echo "directory <${blastdir}/${addname}> already exists" ; exit 0 ; fi
mkdir ${blastdir}/${addname}

## Copy
contigsname=$(basename $contigspath)
cp $fasta ${blastdir}/${addname}/${contigsname}

## Finish 
echo
echo "Done!"
echo
