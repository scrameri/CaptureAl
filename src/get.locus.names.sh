#!/bin/bash

## Usage: get.locus.names.sh <refseqs>

## Value: extracts all fasta headers (w/o leading >)

## Get arguments
refseqs=$1

## Extract
grep '^>' $refseqs | cut -d '>' -f2 > locus.names
