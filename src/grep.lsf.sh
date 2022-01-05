#!/bin/bash

## Usage: grep.lsf.sh <lsf.oXXXXX>

## Details
# looks for bsub jobs that got killed due to TERM_MEMLIMIT or TERM_RUNLIMIT errors etc.

## Get arguments
lsf=$1

grep -E '^Sender:|^TERM_' $lsf -A15 | grep '^TERM_' -B16 -A15 | grep -Ev '^# LSBATCH:|^#!|^Your job looked like:|^---' | grep -E '^Subject:' -A 30 --color

