#!/bin/bash

## Usage: filter.contigs.parallel.sh -s <sample file> -d <directory with exonerate results> -f <exonerate stats file> -m <minimum normalized best exonerate score> -t <number of threads>

## Define arguments
while getopts s:d:f:m:t: opts
do
        case "${opts}"
        in
                s) sfile=${OPTARG};;
                d) dir=${OPTARG};;
                f) stats=${OPTARG};;
                m) minnormbestscore=${OPTARG};;
                t) threads=${OPTARG};;
    	 esac
done

## Check arguments
if [ ! $sfile ] ; then echo "sample file (-s option) not specified, stopping!" ; exit 0 ; fi
if [ ! $dir ] ; then echo "directory with exonerate results (-d option) not specified, stopping!" ; exit 0 ; fi
if [ ! $stats ] ; then echo "exonerate stats file (-f option) not specified, stopping!" ; exit 0 ; fi
if [ ! $minnormbestscore ] ; then echo "minimum normalized best exonerate score not specified (-m option) setting to 2!" ; minnormbestscore=2 ; fi

if [ ! -f $sfile ] ; then echo "file <${sfile}> does not exist, stopping!" ; exit 0 ; fi
if [ ! -d $dir ] ; then echo "directory <${dir}> does not exist, stopping!" ; exit 0 ; fi
if [ ! -f $stats ] ; then echo "file <${stats}> does not exist, stopping!" ; exit 0 ; fi

if [ ! $threads ] ; then echo "number of threads not specified, setting to 4" ; threads=4 ; fi
if [ $threads -gt 30 ] ; then echo "number of threads must not exceed 30, stopping!" ; exit 0 ; fi

## Additional arguments

## Export arguments
export dir=$dir
export stats=$stats
export minnormbestscore=$minnormbestscore

## Generate log file
echo ""                                                             > filter.contigs.log
echo "=== filter.contigs.parallel.sh LOG ==="                       >> filter.contigs.log
echo ""                                                             >> filter.contigs.log
echo "sample file:                                $sfile"           >> filter.contigs.log
echo "input directory:                            $dir"             >> filter.contigs.log
echo "exonerate stats file:                       $stats"           >> filter.contigs.log
echo "minimum normalized best exonerate score:    $minnormbestscore">> filter.contigs.log
echo ""                                                             >> filter.contigs.log
echo "Number of threads:                          $threads"         >> filter.contigs.log
echo "Starting time:                              $(zdump MEC)"     >> filter.contigs.log
echo ""                                                             >> filter.contigs.log
echo "=========================="                                   >> filter.contigs.log

## Define 
do_filter()
{
	sample=$1
	sbase=$(basename $sample)
	echo $sbase
	
	args="${dir}/${sample} ${stats} ${minnormbestscore}"
	
	filter.contigs.R $args >> filter.contigs.log
}

## Execute
export -f do_filter
cat ${sfile} | parallel -j $threads do_filter 2> filter.contigs.err

## Finish
nerr=$(wc -l < filter.contigs.err)
if [ $nerr -eq 0 ] ; then /bin/rm -f filter.contigs.err ; fi
echo
echo "=========================="                                   >> filter.contigs.log
echo "Finish time:                                $(zdump MEC)"     >> filter.contigs.log
echo "All samples processed!"
echo