#!/bin/bash

## Usage: align.multifastas.parallel.sh -d <directory with multifasta files> -m <alignment model> -t <threads>

## Needs:
# mafft

## Define arguments
while getopts d:m:t: opts
do
    case "${opts}"
    in
      d) indir=${OPTARG};;
      m) model=${OPTARG};;
      t) threads=${OPTARG};;
    esac
done

## Check arguments
if [ ! $indir ] ; then echo "input directory not provided (-d option)" ; exit 0 ; fi
if [ ! -d $indir ] ; then echo "input directory $indir not found" ; exit 0 ; fi
if [ ! $threads ] ; then echo "number of threads not specified (-t option), using 4 threads" ; threads=4 ; fi

## Available alignment options (from mafft man page) (-m option)
# - localpair
#   probably most accurate; recommended for <200 sequences; iterative refinement method incorporating local pairwise alignment information
# - globalpair
#   suitable for sequences of similar lengths; recommended for <200 sequences; iterative refinement method incorporating global pairwise alignment information
# - genafpair
#   suitable for sequences containing large unalignable regions; recommended for <200 sequences. The --ep 0 option is recommended to allow large gaps
if [ ! $model ] ; then echo "alignment model not specified (-m option), setting to 'localpair'" ; model='localpair' ; fi
if [ $model != 'localpair' ] && [ $model != 'globalpair' ] && [ $model != 'genafpair' ] ; then echo "alignment model <$model> not supported (choose between 'localpair', 'globalpair' or 'genafpair'" ; exit 0 ; fi

## Additional arguments
suffix=".fasta" # will align all multifasta files with this infix
infix=".aln" # infix of output files
outdir=$(echo $indir | sed -e 's/^multifasta/mafft/g')

## Export
export model=$model
export suffix=$suffix
export infix=$infix
export indir=$indir
export outdir=$outdir

## Create output directory
if [ -d $outdir ] ; then mv $outdir $outdir.bak ; echo ; echo "$outdir exists, renaming it to $outdir.bak" ; echo  ; fi
mkdir $outdir

## Define
do_align()
{
	fas=$1
	fbase=$(basename $fas $suffix)
	echo $fbase
	
	aln="${outdir}/${fbase}${infix}${suffix}"

	# run mafft
	mafft --${model} --maxiterate 1000 --adjustdirection ${fas} > ${aln}

	# post process mafft output (reverse mafft fasta header prefix '_R_' of reverse complements, unwrap)
	sed -i 's/^>_R_/>/g' ${aln}
	awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' ${aln} > ${aln}.tmp
	mv ${aln}.tmp ${aln}

}

## Execute
export -f do_align
ls -1 $indir/*${suffix} | parallel  -j $threads do_align 1> ${outdir}.log 2> ${outdir}.mafft

## Finish
echo
echo "All loci aligned!"
echo
