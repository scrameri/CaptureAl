#!/bin/bash

## Usage: align.multifasta.hyper.sh -d <directory with multifasta files> -L <path to single multifasta file> -a <mafft or muscle> -m <mafft alignment model> -t <threads: parallelize over loci> -h <hyperthreads used by mafft>

## Needs:
# mafft

## Define arguments
while getopts d:L:a:m:t:h: opts
do
    case "${opts}"
    in
      d) indir=${OPTARG};;
      L) infile=${OPTARG};;
      a) method=${OPTARG};;
      m) model=${OPTARG};;
      t) threads=${OPTARG};;
      h) hthreads=${OPTARG};;
    esac
done

## Check arguments
if [ ! $indir ] ; then echo "input directory not provided (-d option)" ; exit 0 ; fi
if [ ! -d $indir ] ; then echo "input directory <$indir> not found" ; exit 0 ; fi
if [ ! $infile ] ; then echo "input multifasta not provided (-L option), aligning all multifasta in <$indir>" ; mode='parallel' ; else echo "input multifasta provided (-L option), aligning single multifasta" ; mode='single' ; fi
if [ $mode == 'single' ] ; then if [ ! -f $infile ] ; then echo "input multifasta <$infile> not found" ; exit 0 ; fi ; fi 

if [ ! $method ] ; then echo "alignment method not specified (-a option), using -a 'muscle'" ; method='muscle' ; fi
if [ $method != 'mafft' ] && [ $method != 'muscle' ] ; then echo "alignment method <$method> not supported (choose between 'mafft' and 'muscle'" ; exit 0 ; fi

if [ ! $threads ] ; then echo "number of threads not specified (-t option), using 1 threads" ; threads=1 ; fi
if [ $method == 'mafft' ] && [ ! $hthreads ] ; then echo "number of hyperthreads for method <$method> not specified (-h option), using --thread 1" ; hthreads=1 ; fi

## Available alignment options (from mafft man page) (-m option)
# - localpair
#   probably most accurate; recommended for <200 sequences; iterative refinement method incorporating local pairwise alignment information
# - globalpair
#   suitable for sequences of similar lengths; recommended for <200 sequences; iterative refinement method incorporating global pairwise alignment information
# - genafpair
#   suitable for sequences containing large unalignable regions; recommended for <200 sequences. The --ep 0 option is recommended to allow large gaps
if [ $method == 'mafft' ] && [ ! $model ] ; then echo "alignment model for method <$method> not specified (-m option), using -m 'localpair'" ; model='localpair' ; fi
if [ $method == 'mafft' ] && [ $model != 'localpair' ] && [ $model != 'globalpair' ] && [ $model != 'genafpair' ] ; then echo "alignment model <$model> not supported for method <$method> (choose between 'localpair', 'globalpair' or 'genafpair'" ; exit 0 ; fi

## Additional arguments
suffix=".fasta" # will align all multifasta files with this infix
infix=".aln" # infix of output files
outdir=$(echo $indir | sed -e "s/^multifasta/$method/g")

## Export
export method=$method
export model=$model
export suffix=$suffix
export infix=$infix
export indir=$indir
export outdir=$outdir
export hthreads=$hthreads

## Create output directory
if [ $mode == 'parallel' ] && [ -d $outdir ] ; then mv $outdir $outdir.bak ; echo ; echo "$outdir exists, renaming it to $outdir.bak" ; echo  ; fi
if [ ! -d $outdir ] ; then mkdir $outdir ; fi

## Define
do_align()
{
	fas=$1
	fbase=$(basename $fas $suffix)
	echo $fbase
	
	aln="${outdir}/${fbase}${infix}${suffix}"

	if [ $method == 'mafft' ]
	then
		# run mafft
		mafft --${model} --maxiterate 1000 --adjustdirection --thread ${hthreads} ${fas} > ${aln}
	fi

	if [ $method == 'muscle' ]
	then
		# run muscle (default -maxiters 16)
		# consider using -stable (keeps input order of sequences rather than grouping them by similarity by default -group)
		# consider using -hours 4 (will allow 4 hours of alignment at most) if time is an issue
		muscle -in ${fas} -out ${aln}
	
		# run muscle (high speed)	
		#muscle -in ${fas} -out ${aln} -maxiters 2
	fi


	# post process mafft output (reverse mafft fasta header prefix '_R_' of reverse complements, unwrap)
	sed -i 's/^>_R_/>/g' ${aln}
	awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' ${aln} > ${aln}.tmp
	mv ${aln}.tmp ${aln}

}

## Execute
export -f do_align

if [ $mode == 'single' ]
then
	# align single multifasta (submit X jobs for each alignment)
	do_align $infile 1>> ${outdir}.log 2>> ${outdir}.${method}
fi

if [ $mode == 'parallel' ]
then
	# parallel alignment (submit 1 job)
	ls -1 $indir/*${suffix} | parallel  -j $threads do_align 1> ${outdir}.log 2> ${outdir}.${method}
	
	echo
	echo "All loci aligned!"
	echo
fi
