#!/bin/bash

## Usage: align.contigs.of.overlapping.loci.sh -l <.blast.filtered.list> -c <multifasta dir> -m <alignment model> -x <seqname prefix> -y <seqname suffix> -t <threads>

## Author: simon.crameri@env.ethz.ch, Apr 2019

## Needs:
# mafft
# awk
# get.overlap.consensus.R

## Define arguments
while getopts l:c:m:x:y:t: opts
do
        case "${opts}"
        in
                l) overlaps=${OPTARG};;
                c) multifastas=${OPTARG};;
                m) model=${OPTARG};;
                x) prefix=${OPTARG};;
                y) suffix=${OPTARG};;
                t) threads=${OPTARG};;
         esac
done

## Check arguments
if [ ! $overlaps ] ; then echo "list of overlapping loci not specified (-l option)" ; exit 0 ; fi
if [ ! $multifastas ] ; then echo "multifastas directory not specified (-c option)" ; exit 0 ; fi

if [ ! -f $overlaps ] ; then echo "list of overlapping loci not <$overlaps> not found" ; exit 0 ; fi
if [ ! -d $multifastas ] ; then echo "multifastas directory <$multifastas> not found" ; exit 0 ; fi

if [ ! $prefix ] ; then echo "prefix of sequence in $overlaps (but not in multifasta names) not specified (-x option), setting to ''" ; prefix='' ; fi
if [ ! $suffix ] ; then echo "suffix of sequence in $overlaps (but not in multifasta names) not specified (-y option), setting to ''" ; suffix='' ; fi

if [ ! $threads ] ; then echo "number of threads not specified, setting to 4" ; threads=4 ; fi
if [ $threads -gt 30 ] ; then echo "number of threads must not exceed 30, stopping" ; exit 0 ; fi

## Additional arguments
noverlaps=$(wc -l < $overlaps)
export outdir="mafft.overlapping.${noverlaps}" # name of output directory
export msuffix=".all.fasta"            # suffix of multifasta files in <multifastas> directory
export ssuffix=".all.fasta"            # suffix of multifasta files in <outdir>/Multicontigfastas directory
export infix=".merged"                 # infix that denotes merged multicontigfastas, alignments, logs and pdfs

## Available alignment options (from mafft man page) (-m option)
# - localpair
#   probably most accurate; recommended for <200 sequences; iterative refinement method incorporating local pairwise alignment information
# - globalpair
#   suitable for sequences of similar lengths; recommended for <200 sequences; iterative refinement method incorporating global pairwise alignment information
# - genafpair
#   suitable for sequences containing large unalignable regions; recommended for <200 sequences. The --ep 0 option is recommended to allow large gaps
if [ ! $model ] ; then echo "alignment model not specified (-m option), setting to 'localpair'" ; model='localpair' ; fi
if [ $model != 'localpair' ] && [ $model != 'globalpair' ] && [ $model != 'genafpair' ] ; then echo "alignment model <$model> not supported (choose between 'localpair', 'globalpair' or 'genafpair'" ; exit 0 ; fi

## Export arguments
export overlaps=$overlaps
export multifastas=$multifastas
export model=$model
export prefix=$prefix
export suffix=$suffix

## Create output directory
if [ -d $outdir ] ; then echo "output directory <$outdir> already exists, moving to <$outdir.bak>" ; mv $outdir $outdir.bak ; fi
mkdir ${outdir}
mkdir ${outdir}/Multicontigfastas
mkdir ${outdir}/AlignmentsN
if [ -f ${outdir}.err ] ; then /bin/rm -f >> ${outdir}.err ; fi

## Define function
do_align_overlapping ()
{
	hit=$1
	
	# test if more than 5 loci are involved (this version cannot currently handle more overlapping loci)
	toomany=$(echo $hit | cut -f6 -d' ')
	toofew=$(echo $hit | cut -f2 -d' ')
	
	third=$(echo $hit | cut -f3 -d' ')
	forth=$(echo $hit | cut -f4 -d' ')
	fifth=$(echo $hit | cut -f5 -d' ')
	
	if [ -z $toofew ] ; then echo "ERROR: less than 2 loci involved, stopping" >> ${outdir}.err ; exit 0 ; fi
	if [ ! -z $toomany ] ; then echo "WARNING: more than 5 loci involved, only aligning first 5 loci in list" >> ${outdir}.err ; fi
	
	# get contig names
	c1=$(echo $hit | cut -f1 -d' ') ; c2=$(echo $hit | cut -f2 -d' ') ; c3=$(echo $hit | cut -f3 -d' ') ; c4=$(echo $hit | cut -f4 -d' ') ; c5=$(echo $hit | cut -f5 -d' ')
	p1=$(echo $c1 | sed -e "s/^${prefix}//g" | sed -e "s/${suffix}$//g") ; p2=$(echo $c2 | sed -e "s/^${prefix}//g" | sed -e "s/${suffix}$//g") ; p3=$(echo $c3 | sed -e "s/^${prefix}//g" | sed -e "s/${suffix}$//g") ; p4=$(echo $c4 | sed -e "s/^${prefix}//g" | sed -e "s/${suffix}$//g") ; p5=$(echo $c5 | sed -e "s/^${prefix}//g" | sed -e "s/${suffix}$//g")
	
	# set new names for output
	string="${p1}${infix}"              # only the name of the first locus is retained but with an infix to distinguish it from the non-merged locus 
	fas="${outdir}/${string}${ssuffix}" # name of multifasta with all overlapping contigs of all taxa to be aligned (mafft input)
	
	# get overlapping multifastas ($fas)
	mf1="${multifastas}/${p1}${msuffix}"
	mf2="${multifastas}/${p2}${msuffix}"
	if [ ! -f $mf1 ] || [ ! -f $mf2 ]
	then 
		echo "ERROR: $mf1 not found, check prefix (-x option) and suffix (-y option)" >> ${outdir}.err
		exit 0
	fi
	
	if [ ! -z $fifth ]
	then
		#string="${p1}_${p2}_${p3}_${p4}_${p5}" # string gets too long. Need to remove $p5 manually from $trimmeddir
		#echo "${p1}_${p2}_${p3}_${p4}_${p5} too long, naming is ${p1}_${p2}_${p3}_${p4}, need to manually replace ${p5}" >> ${outdir}.warn
		
		echo "$p1 <-> $p2 <-> $p3 <-> $p4 <-> $p5" # >> ${outdir}.log to keep track of the merging
		
		#string="${p1}_${p2}_${p3}_${p4}"
		#fas="${outdir}/${string}${ssuffix}"
		
		mf3="${multifastas}/${p3}${msuffix}"
		mf4="${multifastas}/${p4}${msuffix}"
		mf5="${multifastas}/${p5}${msuffix}"
		cat $mf1 $mf2 $mf3 $mf4 $mf5 > $fas
	else
		if [ ! -z $forth ]
		then
			echo "$p1 <-> $p2 <-> $p3 <-> $p4" # >> ${outdir}.log to keep track of the merging
		
			#string="${p1}_${p2}_${p3}_${p4}"
			#fas="${outdir}/${string}${ssuffix}"
			
			mf3="${multifastas}/${p3}${msuffix}"
			mf4="${multifastas}/${p4}${msuffix}"
			cat $mf1 $mf2 $mf3 $mf4 > $fas
		else
			if [ ! -z $third ]
			then
				echo "$p1 <-> $p2 <-> $p3" # >> ${outdir}.log to keep track of the merging
		
				#string="${p1}_${p2}_${p3}"
				#fas="${outdir}/${string}${ssuffix}"
				
				mf3="${multifastas}/${p3}${msuffix}"
				cat $mf1 $mf2 $mf3 > $fas
			else
				echo "$p1 <-> $p2" # >> ${outdir}.log to keep track of the merging
			
				#string="${p1}_${p2}"
				#fas="${outdir}/${string}${ssuffix}"
				
				cat $mf1 $mf2 > $fas
			fi
		fi
	fi
	
	# define output file names ($aln, $mfas $mpdf)
	aln="${outdir}/${string}.aln.fasta"      # name of big alignment (mafft output with all taxa and all contigs of overlapping loci): nseq = ntaxa*noverlappingloci
	mfas="${outdir}/${string}.all.aln.fasta" # name of consensus alignment outputted by get.overlap.consensus.R (all taxa, consensus of overlapping loci): nseq = ntaxa
	mpdf="${outdir}/${string}.all.aln.pdf"   # name of visualization outputted by get.overlap.consensus.R
	
	# sort input order to consecutively represent all loci
	#reorder.fasta.seqs.R ${fas} 'length' # in several tests, no improvement after changing input order to sorted or decreasing length and / or changing mafft alignment algorithm and gap penalty.
		
	# run mafft on 2*N or 3*N or 4*N sequences (try different alogrithms)
 	mafft --${model} --maxiterate 1000 --adjustdirection ${fas} > ${aln}
 	
	# post process mafft output (reverse mafft fasta header prefix '_R_' of reverse complements, unwrap)
	sed -i 's/^>_R_/>/g' ${aln}
	awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' ${aln} > ${aln}.tmp
	mv ${aln}.tmp ${aln}

	# get consensus for each taxon in $aln according to majority rule 
	# (if there is > 1 base for a taxon at a site, take the base that is more frequent at that site)
	get.overlap.consensus.R ${aln} # produces ${mfas} and ${mpdf}
	
	# clean up
	mv ${fas} ${outdir}/Multicontigfastas
	mv ${aln} ${outdir}/AlignmentsN

}

## Execute function
export -f do_align_overlapping
cat $overlaps | parallel -j $threads do_align_overlapping 1> ${outdir}.log 2> ${outdir}.mafft

## Finish
echo
echo "All overlapping loci aligned!"
echo

