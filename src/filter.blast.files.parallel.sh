#!/bin/bash

## Usage: filter.blast.files.parallel.sh -s <sample file> -r <reference .fasta> -d <dir with blast results> -c <contigs fasta>
#										 -e <evalue> 
#                                        -p <min.perc.alignment.length> -a <min.alignment.length> 
#                                        -n <min.subject.length> -x <max.subject.length> -i <min.perc.identity> 
#                                        -t <number of threads>

## Needs: filter-blast-files.R

## Define arguments
while getopts s:r:d:c:e:p:a:n:x:i:t: opts
do
        case "${opts}"
		in
				s) sfile=${OPTARG};;
				r) ref=${OPTARG};;
				d) blastdir=${OPTARG};;
				c) contigs=${OPTARG};;
				e) evalue=${OPTARG};;
				p) percaln=${OPTARG};;
				a) minaln=${OPTARG};;
				n) minlen=${OPTARG};;
				x) maxlen=${OPTARG};;
				i) percid=${OPTARG};;
				t) threads=${OPTARG};;
		
		esac
done

## Check arguments and set defaults
if [ ! $sfile ] ; then echo "sample file not specified (-s option), stopping" ; exit 0 ; fi
if [ ! $ref ] ; then echo "reference fasta not specified (-r option), stopping" ; exit 0 ; fi
if [ ! $blastdir ] ; then echo "directory with blast results not specified (-d option), setting to 'blast.contigs'" ; blastdir=blast.contigs ; fi
if [ ! $contigs ] ; then echo "name of contigs fasta in sample subdirs not specified (-c option), stopping" ; exit 0 ; fi

if [ ! -f $sfile ] ; then echo "sample file <$sfile> not found, stopping" ; exit 0 ; fi
if [ ! -f $ref ] ; then echo "reference fasta <$ref> not found, stopping" ; exit 0 ; fi
if [ ! -d $blastdir ] ; then echo "directory with blast results <$blastdir> not found, stopping" ; exit 0 ; fi

if [ ! $evalue ] ; then echo "evalue threshold (-e option) not set, looking for .blast files generated using evalue threshold of 1e-04" ; evalue=1e-04 ; fi
if [ ! $percaln ] ; then echo "min.perc.alignment.length (-p option) not set, setting to 0" ; percaln=0 ; fi
if [ ! $minaln ] ; then echo "min.alignment.length (-a option) not set, setting to 80" ; minaln=80 ; fi
if [ ! $minlen ] ; then echo "min.subject.length (-n option) not set, setting to 100" ; minlen=100 ; fi
if [ ! $maxlen ] ; then echo "max.subject.length (-x option) not set, setting to 10000" ; maxlen=10000 ; fi
if [ ! $percid ] ; then echo "min.perc.identity (-i option) not set, setting to 80" ; percid=80 ; fi

if [ ! $threads ] ; then echo "number of threads (-t option) not set, setting to 4" ; threads=4 ; fi
if [ $threads -gt 30 ] ; then echo "number of threads must not exceed 30, stopping" ; exit 0 ; fi

## Additional arguments
fstring="${evalue}.${percaln}.${minaln}.${minlen}.${maxlen}.${percid}"  # string used to name output of filter-blast-files.R
outfile="blastres.${fstring}.txt"                                       # output file name

## Export parameters
export ref=$ref
export refbase=$(basename $ref .fasta)
export blastdir=$blastdir
export contigbase=$(basename $contigs .fasta)                           # basename of contig fasta (w/o .fasta extension)

export evalue=$evalue
export percaln=$percaln
export minaln=$minaln
export minlen=$minlen
export maxlen=$maxlen
export percid=$percid

echo
echo "looking for BLAST hits in ${blastdir}/{sample}/${refbase}.on.{sample}.${contigbase}.${evalue}.blast"
echo "minimum percentage of query length in alignment length set to $percaln"
echo "minimum alignment length set to $minaln"
echo "minimum subject length set to $minlen"
echo "maximum subject length set to $maxlen"
echo "minimum percentage of identity between aligned portion of subject and query set to $percid"
echo 

## Define filtering function
doFilterBlast()
{
	sample=$1
	echo $sample
	file="${blastdir}/${sample}/${refbase}.on.${sample}.${contigbase}.${evalue}.blast"
	if [ ! -f $file ] ; then echo "$file not found, stopping" >> filter.blast.err ; exit 0 ; fi

	# Required argument
	# f                           character         path to .blast results file. Must have been generated using blast.fasta.seqs.sh or blast.contigs.parallel.sh

	# Optional arguments (if provided, all must be given, in this order)
	# min.perc.alignment.length   numeric [0, 100]  minimum percentage of query.length in alignment.length [DEFAULT: 0]
	# min.alignment.length        numeric [1, ]     minimum length of aligned portion of subject and query [DEFAULT: 80]
	# min.subject.length          numeric [1, ]     minimum subject length [DEFAULT: 100]
	# max.subject.length          numeric [1, ]     maximum subject length [DEFAULT: 10000]
	# min.perc.identity           numeric [0, 100]  minimum percentage of identity between aligned portion of subject and query [DEFAULT: 80]

	args="$file $percaln $minaln $minlen $maxlen $percid"
	filter.blast.files.R $args >/dev/null

}

## Execute filtering function
export -f doFilterBlast
cat ${sfile} | parallel -j $threads doFilterBlast

## Concatenate results
# header
cat ${blastdir}/$(head -n1 $sfile)/${refbase}.on.$(head -n1 $sfile).${contigbase}.${fstring}.allhits.txt | head -n1 > ${outfile}

# data w/o header
echo
echo "concatenating results..."
for sample in $(cat $sfile)
do
	allhits="${refbase}.on.${sample}.${contigbase}.${fstring}.allhits.txt"
	#passedhits="${refbase}.on.${sample}.${contigbase}.${fstring}.passedhits.txt"
	#besthits="${refbase}.on.${sample}.${contigbase}.${fstring}.besthits.txt"
	cat ${blastdir}/${sample}/${allhits} | tail -n +2 
done >> ${outfile}
echo "done"

## Finish
echo ""
echo "All samples processed!"
echo ""
