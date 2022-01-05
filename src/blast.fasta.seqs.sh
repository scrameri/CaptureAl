#!/bin/bash

## Usage: blast.fasta.seqs.sh -q <query.fasta> -d <db.fasta or nt> -e <max_evalue> -m <max_target_seqs> -t <num_threads>
## Required argument: -q
## Optional arguments: -d, -e, -m, -t

## Define arguments
while getopts q:d:e:m:t: opts
do
        case "${opts}"
        in
                q) q=${OPTARG};;
                d) d=${OPTARG};;
                e) e=${OPTARG};;
                m) m=${OPTARG};;
                t) t=${OPTARG};;
        esac
done

## Default parameters
if [ ! $q ] ; then echo "query sequence(s) [.fasta format] required (-q option)" ; exit 0 ; fi
if [ ! -e $q ] ; then echo "$q not found, stopping." ; exit 0 ; fi
if [ ! $d ] ; then d='nt' ; else if [ ! -e $d ] ; then echo "$d not found, stopping." ; exit 0 ; fi ; fi
if [ ! $e ] ; then e=10 ; fi
if [ ! $m ] ; then m=10000 ; fi
if [ ! $t ] ; then t=1 ; fi

echo ; echo "### BLAST .fasta sequences ###" ; echo 
echo "$q used as BLAST queries"
echo "$d used as BLAST subjects (database to BLAST against)"
echo "evalue threshold set to $e"
echo "max_target_seqs set to $m"
echo "$t threads used"
qname=$(basename ${q} .fasta)
dname=$(basename ${d} .fasta)

## If $q is $d, create BLAST database of $d
if [ -e $d ]
then
	makeblastdb -in ${d} -parse_seqids -dbtype nucl >/dev/null
fi

## BLAST
blastn -db ${d} -query ${q} -evalue ${e} -max_target_seqs ${m} -num_threads ${t} -out ${qname}.tmp -outfmt "7 qseqid qlen sseqid slen length pident qstart qend sstart send evalue bitscore score stitle"

## Sort BLAST output by query.id, subject.id and alignment.length
grep -v '^#' ${qname}.tmp | sort -k1,1 -k3,3 -k5,5nr > ${qname}.tmp2

## Create header and output
echo "$(head -n1 ${qname}.tmp) $(grep '# Fields:' ${qname}.tmp | head -n1 | sed -e 's/% identity/perc.identity/g' | sed -e 's/# //g')" > ${qname}.on.${dname}.blast
cat ${qname}.tmp2 >> ${qname}.on.${dname}.blast

## Clean up
/bin/rm -f ${qname}.tmp*
/bin/rm -f ${d}.n*
echo ; echo "Done!" ; echo
