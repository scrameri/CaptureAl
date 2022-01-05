#!/bin/bash

## USAGE: subset.phylip.sh -p <PHYLIP.phy> -f <indfile with taxon names to be kept or dropped> -d <flag: if turned on, will drop instead of keep>

## VALUE: a subsetted PHYLIP file named <phybase_NewNumberOfTaxa.phy>

## DETAILS: 
# - Only meant to be applied on single-line PHYLIP files:
#   checks whether the number of lines in PHYLIP equals the number of taxa + 1, and exists if not
#
# - Exact taxon string matching:
#   matching of taxa in INDFILE to taxa in PHYLIP includes non-word constituent characters in taxon names (such as '-')
#   under the condition that the taxon string in INDFILE is followed by at least one whitespace in PHYLIP (see EXAMPLE)
#
# - Output control:
#   checks whether the number of kept or dropped taxa actually corresponds to INDFILE

## EXAMPLE: 
# <phylip> has 'IND1 ATC' and 'IND1-1 ATC' and <indfile> has 'IND1' --> 'IND1 ATC' will be grepped but not 'IND1-1 ATC'
# <phylip> has 'IND1ATCG' and 'IND1-1 ATC' and <indfile> has 'IND1' --> nothing will be grepped, with a WARNING

## AUTHOR: simon.crameri@env.ethz.ch, Apr 2019

## Define arguments
drop='false' # default if flag is not turned on 

while getopts p:f:d opts # argument 'd' is a flag and therefore has no ':' succeeding it
do
        case "${opts}"
        in
                p) phylip=${OPTARG};;
                f) indfile=${OPTARG};;
                d) drop='true';;
         esac
done

## Additional arguments
phybase=$(basename $phylip .phy)

## Check arguments
if [ ! $phylip ] ; then echo "phylip file not specified (-p option)" ; exit 0 ; fi
if [ ! $indfile ] ; then echo "taxonfile not specified (-f option)" ; exit 0 ; fi

if [ ! -f $phylip ] ; then echo "phylip file <$phylip> not found" ; exit 0 ; fi
if [ ! -f $indfile ] ; then echo "taxonfile <$indfile> not found" ; exit 0 ; fi

## Check for empty or duplicated individuals in $indfile
NINDfile=$(wc -l < $indfile)
if [ "$NINDfile" -eq 0 ] ; then echo ; echo "no entries in <$indfile>" ; echo ; exit 0 ; fi
NINDfileu=$(cat $indfile | sort | uniq | wc -l)
if [ "$NINDfile" -ne "$NINDfileu" ] ; then echo ; echo "detected $(expr $NINDfile - $NINDfileu) duplicated entries in <$indfile>" ; echo ; exit 0 ; fi

## Check for empty lines in $indfile
NEMPTYfile=$(cat $indfile | grep '^$' | wc -l)
if [ "$NEMPTYfile" -gt 0 ] ; then echo ; echo "detected $NEMPTYfile empty line(s) in <$indfile>" ; echo ; exit 0 ; fi

## Original dimensions
oldheader=$(head -n1 $phylip)
NINDorig=$(echo $oldheader | cut -f1 -d' ')
NALNorig=$(echo $oldheader | cut -f2 -d' ')

## Check that it is a single-line PHYLIP
NLINESoldphy=$(wc -l < $phylip)
NINDoldphy=$(expr "$NLINESoldphy" - 1)
if [ "$NINDorig" -ne "$NINDoldphy" ] ; then echo ; echo "appears to be a multi-line phylip file" ; echo ; exit 0 ; fi

## Create string for grepping (should also work if labels contain non-word-constituent characters such as the hyphen '-')
first=$(head -n1 $indfile)
string=$(echo "\b^$first\b ")
for sample in $(tail -n +2 $indfile) ; do string=$(echo "${string}|\b^${sample}\b ") ; done

## Keep or Drop taxa
if [ "$drop" = true ]
then 

	# drop
	verb="dropped"
	echo ; echo "dropping $NINDfile out of $NINDorig individuals..."
	NINDnew=$(expr "$NINDorig" - "$NINDfile")
	newheader="$NINDnew $NALNorig"
	output="${phybase}_$NINDnew.phy"
	if [ -f $output ] ; then echo ; echo "output file <$output> already exists" ; echo ; exit 0 ; fi

	echo "$newheader" > $output
	string="^$oldheader|$string"
	grep -vE "$string" $phylip >> $output

else

	# keep
	verb="kept"
	echo ; echo "keeping $NINDfile out of $NINDorig individuals..."
	NINDnew=$NINDfile
	newheader="$NINDnew $NALNorig"
	output="${phybase}_$NINDnew.phy"
	if [ -f $output ] ; then echo ; echo "output file <$output> already exists" ; echo ; exit 0 ; fi

	echo "$newheader" > $output
	grep -E "$string" $phylip >> $output

fi

## Check consistency of $indfile and $output
NLINESnewphy=$(wc -l < $output)
NINDnewphy=$(expr "$NLINESnewphy" - 1)

if [ "$NINDnewphy" -ne "$NINDnew" ]
then
	# adjust output header
	sed -i "1 s/^$NINDnew /$NINDnewphy /g" $output
	
	# warn	
	echo ; echo "WARNING: Unexpected number of kept taxa ($NINDnewphy)!"
	echo ; echo "<$indfile> appears to contain taxa that do not match any taxon in <$phylip>!"
	
	# find string(s) that are not present in $phylip
	grep -F -f $indfile $phylip | cut -f1 -d' ' > $indfile.log
	grep -v -f $indfile.log $indfile > $indfile.err
	echo ; echo "not found: $(grep -v -f $indfile.log $indfile)"
fi

## Be verbose
echo ; echo ; echo "Input header:  $oldheader"
echo "Number of taxa to be $verb: $NINDfile"
echo "Output header: $(head -n1 $output)"
echo
