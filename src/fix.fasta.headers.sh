#!/bin/bash

## Usage: fix.fasta.headers.sh <.fasta> <i[nteractive]>

## Value: replaces any '|' or ' ' character in FASTA headers with '_', preventing any problems in downstream analyses.

## Get arguments
ref=$1
refbase=$(basename $ref .fasta)
refbase=$(basename $refbase .fas)
interactive=$2

## Check fasta headers
grep '^>' ${ref} > tmp.headers
if [ -f tmp.badchar ] ; then /bin/rm -f tmp.badchar ; fi
touch tmp.badchar

while read i
do 
        t1=$(echo $i | cut -d " " -f1)
        t2=$(echo $i | cut -d "|" -f1)
        if [[ ( $t1 != $i ) || ( $t2 != $i ) ]] ; then echo "stop" >> tmp.badchar ; fi
done < tmp.headers

ncheck=$(wc -l tmp.badchar | cut -d " " -f1)

if [ $ncheck != 0 ]
then	
	if [[ $interactive == "i" ]]
	then
		echo "Found ' ' and/or '|' characters in headers of <$(basename $ref)>"
        	read -r -p "Fix headers in new file? (replace ' ' and/or '|' with '_') [y/n] " response
        	case $response in

		[yY] | [yY][eE][sS] )
           
               	 	# Save original file as *.badheaders.fasta
              	 	mv ${ref} ${refbase}.badheaders.fasta

			# Change fasta headers
              	 	cat ${refbase}.badheaders.fasta |sed -e 's/ /_/g' |sed -e 's/|/_/g'> ${ref}

                        /bin/rm -f tmp.badchar
                        /bin/rm -f tmp.headers

			;;

		[nN] | [nN][oO] )

			/bin/rm -f tmp.badchar
			/bin/rm -f tmp.headers
			
			;;

		esac

	else
		# Save original file as *.badheaders.fasta
                mv ${ref} ${refbase}.badheaders.fasta

                # Change fasta headers
                cat ${refbase}.badheaders.fasta |sed -e 's/ /_/g' |sed -e 's/|/_/g'> ${ref}

                /bin/rm -f tmp.badchar
                /bin/rm -f tmp.headers

		echo "Fixed ' ' and/or '|' characters (replaced with '_') in <$(basename $ref)>"

	fi
else
	echo "All headers passed test!"
fi
