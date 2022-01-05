#!/bin/bash

## Usage: setup.CaptureAl.sh -d <directory with bash / R / python scripts> 
#                            -b <absolute path to bash executable> 
#                            -r <absolute path to Rscript executable> 
#                            -p <absolute path to python executable>

## Details:
# makes all .sh / .R / .py scripts in ${dir} executable
# sets the first line in each .sh / .R / .py script in ${dir} to #!${bashpath} or #!${rpath} or #!${pypath}, respectively
# default paths are /bin/bash for bash, /usr/bin/Rscript for Rscript and /usr/bin/python for python
# only changes the first line in a script if it is empty or begins with #!, writes any non-conforming script to setup.err

## Author: simon.crameri@env.ethz.ch, Apr 2019

## Get arguments (add more here if you need)
while getopts 'd:b:r:p:' flag; 
do 
  case "${flag}" 
  in
    d) dir="${OPTARG}" ;;
    b) bashpath="${OPTARG}" ;;
    r) rpath="${OPTARG}" ;;
    p) pypath="${OPTARG}" ;;
  esac
done

## Default arguments
bashdefault='/bin/bash'
rdefault='/usr/bin/Rscript'
pydefault='/usr/bin/python'

## Define functions
check_input()
{
	input=$1 ; language=$2 ; option=$3
	
	if [ ! -f $input ] ; then echo "$language executable <$input> (-$option option) not found" ; exit 0 ; fi
	if [ ! -x $input ] ; then echo "$input (-$option option) is not executable" ; exit 0 ; fi
}
export -f check_input

set_executable()
{
	scriptdir=$1 ; extension=$2 ; execpath=$3

	count=$(ls -1 ${scriptdir}/*${extension} 2>/dev/null | wc -l)
	
	if [ $count -gt 0 ] 
	then
		line="#!${execpath}"
		for script in $(ls -1 ${scriptdir}/*${extension})
		do
			chmod -u+rwx $script
		
			header=$(head -n1 $script)
			if [[ $header =~ ^#!* ]] || [[ $header =~ ^$ ]]
			then
				sed -i "1 s@.*@$line@g" $script
			else
				echo "$(basename $script)" >> setup.err
			fi
		done
	fi
}
export -f set_executable

## Check arguments
if [ ! $dir ] ; then echo "scripts directory not specified (-d option)" ; exit 0 ; fi
if [ ! -d $dir ] ; then echo "directory <$dir> not found" ; exit 0 ; fi

if [ ! $bashpath ] ; then  bashpath=$bashdefault ; fi
if [ ! $rpath ] ; then rpath=$rdefault ; fi
if [ ! $pypath ] ; then pypath=$pydefault ; fi

check_input $bashpath 'bash' 'b' $bashdefault
check_input $rpath 'Rscript' 'r' $rdefault
check_input $pypath 'python' 'p' $pydefault
# add more here if you need

## Execute
set_executable $dir '.sh' $bashpath
set_executable $dir '.R' $rpath ; set_executable $dir '.r' $rpath
set_executable $dir '.py' $pypath
# add more here if you need

## Finish
echo 
echo "Setup completed!"
echo
