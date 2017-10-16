#!/bin/bash
#Wrapper to blast an assembly against a custom database, and positively select for a pre-defined set of reference sequences.

#MIT liscense
#Copyright 2017 Kyle Fletcher
#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, 
#including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, 
#subject to the following conditions:

#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
#IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE 
#OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

usage="$basename "$0") [-h] [-i] [-p] [-d] [-t] [-l] [-T] -- A wrapper to BLAST an assembly and filter for similarity to pre-determined reference sequences
	A new assembly file, filtered for similarity to the pre-defined reference sequences will be created in the current working directory called PREFIX.filt.fasta .
	We will also retain the blastn results in a directory called blastn, in case the user wishes to inspect them. Otherwise this directory can be safely deleted after run completeion.

	Please note, any pipes '|' present in the sequence headers will be removed and replaced with underscores '_'

Options:
	-h show this help message
	
Required Arguments:
	-i Assembly file to be filtered - can be gzipped.
	-p Output prefix
	-d BLAST formatted nucleotide reference database to query sequnces against
	-t BLASTdb type:
		1 = NCBI formatted db
		2 = Custome formatted db
	   We need this information because it will tell us how to filter the blast results.
	   NCBI formatted databases have gi information we want to strip away.
	-l List of reference sequences to positively filter assembly for
	
Optional Arguments:
	-T Number of threads to use [Default = 1]
	-W Word size for Blast [Default = 24]
		Smaller word size have a higher specificity at a trade-off with cpu time.
		If unsure compare results of 24 with 16 and optimize from there.

"

while getopts ':hi:p:d:l:t:T;' option; do
	case "$option" in
		h)  echo "$usage"
			 exit
			 ;;
		p)  Prefix=$OPTARG
			 ;;
		i)  Query=$OPTARG
			 ;;
		d)  DB=$OPTARG
			 ;;
		t)  RunType=$OPTARG
			 ;;
		l)  Ref=$OPTARG
			 ;;
		T)  Threads=$OPTARG
			 ;;
		W)  WS=$OPTARG
			 ;;
		\?) printf "illegal option: -%s\n\n" "$OPTARG" >&2
		    echo "$usage"
		    exit 1
		    	 ;;
	esac
done

#Check Arguments provided
if [[ -z $Prefix || -z $Query || -z $DB || -z $Ref ]]; then
	echo "
	ERROR: Some Variables are not set
"
	echo "$usage"
	exit 1
fi
if [[ $RunType != [1/2] ]]; then
echo "
           Error RunType set incorrectly. Please set to either:
                1 = NCBI formatted db
                2 = Custome formatted db
           We need this information because it will tell us how to filter the blast results.
           NCBI formatted databases have gi information we want to strip away.
"
exit 1
fi
#Set Default threads if not provided
if [[ -z $Threads ]];then
echo "Threads not specified, will proceed with default [1]
"
Threads=1
fi

#Set Default word size if not provided
if [[ -z $WS ]]; then
echo "Word Size not specified, will proceed with default [24]
"
WS=24
fi

#Check for existence of files to prevent over-writing
if [[ -e blastn/$Prefix.blastn ]]; then
echo "
$Prefix.blastn already exists. We don't want to delete something valuable. Please either:
	1. Choose a different prefix
	2. Move the blast results somewhere else
	3. Manually delete the blast results using:
		rm blastn/$Prefix.blastn
	Thanks
" >&2
exit 1
fi
if [[ -e $Prefix.filt.fasta ]]; then
echo "
$Prefix.filt.fasta already exists. We don't want to delete something valuable. Please either:
        1. Choose a different prefix
        2. Move the blast results somewhere else
        3. Manually delete the blast results using:
                rm $Prefix.filt.fasta
        Thanks
" >&2
exit 1
fi
#Check Dependencies
SAM=$(samtools 2> SAMtmp ; grep 'Version:' SAMtmp | awk '$2 >= 1.3 {print "Samtools version", $2, "detected, ok to proceed"} ; $2 < 1.3 {print "Samtools version", $2, "detected. We required Samtools version >= 1.3"}')
#clean up
rm SAMtmp
if [[ $SAM =~ 1.3$ ]]; then
	echo "$SAM"
	exit 1
else
	echo "$SAM"
fi

BLA=$(command -v blastn >/dev/null 2>&1 || { echo >&2 "blastn required"; })
if [[ $BLA != "" ]]; then
	echo "$BLA"
	exit 1
else
	echo "Blast detected"
fi
#BDB=$(command -v makeblastdb >/dev/null 2>&1 || { echo >&2 "makeblastdb"; })
#if [[ $BDB != "" ]]; then
#	echo "makeblastdb command not detected, if not pre-formatted then command will fail"
#else
#	echo "makeblastdb command detected, we may not need it, but good to know it's there!"
#fi

echo "Everything seems in order, ready to proceed
"

#Check if blastdb present
#DBN=$DB.nin
#if [ -e $DBN ] ; then
#	echo "Suitable nucleotide blast database detected, we will use this
#"
#else
#	echo "No blastdb found
#"
#	if [[ $BDB == "" ]]; then
#	echo "Good thing we have makeblastdb command available!"
#	makeblastdb -dbtype nucl -in $DB
#	else "Uh Oh! We require either the already formatted nucleotide database or makeblastdb to be in your PATH."
#	exit 1
#	fi
#fi

#Run with it
#Gunzip file - Not a common input, but good to have.
if [[ $Query =~ .gz$ ]]; then
gunzip $Query
Query2=$(echo $Query | sed 's/.gz$//')
else
Query2=$Query
fi

#Remove any pipes that may be lurking
sed 's/|/_/g' $Query2 > $Prefix.input.fa

#Run blast
mkdir -p blastn
blastn -query $Prefix.input.fa -db $DB -max_target_seqs 1 -outfmt 6 -word_size $WS -num_threads $Threads -out blastn/$Prefix.blastn
if [[$RunType == 1]]; then
awk -v FS='|' '{print $4}' blastn/$Prefix.blastn | sed 's/\..*//' | sort -u | comm -12 - <(sort $Ref | sed 's/\..*//') | join -2 3 - <(sed 's/|/ /3' blastn/$Prefix.blastn | sed 's/\..*//' | sort -k3,3) | awk '{print $2}' | sort -u > $Prefix.filt.h
fi
if [[$RunType == 2]]; then
awk '{print $2}' blastn/$Prefix.blastn | sort -u | comm -12 - <(sort $Ref) | join -2 2 - <(sort -k2,2 blastn/$Prefix.blastn) | awk '{print $2}' | sort -u > $Prefix.filt.h
fi
xargs samtools faidx $Prefix.input.fa < $Prefix.filt.h > $Prefix.filt.fasta
rm $Prefix.input.fa
rm $Prefix.filt.h

#Restore reference to gzipped state, if provided that way by the user
if [[ $Query =~ .gz$ ]]; then
gzip $Query2
fi

if [ -e $Prefix.filt.fasta ]; then
echo "Congratulations, it looks like everything went well. Please check the output file $Prefix.filt.fasta for your results.

If this file is empty then please ensure that you were using a NCBI formatted database
"
else 
echo "Oh dear! Something went wrong and the output file was not produced. Please check:
1. That all input files are correct
2. The the reference list contains the NCBI accessions that you wish to positively filter for
"
fi

exit
