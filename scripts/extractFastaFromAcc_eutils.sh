#!/bin/bash

if [ $# -lt 2 ] ; then
	echo -e "\nERROR: Provide list with Accession Numbers (one per line)\n"
    echo -e "Usage: $0 acc.list outFileName\n"
    exit
fi

die () {
    echo "ERROR: $*. Aborting." >&2
    exit 1
}

acc_list=$1
outFile=$2

if [[ -z "$outFile" ]] ;then
	die "Please provide argv[2] output file name";
fi

accEutilsList=""

while read ACC
do
	accEutilsList="$accEutilsList$ACC,"
   #echo -n -e "$ACC\t"
   #curl -s "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${ACC}&rettype=fasta&retmode=xml" |\
   #grep TSeq_taxid |\
   #cut -d '>' -f 2 |\
   #cut -d '<' -f 1 |\
   #tr -d "\n"
   #echo

done <$acc_list

echo -e "\nRetrieving Fasta Sequences for accession list:\n$accEutilsList"
# curl -s "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=$accEutilsList&rettype=fasta&retmode=text" > $outFile
curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=$accEutilsList&rettype=fasta&retmode=text" > $outFile
echo -e "Done! Sequences have been written to $outFile"
