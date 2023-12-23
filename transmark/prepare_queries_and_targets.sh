#!/bin/bash

if [ "$#" -ne 9 ]; then
	echo "useage: prepare_queries_and_targets.sh <basename> <benchmark directory> <blast directory> <last directory> <mmseqs directory> <hmmer directory> <bathsearch directory> <diamond directory> <cpus available>"
        echo "e.g. prepare_queries_and_targets.sh transmark-00 /home/user/transmark_files /home/user/ncbi-blast-2.13.0+/bin /home/user/last/bin /home/user/MMseqs2/build/src /home/user/hmmer/src /home/user/BATH/src /home/user/diamond 16"
	exit 1;
fi


echo "benchmark basename " $1
echo "benchmark directory " $2
echo "blast directory " $3
echo "last directory " $4
echo "mmseqs diectory " $5
echo "hmmer directory " $6
echo "bath directory " $7 
echo "diamons directory " $8
echo "cpus available " $9

#cd into benchmark directory
cd $2

#BATH
# make bath query hmms
$7/bathbuild --cpu $9 $1.AA.bhmm $1.AA.msa

#TBLASTN
#make blast db from taget sequences
~/ncbi-blast-2.13.0+/bin/makeblastdb -dbtype nucl -in $1.fa

#create unlaigned sequence file from the split Amino Acid query MSA
$7/../easel/miniapps/esl-reformat fasta $1.AA.msa > $1.AA.fa

#MMSEQS
# make mmseqs db from taget sequences
$5/mmseqs createdb $1.fa $1.mmseqs

# make mmseqs profile from query sequences
$5/mmseqs convertmsa --identifier-field 0 $1.AA.msa $1.AA.afa
$5/mmseqs msa2profile $1.AA.afa $1.AA.profile --match-mode 1

#LAST
# make last train file file using full query set
$4/lastdb -P$9 -q -c $1.AA.last $1.AA.fa
$4/last-train -P$9 --codon $1.AA.last $1.fa > $1.AA.train

#NHMMER
# make hmmer DNA hmms
$6/hmmbuild --cpu $9 $1.DNA.hmm $1.DNA.msa

#DIAMOND
# make diamond database from query file
$8/diamond makedb --in $1.AA.fa -d $1.AA.diamond

