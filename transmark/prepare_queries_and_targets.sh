#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --partition=standard
#SBATCH --account=twheeler
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16

if [ "$#" -ne 7 ]; then
	echo "useage: prepare_queries_and_targets.sh <basename> <benchmark directory> <blast directory> <last directory> <hmmer directory> <bathsearch directory> <split number>"
        echo "e.g. prepare_queries_and_targets.sh transmark-00 /home/user/transmark_files /home/user/ncbi-blast-2.13.0+/bin /home/user/last/bin /home/user/hmmer/src /home/user/BATH/src 100"
	exit 1;
fi


echo "benchmark basename " $1
echo "benchmark directory " $2
echo "blast directory " $3
echo "last directory " $4
echo "hmmer directory " $5
echo "bathsearch directory " $6 
echo "split number " $7

#cd into benchmark directory
cd $2

#make blast db from taget sequences
~/ncbi-blast-2.13.0+/bin/makeblastdb -dbtype nucl -in $1.fa

# make last train file file using full query set
$6/../easel/miniapps/esl-reformat fasta $1.AA.msa > $1.AA.fa
$4/lastdb -P16 -q -c $1.AA.fa.db $1.AA.fa
$4/last-train -P16 --codon $1.AA.fa.db $1.fa > $1.train

#index query msas
$6/../easel/miniapps/esl-afetch --index $1.AA.msa
$6/../easel/miniapps/esl-afetch --index $1.DNA.msa

#create direcory for split queries
mkdir queries

# get names of all query families
$6/../easel/miniapps/esl-alistat -1 $1.DNA.msa | grep -v '#' | awk {'print $2'}  > names.txt

readarray -t names < names.txt

#split query msas in "split number" serperate files and build the necessary query files for each one
for (( i=1; i<=$7; i++ ))
do
   echo "file " $i
   rm split_names.txt
   j=$(($i-1))
   
	while [ $j -lt ${#names[@]} ]
	do
             echo ${names[$j]} >> split_names.txt
             j=$(($j + $7))
	done

    #fetch query subsetd from msa
    $6/../easel/miniapps/esl-afetch -o queries/tbl$i.AA.msa -f $1.AA.msa split_names.txt
    $6/../easel/miniapps/esl-afetch -o queries/tbl$i.DNA.msa -f $1.DNA.msa split_names.txt

    #create unlaigned sequence file from the split Amino Acid query MSA
    $6/../easel/miniapps/esl-reformat fasta queries/tbl$i.AA.msa > queries/tbl$i.AA.fa

    # build last db file
    $4/lastdb -P16 -q -c queries/tbl$i.AA.fa.db queries/tbl$i.AA.fa

    #build DNA hmms for nhmmer
    $5/hmmbuild --cpu 16 queries/tbl$i.DNA.hmm queries/tbl$i.DNA.msa

    #build frameshift aware hmms for bathsearch
    $6/bathbuild --cpu 16 queries/tbl$i.AA.bhmm queries/tbl$i.AA.msa

    #build Amino Acid hmms to create consensus sequences
    $5/hmmbuild --cpu 16 queries/tbl$i.AA.hmm queries/tbl$i.AA.msa

    #emit a consensus sequence for each family
    $5/hmmemit -o queries/tbl$i.AA.cons.fa -c queries/tbl$i.AA.hmm

    #build a bathsearch hmm from the consensus sequence
    $6/bathbuild --cpu 16 --unali queries/tbl$i.AA.cons.bhmm queries/tbl$i.AA.cons.fa

done

