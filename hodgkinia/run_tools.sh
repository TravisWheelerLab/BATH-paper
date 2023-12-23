#!/bin/bash


# This Program runs three tools (bathsearch, bathsearch --nofs, tbalstn) to annotate the DNA of the Candidatus Hodgkinia genommes from 7 species of Magicicada using a set of protiens from Candidatus Hodgkinia in other cicada species

#stop the script if it attempts to use any unset variables
set -o nounset
#stop the script if any command fails (returns non zero)
set -o errexit


if [ $# -lt 3 ]; then
    echo "ERROR: Your command line contains less than four arguments ($bathsearch_dir, $tblastn_dir, $easel_dir),"
    exit 1
fi


# Prepare query and target files

# Make BATH hmm file from protien MSA
$1/bathbuild --ct 4 query_files/query_hmms.bhmm query_files/query_msas.stk

# Make sequences from protien MSA for tblastn
$3/miniapps/esl-reformat -o query_files/query_seqs.fa fasta query_files/query_msas.stk

# Make tblast tarrget database
cd target_files
$2/makeblastdb -dbtype nucl -in MAGI_seqs
cd ..

# Create results directory to hold results files
if ! [ -d "results_files" ]
then
mkdir results_files
fi

# Run the tools 

# Run bathsearch
$1/bathsearch --ct 4 -o results_files/magi-bath.out --tblout results_files/magi-bath.tbl query_files/query_hmms.bhmm target_files/MAGI_seqs

# Run bathsearch
$1/bathsearch --nofs --ct 4 -o results_files/magi-bath-nofs.out --tblout results_files/magi-bath-nofs.tbl query_files/query_hmms.bhmm target_files/MAGI_seqs

# Run tblastn
$2/tblastn -db_gencode 4 -outfmt 7 -out results_files/magi-tblastn.tbl -db target_files/MAGI_seqs -query query_files/query_seqs.fa


