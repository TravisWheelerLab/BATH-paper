#!/bin/bash


# This Program runs three tools (bathsearch, tbalstn, and last) to annotate the DNA of the Candidatus Hodgkinia genommes from 7 species of Magicicada using a set of protiens from Candidatus Hodgkinia in other cicada species

#stop the script if it attempts to use any unset variables
set -o nounset
#stop the script if any command fails (returns non zero)
set -o errexit


if [ $# -lt 4 ]; then
    echo "ERROR: Your command line contains less than four arguments ($bathsearch_dir, $tblastn_dir, $lastal_dir, $easel_dir),"
    exit 1
fi


# Prepare query and target files

# Make BATH hmm file from protien MSA
$1/bathbuild --ct 4 query_files/query_hmms.bhmm query_files/query_msas.stk

# Make sequences from protien MSA for tblastn
$4/miniapps/esl-reformat -o query_files/query_seqs.fa fasta query_files/query_msas.stk

# Get list of query sequence lengths for use in full length hit anaylsis
$4/miniapps/esl-seqstat -a query_files/query_seqs.fa > query_files/query_seqs.stat 

# Make tblastn target database files
$2/makeblastdb -dbtype nucl -in target_files/MAGI_seqs -out target_files/MAGI_seqs

# Make last query database files
$3/lastdb -P16 -q -c query_files/query_seqs.fa.db query_files/query_seqs.fa

# Make last codon training file
$3/last-train -P16 --codon query_files/query_seqs.fa.db target_files/MAGI_seqs  > target_files/MAGI_seqs.train

# Create results directory to hold results files
if ! [ -d "results_files" ]
then
mkdir results_files
fi

# Run the tools 

# Run bathsearch
$1/bathsearch --ct 4 -o results_files/magi-bath.out --tblout results_files/magi-bath.tbl query_files/query_hmms.bhmm target_files/MAGI_seqs

# Run tblastn
$2/tblastn -db_gencode 4 -outfmt 7 -out results_files/magi-tblastn.tbl -db target_files/MAGI_seqs -query query_files/query_seqs.fa

# Run lastal
$3/lastal -p target_files/MAGI_seqs.train -K1 -fTAB query_files/query_seqs.fa.db target_files/MAGI_seqs > results_files/magi-last.tbl 

