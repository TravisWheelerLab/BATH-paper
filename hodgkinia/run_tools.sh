#!/bin/bash


# This Program runs three tools (frahmmer, hmmsearcht, and tblastn) to annotate the DNA of the Candidatus Hodgkinia genommes from 7 species of Magicicada using a set of protiens from Candidatus Hodgkinia in other cicada species

#stop the script if it attempts to use any unset variables
set -o nounset
#stop the script if any command fails (returns non zero)
set -o errexit


if [ $# -lt 4 ]; then
    echo "ERROR: Your command line contains less than four arguments ($frahmmerdir, $hmmsearchtdir, $tblastndir, $easeldir),"
    exit 1
fi


# Prepare query and target files

# Make frahmmer hmm file from protien MSA
$1/frahmmbuild --ct 4 query_files/query_hmms.fhmm query_files/query_msas.stk

# Make hmmer hmm file from protien MSA
$2/hmmbuild query_files/query_hmms.hmm query_files/query_msas.stk

# Make sequences from protien MSA for tblastn
$4/miniapps/esl-reformat -o query_files/query_seqs.fa fasta query_files/query_msas.stk

# Get list of query sequence lengths for use in full length hit anaylsis
$4/miniapps/esl-seqstat -a query_files/query_seqs.fa > query_files/query_seqs.stat 

# Make tbalstn target database files
$3/makeblastdb -dbtype nucl -in target_files/MAGI_seqs -out target_files/MAGI_seqs


# Create results directory to hold results files
if ! [ -d "results_files" ]
then
mkdir results_files
fi

# Run the tools 

# Run frahmmer
$1/frahmmer --ct 4 -o results_files/magi-frahmmer.out --tblout results_files/magi-frahmmer.tbl query_files/query_hmms.fhmm target_files/MAGI_seqs

# Run hmmersearcht
$2/hmmsearcht -c 4 -o results_files/magi-hmmsearcht.out --domtblout results_files/magi-hmmsearcht.tbl query_files/query_hmms.hmm target_files/MAGI_seqs

# Run tblastn
$3/tblastn -db_gencode 4 -outfmt 7 -out results_files/magi-tblastn.tbl -db target_files/MAGI_seqs -query query_files/query_seqs.fa


