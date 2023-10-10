# Transmark

Transmark is a benchmark for testing the speed and accuracy of a series of sequence alignment tools in the annotation of protein-coding DNA with or without frameshift errors.

This directory includes:
1. The seed-files directory which includes zipped files containing all pre-made test and train sets
2. All code necessary to build a transmark benchmark
3. All code necessary to create the required query and target formats for the following sequence alignment software:
     1. bathsearch
     2. tblastn
     3. tfasty
     4. lastal
     5. nhmmer
4. Code to run and analyze the results from the above mentioned software

This directory DOES NOT include the following:
1. A protein sequence dataset from which to create decoy sequences for the benchmark.  We recommend using PFAM seed alignments available at: https://www.ebi.ac.uk/interpro/download/Pfam/
2. The actual sequence alignment software. Transmark is designed to run all sequence alignment tools from local copies of source code or binaries.  It accepts the full local paths to all executables so you do not need to have installed them via root access. The tools tested by transmark can be cloned or downloaded from the following locations. Be sure to follow directions to compile and/or make.
    1. BATH - https://github.com/TravisWheelerLab/BATH
    2. BLAST - https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
    3. FASTA - https://github.com/wrpearson/fasta36
    4. LAST - https://gitlab.com/mcfrith/last
    5. HMMER - https://github.com/EddyRivasLab/hmmer
  
Steps to create a transmark benchamrk
1. Clone the BATH-paper repository
2. Ensure you have a local copy of the HMMER3 repository - including EASLE - as well as the BLAST software (see above)
3. Ensure you have a dataset of protein sequences to use in the construction of decoys (see above)
4. cd into the transmark directory
5. Open the MAKEFILE and edit the ESLDIR and HMMERDIR variable to include the path to your local copies
6. run make
7. cd into the seed-files directory and unzip both the train-test.AA.tfile.zip and train-test.DNA.tfile.zip files
8. Run transmark-create.  To create a benchmark similar to the ones used in the BATH manuscript use the following set of command line arguments:



transmark-create --pre -N -L 100000000 -R 10 -E 10 --maxtrain 30 --maxtest 20 -D <decoy file> <benchmark name> <DNA seed file> <AA seed file> <blast path> <easle path> <background hmm file> <frameshift rate>

