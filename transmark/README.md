# Transmark

Transmark is a benchmark for testing the speed and accuracy of a series of sequence alignment tools in the annotation of protein-coding DNA with or without frameshift errors.

This directory includes:
1. The seed-files directory which includes zipped files containing all pre-made test and train sets
2. All code necessary to build a transmark benchmark
3. All code necessary to create the required query and target formats for the following sequence alignment software:
     1. bathsearch
     2. tblastn
     3. mmseqs
     4. lastal
     5. nhmmer
     6. diamond
4. Code to run and analyze the results from the above mentioned software

This directory DOES NOT include the following:
1. A protein sequence dataset from which to create decoy sequences for the benchmark.  We recommend using PFAM seed alignments available at: https://www.ebi.ac.uk/interpro/download/Pfam/
2. The actual sequence alignment software. Transmark is designed to run all sequence alignment tools from local copies of source code or binaries.  It accepts the full local paths to all executables so you do not need to have installed them via root access. The tools tested by transmark can be cloned or downloaded from the following locations. Be sure to follow directions to compile and/or make.
    1. BATH    - https://github.com/TravisWheelerLab/BATH
    2. BLAST   - https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
    3. MMSEQS2 - https://github.com/soedinglab/MMseqs2
    4. LAST    - https://gitlab.com/mcfrith/last
    5. HMMER3  - https://github.com/EddyRivasLab/hmmer
    6. DIAMOND - https://github.com/bbuchfink/diamond

Steps to create a transmark benchmark
1. Clone the BATH-paper repository
2. Ensure you have a local copy of the HMMER3 repository - including EASLE - as well as the BLAST software (see above)
3. Ensure you have a dataset of protein sequences to use in the construction of decoys (see above)
4. cd into the transmark directory
5. Open the MAKEFILE and edit the ESLDIR and HMMERDIR variable to include the path to your local copies
6. run make
7. cd into the seed-files directory and unzip both the train-test.AA.tfile.zip and train-test.DNA.tfile.zip files
8. Run transmark-create.  To create a benchmark similar to the ones used in the BATH manuscript use the following set of command line arguments:
```
./transmark-create --pre -N 10 -L 100000000 -R 10 -E 10 --maxtrain 30 --maxtest 20 --maxfams 1500 -D <decoy file> <benchmark name> <DNA seed file> <AA seed file> <blast path> <easel path> <background hmm file> <frameshift rate>
```
The \<DNA seed file\>, \<AA seed file\>, and \<background hmm file\> are all files included in this directory.
The \<decoy file\>, \<blast path\>, and \<easel path\> are the paths to your local copies of the executables or data sets
The \<benchmark name\> is whatever you want the base names of the output files to be
The \<frameshift rate\> is your desired rate of simulated frameshifts in the embedded test sequences of this benchmark

For example:
```
./transmark-create --pre -N 10 -L 100000000 -R 10 -E 10 --maxtrain 30 --maxtest 20 --maxfams 1500 -D ~/PFAM/Pfam-A.seed tanskmark00 seed-files/train-test.DNA.tfile seed-files/train-test.AA.tfile  ~/ncbi-blast-2.14.1+/bin/ ~/hmmer/easel/  0.0
```
To create an overextension benchmark like the ones used in the BATH manuscript (which embeds only the middle 50% of each test sequence) simply add the --over flag. 

Steps to run a transmark benchmark
1. Ensure you have local copies of all 6 tested software packages (see above)
2. Build a transmark benchmark (see above)
3. Run prepare_queries_and_targets.sh to create all the necessary file formats
```
bash prepare_queries_and_targets.sh <benchmark name> <benchmark directory> <blast directory> <last directory> <mmseqs directory> <hmmer directory> <bath directory> <diamond directory> <cpus available>
```
For example:
```
bash prepare_queries_and_targets.sh  tanskmark00 ~/BATH-paper/transmark/ ~/ncbi-blast-2.14.1+/bin/ ~/last/bin/ ~/MMseqs2/build/src ~/hmmer/src/ ~/BATH/src ~/diamond 16
```
The /<cpus available/> is what it sounds like - the number of CPUs available for paralyzing on your system. 

4. You should now have all the files you need to run the tools on your benchmark. To do this you will run each of the files starting with "x-" followed by the toolname. Each of these files takes the same set of command line arguments (with the exception of x-last):
```
x-toolname <tool directory> <benchmark directory> <results directory> <results filename>  <queryfile> <posfile> <targetfile> <cpus> <e-max>
```
For example:
```
x-bathsearch ~/BATH/src ~/BATH-paper/transmark/ ~/transmark00-results/BATH/ bath-results  transmark00-AA.bhmm transamark00.pos transmark00.fa 16 100
```

To run x-last you need one additional argument at the end of the command giving the average length of your target sequences. This can be determined by running esl-seqstat on the target file.

5. Once you have run a tool on a transmark benchmark you can analyze the results by running transmark-pp.sh with the following command line arguments.
```
bash transmark-pp.sh <benchmark directory> <results filename> <results directory> <analysis filename> <toolname> 
```
For example:
```
bash transmark-pp.sh ~/BATH-paper/transmark bath-results ~/transmark00-results/BATH/ bath-analysis bathsearch
``` 
This will produce the following files:
	1. bath-analysis.xy  - contains the x and y coordinates for creating ROC plots
	2. bath-analysis.mer - contains important per hit and per family information
        3. bath-analysis.tpr - contains the per-family data needed to make RBFP ridgeline plots
	4. bath-analysis.cover - contains coverage and overextension data for all true positives
        5. bath-analysis.time - contains run time information
 
