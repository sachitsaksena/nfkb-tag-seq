#!/bin/bash
./tagseq_trim_launch.pl '\.fastq$' > clean
## open screen, run
## . clean
## ctrl + A + D to exit screen while still runnning process
## "exit" to exit screen and kill process

## http://weizhongli-lab.org/RNA-seq/Data/reference-genomes/
bowtie2-build Mus_musculus.GRCm38.cdna.all.fa Mus_musculus.GRCm38.cdna.all.fa


./tagseq_bowtie2map.pl "trim$" Mus_musculus.GRCm38.cdna.all.fa > maps
## open screen, run
## . maps
## ctrl + A + D to exit screen while still runnning process
## "exit" to exit screen and kill process

grep '^>' Mus_musculus.GRCm38.cdna.all.fa |\
 perl -pe 's/>(\S+)\s.*gene:(\S+).*/\1\t\2/g' >\
 Mus_musculus.GRCm38.cdna.all_seq2iso.tab


./samcount_launch_bt2.pl '\.sam' Mus_musculus.GRCm38.cdna.all_seq2iso.tab > sc
## open screen, run
## . sc
## ctrl + A + D to exit screen while still runnning process
## "exit" to exit screen and kill process


./expression_compiler.pl *.sam.counts > allcounts.txt
