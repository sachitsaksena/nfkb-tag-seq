# nfkb-tag-seq
The following readme is the code used for processing and analysis of TAGseq data in "The Effects of IKKβ Inhibition on Early NF-κB Activation and Transcription of Downstream Genes". The link to the article is provided here: https://www.sciencedirect.com/science/article/pii/S0898656818303061. In addition, the abstract is provided below. 

### Abstract
Small molecule approaches targeting the nuclear factor kappa B (NF-kB) pathway, a regulator of inflammation, have thus far proven unsuccessful in the clinic in part due to the complex pleiotropic nature of this network. Downstream effects depend on multiple factors including stimulus-specific temporal patterns of NF-kB activity. Despite considerable advances, genome-level impact of changes in temporal NF-kB activity caused by inhibitors and their stimulus dependency remains unexplored. This study evaluates the effects of pathway inhibitors on early NF-κB activity and downstream gene transcription. 3T3 fibroblasts were treated with SC-514, an inhibitor targeted to the NF-kB pathway, prior to stimulation with interleukin 1 beta (IL-1β) or tumor necrosis factor alpha (TNF-α). Stimulus induced NF-κB activation was quantified using immunofluorescence imaging over 90-minutes and gene expression tracked over 6-hours using mRNA TagSeq. When stimulated with IL-1β or TNF-α, significant differences (P < 0.05, two-way ANOVA), were observed in the temporal profiles of NF-κB activation between treated and untreated cells. Increasing numbers of differentially expressed genes (P < 0.01) were observed at higher inhibitor concentrations. Individual gene expression profiles varied in an inhibitor concentration and stimulus-dependent manner. The results in this study demonstrate small molecule inhibitors acting on pleiotropic pathway components can alter signal dynamics in a stimulus-dependent manner and affect gene response in complex ways.

This code relies on the code used for the original TagSeq analysis paper located at https://github.com/z0on/tag-based_RNAseq. The original paper can be found here: 

### Dependencies

#### Processing 
Command line interface. Perl. 

#### Analysis
Software | Version
------------ | -------------
R version 3.5.1 (2018-07-02)| 
tidyverse |
dplyr |
tibble |
gplots |
RColorBrewer |
genefilter |
calibrate |



library(RColorBrewer)
library(gplots)
library(genefilter)
library(calibrate)
library(DESeq2)
library(data.table)
library(wesanderson)

### Sequencing data format
This code is also generated using the file format for reads downloaded from Illumina's BaseSpace Sequence Hub: https://login.illumina.com/platform-services-manager/?rURL=https://basespace.illumina.com&clientId=basespace&clientVars=aHR0cHM6Ly9iYXNlc3BhY2UuaWxsdW1pbmEuY29tL2Rhc2hib2FyZA&redirectMethod=GET#/. The original data can be downloaded here: https://www.ncbi.nlm.nih.gov/sra?term=SRP163157. This repository will be updated with scripts to analyze the data in SRA format after installing sra-toolkit (https://ncbi.github.io/sra-tools/install_config.html) as follows:
```
$ fastq-dump SRP163157
```
For now, data in BaseSpace format can be requested from the original authors if needed (meghanbloom@utexas.edu). 

### Processing TagSeq data
First, clone this repository: 
```
$ git clone https://github.com/sachitsaksena/nfkb-tag-seq.git
$ cd nfkb-tag-seq
```
Move the BaseSpace downloaded TagSeq data to this repository: 
```
$ mv /path/to/Tag/Seq/data . 
```
Then run the following script while still in the nfkb-tag-seq directory to organize reads by experimental condition and concatenate lanes of each RNAseq run: 
```
$ bash concatenate_lanes.sh
```
Then clone the original repository while stil in the same directory for access to processing and clipping scripts scripts:
```
$ git clone https://github.com/z0on/tag-based_RNAseq
```
Now, to process simply run the following set of commands within the same directory:

First generate a set of commands that trims reads and run the process in "screen," which will take some time:
```
$ ./tagseq_trim_launch.pl '\.fastq$' > clean
# open a screen session
$ screen
# now in screen run
$ . clean
# ctrl + A + D to exit screen while still runnning process
# exit screen
$ exit
```
Then generate a set of commands to map to a reference genome, which will take much more time
```
$ ./tagseq_bowtie2map.pl "trim$" Mus_musculus.GRCm38.cdna.all.fa > maps
# open screen
$ screen
# now in screen
$ . maps
# ctrl + A + D to exit screen while still runnning process
# exit screen
$ exit
```
Now, run the following command:
```
$ grep '^>' Mus_musculus.GRCm38.cdna.all.fa |\
 perl -pe 's/>(\S+)\s.*gene:(\S+).*/\1\t\2/g' >\
 Mus_musculus.GRCm38.cdna.all_seq2iso.tab
```
Next, run samcounts via a perl script in the original author's repository:
```
$ ./samcount_launch_bt2.pl '\.sam' Mus_musculus.GRCm38.cdna.all_seq2iso.tab > sc
$ screen
$ . sc
## ctrl + A + D to exit screen while still runnning process
## "exit" to exit screen and kill process
$ exit
```

Finally, compile expression measurements from all the disparate samples: 
```
$ ./expression_compiler.pl *.sam.counts > allcounts.txt
```
A special thanks to Dennis Wylie at UT Austin for help with the command generating scripts.

### Analysis reproduction 
R scripts can be run from start to finish on the data produced from processing above. 
**check_replicates.R** is quality control code for ensuring high quality data by evaluating clustering of biological replicates. 

**diff_exp_analysis.R** is differential expression analysis used in the paper. 
