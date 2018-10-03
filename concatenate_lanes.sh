#!/bin/bash

# set environment variables

DIR="./fq"
# get names of each replicate for each condition
# make directory for concatenated file to go
get_names() {
  FILE_STRING="$1/*"
  for file in $FILE_STRING; do
    filename=$(basename -- "$file")
    name=${filename:10:3}
    if [[ ! -e $name ]]; then
      mkdir $name
    fi
  done
}

# # extract all fastqs per replicate
# # move all to the corresponding
extract() {
  FILE_STRING="$1/*"
  for file in $FILE_STRING;do
    filename=$(basename -- "$file")
    name=${filename:10:3}
    mv $file/* ./$name
  done
}
# concatenate
concatenate() {
  for directory in */;do
    dir=$(basename -- "$directory")
    cd $directory
    cat `find . -maxdepth 1 | egrep '1.'` > $dir.fastq
    mv $dir.fastq ..
    cd ..
  done
}

####
# EXECUTE COMMANDS
####

get_names $DIR
extract $DIR
concatenate
