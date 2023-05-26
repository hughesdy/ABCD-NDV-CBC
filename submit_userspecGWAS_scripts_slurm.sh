#!/bin/bash

while getopts d: flag
do
    case "${flag}" in
        d) =${OPTARG};;
    esac
done

cd $d

for file in *
do
	
	printf "jobsubmit -A hughesdy -m 7G -p basic -t 0-16:00:00 /usr/pubsw/packages/R/4.0.2/bin/Rscript ${file}\n"

done
