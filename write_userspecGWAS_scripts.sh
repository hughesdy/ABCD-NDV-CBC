#!/bin/bash

nums=(1  123456 246911  370366  493821  617276  740731  864186  987641 1111096 1234551 1358006 1481461 1604916 1728371 1851826 1975281 2098736 2222191 2345646 2469101 2592556 2716011 2839466 2962921 3086376 3209831 3333286 3456741 3580196 3703651 3827106 3950561 4074016 4197471 4320926 4444381 4567836 4691291 4814746 4938201 5061656 5185111 5308566 5432021 5555476 5678931 5802386 5925841)

while getopts l:s:o: flag
do
    case "${flag}" in
        l) path2ldsc=${OPTARG};;
        s) path2sumstats=${OPTARG};;
        o) outputdir=${OPTARG};;
    esac
done

## This loop will chunk the summary statistics and create R scripts to generate factor summary statistics.

for i in ${nums[@]}
do
printf "library(dplyr)

library(GenomicSEM)

ldsc.all <- readRDS('$path2ldsc')

sumstats <- readRDS('$path2sumstats')

sumstats <- pfac.sumstats[c($i:$((i+123454))),]

model <- 'F1 =~ ADHD + ASD + MDD + TS
F2 =~ AN + OCD + TS
F1~~F2
F3 =~ BIP + SCZ + MDD
F2~~F3
F1~~F3
F1~SNP
F2~SNP
F3~SNP'

correlatedfactors <- userGWAS(covstruc = ldsc.all, SNPs = sumstats, model = model, sub =c('F1~SNP', 'F2~SNP', 'F3~SNP'), parallel = FALSE)


f1 <- correlatedfactors[[1]]
f2 <- correlatedfactors[[2]]
f3 <- correlatedfactors[[3]]

saveRDS(f1, '$outputdir/f1/userspecgwas_$i')
saveRDS(f2, '$outputdir/f2/userspecgwas_$i')
saveRDS(f3, '$outputdir/f3/userspecgwas_$i')" >> $i.userspecgwas.r
done

