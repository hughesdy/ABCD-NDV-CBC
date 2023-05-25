## This script produces correlated factor summary statistics (using userGWAS function of gSEM); please refer to https://github.com/GenomicSEM/GenomicSEM/wiki/5.-User-Specified-Models-with-SNP-Effects for a terrific walkthrough

# the first part of this (munging) is the same as commonfactorGWAS_NDV.r

library(devtools)
library(Matrix)
library(stats)

devtools::install_github('MichelNivard/GenomicSEM')
library(GenomicSEM)
# ------------------------------------------------

# First we must munge! List your summary statistics here

## Set the sample sizes for each GWAS
adhd = 53293
an = 72517
asd = 46350
bip = 51710
mdd = 173005
ocd = 9725
scz = 105318
ts = 14307

# set the file paths to each set of summary statistics
adhd.path = '/path/to/ADHDsumstats'
an.path = '/path/to/ANsumstats'
asd.path = '/path/to/ASDsumstats'
bip.path = '/path/to/BIPsumstats'
mdd.path = '/path/to/MDDsumstats'
ocd.path = '/path/to/OCDsumstats'
scz.path = '/path/to/SCZsumstats'
ts.path = '/path/to/TSsumstats'

# Munge should output/save munged summary statistics in the format <trait.name>.sumstats.gz
munge(files = c(adhd.path, an.path, asd.path, bip.path, mdd.path, ocd.path, scz.path, ts.path, trait.names=c('adhd','an','asd','bip','mdd','ocd','scz','ts'),N = c(adhd, an, asd, bip, mdd, ocd, scz, ts)))
# ---------------------------------------------------------------------------------
# ------------------------------- running LDSC --------------------------------------------

# First initialize some arguments

# set working directory to where your munged summary statistics live
mungedLocation = '/path/to/mungedSumStats'

setwd(mungedLocation)

# **** This is where we diverge from commonfactorGWAS_NDV.r
#----------                           ------------                  ---------           -

## Specify munged summary statistics names
traits <- c('adhd.sumstats.gz','an.sumstats.gz','asd.sumstats.gz','bip.sumstats.gz','mdd.sumstats.gz','ocd.sumstats.gz','scz.sumstats.gz', 'ts.sumstats.gz')

# Prevalence of disorder in GWAS sample
sample.prev <- c(.36,.24,.40,.39,.35,.28,.25,.34)

# Global prevalence of disorder
population.prev <- c(.05,.01,.01,.01,.15,.025,.01,.008)

## Point to reference files; can be downloaded here: https://utexas.box.com/s/vkd36n197m8klbaio3yzoxsee6sxo11v
ld <- "eur_w_ld_chr/"
wld <- "eur_w_ld_chr/"

## Name traits
trait.names<-c('ADHD','AN','ASD','BIP','MDD','OCD','SCZ','TS')

# Specify name of log file
myLogFile = 'i_am_not_natural'

# Run LDSC on all traits
LDSCoutput <- ldsc(traits, sample.prev, population.prev, ld, wld, trait.names, ldsc.log = myLogFile)

# Save LDSC output as .rds so that you can access it later and not have to run ldsc again (optional)
saveRDS(LDSCoutput, 'all_traits_LDSC')

## -----------------------------------------------------------------------------------------
##---------------------------  prepare summary statistics

files.all = c('adhd.summarystatistics.txt','an.summarystatistics.txt','asd.summarystatistics.txt',
              'bip.summarystatistics.txt','mdd.summarystatistics.txt','odc.summarystatistics.txt',
              'scz.summarystatistics.txt','ts.summarystatistics.txt')

trait.names.all <- c('ADHD','AN','ASD','BIP','MDD','OCD','SCZ','TS')

se.logit.all <- c(T,T,T,T,T,T,T,T)

ref='PGC_Files/reference.1000G.maf.0.005.txt'

# run sumstats on all files
alltraits_sumstats = sumstats(files = files.all, ref = ref, trait.names = trait.names.all, se.logit = se.logit.all)

# save output in rds for easy acces (optional)
saveRDS(alltraits_sumstats, 'all_traits_sumstats')

# -----------------------
## -----------------------------------------------------------------------------------------
# So, now that you've run all of that, you can either continue on to run it all on one machine, which is certainly doable. Might be worth parellelizing it if that's the case though. Alternatively, you can split the summary statistics into chunks, run gSEM on the chunks and then combine them back after. 
## -----------------------------------------------------------------------------------------
# -----------------------


# Specify model formula
model <- "F1 =~ ADHD + ASD + MDD + TS
F2 =~ AN + OCD + TS
F1~~F2
F3 =~ BIP + SCZ + MDD
F2~~F3
F1~~F3
F1~SNP
F2~SNP
F3~SNP"

correlatedfactors <- userGWAS(covstruc = LDSCoutput, SNPs = alltraits_sumstats, model = model, sub =c('F1~SNP', 'F2~SNP', 'F3~SNP'), parallel = F)

## NDV sumstats
f1 <- correlatedfactors[[1]]

## COMP sumstats
f2 <- correlatedfactors[[2]]

## MP sumstats
f3 <- correlatedfactors[[3]]

