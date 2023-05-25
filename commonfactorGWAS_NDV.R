## This script uses Genomic SEM to prepare factor summary statistics. For a way better explanation of how to run this and what it does, see https://github.com/GenomicSEM/GenomicSEM

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

## NDV
traits.snp.f1 <- c('adhd.sumstats.gz','asd.sumstats.gz','mdd.sumstats.gz','ts.sumstats.gz') # names of munged sum stats
sampleprev.snp.f1 <- c(.36,.40,.36,.34) # prevalence of disorder within sum stats
popprev.snp.f1 <- c(.05,.01,.15,.008) # global rate of disorder
trait.names.snp.f1 <- c('ADHD','ASD','MDD','TS') # names of traits

## COMP
traits.snp.f2 <- c('ocd.sumstats.gz','ts.sumstats.gz','an.sumstats.gz')
sampleprev.snp.f2 <- c(.28,.34,.23)
popprev.snp.f2 <- c(.025,.008,.01)
trait.names.snp.f2 <- c('OCD','TS','AN')

## MP
traits.snp.f3 <- c('bip.sumstats.gz','mdd.sumstats.gz','scz2018.sumstats.gz')
sampleprev.snp.f3 <- c(.39,.35,.39)
popprev.snp.f3 <- c(.01,.15,.01)
trait.names.snp.f3 <- c('BIP','MDD','SCZ')

ld <- 'eur_w_ld_chr/'
wld <- 'eur_w_ld_chr/'

## Now run LDSC

## NDV
ldsc.f1 <- ldsc(traits.snp.f1, sampleprev.snp.f1, popprev.snp.f1, ld, wld, trait.names.snp.f1, ldsc.log = 'OUR.F1.ldsc.log')

##F2
ldsc.f2 <- ldsc(traits.snp.f2, sampleprev.snp.f2, popprev.snp.f2, ld, wld, trait.names.snp.f2, ldsc.log = 'OUR.F2.ldsc.log')

##F3
ldsc.f3 <- ldsc(traits.snp.f3, sampleprev.snp.f3, popprev.snp.f3, ld, wld, trait.names.snp.f3, ldsc.log = 'OUR.F3.ldsc.log')

## ----------------------------------------------------------------------------
# ------------------------ Prep summary statistics --------------------------------

# set path to your (zipped) raw summary statistics
path2rawsumstats <- 'path/to/sumstats'

# Set path to ref data (downloaded here: https://utexas.box.com/s/vkd36n197m8klbaio3yzoxsee6sxo11v)
ref='reference.1000G.maf.0.005.txt'

##NDV
files.f1 = c('adhd.summarystatistics.txt', 'asd.summarystatistics.txt','mdd.summarystatistics','ts.summarystatistics.txt')
trait.names.p.f1 <- c('ADHD','ASD','MDD','TS')
se.logit.f1=c(T,T,T,T)

f1.sumstats <- sumstats(files=files.f1, ref=ref, trait.names = trait.names.p.f1, se.logit=se.logit.f1)

##COMP
files.f2 = c('ocd.summarystatistics.txt', 'ts.summarystatistics.txt','an.summarystatistics.txt')
trait.names.p.f2 <- c('OCD','TS','AN')
se.logit.f2=c(T,T,T)


f2.sumstats <- sumstats(files=files.f2, ref = ref, trait.names = trait.names.p.f2, se.logit = se.logit.f2)

##MP
files.f3 = c('bip.summarystatistics.txt', 'mdd.summarystatistics.txt','scz.summarystatistics.txt')
trait.names.p.f3 <- c('BIP','MDD','SCZ')
se.logit.f3=c(T,T,T)

f3.sumstats <- sumstats(files = files.f3, ref = ref, trait.names = trait.names.p.f3, se.logit = se.logit.f3)

# --------------------------------------------------------------------------
# ----------------- Now to generate summary statistics for each cluster of disorders ---------

f1_gwas <- commonfactorGWAS(covstruc = ldsc.f1, SNPs = f1.sumstats, parallel = TRUE, cores=38, toler = 1e-30)

f2_gwas <- commonfactorGWAS(covstruc = ldsc.f2, SNPs = f2.sumstats, parallel = TRUE, cores = 24, toler = 1e-30)

f3_gwas <- commonfactorGWAS(covstruc = ldsc.f3, SNPs = f3.sumstats, parallel = TRUE, cores = 20, toler = 1e-30)
    
