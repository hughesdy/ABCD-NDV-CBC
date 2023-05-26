# ABCD-NDV-CBC

Hey there, thanks for stopping by. 
--------

This is a page housing the code relevant to our work: Hughes, D.E., Kunitoki, K., Elyounssi, S. et al. Genetic patterning for child psychopathology is distinct from that for adults and implicates fetal cerebellar development. Nat Neurosci (2023). https://doi.org/10.1038/s41593-023-01321-8

Link to full text: https://rdcu.be/dcJWk
#### ----------------------------------------------------

## partitioningSumStats.r
Hopefully, the comments within the script will guide you. In brief, this is the script that generates partitioned NDV summary statistics. For example, this is how the NDV scores representing genes prenatally, postnatally, and continuously expressed in the cerebellum were created. 

The data referenced in partitioningSumStats.r can be found in the "partitioning" directory within this repository. 

Side note: I think this is a pretty cool tool. By using the annotate.gene2snps function, you can really feed it any set of genes and thus partition your summary statistics based on any arbitrary (hypothesis- / data-driven) set of genes. This version uses an annotation file that was created via MAGMA (log file is available in the "partitioning" directory) and annotates SNPs to genes within the gene range and also 35 kb upstream and 10 kb downstream of the gene. Changing the window size will change the results of course. You can view the log file in the "partitioning" directory to look at the commands to use. More MAGMA documentation can be found here: https://ctg.cncr.nl/software/magma 
#### ----------------------------------------------------

## makingHeatmaps.r
This script shows how the main heatmaps were created (see Figure 2a for an example). The example shown here is for baseline data. 
#### ----------------------------------------------------

## commonfactorGWAS_NDV.R
This script generates summary statistics for the commonfactor NDV. The GenomicSEM site has a wonderful explanation + tutorial for more information: https://github.com/GenomicSEM/GenomicSEM
#### ----------------------------------------------------

## userspecGWAS_NDV.R
This script generates summary statistics for the userspec NDV (aka correlated factor model).
#### ----------------------------------------------------

## write_userspecGWAS_scripts.sh
This is a bash script that chunks the summary statistics into even parts so that you can submit each to a compute cluster; an example code of how you would submit the outputs of "write_userspecGWAS_scripts.sh" can be found in "submit_userspecGWAS_scripts_slurm.sh". This script takes 3 arguments: l, s, o. 

-l and -s will specify the path to the ldsc and sumstats (respectively) output from 'userspecGWAS_NDV.R' i.e., the 2 instances where saveRDS is called. -o will specify the output directory for the final chunks of the summary statistics. IMPORTANTLY, the output directory should already be created AND should have three folders within names 'f1', 'f2', and 'f3'. If you are running your own version of this with more or less than 3 factors, you'll have to change that part of this script.

An example of a call would be:
      bash write_userspecGWAS_scripts.sh -l path/to/ldsc -s path/to/sumstats -o path/to/where/final/rscripts/will/live
      
Of note, this same logic can be applied to the commonfactor model to make generation of summary statistics more efficient/quick
#### ----------------------------------------------------

## submit_userspecGWAS_scripts_slurm.sh
As noted above, this script is just an example of how you could submit your auto-generated Rscripts from "write_userspecGWAS_scripts.sh". If you happen to be from Martinos (shout out Martinos), this example is formatted for slurm so you'll need minimum editing :)
#### ----------------------------------------------------

I am always in awe of reading others' scripts when they are well-commented. Comments are a pillar of reproducible science! I am always trying to improve, so if you have suggestions please share. Feel free to reach out to me with any questions about the scripts: hughesdy@g.ucla.edu 

