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

I am always in awe of reading others' scripts when they are well-commented. Comments are the bedrock of reproducible science! I am always trying to improve, so if you have suggestions please share. Feel free to reach out to me with any questions about the scripts: hughesdy@g.ucla.edu 

