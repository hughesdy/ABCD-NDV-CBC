## by: Dylan Hughes; this script tests differences in expression across brain structures available in brainspan, "reverse" annotates SNPs to genes (i.e. translate genes to SNPs), and then partitions summary statistics by genes showing unique expression

library(dplyr)

working.directory <- 'path/to/data'

setwd(working.directory)

#-----------------------------------------------------------
## Read in brainspan data
rows.cel <- read.csv('brainspan_expression/rows_metadata.csv', header=TRUE)

cols.cel <- read.csv('brainspan_expression/columns_metadata.csv', header=TRUE)

cel <- read.csv('brainspan_expression/expression_matrix.csv', header=FALSE)

cel.df <- cel[,-1]

# Pull unique structures out to easily reference later on
structures = unique(cols.cel$structure_acronym)

genenames <- dplyr::select(rows.cel, row_num, gene_symbol)

## Read in a gene "dictionary"
genes.dic <- read.table('gene_loc/NCBI37.3.gene.loc', header=FALSE)%>%
  dplyr::select(V1, V6)%>%
  dplyr::rename(GENE = V1, Gene_name = V6)


fornames <- arrange(cols.cel, structure_name)%>%
  filter(!duplicated(structure_name))%>%
  filter(!grepl('eminence', structure_name))%>%
  dplyr::select(structure_acronym, structure_name)


# ------------------------------------------------

#---- Before vs after birth (all rois, all genes) ---

# The first 239 columns correspond to expression before birt
b4birth <- c(1:239)

# The remainder represent expression after birth
afterbirth <- c(240:nrow(cols.cel))

# Select columns of cel corresponding to prenatal expression and flip it so that you have one row per timepoint/structure and one column per gene; create a variable called time and set it to string "pre"
b4birth.genes <- as.data.frame(t(cel.df[,b4birth]))%>%
  mutate(time = 'pre')

# Do the same but for postnatal expression
afterbirth.genes <- as.data.frame(t(cel.df[,afterbirth]))%>%
  mutate(time = 'post')

# Row bind the two above dataframes
recomb.b4after <- rbind(b4birth.genes, afterbirth.genes)

# -----------------------------------------------
# This function will compare pre to postnatal expression of genes within specific ROIs/structures. The user feeds the function a list of indices corresponding to expression values for an ROI and the function will automatically separate by pre and postnatal expression, run t-tests to determine significant overexpression, and spit out p-values and estimates.

compare.expression <- function(roicols) {
  firstcols = roicols[which(roicols<=239)]
  secondcols = roicols[which(roicols>239)]
  one <- recomb.b4after[firstcols,]%>%
    mutate(time = 'pre')
  two <- recomb.b4after[secondcols,]%>%
    mutate(time = 'post')
  both <- rbind(one, two)

  ncols <- ncol(both)

  peas = c(1:(ncols-1))
  bees = c(1:(ncols-1))
  pre=c()
  post=c()
  greater = c()

  for (i in c(1:(ncols-1))) {
    mod <- t.test(both[,i] ~ time, both)

    peas[i] = mod$p.value

    post[i] = mod$estimate[[1]]

    pre[i] = mod$estimate[[2]]

    bees[i] = post[i] - pre[i]

    greater[i] = ifelse(pre[i]>post[i], 'pre', 'post')

    print(i)
  }

  b4afterbirth.res <- data.frame('gene' = rows.cel$gene_symbol, 'peas' = peas, 'prenatal' = pre, 'postnatal' = post, 'bees'=bees, 'which.greater' = greater)

  b4afterbirth.res$fdr = p.adjust(b4afterbirth.res$peas, method = 'fdr')

  return(b4afterbirth.res)
}
# -------------------------------------

## Initialize some variables for each of the 6 structures analyzed

##Get cerebellum rows
cbc <- which(cols.cel$structure_acronym %in% c('CB', 'CBC','URL'))

# Amygdala
amyg <- which(cols.cel$structure_acronym%in%('AMY'))

# Hippocampus
hippo <- which(cols.cel$structure_acronym == 'HIP')

# medialdorsal nucleus of the thalamus
mdth <- which(cols.cel$structure_acronym%in%c('MD','DTH'))

# striatum
str <- which(cols.cel$structure_acronym == 'STR')

# neocortex
ncx <- which(cols.cel$structure_acronym%in%c('A1C','DFC','IPC','ITC','M1C','MFC','OFC','S1C','STC','V1C','VFC','Ocx','PCx','TCx','M1C-S1C'))
# ----------------------------------------------------------------------


# Run the above defined function for cerebellar expression
cbc.compare <- compare.expression(cbc)

# Pull genes that show significant prenatal expression
cbc.prenatal.genes <- cbc.compare$gene[which(cbc.compare$fdr<0.05 & cbc.compare$which.greater=='pre')]

# Pull genes that show significant postnatal expression
cbc.postnatal.genes <- cbc.compare$gene[which(cbc.compare$fdr<0.05 & cbc.compare$which.greater=='post')]

# Pull genes that show "continuous" (written as "constitutive" here) expression i.e., no significant difference in expression between pre and postnatal time points
cbc.constitutive.genes <- cbc.compare$gene[which(cbc.compare$fdr >= 0.05)]

# ----------------------
##Amygdala
amyg.compare <- compare.expression(amyg)

amyg.prenatal.genes <- amyg.compare$gene[which(amyg.compare$fdr<0.05 & amyg.compare$which.greater=='pre')]

amyg.postnatal.genes <- amyg.compare$gene[which(amyg.compare$fdr<0.05 & amyg.compare$which.greater=='post')]

amyg.constit.genes <- amyg.compare$gene[which(amyg.compare$fdr>=0.05)]

#----------------------------
#Hippocampus
hippo.compare <- compare.expression(hippo)

hippo.prenatal.genes <- hippo.compare$gene[which(hippo.compare$fdr < 0.05 & hippo.compare$which.greater == 'pre')]

hippo.postnatal.genes <- hippo.compare$gene[which(hippo.compare$fdr < 0.05 & hippo.compare$which.greater == 'post')]

hippo.constit.genes <- hippo.compare$gene[which(hippo.compare$fdr >= 0.05)]

# ------------------------------------
##mdthal
mdthal.compare <- compare.expression(mdth)

mdthal.prenatal.genes <- mdthal.compare$gene[which(mdthal.compare$fdr<0.05 & hippo.compare$which.greater == 'pre')]

mdthal.postnatal.genes <- mdthal.compare$gene[which(mdthal.compare$fdr<0.05 & mdthal.compare$which.greater == 'post')]

mdthal.constit.genes <- mdthal.compare$gene[which(mdthal.compare$fdr>=0.05)]

# ------------------------------------------
##striatum
str.compare <- compare.expression(str)

str.prenatal.genes <- str.compare$gene[which(str.compare$fdr < 0.05 & str.compare$which.greater == 'pre')]

str.postnatal.genes <- str.compare$gene[which(str.compare$fdr < 0.05 & str.compare$which.greater == 'post')]

str.constit.genes <- str.compare$gene[which(str.compare$fdr >= 0.05)]

# ---------------------------------------------
##neocortex
ncx.compare <- compare.expression(ncx)

ncx.prenatal.genes <- ncx.compare$gene[which(ncx.compare$fdr < 0.05 & ncx.compare$which.greater == 'pre')]

ncx.postnatal.genes <- ncx.compare$gene[which(ncx.compare$fdr < 0.05 & ncx.compare$which.greater == 'post')]

ncx.constit.genes <- ncx.compare$gene[which(ncx.compare$fdr >= 0.05)]
# -----------------------------------
# ---------------
# ------
# ------ Conceptual break point. Now we will move onto converting genes to SNPs
# ---------------
# ------------------------------------


# Read in the gene annotation file from MAGMA
annot <- read.table('annot/annot_35up10down.genes.annot', header=F, fill = T)


## Create tracker because: genes are in the first column followed by a list of their associated snps. Sometimes there are so many snps that they bleed into the next row. Thus we go row by row and categorize them as such: 0 the first column represents a gene and the following columns represent the associated snps. 1 means the first column represents a SNP meaning that it is a continuation of the previous line. Then we go back through the tracker vector and set it to a 2 if it is the beginning of a multi-row list of snps. Then we attach the vector to the annotation file such that we can quickly identify which SNPs belong to which genes. There's probably some more efficient way of doing this, but... oh well? I wrote this a while ago :)

# in summary:
## 0 says "this is a gene! there are some associated SNPs in this row, BUT not enough to fill the entire row"
## 1 says "this is NOT a gene! it is instead a SNP and is a continuation of the previous line. in other words, the SNPs in this row all belong to the gene above"
## 2 says "this is a gene! AND it has so many SNPs annotated to it that the list bleeds into the next row - watch out!

tracker = c(1:nrow(annot))

count = c()
for (i in c(1:nrow(annot))) {
  split = strsplit(annot[i,1], split = '')[[1]]
  first = paste(split[c(1,2)], sep = '', collapse = '')
  if (first == 'rs') {
    tracker[i] = 1
    next
  }
  if (first != 'rs') {
    tracker[i] = 0
    next
  }
}
for (i in c(1:(length(tracker)-1))) {
  if (tracker[i]==0 & tracker[i+1]==1) {
    tracker[i] = 2
    next
  }
}
for (i in c(1:nrow(annot))) {
  if (annot[i,1]%in%genes.dic$GENE) {
    annot[i,1]=genes.dic$Gene_name[which(annot[i,1]==genes.dic$GENE)]
    print(i)
    next
  }
  if (!annot[i,1]%in%genes.dic$GENE) {
    print(i)
    next
  }
}

annot$tracker = tracker


df <- data.frame('gene' = annot$V1, 'tracker' = annot$tracker)

# ----------------------------------------------------
## This function is used in the annotate.gene2snps function. When the loop reaches a gene with a long SNP list (i.e., a 2 in the tracker; or the start of a gene annotation with > 1 rows of SNPs), it iterates through the next set of rows and keeps track of the 1s - which again indicate rows (of SNPs) that belong/are annotated to the gene above. This effectively tells the function where all the SNPs are that are annotated to the gene of interest in the loop. This will hopefully become clearer when we get to the annotate.gene2snps function.

find1s <- function(x, start) {
  ones <- c()
  for (i in c((start+1):length(x))) {
    if (x[i]==1) {
      ones = append(ones, i)
      next
    }
    if (x[i]!=1) {
      break
    }
  }
  return(ones)
}
# -----------------------

# Converts a dataframe to a single vector..there's probably a built in function for this but I was clearly going through a loop phase....

dftolist <- function(df) {
  new <- c()
  for (i in c(1:ncol(df))) {
    new = append(new, df[,i])
  }
  return(new)
}
# ---------------------------

# This function converts genes to snps. I'll explain this function directly within the function (below). Oh right, it accepts a list of genes as an argument.

annotate.gene2snps <- function(genes) {
  snps = c()
  eye = 1
  for (i in c(1:nrow(annot))) {  # Iterate through annotation dataframe

    ## below: if the i-th row contains a gene that exists in the gene list provided to the function AND the tracker denotes that the first element of the row is a gene with its associated SNPs taking up less than a full row (i.e., marked by a 0 in the tracker), pull all of the SNPs from that row into the snps variable for safekeeping.
    if (annot[i,1]%in%genes & annot$tracker[i]==0) {
      dist = length(which(annot[i,c(3:(ncol(annot)-1))]!='')) # count the non-empty cells
      snps[c(eye:(eye+dist-1))] = as.character(annot[i, c(3:(1+(dist+1)))]) # pull the non empty cells
      eye = eye+dist # update eye
      print(nrow(annot)-i)
      next
    }

    # below: if the gene is in the provided gene list AND the tracker denotes the gene as a 2 or the start of a multi-row SNP list, pull and save those SNPs
    if (annot[i,1]%in%genes & annot$tracker[i]==2) {
      mainsnps <- as.character(annot[i,c(3:(ncol(annot)-1))])
      snps[c(eye:(eye+length(mainsnps)-1))] = mainsnps
      eye = eye+length(mainsnps)

      ones = find1s(annot$tracker, i) # here's the find1s function from above. it goes through the next set of rows in annot and determines how many of the following rows house SNPs that are annotated to the current gene

      one.snps = dftolist(annot[ones, c(1:(ncol(annot)-1))])
      one.snps = one.snps[which(one.snps!='')]
      snps[c(eye:(eye+length(one.snps)-1))] = one.snps # pull the SNPs and get rid of the elements that are blank

      eye = eye+length(one.snps)
      print(nrow(annot)-i)
      next
    }
    if (!annot[i,1]%in%genes) { # if gene not in gene list, skip iteration
      print(nrow(annot)-i)
      next
    }
  }
  return(snps)
}
# ---------------------------

# Okay, now we're going to use the annotate function to convert our pre and postnatally expressed genes into SNPs that we can use to filter our summary statistics. woo!

###Cerebellum
prenatal.snps.2 <- annotate.gene2snps(cbc.prenatal.genes)

postnatal.snps.2 <- annotate.gene2snps(cbc.postnatal.genes)

neither.snps.2 <- annotate.gene2snps(cbc.constitutive.genes)

##Amygdala
amyg.prenatal.snps <- annotate.gene2snps(amyg.prenatal.genes)

amyg.postnatal.snps <- annotate.gene2snps(amyg.postnatal.genes)

amyg.constit.snps <- annotate.gene2snps(amyg.constit.genes)

##Hippocampus
hippo.prenatal.snps <- annotate.gene2snps(hippo.prenatal.genes)

hippo.postnatal.snps <- annotate.gene2snps(hippo.postnatal.genes)

hippo.constit.snps <- annotate.gene2snps(hippo.constit.genes)

##mdthal
mdthal.prenatal.snps <- annotate.gene2snps(mdthal.prenatal.genes)

mdthal.postnatal.snps <- annotate.gene2snps(mdthal.postnatal.genes)

mdthal.constit.snps <- annotate.gene2snps(mdthal.constit.genes)

##striatum
str.prenatal.snps <- annotate.gene2snps(str.prenatal.genes)

str.postnatal.snps <- annotate.gene2snps(str.postnatal.genes)

str.constit.snps <- annotate.gene2snps(str.constit.genes)

##ncx
ncx.prenatal.snps <- annotate.gene2snps(ncx.prenatal.genes)

ncx.postnatal.snps <- annotate.gene2snps(ncx.postnatal.genes)

ncx.constit.snps <- annotate.gene2snps(ncx.constit.genes)
# ---------------------------------
# ---------------------
# Now that we've annotated genes to snps, we want to partition (filter) summary statistics by those SNPs so that we can generate partitioned polygenic scores.

## Load neurodev sumstats (or whatever sumstats you want to partition)
neurodev <- read.table('neurodev.for.clump.mich.txt', header=T)


## If your SNPs are in rs format, you don't need this "translation" part, but if they are in chrbp format, then you'll have to translate from rs format (which is how the genes were translated to SNPs) to chrbp format.

## Load translation
trans.snp <- read.table('g1000_eur.bim', header=F)

snp.dictionary <- mutate(trans.snp, chrome = paste(V1, V4, V5, V6, sep = ':'))%>%
  mutate(chrome2 = paste(V1,V4,V6,V5, sep = ':'))


# ----------------------------------
# This function takes sumstats, partitions them, and saves the files in specified directories
#1. sumstats: full summary statistics
#2. presnps: list of prenatal snps
#3. postsnps: list of postnatal snps
#4. insig: list of continuous/constitutive/insig snps
#5. topdirectory: directory to which you would like to save your summary statistics (this should already exist)
#6. roiname: character argument to append to the beginning of your file name that specifies which ROI (within BrainSpan) these summary statistics represent; eg. 'CBC'
#7. sumstatname: character argument to append to the beginning of your filename that specifiies what summary statistics were partitioned; eg. 'NDV'
    # ---> an example of the final file path would be: 'topdirectory/sumstatname.roiname.prenatal.sumstats.allFDR.txt'
#8. snpformat: default is 'chrbp'; you can also change it to 'rs' if your snps are in rs format

write.sumstats <- function(sumstats, presnps, postsnps, insig, topdirectory, roiname, sumstatname, snpformat='chrbp') {
  if (snpformat == 'chrbp') {
    prenatal.snp.translated <- snp.dictionary$chrome2[which(snp.dictionary$V2%in%presnps)]

    prenatal.snp.translated.2 <- snp.dictionary$chrome[which(snp.dictionary$V2%in%presnps)]

    prenatal.snp.trans.comb <- c(prenatal.snp.translated, prenatal.snp.translated.2)

    sumstats.prenatal <- filter(sumstats, SNP%in%prenatal.snp.trans.comb)

    write.table(sumstats.prenatal, paste(topdirectory,sumstatname,'.', roiname,'.prenatal.sumstats.allFDR.txt', sep = ''), col.names=T, row.names=F, quote = F)


    postnatal.snp.translated <- snp.dictionary$chrome2[which(snp.dictionary$V2%in%postsnps)]

    postnatal.snp.translated.2 <- snp.dictionary$chrome[which(snp.dictionary$V2%in%postsnps)]

    postnatal.snp.trans.comb <- c(postnatal.snp.translated, postnatal.snp.translated.2)

    sumstats.postnatal <- filter(sumstats, SNPID%in%postnatal.snp.trans.comb)

    write.table(sumstats.postnatal, paste(topdirectory,sumstatname, '.', roiname,'.postnatal.sumstats.allFDR.txt', sep = ''), col.names=T, row.names=F, quote = F)


    insig.snp.translated <- snp.dictionary$chrome2[which(snp.dictionary$V2%in%insig)]

    insig.snp.translated.2 <- snp.dictionary$chrome[which(snp.dictionary$V2%in%insig)]

    insig.snp.trans.comb <- c(insig.snp.translated, insig.snp.translated.2)

    sumstats.insig <- filter(sumstats, SNPID%in%insig.snp.trans.comb)

    write.table(sumstats.insig, paste(topdirectory, sumstatname, '.', roiname, '.insig.sumstats.allFDR.txt', sep = ''), col.names=T, row.names=F, quote=F)
  } else if (snpformat == 'rs') {

    sumstats.prenatal <- filter(sumstats, SNP%in%presnps)

    write.table(sumstats.prenatal, paste(topdirectory,sumstatname,'.', roiname,'.prenatal.sumstats.allFDR.rsformat.txt', sep = ''), col.names=T, row.names=F, quote = F)


    sumstats.postnatal <- filter(sumstats, SNP%in%postsnps)

    write.table(sumstats.postnatal, paste(topdirectory,sumstatname, '.', roiname,'.postnatal.sumstats.allFDR.rsformat.txt', sep = ''), col.names=T, row.names=F, quote = F)


    sumstats.insig <- filter(sumstats, SNP%in%insig)

    write.table(sumstats.insig, paste(topdirectory, sumstatname, '.', roiname, '.insig.sumstats.allFDR.rsformat.txt', sep = ''), col.names=T, row.names=F, quote=F)
  }

}
# -------------
# Now to use the write.sumstats function. The examples below save to the working.directory variable specified at the beginning of the script, but that can be changed.

##cbc
write.sumstats(sumstats = neurodev, presnps = prenatal.snps.2, postsnps = postnatal.snps.2, insig = neither.snps.2, topdirectory = working.directory, roiname = 'cbc', sumstatname = 'ndv')

##amyg
write.sumstats(sumstats = neurodev, presnps = amyg.prenatal.snps, postsnps = amyg.postnatal.snps, insig = amyg.constit.snps, topdirectory = working.directory, roiname = 'amyg', sumstatname = 'ndv')

##hip
write.sumstats(sumstats = neurodev, presnps = hippo.prenatal.snps, postsnps = hippo.postnatal.snps, insig = hippo.constit.snps, topdirectory = working.directory, roiname = 'hippo', sumstatname='ndv')

##mdthal
write.sumstats(sumstats = neurodev, presnps = mdthal.prenatal.snps, postsnps = mdthal.postnatal.snps, insig = mdthal.constit.snps, topdirectory = working.directory, roiname ='mdthal', sumstatname = 'ndv')

##str
write.sumstats(sumstats = neurodev, presnps = str.prenatal.snps,  postsnps = str.postnatal.snps, insig = str.constit.snps, topdirectory = working.directory, roiname = 'striatum', sumstatname = 'ndv')

##neocortex
write.sumstats(sumstats = neurodev, presnps=ncx.prenatal.snps, postsnps = ncx.postnatal.snps, insig = ncx.constit.snps, topdirectory = working.directory, name = 'neocortex', sumstatname = 'ndv')

