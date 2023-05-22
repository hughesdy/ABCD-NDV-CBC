#THIS SCRIPT REGRESSES PHENOTYPE/PSYCHOPATHOLOGY ON PRS AND CREATES A HEATMAP TO VISUALIZE THE RESULTS
#By: Dylan Hughes


# Load relevant libraries
library(dplyr)
library(pheatmap)
conflict_prefer("select", "dplyr")

## Load dylanfuncs
library(devtools)
document('path/to/dylanfuncs')

## Define home directory (e.g., where your data live)
homedir <- '/Volumes/My Passport/'

# Read in data
master <- read.csv(paste0(homedir, 'dataFileName.csv'))

#----------------------------------------------------------------
# This function below converts the p-values to n different categories for the heatmap function. Args: x is a matrix of p-values, levels are the thresholds. Ex) from the below initialization of the levels variable, for any p-value greater than 0.05, the output would be 0.

levels = c(0.05, 0.005, 5e-4, 5e-5, 5e-6, 5e-7, 5e-8)

adjustmatrix.any <- function(x, levels) {  
  
  new=as.data.frame(matrix(ncol=ncol(x), nrow=nrow(x)))
  
  for (i in c(1:ncol(x))) {
    for (j in c(1:nrow(x))) {
      for (a in c(1:length(levels))) {
        if (a == 1) {
          if (x[j,i] > levels[a]) {
            new[j,i] = 0
          } else {
            next
          }
        } else {
          if (x[j,i] <= levels[a-1] & x[j,i] > levels[a]) {
            new[j,i] = a-1
          } else if (a==length(levels)) {
            if (x[j,i]<=levels[a]) {
              new[j,i] = a
            }
          }
        } 
      }  
    }
  }
  return(new)
}

# ----------------------------------------------
# This function replaces P < 0.05 with an asterisk (*) to add to the heatmap if the FDR adjusted P value is less than 0.05. Args: x is the matrix of significance values

updatefdr <- function(x){   
  peasfhm <- matrix(ncol=ncol(x),nrow=nrow(x))
  for (a in c(1:(ncol(x)))) {
    for (i in c(1:nrow(x))) {
      if (x[i,a] < 0.05) {
        peasfhm[i,a] = '*'
      } 
      else {
        peasfhm[i,a] = ' '
      }
    }
  }
  return(peasfhm)
}

#----------------------------------------
# This function determines the direction of the effect and adjusts mapping parameters accordingly. The two matrices (x and y) must have the same dimensions and each corresponding location (e.g., x[i,j] and y[i,j]) must correspond to the same model. Specifically, y has the effect size of the term of interest from the model, and x has a number that represents the significance of the term (calculated from adjustmatrix.any function).

direction <- function(x,y) {   
  for (i in c(1:ncol(x))) {
    for (j in c(1:nrow(x))) {
      if (y[j,i] < 0) {
        x[j,i] = x[j,i]*-1
      }
    }
  }
  return(x)
}

## --------------------------------------------
# This function will iterate and permute through your specified dependent (dvs) and independent (ivs) variables and generate effect sizes, significance values, and standard errors from each of the models. By default the models are lmer (to include random effects for site). Args: 
# 1 and 2. dvs and ivs are the lists of dependent and independent variables respectively. 
# 3. pcs: names of the variables in your dataframe that contain genetic principal components. 
# 4. age: specifies the name of the age variable. sex the sex variable. x is the dataframe. 
#5. breaks: indicate the breaks at which to threshold your heatmap (these are just passed through this function to the mapping function). 
# 6. sex.interaction: a boolean value and allows you to add a term for sex interaction (not included in Hughes et al. 2023 results)

cbclheatmap <- function(dvs, ivs, pcs, age, sex, x, breaks, sex.interaction=F) {
  require(lme4)
  require(pheatmap)
  
  name = colnames(x)
  cbcl=c(1:length(dvs))
  for (i in cbcl) {
    cbcl[i] = grep(dvs[i],name)
  }
  prs=c(1:length(ivs))
  for (i in prs) {
    prs[i] = grep(ivs[i], name)
  }
  
  peas = c()
  bees = c()
  sees = c()
  int.p = c()
  int.b = c()
  nobs=9999
  
  pea.df = as.data.frame(matrix(ncol=length(prs), nrow=length(cbcl)))
  bee.df = as.data.frame(matrix(ncol=length(prs), nrow=length(cbcl)))
  see.df = as.data.frame(matrix(ncol=length(prs), nrow=length(cbcl)))
  
  if (sex.interaction == F) {
    for (i in prs) {
      for (a in cbcl) {
        model = paste("scale(",name[a],')~scale(',name[i],')+
                    scale(',pcs[1],')+scale(',pcs[2],')+scale(',pcs[3],')+scale(',pcs[4],')+scale(',pcs[5],')+scale(',age,')+(1|SITE_ID)+scale(',sex,')')
        mod=lmer(model,x)
        summ = data.frame(lmer.interpret(mod))
        peas = append(peas, summ[2,4])
        bees = append(bees, round(summ[2,1], 3))
        sees = append(sees, round(summ[2,2], 3))
        if (nobs(mod)<nobs) {
          nobs=nobs(mod)
        } else {
          next
        }
        print(nobs)
      }
    }
    
    for (i in c(1:ncol(pea.df))) {
      pea.df[,i]=peas[c((((i-1)*length(cbcl))+1):(i*length(cbcl)))]
      bee.df[,i]=bees[c((((i-1)*length(cbcl))+1):(i*length(cbcl)))]
      see.df[,i]=sees[c((((i-1)*length(cbcl))+1):(i*length(cbcl)))]
      colnames(pea.df)[i]=paste(ivs[i],'P',sep='.')
      colnames(bee.df)[i]=paste(ivs[i],'Beta',sep='.')
      colnames(see.df)[i]=paste(ivs[i],'SE',sep='.')
    }
    
    fhm.fdr.p <- p.adjust(peas, method = 'fdr')
    
  } else if (sex.interaction == T) {
    for (i in prs) {
      for (a in cbcl) {
        model = paste("scale(",name[a],')~scale(',name[i],')+
                    scale(',pcs[1],')+scale(',pcs[2],')+scale(',pcs[3],')+scale(',pcs[4],')+scale(',pcs[5],')+scale(',age,')+(1|SITE_ID)+scale(',sex,')+scale(',sex,'):scale(',name[i],')')
        mod=lmer(model,x)
        summ = data.frame(lmer.interpret(mod))
        peas = append(peas, summ[2,4])
        bees = append(bees, round(summ[2,1], 3))
        int.p = append(int.p, summ[10,4])
        int.b = append(int.b, summ[10,1])
        if (nobs(mod)<nobs) {
          nobs=nobs(mod)
        } else {
          next
        }
        print(nobs)
      }
    }
    
    for (i in c(1:ncol(pea.df))) {
      pea.df[,i]=int.p[c((((i-1)*length(cbcl))+1):(i*length(cbcl)))]
      bee.df[,i]=int.b[c((((i-1)*length(cbcl))+1):(i*length(cbcl)))]
      colnames(pea.df)[i]=paste(ivs[i],'P',sep='.')
      colnames(bee.df)[i]=paste(ivs[i],'Beta',sep='.')
    }
    
    fhm.fdr.p <- p.adjust(int.p, method = 'fdr')
  }
  
  
  fhm.fdr.mat <- matrix(fhm.fdr.p,ncol=ncol(pea.df),nrow=nrow(pea.df))
  
  fhm.fdr <- updatefdr(fhm.fdr.mat)
  
  pea.df.adj <- adjustmatrix.any(pea.df, breaks)
  pea.df.adj.dir <- direction(pea.df.adj, bee.df)
  
  rownames(pea.df.adj.dir) = dvs
  
  peasfhm = matrix(nrow = nrow(pea.df),ncol=ncol(pea.df))
  for (a in c(1:ncol(peasfhm))) {
    for (i in c(1:nrow(peasfhm))) {
      if (pea.df[i,a]>0.0001) {
        peasfhm[i,a] = sprintf("%.4f", pea.df[i,a])
      }
      else {
        peasfhm[i,a] = sprintf("%.1e", pea.df[i,a])
      }
    }
  }
  
  list <- list("adjusted" = pea.df.adj.dir, "fdr" = fhm.fdr, 'raw'=df, 'peas' = pea.df, 'breaks' = breaks, 'betas' = bee.df, 'min.nobs'=nobs, 'sees' = see.df)
  return(list)
}

# -----------------------------------------
# This function below maps the results from the output of the above function. Args:
# 1. fhm: output from cbclheatmap function
# 2. col.names defaults to setting the column names as the names of your IVs. You can pass a list of reader friendly names instead of length(ivs)
# 3. row.names defaults setting the row names as the names of your DVs. You can pass a list of reader friendly names instead of length(dvs)
# 4. title: character vector that will appear on the heatmap as the title
# 5. gaps: just adds gaps to delineate disorder-specific from cross-disorder PGS as well as categories of CBCL subscales

fhmtopheatmap <- function(fhm, col.names = NULL, row.names = NULL, title, gaps=T) {
  
  if (!is.null(col.names)) {
    colnames(fhm$adjusted) = col.names  
  }
  
  if (!is.null(row.names)) {
    rownames(fhm$adjusted) = row.names  
  }
  
  
  #c('darkblue','turquoise4','cornflowerblue','darkturquoise','turquoise','turquoise2','turquoise1','snow','wheat','lightgoldenrod1','orange1','darkorange2','brown2','red3','red4')
  
  # ** Depending on your data/analyses, you might get some funky color schemes/choices. If your heatmap looks wonky (i.e., the color scale is shifted so that the colors don't represent the significance correctly), you might be having some trouble with this part below. I tried to make it as generalizable as I could. **
  
  legendbreaks = c(range(fhm$adjusted)[1]:range(fhm$adjusted)[2])
  #legendbreaks[which(legendbreaks>0)]
  if (min(legendbreaks)==-2) {
    min='neg'
    legendlabels = c('<0.005','<0.05','>0.05', rep(NA, length(fhm$breaks)))
    start = 4
  } else if (min(legendbreaks)==-1) {
    min='neg'
    legendlabels = c('<0.05','>0.05', rep(NA, length(fhm$breaks)))
    start = 3
  } else if (min(legendbreaks)==0) {
    min='zero'
    legendlabels = c('>0.05', rep(NA, length(fhm$breaks)))
    start = 2
  }
  
  
  for (i in c(start:length(legendlabels))) {
    legendlabels[i] = paste('<', fhm$breaks[i-(start-1)], sep='')
  }
  
  if (min == 'neg') {
    legendlabels = legendlabels[c(1,2,3:((start-1)+length(which(legendbreaks>0))))]
  } else if (min == 'zero') {
    legendlabels = legendlabels[c(1,2:(1+length(which(legendbreaks>0))))]
  }
  
  
  colors = c(rgb(0,1,1),rgb(0,1,1),'white',rgb(0.6, 0, 0), rgb(0.8, 0, 0), rgb(1, 0, 0), rgb(1, 0.3, 0), rgb(1, 0.6, 0), rgb(1, 0.8, 0), rgb(1, 1, 0))
  colorspond <- c(-2,-1,0,1,2,3,4,5,6,7)
  
  if (min == 'neg') {
    colors = colors[c((which(colorspond==min(legendbreaks))):(which(colorspond==max(legendbreaks))))]
  } else if (min == 'zero') {
    colors = colors[c(3:which(colors==colors[which(colorspond==max(legendbreaks))]))]
  } 
  
  
  if (gaps == T) {
    map <- pheatmap(fhm$adjusted, cluster_rows = FALSE, cluster_cols = FALSE,   
                    color = colors,
                    legend_breaks = legendbreaks,
                    legend_labels = legendlabels,
                    fontsize_col = 10,
                    fontsize_row = 10,
                    fontsize = 10,
                    display_numbers = fhm$fdr,
                    fontsize_number = 20,
                    number_color = "black",
                    main = title,
                    gaps_row = c(8,11),
                    gaps_col = c(8,9))
  } else {
    map <- pheatmap(fhm$adjusted, cluster_rows = FALSE, cluster_cols = FALSE,   
                    color = colors,
                    legend_breaks = legendbreaks,
                    legend_labels = legendlabels,
                    fontsize_col = 10,
                    fontsize_row = 10,
                    fontsize = 10,
                    display_numbers = fhm$fdr,
                    fontsize_number = 20,
                    number_color = "black",
                    main = title)
  }
  return(map)
}

#---------------------------------------------------------------
## Define your dependent variables
dvs = c('BL.CBCL_ANXDEP_T','BL.CBCL_WITHDEP_T','BL.CBCL_SOMATIC_T','BL.CBCL_SOCIAL_T','BL.CBCL_THOUGHT_T','BL.CBCL_ATTENTION_T','BL.CBCL_RULEBREAK_T','BL.CBCL_AGGRESSIVE_T','BL.CBCL_INTERNAL_T','BL.CBCL_EXTERNAL_T','BL.CBCL_TOTPROB_T','BL.DSCORE')

# Define your independent variables
ivs = c('z.an','z.ocd','z.ts','z.adhd','z.asd','z.mdd','z.bip','z.scz.new','z.cross','z.compuls','z.moodpsych','z.neurodev')

# Define age variable 
age='INTERVIEW_AGE'

# Define sex variable 
sex = 'Gender.0M.1F'

# Define principal component variables
pcs=c('PC1.3','PC2.3','PC3.3','PC4.3','PC5.3')

# Define your breaks in the legend. These correspond to p-value ranges (i.e., 0.05 means p=0.05) 
breaks = c(0.05, 0.005, 5e-4, 5e-5, 5e-6, 5e-7, 5e-8)


#---- Run models to generate effect sizes and significance values ----
bl.fhm.lme4 <- cbclheatmap(ivs=ivs, dvs=dvs, x=master, age=age, sex=sex, pcs=pcs, breaks = breaks)
# -----------------------------------------

# Set rownames of beta matrix
rownames(bl.fhm.lme4$betas) = dvs

# Set rownames of FDR adjusted matrix
rownames(bl.fhm.lme4$fdr) = dvs

# Set rownames of raw p-value matrix
rownames(bl.fhm.lme4$peas) = dvs

# Set rownmaes of standard error matrix
rownames(bl.fhm.lme4$sees) = dvs

#---------------------------------
##------ Set mapping parameters

# Set readable column names (e.g., abbreviations for PGS)
fhmcolnames =  c('AN','OCD','TS','ADHD','ASD','MDD','BIP','SZ','CROSS','COMP','MP','NDV')

# Set readable rownames (e.g., CBCL subscales + PQB)
fhmrownames = c('Anx/Dep','With/Dep','Somatic','Social','Thought','Attention','Rulebreak','Aggressive','Internal','External','Total','Psychosis Spectrum')

## Run mapping function
fhmtopheatmap(fhm=bl.fhm.lme4, col.names = fhmcolnames, row.names = fhmrownames, title = 'Baseline psychopathology regressed on PGS', gaps = T)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@