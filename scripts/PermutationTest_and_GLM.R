### Load Libraries
library(ggplot2)
library(corrplot)
require(tidyverse)
require(rcompanion)


### Load Data 
infant_dat <- read.csv("~/table/Lou2021_Persister_Metadata_for_TermPreterm_Regression_Analysis.csv")


### Create Functions for Association Analysis and Permutation Testing

# Function 1 - permutation.test
permutation.test <- function(treatment, outcome, n){
  
  distribution <- c()
  result <- 0
  
  for(i in 1:n){
    distribution[i] <- diff(by(outcome, sample(treatment, length(treatment), FALSE), mean))
  }
  
  result <- sum(abs(distribution) >= abs(original))/(n)
  
  return(list(result, distribution))
}

# Function 2 - mixed.assoc
# Calculate a pairwise association between all variables in a data-frame. In particular nominal vs nominal with Chi-square, numeric vs numeric with Pearson correlation, and nominal vs numeric with ANOVA.
# Adopted from https://stackoverflow.com/a/52557631/590437
mixed.assoc <- function(df, cor_method="spearman", adjust_cramersv_bias=TRUE){
  df_comb <- expand.grid(names(df), names(df),  stringsAsFactors = F) %>% set_names("X1", "X2")
  
  is_nominal <- function(x) class(x) %in% c("factor", "character")
  is_numeric <- function(x) { is.integer(x) || is_double(x)}
  
  f = function(xName,yName) {
    x =  pull(df, xName)
    y =  pull(df, yName)
    
    result = if(is_nominal(x) && is_nominal(y)){
      # use bias corrected cramersV as described in https://rdrr.io/cran/rcompanion/man/cramerV.html
      cv = cramerV(as.character(x), as.character(y), bias.correct = adjust_cramersv_bias)
      data.frame(xName, yName, assoc=cv, type="cramersV")
      
    }else if(is_numeric(x) && is_numeric(y)){
      correlation = cor(x, y, method=cor_method, use="complete.obs")
      data.frame(xName, yName, assoc=correlation, type="correlation")
      
    }else if(is_numeric(x) && is_nominal(y)){
      # from https://stats.stackexchange.com/questions/119835/correlation-between-a-nominal-iv-and-a-continuous-dv-variable/124618#124618
      r_squared = summary(lm(x ~ y))$r.squared
      data.frame(xName, yName, assoc=sqrt(r_squared), type="anova")
      
    }else if(is_nominal(x) && is_numeric(y)){
      r_squared = summary(lm(y ~x))$r.squared
      data.frame(xName, yName, assoc=sqrt(r_squared), type="anova")
      
    }else {
      warning(paste("unmatched column type combination: ", class(x), class(y)))
    }
    
    # finally add complete obs number and ratio to table
    result %>% mutate(complete_obs_pairs=sum(!is.na(x) & !is.na(y)), complete_obs_ratio=complete_obs_pairs/length(x)) %>% rename(x=xName, y=yName)
  }
  
  # apply function to eachw variable combination
  map2_df(df_comb$X1, df_comb$X2, f)
}


### Step 1: Calculate Bootstrapped Confidence Intervals for Percent of Persisters by Term
pp_95ci <- groupwiseMean(perc_persisters ~ Term, data = infant_dat, R = 9999, boot = T, normal = T, bca = T)


### Step 2: Plot Percent Persisters Data with 95% Bootstraped CI
theme_set(theme_bw())
ggplot(infant_dat, aes(x = Term, y = perc_persisters)) +
  geom_jitter(shape = 21, width = 0.1, size = 3, fill = "steelblue3") +
  geom_segment(data = pp_95ci, aes(x = Term, xend = Term, y = Bca.lower, yend = Bca.upper)) +
  geom_point(data = pp_95ci, aes(x = Term, y = Mean), shape = 22, size = 4, fill = "firebrick") 


### Step 3: Run Permutation Test to Compare Persister Percentages in Term vs Pre-Term

## Calculate the difference in the means between the two groups
original <- diff(tapply(infant_dat$perc_persisters, infant_dat$Term, mean))

## Generate a large number mean differences calculated for randomly permuted groups
set.seed(123)
perm_out <- permutation.test(infant_dat$Term, infant_dat$perc_persisters, 9999)

perm_pval <- perm_out[[1]] ## This is the 2-sided permutation p-value

## This is a histogram showing actual mean difference vs permuted mean difference
hist(perm_out[[2]], breaks=100, col='grey', main="Permutation Distribution", las=1, xlab='')
abline(v=original, lwd=3, col="red")



### Step 4: Calculate Associations Between Variables Using Mixed Association Methods

## Output here is a table that gives each method used ect 
var_associations <- mixed.assoc(infant_dat[,4:ncol(infant_dat)])

## Reshape data into format for association matrix 
ass_mat <- reshape(var_associations[,1:3], direction="wide", idvar="x", timevar="y")
names <- ass_mat$x
ass_mat <- ass_mat[,-1] 
rownames(ass_mat) <- names
colnames(ass_mat) <- names

## Plot correlation matrix and clustered using hierarchical clustering 
corrplot(as.matrix(ass_mat), order = "hclust",
         hclust.method = "ward.D2", type = "upper")

## Plot correlation matrix and clustered using first principal component ordering
corrplot(as.matrix(ass_mat), order = "FPC", type = "upper")



### Step 5: Run GLM to assess multivariate impacts on number of persisters in each babies with highly correlated variables removed

## Modeled with Solid Food DOL in
m <- glm(num_persisters ~ Term + Race + Gender + Solid_Food_DOL + Cessation_BRM + Delivery + Feeds + Mid_Abx_Courses + Late_Abx_Courses + offset(log(num_early_colonizers)), family = 'poisson', data = infant_dat)
summary(m)
plot(m)
