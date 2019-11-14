#### This file uses a Cox Proportional-Hazards model (no regularization or penalty) to predict patient survival (using survival analysis) using 165 drug gene signatures (i.e., RNA expression of corresponding genes in BRCA tumors)

## Output file contains six fields for each of the 165 drugs from the small molecule suite based on breast cancer subtypes of patient clinical data: ER+, HER2, Triple Negative
# a) Drug ID - index of the drug in HMS LINCS
# b) "test.cindex" - mean concordance index (as a measure of model performance) after prediction on test data after 10-fold cross-validation [main measure of model performance]
# c) "train.cindex" - mean concordance (as a measure of model performance) on training data after 10-fold cross-validation [this will always be higher than test.cindex]
# d) "testindex.sd" - standard deviation of concordance index after prediction on test data after 10-fold cross-validation
# e) "trainindex.sd" - standard deviation of concordance index after prediction on train data after 10-fold cross-validation
# f) "genecount" - number of genes used as features in the cox-ph regression model; also referred to as gene signature

# install.packages("caret", dependencies = c("Depends", "Suggests"))
# source("https://bioconductor.org/biocLite.R")
# biocLite("survcomp") # for concordance.index
require(data.table) # for fread
library(dplyr)
library(ggplot2)
library(magrittr)
library(tidyr)
library(survcomp) # for concordance.index function
library(survival) # for survival analysis
library(caret)
library(modelr) # crossv_kfold for cross validation
library(purrr) # for map, map2 functions
library(broom) # augment
library(randomForestSRC) # Random Survival Forest package for random forest approach for prediction
# library(ROCR)
# library(pROC)
# library(limma) # quantiles normalization
# library(glmnet) # for coxnet to regularize coxph model to select relevant covariates when data is sparse

## set directory to where generated data is stored (refer to file dr_generate_corefiles.r)
sourcedir <- './generateddata/'

## Load expression data, clinical characteristics (event - death, time - days), and drug-gene signatures
## These dataframes are loaded ->
# 1) brca_cancer - BRCA expression data with mRNA expression of 1097 cancer patient samples + 1 column for gene_IDs
# 2) brca_days - clinical variables (vital status, time in days) of 1097 patients for survival analysis
# 3) gene_array - contains list of genes stored as list object for 165 drugs
# 4) gene_index - contains the number associated with drug in HMS LINCS to link back during analysis
# 5) subtype_cancer - contains molecular subtype (ER+, HER2+, PR+, Triplenegative) for 1097 brca patients
load(paste0(sourcedir,'drugrepoData.RData'))

# SETUP MODEL DATA: transpose data table - samples as rows, genes as columns
brca_exp <- select(brca_cancer, -one_of('gene_name'))

# remove genes where mean expression of any gene is less that 5.  Doing this avoids issues with prediction later on as we will end up with genes that cannot be tested.
brca_low_rowsum <- which(rowMeans(brca_exp) < 5) # 3801 of the 20533 genes
brca_exp <- brca_exp[-brca_low_rowsum,]
brca_cancer <- brca_cancer[-brca_low_rowsum, ] # only 41 of the genes associated with any of the drugs are part of this 3801 genes

# setup model data to contain transpose of the brca data frame with patient samples as rows and genes in columns
model_days_new <- data.table( t(brca_exp) ) # transpose only expression values, drop gene_name column
colnames(model_days_new) <- brca_cancer$gene_name  # re-assign gene names to columns (lost during transposition)

# add clinical variables to the data table
model_days_new <- model_days_new %>% mutate(days=brca_days$days, class = brca_days$status)

## Remove all patients where days count is 0 as COXNET model does not handle 0 values
# model_days_new <- model_days_new %>% filter(days != 0)

rownames(model_days_new) <- colnames(brca_exp)  # re-assign patient IDs to data table

# Separate brca data by subtype for analysis
er_patients <- subtype_cancer %>% filter(er_status_by_ihc=="Positive") %>% select(bcr_patient_barcode)
# for Her2 patients include those that are categorized as "equivocal" too: https://www.cancer.org/cancer/breast-cancer/understanding-a-breast-cancer-diagnosis/breast-cancer-her2-status.html
her2_patients <- subtype_cancer %>% filter(her2_status_by_ihc=="Positive" | her2_status_by_ihc=="Equivocal") %>% select(bcr_patient_barcode)
tripneg_patients <- subtype_cancer %>% filter(mol_type=="TripleNegative") %>% select(bcr_patient_barcode)


## Setup function to get concordance index (performance metric) for random gene signatures
get_index <- function(model_matrix, num_folds){
  ## use crossv_kfold from modelr package to run n-fold cross validation
  samp_model <- crossv_kfold(model_matrix, k=num_folds, id="fold")
  
  ## train and test model on n-folds; extract average estimate of final c.index score
  # train - coxph
  samp_model <- samp_model %>% mutate(model = map(train, ~ rfsrc(Surv(as.numeric(days),class) ~ ., 
                                                                 data= as.data.frame(.), ntree = 100)))
  # train - save training model concordance
  samp_model <- samp_model %>% mutate(train.index = map(model, ~ 1-.x$err.rate[.x$ntree])) # concordance index is 1 - err.rate for random survival forest
  mean_train_cindex <- mean(unlist(samp_model$train.index), na.rm=T)
  sd_train_cindex <- sd(unlist(samp_model$train.index), na.rm=T)
  
  # test - predict
  samp_model <- samp_model %>% mutate(predicted = map2(model, test, ~ predict(.x, newdata=as.data.frame(.y))))
  # get concordance index based on predicted values
  samp_model <- samp_model %>% mutate(test.index = map(predicted, ~ 1-.x$err.rate[.x$ntree])) # concordance index is 1 - err.rate for random survival forest
  # get mean c.index
  mean_test_cindex <- mean(unlist(samp_model$test.index), na.rm=T)
  sd_test_cindex <- sd(unlist(samp_model$test.index), na.rm=T)
  indices <- list(train_index = mean_train_cindex, sd_train = sd_train_cindex,
                  test_index = mean_test_cindex, sd_test = sd_test_cindex)
  return(indices)
}


## Train model on subset of expression data and predict on testing data using gene signatures of 165 drugs
## Train model for each subtype in separate iterations and save the results in 3 separate files

subtype_objects <- c(er_patients, her2_patients, tripneg_patients)
names(subtype_objects) <- c("er", "her2", "tneg")
drug_count <- length(gene_array)

for (obj in 1:length(subtype_objects)){
  
  ## Final model data to use in training and prediction
  model_data_subtype <- model_days_new[which( substr(colnames(brca_cancer),1,12) %in% subtype_objects[[obj]] ),]
  
  print(paste("Subtype object",obj))
  
  #### For 165 drug gene signatures
  test_cindex <- vector()
  testindex_sd <- vector()
  train_cindex <- vector()
  trainindex_sd <- vector()
  drug_num <- vector()
  num_genes <- vector()
  
  ## Iterate over each of the 165 gene signatures corresponding to 165 drugs
  for (x in 1:drug_count){
    ## Set seed for reproducibility
    set.seed(x+523)
    
    ## Iterate over the 165 genes
    print(paste("Iteration number",x))
    gene_list <- gene_array[[x]] # get gene signature for drug x
    
    
    # SETUP PREDICTORS: select sample predictor genes based on target gene set for specific drugs
    model_genes <- grep(paste0("\\b",gene_list,"\\b",collapse="|"), colnames(model_data_subtype), ignore.case = T)
    model_matrix <- model_data_subtype %>% select(c(model_genes, "days", "class"))
    
    
    # remove gene columns where majority (>90%) of values are 0 - this pre-processing steps avoids warnings from cross validation if a large number of expression values are all 0
    # zero_cols <- which( colSums(model_matrix==0) > 0.9 * nrow(model_matrix) )
    # ifelse(length(zero_cols) > 0, model_matrix <- model_matrix[,-zero_cols], model_matrix <- model_matrix)
    predictor_matrix <- model_matrix
    
    
    ### Use 10 fold cross validation for training and predicting survival to get concordance index
    ## Occasionally this statement errors with the following issue: NA/NaN/Inf in foreign function call (arg 6)
    ## To prevent the error in a single iteration to break an entire run, we catch the error using try
    # use 10 fold cv for ER subtype, and 5 fold cv for her2 and tneg because of smaller sample sizes for latter
    try_err <- try (
      if (obj==1){model_cox <- get_index(predictor_matrix, 10)}
      else{model_cox <- get_index(predictor_matrix, 10)}
      , silent = T)
    
    # skip the current iteration if there is an error, else continues
    if (class(try_err) == "try-error"){
      print("error")
      next
    }
    
    ## Store train and test cindex
    test_cindex <- c(test_cindex, model_cox$test_index)
    train_cindex <- c(train_cindex, model_cox$train_index)
    testindex_sd <- c(testindex_sd, model_cox$sd_test)
    trainindex_sd <- c(trainindex_sd, model_cox$sd_train)
    drug_num <- c(drug_num, gene_index[x])
    num_genes <- c(num_genes, ncol(predictor_matrix)-2)
    
  }
  
  ## Save the three performance related variables to a data frame and write to file
  
  cox_drug_perf_subtype_df <- data.table(drug=drug_num, test.cindex=test_cindex, train.cindex=train_cindex, 
                                         testindex.sd = testindex_sd, trainindex.sd = trainindex_sd, genecount=num_genes)

  # write to file
  dest_folder <- './generateddata/cox_rsf/subtype/'
  cox_drug_subtype_file <- paste0(dest_folder,"cox_drug_results_subtype_",names(subtype_objects)[obj],"_rsf.txt")
  write.table(cox_drug_perf_subtype_df, file = cox_drug_subtype_file, row.names = F)
  
}



