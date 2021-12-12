cat("\014")  

pkg <- c("MASS", "ISLR", "class", "data.table", "caret", "Metrics","tidyverse", 
         "tidymodels", "mltools", "nortest", "devtools", "ggbiplot", "glmnet")
new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
if (length(new.pkg)) {install.packages(new.pkg, dependencies = TRUE)}
sapply(pkg, require, character.only = TRUE)

####### ####### ####### ####### ####### ####### #######
################ IMPORTING LIBRARIES  #################
####### ####### ####### ####### ####### ####### #######

library(tidyverse)
library(tidymodels)
library(caret)
library(Metrics)
library(mltools)
library(nortest)
library(usdm)
library(devtools)
install_github("vqv/ggbiplot")
library(ggbiplot)
library(glmnet)


####### ####### ####### ####### ####### ####### ####### #######
################ DEFINING ALL THE FUNCTIONS  ##################
####### ####### ####### ####### ####### ####### ####### #######

########### Co-Linearity ###########
computing_colinearity <- function(x) {
  colinearity <- vifcor(x)
  colinearity
}

########### Correlation Matrix ###########
computing_correlation_matrix <- function(x, t) {
  cor_mat <- cor(t)
  index <- findCorrelation(cor_mat, .75)
  to_be_removed <- colnames(cor_mat)[index] 
  x <- x[!names(x) %in% to_be_removed] 
  #print("After computing correlation")
  return(x)
}

########### Lasso Regression ###########
# Fit the LASSO model (Lasso: Alpha = 1)
plotting_lasso <- function(x, y, t){
  cv.lasso <- cv.glmnet(x, y, family='binomial', alpha=1, parallel=TRUE, standardize=TRUE, type.measure='auc')
  # Results
  #plot(cv.lasso)
  #plot(cv.lasso$glmnet.fit, xvar="lambda", label=TRUE)
  cat('Min Lambda: ', cv.lasso$lambda.min, '\n 1Sd Lambda: ', cv.lasso$lambda.1se)
  df_coef <- round(as.matrix(coef(cv.lasso, s=cv.lasso$lambda.min)), 2)
  # See all contributing variables
  df_coef <- df_coef[df_coef[, 1] == 0,]
  #names(t)
  t <- t[!names(t) %in% names(df_coef)] 
  print("After lasso feature selection")
  return(t)
}

########### Training the model and printing evaluation score ###########

ctrl <- trainControl(method="cv", number = 10) 
######### FITTING K NEAREST NEIGHBOUR #########
fitting_knn_model <- function(x){
  knn_fit <- train(Labels ~ ., data = x, method="knn", trControl = ctrl, 
                   preProcess = c("center", "scale"))
  return(knn_fit)
}

######### FITTING LINEAR DISCRIMINATIVE ANALYSIS #########
fitting_lda_model <- function(x){
  lda_fit <- train(Labels ~ ., data = x, method="lda", trControl = ctrl, 
                   preProcess = c("center", "scale"))
  return(lda_fit)
}

######### FITTING LOGISTIC REGRESSION #########
fitting_lr_model <- function(x){
  lr_fit <- train(Labels ~ ., data = x, method="glm", trControl = ctrl, 
                  preProcess = c("center", "scale"))
  return(lr_fit)
}
######### FITTING SUPPORT VECTOR MACHINE #########
fitting_svm_model <- function(x){
  svm_fit <- train(Labels ~ ., data = x, method="svmLinear", 
                   trControl = ctrl, preProcess = c("center", "scale"))
  return(svm_fit)
}

######### FITTING QUADRATIC DISCRIMINATIVE ANALYSIS #########
fitting_qda_model <- function(x){
  qda_fit <- train(Labels ~ ., data = x, method="qda", trControl = ctrl, 
                   preProcess = c("center", "scale"))
  return(qda_fit)
}

######### PREDICTING ON A DATASET #########
prediciton_using_model <- function(x_fit, x){
  model_result <- predict(x_fit, x)
  return(model_result)
}

######### COMPUTING AUPR SCORE #########
computing_aupr_score <- function(x, x_pred){
  print("auPR score")
  auc_score <- auc(x$Labels, x_pred)
  print(auc_score)
}
######### COMPUTING MCC SCORE #########
computing_mcc_score <- function(x, x_pred){
  print("MCC score")
  mccr_score <- mcc(x$Labels, x_pred)
  print(mccr_score)
}


####### ####### ####### ####### ####### ####### #######
############ IMPORTING THE ADCN DATASET ###############
####### ####### ####### ####### ####### ####### #######
set.seed(400)

print("Working on the ADCN DATASET....")

# Import the data 
train_full_adcn <- read.csv("ADCNtrain.csv")
test_data_adcn <- read.csv("ADCNtest.csv")
#dim(train_full_adcn)

######################## PRE PROCESSING #######################

########### Preparing the Dataset ###########
best_model_pred_adcn <- subset(test_data_adcn, select = c(1))
train_full_adcn <- subset(train_full_adcn, select = -c(1))
train_data_adcn <- subset(train_full_adcn, select = -c(Labels))
test_data_adcn <- subset(test_data_adcn, select  = -c(1))

#names(train_data_adcn)
########### Encoding Label Column ###########
train_full_adcn$Labels <- ifelse(train_full_adcn$Labels == "AD",1,0)
train_full_adcn$Labels <- as.factor(train_full_adcn$Labels)
#train_full_adcn$Labels

###################### FEATURE SELECTION ######################

########### computing colinearity and removing problematic features from ADCN ###########
# computing_colinearity(train_data_adcn) 
train_full_adcn <- within(train_full_adcn, rm("SLC6A8", "SELENBP1", "SLC6A10P....SLC6A8....SLC6A10PB*","TRIM58","TSPAN5","FECH","SLC25A39","PIM1", "GMPR", "S_Anterior_Rostral.1.R", "KLF1", "RNF10", "PIP4K2A", "STRADB", "DMTN", "BCL2L1", "N_Caudate.3.R", "G_subcallosal.1.R", "TNS1", "PITHD1", "N_Putamen.2.R", "TMOD1", "SNCA", "G_Paracentral_Lobule.3.R", "DCAF12", "GLRX5", "G_Fusiform.3.R", "S_Sup_Frontal.1.R","YBX3", "IFIT1B", "G_Hippocampus.2.R", "S_Sup_Frontal.2.L", "N_Thalamus.9.R" ))
#dim(train_full_adcn)

########### computing correlation matrix and removing correlated features from ADCN ###########
train_full_adcn <- computing_correlation_matrix(train_full_adcn, train_data_adcn)
#dim(train_full_adcn)

########### studying the principal components of ADCN features ###########
# pca_adcn <- prcomp(train_data_adcn, scale. = TRUE)
# str(pca_adcn)
# pca_adcn$rotation
# ggbiplot(pca_adcn, choices=c(1,2), groups = c("AD", "CN"), ellipse = TRUE) 
# ggbiplot(pca_adcn, choices=c(1,3), groups = c("AD", "CN"), ellipse = TRUE) 
# ggbiplot(pca_adcn, choices=c(1,4), groups = c("AD", "CN"), ellipse = TRUE) 
# ggbiplot(pca_adcn, choices=c(2,3), groups = c("AD", "CN"), ellipse = TRUE) 

# Lasso Regression
train_full_adcn <- plotting_lasso(as.matrix(train_data_adcn), train_full_adcn$Labels, train_full_adcn)

print("Final dimension of ADCN")
dim(train_full_adcn)

X_full <- train_full_adcn

# ########### K Nearest Neighbour ###########
# knn_fit_adcn <- fitting_knn_model(X_full)
# knn_fit_adcn
# knn_result_train_adcn <- prediciton_using_model(knn_fit_adcn, X_full)
# computing_mcc_score(X_full, knn_result_train_adcn)
# computing_aupr_score(X_full, knn_result_train_adcn)
# 
# ########### LDA ###########
# lda_fit_adcn <- fitting_lda_model(X_full)
# lda_fit_adcn
# lda_result_train_adcn <- prediciton_using_model(lda_fit_adcn, X_full)
# computing_mcc_score(X_full, lda_result_train_adcn)
# computing_aupr_score(X_full, lda_result_train_adcn)
# 
# ########### LR ###########
# lr_fit_adcn <- fitting_lr_model(X_full)
# lr_fit_adcn
# lr_result_train_adcn <- prediciton_using_model(lr_fit_adcn, X_full)
# computing_mcc_score(X_full, lr_result_train_adcn)
# computing_aupr_score(X_full, lr_result_train_adcn)

########### SVM ###########
svm_fit_adcn <- fitting_svm_model(X_full)
svm_fit_adcn
svm_result_train_adcn <- prediciton_using_model(svm_fit_adcn, X_full)
#svm_result_train_adcn
computing_mcc_score(X_full, svm_result_train_adcn)
computing_aupr_score(X_full, svm_result_train_adcn)

########### Predicting on ADCN DATASET ###########
best_model_test_adcn <- prediciton_using_model(svm_fit_adcn, test_data_adcn)
temp_adcn <- ifelse(best_model_test_adcn == "1","AD","CN")
best_model_pred_adcn$Labels <- temp_adcn

best_model_features_adcn <- X_full[, !names(X_full) %in% c("Labels")] 
best_model_features_adcn <- names(best_model_features_adcn)
#best_model_features_adcn
########### Saving variables in RData format ###########
save(best_model_pred_adcn, file = "0063769_Adilina_challenge2_ADCNres.RData")
save(best_model_features_adcn, file = "0063769_Adilina_challenge2_ADCNfeat.RData")


####### ####### ####### ####### ####### ####### #######
########### IMPORTING THE ADMCI DATASET ###############
####### ####### ####### ####### ####### ####### #######
set.seed(400)
print("Working on the ADMCI DATASET....")

# Import the data 
train_full_admci <- read.csv("ADMCItrain.csv")
test_data_admci <- read.csv("ADMCItest.csv")
#dim(train_full_admci)

######################## PRE PROCESSING #######################

########### Preparing the Dataset ###########
best_model_pred_admci <- subset(test_data_admci, select = c(1))
train_full_admci <- subset(train_full_admci, select = -c(1))
train_data_admci <- subset(train_full_admci, select = -c(Labels))
test_data_admci <- subset(test_data_admci, select  = -c(1))
#dim(train_data_admci)

########### Encoding Label Column ###########
train_full_admci$Labels <- ifelse(train_full_admci$Labels == "AD",1,0)
train_full_admci$Labels <- as.factor(train_full_admci$Labels)
#train_full_admci$Labels

###################### FEATURE SELECTION ######################

########### computing colinearity and removing problematic features from ADMCI ###########
# computing_colinearity(train_data_admci) 
train_full_admci <- within(train_full_admci, rm("SLC6A8","UBXN6","GMPR","STRADB","TSPAN5","SLC6A10P....SLC6A8....SLC6A10PB.","TMOD1","SNCA","N_Thalamus.9.R","RNF10","SELENBP1","FECH","PIM1","G_Paracentral_Lobule.3.R","G_Frontal_Med_Orb.1.R","SLC25A39","TRIM58","G_subcallosal.1.R","GUK1","GLRX5","DMTN","DCAF12","TNS1","G_Frontal_Sup_Orb.1.R","N_Putamen.2.R","N_Caudate.3.R","PIP4K2A","S_Anterior_Rostral.1.R","S_Orbital.1.R","S_Olfactory.1.L","G_Frontal_Sup_Orb.1.L","S_Sup_Frontal.2.L","PITHD1","TESC","YBX3","G_Fusiform.1.R","G_Fusiform.3.R","BCL2L1","G_ParaHippocampal.1.L","DPM2"))   
dim(train_full_admci)
########### computing correlation matrix and removing correlated features from ADMCI ###########
train_full_admci <- computing_correlation_matrix(train_full_admci, train_data_admci)
#dim(train_full_admci)

########### studying the principal components of ADCN features ###########
# pca_admci <- prcomp(train_data_admci, scale. = TRUE)
# str(pca_admci)
# pca_admci$rotation
# ggbiplot(pca_admci, choices=c(1,2), groups = c("AD", "MCI"), ellipse = TRUE)
# ggbiplot(pca_admci, choices=c(1,3), groups = c("AD", "MCI"), ellipse = TRUE)
# ggbiplot(pca_admci, choices=c(1,4), groups = c("AD", "MCI"), ellipse = TRUE)
# ggbiplot(pca_admci, choices=c(2,3), groups = c("AD", "MCI"), ellipse = TRUE)

# Lasso Regression
train_full_admci <- plotting_lasso(as.matrix(train_data_admci), train_full_admci$Labels, train_full_admci)
print("Final dimension of ADMCI")
dim(train_full_admci)

X_full <- train_full_admci

# ########### K Nearest Neighbour ###########
# knn_fit_admci <- fitting_knn_model(X_full)
# knn_fit_admci
# knn_result_train_admci <- prediciton_using_model(knn_fit_admci, X_full)
# computing_mcc_score(X_full, knn_result_train_admci)
# computing_aupr_score(X_full, knn_result_train_admci)

########### LDA ###########
lda_fit_admci <- fitting_lda_model(X_full)
lda_fit_admci
lda_result_train_admci <- prediciton_using_model(lda_fit_admci, X_full)
computing_mcc_score(X_full, lda_result_train_admci)
computing_aupr_score(X_full, lda_result_train_admci)

# ########### LR ###########
# lr_fit_admci <- fitting_lr_model(X_full)
# lr_fit_admci
# lr_result_train_admci <- prediciton_using_model(lr_fit_admci, X_full)
# computing_mcc_score(X_full, lr_result_train_admci)
# computing_aupr_score(X_full, lr_result_train_admci)
# 
# ########### SVM ###########
# svm_fit_admci <- fitting_svm_model(X_full)
# svm_fit_admci
# svm_result_train_admci <- prediciton_using_model(svm_fit_admci, X_full)
# #svm_result_train_adcn
# computing_mcc_score(X_full, svm_result_train_admci)
# computing_aupr_score(X_full, svm_result_train_admci)

########### Predicting on ADMCI DATASET ###########
best_model_test_admci <- prediciton_using_model(lda_fit_admci, test_data_admci)
temp_admci <- ifelse(best_model_test_admci == "1","AD","MCI")
best_model_pred_admci$Labels <- temp_admci

best_model_features_admci <- X_full[, !names(X_full) %in% c("Labels")] 
best_model_features_admci <- names(best_model_features_admci)
#best_model_features_admci
########### Saving variables in RData format ###########
save(best_model_pred_admci, file = "0063769_Adilina_challenge2_ADMCIres.RData")
save(best_model_features_admci, file = "0063769_Adilina_challenge2_ADMCIfeat.RData")


####### ####### ####### ####### ####### ####### #######
########### IMPORTING THE MCICN DATASET ###############
####### ####### ####### ####### ####### ####### #######
set.seed(400)
print("Working on the MCICN DATASET....")

###########  Import the data ###########
train_full_mcicn <- read.csv("MCICNtrain.csv")
test_data_mcicn <- read.csv("MCICNtest.csv")
#dim(train_full_mcicn)

######################## PRE PROCESSING #######################

########### Preparing the Dataset ###########
best_model_pred_mcicn <- subset(test_data_mcicn, select = c(1))
train_full_mcicn <- subset(train_full_mcicn, select = -c(1))
train_data_mcicn <- subset(train_full_mcicn, select = -c(Labels))
test_data_mcicn <- subset(test_data_mcicn, select  = -c(1))
#dim(train_data_mcicn)

########### Encoding Label Column ###########
train_full_mcicn$Labels <- ifelse(train_full_mcicn$Labels == "MCI",1,0)
train_full_mcicn$Labels <- as.factor(train_full_mcicn$Labels)
#dim(train_data_mcicn)

###################### FEATURE SELECTION ######################

########### computing colinearity and removing problematic features from MCICN ###########
# computing_colinearity(train_data_mcicn) 
train_full_mcicn <- within(train_full_mcicn, rm("Amygdala_R", "IGHG2", "ENSG00000211896....ENSG00000211897....ENSG00000233855", "ENSG00000211896", "ENSG00000211890....ENSG00000211895", "IGHV4.31", "ENSG00000211893....ENSG00000211896....ENSG00000211897....ENSG00", "ENSG00000211890", "KDM5D", "XIST", "IGHG3", "EIF1AY", "ORM1", "ENSG00000211895", "KIR2DS2....LOC100996743....KIR2DL3....KIR2DS1....KIR2DS4", "DDX3Y", "APOBEC3B", "SLC6A8", "IGHA1", "TMEM176B", "HBG1....HBG2", "LAIR2", "PRKY", "ENSG00000231486....ENSG00000239975....ENSG00000242076", "EPB42", "LCN2", "KLRC1....KLRC2", "ENSG00000211625....ENSG00000239951", "RSAD2", "IFIT1", "Cerebelum_9_L", "IFI44", "Frontal_Med_Orb_R", "HBD", "GMPR", "Cerebelum_8_R", "FAM46C", "Cerebelum_7b_L", "Rectus_L", "RPS4Y1", "BPGM", "Cerebelum_6_L", "Frontal_Sup_2_R", "SLC4A1", "GSTM1....GSTM2....GSTM4....GSTM5....GSTM2P1", "ISG15", "ANK1", "HERC5", "SELENBP1", "Hippocampus_R", "OAS3", "Putamen_L", "Lingual_L", "Insula_R", "TMOD1", "Precuneus_R", "NEK7", "HLA.DQB1"))
#dim(train_full_mcicn)

########### computing correlation matrix and removing correlated features from MCICN ###########
train_full_mcicn <- computing_correlation_matrix(train_full_mcicn, train_data_mcicn)
#dim(train_full_mcicn)

########### studying the principal components of ADCN features ###########
# pca_mcicn <- prcomp(train_data_mcicn, scale. = TRUE)
# str(pca_mcicn)
# pca_mcicn$rotation
# ggbiplot(pca_mcicn, choices=c(1,2), groups = c("CN", "MCI"), ellipse = TRUE)
# ggbiplot(pca_mcicn, choices=c(1,3), groups = c("CN", "MCI"), ellipse = TRUE)
# ggbiplot(pca_mcicn, choices=c(1,4), groups = c("CN", "MCI"), ellipse = TRUE)
# ggbiplot(pca_mcicn, choices=c(2,3), groups = c("CN", "MCI"), ellipse = TRUE)

# Lasso Regression
# train_full_mcicn <- plotting_lasso(as.matrix(train_data_mcicn), train_full_mcicn$Labels, train_full_mcicn)

print("Final dimension of MCICN")
dim(train_full_mcicn)

X_full <- train_full_mcicn

# ########### K Nearest Neighbour ###########
# knn_fit_mcicn <- fitting_knn_model(X_full)
# knn_fit_mcicn
# knn_result_train_mcicn <- prediciton_using_model(knn_fit_mcicn, X_full)
# computing_mcc_score(X_full, knn_result_train_mcicn)
# computing_aupr_score(X_full, knn_result_train_mcicn)
# 
# ########### LDA ###########
# lda_fit_mcicn <- fitting_lda_model(X_full)
# lda_fit_mcicn
# lda_result_train_mcicn <- prediciton_using_model(lda_fit_mcicn, X_full)
# computing_mcc_score(X_full, lda_result_train_mcicn)
# computing_aupr_score(X_full, lda_result_train_mcicn)
# 
# ########### LR ###########
# lr_fit_mcicn <- fitting_lr_model(X_full)
# lr_fit_mcicn
# lr_result_train_mcicn <- prediciton_using_model(lr_fit_mcicn, X_full)
# computing_mcc_score(X_full, lr_result_train_mcicn)
# computing_aupr_score(X_full, lr_result_train_mcicn)

########### SVM ###########
svm_fit_mcicn <- fitting_svm_model(X_full)
svm_fit_mcicn
svm_result_train_mcicn <- prediciton_using_model(svm_fit_mcicn, X_full)
#svm_result_train_mcicn
computing_mcc_score(X_full, svm_result_train_mcicn)
computing_aupr_score(X_full, svm_result_train_mcicn)

########### Predicting on ADMCI DATASET ###########
best_model_test_mcicn <- prediciton_using_model(svm_fit_mcicn, test_data_mcicn)
temp_mcicn <- ifelse(best_model_test_mcicn == "1","MCI", "CN")
best_model_pred_mcicn$Labels <- temp_mcicn

best_model_features_mcicn <- X_full[, !names(X_full) %in% c("Labels")] 
best_model_features_mcicn <- names(best_model_features_mcicn)
#best_model_features_mcicn
########### Saving variables in RData format ###########
save(best_model_pred_mcicn, file = "0063769_Adilina_challenge2_MCICNres.RData")
save(best_model_features_mcicn, file = "0063769_Adilina_challenge2_MCICNfeat.RData")

