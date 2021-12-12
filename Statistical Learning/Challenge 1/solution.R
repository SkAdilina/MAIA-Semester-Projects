
pkg <- c("MASS", "ISLR", "class", "data.table", "caret", "Metrics","tidyverse", "tidymodels")
new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
if (length(new.pkg)) {install.packages(new.pkg, dependencies = TRUE)}
sapply(pkg, require, character.only = TRUE)

####### ####### ####### ####### ####### ####### #######
########## IMPORTING LIBRARIES AND DATASETS ###########
####### ####### ####### ####### ####### ####### #######

library(tidyverse)
library(tidymodels)
library(caret)
library(Metrics)

set.seed(400)

# Import the data and look at the first six rows
train_full <- read.csv("train_ch.csv")
test_data <- read.csv("test_ch.csv")

########### Plotting the training and testing datasets ###########
#plot(train_full, main="Training Data")
#plot(test_data, main="Test Set")
#boxplot(train_full)$out
#boxplot(test_data)$out


####### ####### ####### ####### #######
############ PRE PROCESSING ###########
####### ####### ####### ####### #######

########### Delete outliers ###########
#boxplot(train_full)$out
outliers <- boxplot(train_full$Y, plot=FALSE)$out
train_full <- train_full[-which(train_full$Y %in% outliers),]
# print(is.data.frame(train_data))

########### Deleting the "X" column ###########
train_full <- subset(train_full, select = -c(1))
test_data <- subset(test_data, select  = -c(1))

########### Centering and Scaling ###########
preProcess(x=train_full, method = c("center", "scale"))
preProcess(x=test_data, method = c("center", "scale"))

########### Plotting the boxplot of the training and testing datasets ###########
#boxplot(train_full, main="Preprocessed Training Dataset")$out
#boxplot(test_data, main="Preprocessed Testing Dataset")$out

####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### #######

####### ####### ####### ####### #######
########## FEATURE SELECTION ##########
####### ####### ####### ####### #######

########### Correlation Matrix ###########
#correlation_matrix = cor(train_full, method = c("pearson", "kendall", "spearman"))
#correlation_matrix
#plot(correlation_matrix)

########### Forward Feature Selection ###########
train.control <- trainControl(method = "cv", number = 10)
step.model <- train(Y ~., data = train_full,
                    method = "leapForward",
                    trControl = train.control)

#step.model$results
#step.model$bestTune
summary(step.model$finalModel)

####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### #######

####### ####### ####### ####### #######
########## LINEAR REGRESSION ##########
####### ####### ####### ####### #######

########### Fitting the model ###########
ctrl <- trainControl(method="cv", number = 10) 
lm_fit <- lm(Y ~ v1*v2*v3*v4*v5*v6*v7*v8*v9, data = train_full, trControl = ctrl)
summary(lm_fit)

# Predicting on the target values
lmPredict <- predict(lm_fit, train_full)
lm_pred <- predict(lm_fit, test_data)

# Calculating RMSE
RMSE_lr = RMSE(lmPredict, train_full$Y)
RMSE_lr 


####### ####### ####### ####### #######
######### K NEAREST NEIGHBOUR #########
####### ####### ####### ####### #######

########### Fitting the model ###########
ctrl <- trainControl(method="cv", number = 10) 
knn_fit <- train(Y ~ v1+v2+v3, data = train_full, method="knn", trControl = ctrl)
knn_fit

# Predicting the target values
knnPredict <- predict(knn_fit, train_full)
knn_pred <- predict(knn_fit, test_data)

# Calculating RMSE
RMSE_knn_train = RMSE(knnPredict, train_full$Y)
RMSE_knn_train


####### ####### ####### ####### #######
######## SAVING THE VARIABLES #########
####### ####### ####### ####### #######

########### Saving variables in RData format ###########
save(knn_fit, knn_pred, lm_fit, lm_pred, file = "Adilina_challenge1.RData")

