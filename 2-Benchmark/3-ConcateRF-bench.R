##----------------------------------
##  Code Version 1.0
##  This real data code for Random Forest embeded in the concatenation-based framework.
##  Concate_RF
##  Created by Zi-Yi Yang 
##  Modified by Zi-Yi Yang on July 23, 2019
##  Concact: yangziyi091100@163.com
##----------------------------------

## load libraries
library(Matrix)
library(tseries)
library(randomForest)
library(glmnet)
library(ROCR)

setwd("D:/Ziyi/School/PMO/12.Multi_omics/7.MSPL-code/2-Benchmark")
source('Function_performance.R')

## Load data
load("D:/Ziyi/School/PMO/12.Multi_omics/7.MSPL-code/2-Benchmark/1-data/SNFdatasets.RDATA")
  
#----------------
## 1. Colon dataset
#----------------
# Preparing data
data_colon <- snf_data$colon
label_colon <- snf_group$colon

# random select samples and Setting training and testing
randidx_colon = sample(c(1:length(label_colon)),size=length(label_colon))
splits_colon = vector(mode = "numeric", length = length(label_colon))
splits_colon_trainidx = randidx_colon[1:(0.7*length(randidx_colon))]
splits_colon_testidx = randidx_colon[(0.7*length(randidx_colon)+1):length(randidx_colon)]

splits_colon[splits_colon_trainidx] = 0
splits_colon[splits_colon_testidx] = 1
splits_colon = as.matrix(splits_colon)

trainidx_colon = which(splits_colon[,1]==0)
testidx_colon = which(splits_colon[,1]==1)

train_colon_mrna<-data_colon$mrna[trainidx_colon,]
train_colon_mirna<-data_colon$mirna[trainidx_colon,]
train_colon_cpg<-data_colon$cpg[trainidx_colon,]

test_colon_mrna<-data_colon$mrna[testidx_colon,]
test_colon_mirna<-data_colon$mirna[testidx_colon,]
test_colon_cpg<-data_colon$cpg[testidx_colon,]

label_colon[which(label_colon=="high")] <- 1
label_colon[which(label_colon=="low")] <- 0

train_colon_label <- label_colon[trainidx_colon]
test_colon_label <- label_colon[testidx_colon]

## Generate Colon data
colon_train_data = list(mrna=train_colon_mrna, mirna=train_colon_mirna, cpg=train_colon_cpg)
colon_train_lable = train_colon_label
colon_test_data = list(mrna=test_colon_mrna, mirna=test_colon_mirna, cpg=test_colon_cpg)
colon_test_lable = test_colon_label

#----------------
## 1. GBM dataset
#----------------
# Preparing data
data_gbm <- snf_data$gbm
label_gbm <- snf_group$gbm

# random select samples and Setting training and testing
randidx_gbm = sample(c(1:length(label_gbm)),size=length(label_gbm))
splits_gbm = vector(mode = "numeric", length = length(label_gbm))
splits_gbm_trainidx = randidx_gbm[1:(0.7*length(randidx_gbm))]
splits_gbm_testidx = randidx_gbm[(0.7*length(randidx_gbm)+1):length(randidx_gbm)]

splits_gbm[splits_gbm_trainidx] = 0
splits_gbm[splits_gbm_testidx] = 1
splits_gbm = as.matrix(splits_gbm)

trainidx_gbm = which(splits_gbm[,1]==0)
testidx_gbm = which(splits_gbm[,1]==1)

train_gbm_mrna<-data_gbm$mrna[trainidx_gbm,]
train_gbm_mirna<-data_gbm$mirna[trainidx_gbm,]
train_gbm_cpg<-data_gbm$cpg[trainidx_gbm,]

test_gbm_mrna<-data_gbm$mrna[testidx_gbm,]
test_gbm_mirna<-data_gbm$mirna[testidx_gbm,]
test_gbm_cpg<-data_gbm$cpg[testidx_gbm,]

label_gbm[which(label_gbm=="high")] <- 1
label_gbm[which(label_gbm=="low")] <- 0

train_gbm_label <- label_gbm[trainidx_gbm]
test_gbm_label <- label_gbm[testidx_gbm]

## Generate gbm data
gbm_train_data = list(mrna=train_gbm_mrna, mirna=train_gbm_mirna, cpg=train_gbm_cpg)
gbm_train_lable = train_gbm_label
gbm_test_data = list(mrna=test_gbm_mrna, mirna=test_gbm_mirna, cpg=test_gbm_cpg)
gbm_test_lable = test_gbm_label

#----------------
## 1. Kidney dataset
#----------------
# Preparing data
data_kidney <- snf_data$kidney
label_kidney <- snf_group$kidney

# random select samples and Setting training and testing
randidx_kidney = sample(c(1:length(label_kidney)),size=length(label_kidney))
splits_kidney = vector(mode = "numeric", length = length(label_kidney))
splits_kidney_trainidx = randidx_kidney[1:(0.7*length(randidx_kidney))]
splits_kidney_testidx = randidx_kidney[(0.7*length(randidx_kidney)+1):length(randidx_kidney)]

splits_kidney[splits_kidney_trainidx] = 0
splits_kidney[splits_kidney_testidx] = 1
splits_kidney = as.matrix(splits_kidney)

trainidx_kidney = which(splits_kidney[,1]==0)
testidx_kidney = which(splits_kidney[,1]==1)

train_kidney_mrna<-data_kidney$mrna[trainidx_kidney,]
train_kidney_mirna<-data_kidney$mirna[trainidx_kidney,]
train_kidney_cpg<-data_kidney$cpg[trainidx_kidney,]

test_kidney_mrna<-data_kidney$mrna[testidx_kidney,]
test_kidney_mirna<-data_kidney$mirna[testidx_kidney,]
test_kidney_cpg<-data_kidney$cpg[testidx_kidney,]

label_kidney[which(label_kidney=="high")] <- 1
label_kidney[which(label_kidney=="low")] <- 0

train_kidney_label <- label_kidney[trainidx_kidney]
test_kidney_label <- label_kidney[testidx_kidney]

## Generate kidney data
kidney_train_data = list(mrna=train_kidney_mrna, mirna=train_kidney_mirna, cpg=train_kidney_cpg)
kidney_train_lable = train_kidney_label
kidney_test_data = list(mrna=test_kidney_mrna, mirna=test_kidney_mirna, cpg=test_kidney_cpg)
kidney_test_lable = test_kidney_label


#----------------
## 1. Lung dataset
#----------------
# Preparing data
data_lung <- snf_data$lung
label_lung <- snf_group$lung

# random select samples and Setting training and testing
randidx_lung = sample(c(1:length(label_lung)),size=length(label_lung))
splits_lung = vector(mode = "numeric", length = length(label_lung))
splits_lung_trainidx = randidx_lung[1:(0.7*length(randidx_lung))]
splits_lung_testidx = randidx_lung[(0.7*length(randidx_lung)+1):length(randidx_lung)]

splits_lung[splits_lung_trainidx] = 0
splits_lung[splits_lung_testidx] = 1
splits_lung = as.matrix(splits_lung)

trainidx_lung = which(splits_lung[,1]==0)
testidx_lung = which(splits_lung[,1]==1)

train_lung_mrna<-data_lung$mrna[trainidx_lung,]
train_lung_mirna<-data_lung$mirna[trainidx_lung,]
train_lung_cpg<-data_lung$cpg[trainidx_lung,]

test_lung_mrna<-data_lung$mrna[testidx_lung,]
test_lung_mirna<-data_lung$mirna[testidx_lung,]
test_lung_cpg<-data_lung$cpg[testidx_lung,]

label_lung[which(label_lung=="high")] <- 1
label_lung[which(label_lung=="low")] <- 0

train_lung_label <- label_lung[trainidx_lung]
test_lung_label <- label_lung[testidx_lung]

## Generate lung data
lung_train_data = list(mrna=train_lung_mrna, mirna=train_lung_mirna, cpg=train_lung_cpg)
lung_train_lable = train_lung_label
lung_test_data = list(mrna=test_lung_mrna, mirna=test_lung_mirna, cpg=test_lung_cpg)
lung_test_lable = test_lung_label


##---------------------
## 2. RandomForest
##---------------------
##-----2.1 colon dataset------
combined_train_colon <- do.call(cbind, colon_train_data)
combined_test_colon <- do.call(cbind, colon_test_data)

net.colon <- randomForest(x = combined_train_colon, y = factor(colon_train_lable), 
                          importance = TRUE, ntree = 10)
pred_train_colon<- predict(net.colon, combined_train_colon)
pred_test_colon<- predict(net.colon, combined_test_colon)

##-----2.2 gbm dataset-----
combined_train_gbm <- do.call(cbind, gbm_train_data)
combined_test_gbm <- do.call(cbind, gbm_test_data)

net.gbm <- randomForest(x = combined_train_gbm, y = factor(gbm_train_lable), 
                        importance = TRUE, ntree = 10)
pred_train_gbm<- predict(net.gbm, combined_train_gbm)
pred_test_gbm<- predict(net.gbm, combined_test_gbm)

##-----2.3 kidney dataset-----
combined_train_kidney <- do.call(cbind, kidney_train_data)
combined_test_kidney <- do.call(cbind, kidney_test_data)

net.kidney <- randomForest(x = combined_train_kidney, y = factor(kidney_train_lable), 
                           importance = TRUE, ntree = 10)
pred_train_kidney<- predict(net.kidney, combined_train_kidney)
pred_test_kidney<- predict(net.kidney, combined_test_kidney)


##-----2.4 lung dataset-----
combined_train_lung <- do.call(cbind, lung_train_data)
combined_test_lung <- do.call(cbind, lung_test_data)

net.lung <- randomForest(x = combined_train_lung, y = factor(lung_train_lable), 
                         importance = TRUE, ntree = 10)
pred_train_lung<- predict(net.lung, combined_train_lung)
pred_test_lung<- predict(net.lung, combined_test_lung)

##---------------------
## 3. Performance
##---------------------
##-----3.1 colon dataset------
confusion.matrix.train.colon <- table(observed=colon_train_lable,predicted=pred_train_colon)
confusion.matrix.test.colon <- table(observed=colon_test_lable,predicted=pred_test_colon)

select.feature.colon.idx <- which(net.colon$importance[,4]!=0)  ## MeanDecreaseGini
select.feature.colon.vaule <- net.colon$importance[which(net.colon$importance[,4]!=0),4] ## MeanDecreaseGini
select.fearure.colon.name <- colnames(combined_train_colon)[select.feature.colon.idx]

pref.train.colon <- evaluate.randforest.performance(confusion.mat = confusion.matrix.train.colon, 
                                                    true_lable = colon_train_lable, predict_label = pred_train_colon)
pref.test.colon <- evaluate.randforest.performance(confusion.mat = confusion.matrix.test.colon, 
                                                   true_lable = colon_test_lable, predict_label = pred_test_colon)

perf.colon <- list("pref.train" = pref.train.colon, "perf.test" = pref.test.colon)

results.colon <- list("net" = net.colon, "feature.idx" = select.feature.colon.idx, 
                      "feature.value" = select.feature.colon.vaule,
                      "feature.name" = select.fearure.colon.name)


##-----3.2 gbm dataset------
confusion.matrix.train.gbm <- table(observed=gbm_train_lable,predicted=pred_train_gbm)
confusion.matrix.test.gbm <- table(observed=gbm_test_lable,predicted=pred_test_gbm)

select.feature.gbm.idx <- which(net.gbm$importance[,4]!=0)  ## MeanDecreaseGini
select.feature.gbm.vaule <- net.gbm$importance[which(net.gbm$importance[,4]!=0),4] ## MeanDecreaseGini
select.fearure.gbm.name <- colnames(combined_train_gbm)[select.feature.gbm.idx]

pref.train.gbm <- evaluate.randforest.performance(confusion.mat = confusion.matrix.train.gbm, 
                                                  true_lable = gbm_train_lable, predict_label = pred_train_gbm)
pref.test.gbm <- evaluate.randforest.performance(confusion.mat = confusion.matrix.test.gbm, 
                                                 true_lable = gbm_test_lable, predict_label = pred_test_gbm)

perf.gbm <- list("pref.train" = pref.train.gbm, "perf.test" = pref.test.gbm)

results.gbm <- list("net" = net.gbm, "feature.idx" = select.feature.gbm.idx, 
                    "feature.value" = select.feature.gbm.vaule,
                    "feature.name" = select.fearure.gbm.name)

##-----3.3 kidney dataset------
confusion.matrix.train.kidney <- table(observed=kidney_train_lable,predicted=pred_train_kidney)
confusion.matrix.test.kidney <- table(observed=kidney_test_lable,predicted=pred_test_kidney)

select.feature.kidney.idx <- which(net.kidney$importance[,4]!=0)  ## MeanDecreaseGini
select.feature.kidney.vaule <- net.kidney$importance[which(net.kidney$importance[,4]!=0),4] ## MeanDecreaseGini
select.fearure.kidney.name <- colnames(combined_train_kidney)[select.feature.kidney.idx]

pref.train.kidney <- evaluate.randforest.performance(confusion.mat = confusion.matrix.train.kidney, 
                                                     true_lable = kidney_train_lable, predict_label = pred_train_kidney)
pref.test.kidney <- evaluate.randforest.performance(confusion.mat = confusion.matrix.test.kidney, 
                                                    true_lable = kidney_test_lable, predict_label = pred_test_kidney)

perf.kidney <- list("pref.train" = pref.train.kidney, "perf.test" = pref.test.kidney)

results.kidney <- list("net" = net.kidney, "feature.idx" = select.feature.kidney.idx, 
                       "feature.value" = select.feature.kidney.vaule,
                       "feature.name" = select.fearure.kidney.name)

##-----3.4 lung dataset------
confusion.matrix.train.lung <- table(observed=lung_train_lable,predicted=pred_train_lung)
confusion.matrix.test.lung <- table(observed=lung_test_lable,predicted=pred_test_lung)

select.feature.lung.idx <- which(net.lung$importance[,4]!=0)  ## MeanDecreaseGini
select.feature.lung.vaule <- net.lung$importance[which(net.lung$importance[,4]!=0),4] ## MeanDecreaseGini
select.fearure.lung.name <- colnames(combined_train_lung)[select.feature.lung.idx]

pref.train.lung <- evaluate.randforest.performance(confusion.mat = confusion.matrix.train.lung, 
                                                   true_lable = lung_train_lable, predict_label = pred_train_lung)
pref.test.lung <- evaluate.randforest.performance(confusion.mat = confusion.matrix.test.lung, 
                                                  true_lable = lung_test_lable, predict_label = pred_test_lung)

perf.lung <- list("pref.train" = pref.train.lung, "perf.test" = pref.test.lung)

results.lung <- list("net" = net.lung, "feature.idx" = select.feature.lung.idx, 
                     "feature.value" = select.feature.lung.vaule,
                     "feature.name" = select.fearure.lung.name)

##---------------------
## 4. Output
##---------------------
results.benchmark <- list("results.colon" = results.colon,
                          "results.gbm" = results.gbm,
                          "results.kidney" = results.kidney,
                          "results.lung" = results.lung)
performance.benchmark <- list("perf.colon" = perf.colon,
                              "perf.gbm" = perf.gbm,
                              "perf.kidney" = perf.kidney,
                              "perf.lung" = perf.lung)

