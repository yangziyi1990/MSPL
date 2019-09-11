setwd("D:/Ziyi/School/PMO/12.Multi_omics/7.MSPL-code/1-Simulation")

##----------------------------------
##  Code Version 1.0
##  This simulation code for RandomForest embeded in the concatenation-based framework.
##  Concate_RF
##  Created by Zi-Yi Yang 
##  Modified by Zi-Yi Yang on July 23, 2019
##  Concact: yangziyi091100@163.com
##----------------------------------

##---------------------
## load libraries
##---------------------
library(Matrix)
library(tseries)
library(randomForest)
library(glmnet)
library(ROCR)

source("Function_simu_performance.R")


##--------------------------------
## 1. Generate simulation data
##--------------------------------
J <- 3                      # number of omics
p <- c(2000,500,1500)       # number of features in each omics
s <- c(10, 9, 8)            # number of nonzero coefficients
sam_num <- 1500
train_num <- c(100,150,200) # sample numbers of training
test_num <- 100             # sample numbers of testing
para_correlation <- c(0)    # c(0,0.2,0.4,0.6,0.8)
para_noise <- c(0,0.4,0.8)
S <- 3                      # index of train number  (100,150,200)
N <- 3                      # index of noise number (0,0.4,0.8)

## Initalize ##
beta.matrix <- list()
data.train <- list()
data.test <- list()
sample.train <- list()
sample.test <- list()
beta <- list()

y_train <- NULL
y_test <- NULL
label.matrix.train <- matrix(0, nrow = sam_num, ncol = J)
label.matrix.test <- matrix(0, nrow = sam_num, ncol = J)

for(i in 1:J){
  data.train[[i]] <- list()
  data.test[[i]] <- list()
  sample.train[[i]] <- list()
  sample.test[[i]] <- list()
  ## Initial beta
  idx <- sample(seq(from = 1, to = p[i], by = 1), size = s[i], replace = FALSE)
  beta[[i]] <- rep(0,p[i])
  beta[[i]][idx] <- runif(idx, min = -2.5, max = 2.5)
}

beta.matrix <- list("beta1" = beta[[1]],
                    "beta2" = beta[[2]],
                    "beta3" = beta[[3]])

## Inital pool of X ##
for(i in 1:J){
  sample.train[[i]] <- matrix(rnorm(sam_num * p[i], mean = 0, sd = 1), 
                              sam_num, p[i])
  sample.test[[i]] <- matrix(rnorm(sam_num * p[i], mean = 0, sd = 1), 
                             sam_num, p[i])
}

## Select sample with the sample label in all omics ##
for(i in 1:J){
  ## Pool of training
  noise_rnorm <- rnorm(sam_num,mean = 0, sd = 4)
  l.train <- sample.train[[i]] %*% beta.matrix[[i]] + para_noise[N] * noise_rnorm
  prob.train <- exp(l.train)/(1+exp(l.train))
  prob.train[which(prob.train>0.5)] <- 1
  prob.train[which(prob.train<=0.5)] <- 0
  label.matrix.train[,i] <- prob.train
  ## Pool of testing
  l.test <- sample.test[[i]] %*% beta.matrix[[i]]
  prob.test <- exp(l.test)/(1+exp(l.test))
  prob.test[which(prob.test>0.5)] <- 1
  prob.test[which(prob.test<=0.5)] <- 0
  label.matrix.test[,i] <- prob.test
  
}

sam0_train.idx <- which(apply(label.matrix.train, 1, sum)==0)
sam1_train.idx <- which(apply(label.matrix.train, 1, sum)==3)
sam0_test.idx <- which(apply(label.matrix.test, 1, sum)==0)
sam1_test.idx <- which(apply(label.matrix.test, 1, sum)==3)

select0.train.idx <- sample(x = sam0_train.idx, size = train_num[S]/2)
select1.train.idx <- sample(x = sam1_train.idx, size = train_num[S]/2)
select0.test.idx <- sample(x = sam0_test.idx, size = test_num/2)
select1.test.idx <- sample(x = sam1_test.idx, size = test_num/2)

select.train.idx <- c(select0.train.idx,select1.train.idx)
select.test.idx <- c(select0.test.idx,select1.test.idx)

##------------------------------
## generate simulation data
##------------------------------

for(i in 1:J){
  data.train[[i]] <- sample.train[[i]][select.train.idx,]
  y_train <- c(rep(0,train_num[S]/2),rep(1,train_num[S]/2)) 
  data.test[[i]] <- sample.test[[i]][select.test.idx,]
  y_test <- c(rep(0,test_num/2),rep(1,test_num/2)) 
}


##--------------------------------
## 2. Concatenation RandomForest
##--------------------------------
combined_train <- do.call(cbind, data.train)
combined_test <- do.call(cbind, data.test)

net <- randomForest(x = combined_train, y = factor(y_train), 
                    importance = TRUE, ntree = 30)
pred_train <- predict(net, combined_train)
pred_test <- predict(net, combined_test)

##---------------------
## 3. Performance
##---------------------

confusion.matrix.train <- table(observed = y_train, predicted = pred_train)
confusion.matrix.test <- table(observed = y_test, predicted = pred_test)

select.feature.idx <- which(net$importance[,4]!=0)  ## MeanDecreaseGini
select.feature.vaule <- net$importance[which(net$importance[,4]!=0),4] ## MeanDecreaseGini
select.fearure.name <- colnames(combined_train)[select.feature.idx]

pref.train <- evaluate.concatenationRF.performance(confusion.mat = confusion.matrix.train, 
                                                   true_lable = y_train, predict_label = pred_train)
pref.test <- evaluate.concatenationRF.performance(confusion.mat = confusion.matrix.test, 
                                                  true_lable = y_test, predict_label = pred_test)
pref.beta <- evaluate.beta.concatenationRF.performance(coef = select.feature.idx, beta = beta.matrix)

perf <- list("perf.train" = pref.train, "perf.test" = pref.test, "perf.beta" = pref.beta)

results <- list("net" = net, "feature.idx" = select.feature.idx, 
                "feature.value" = select.feature.vaule,
                "feature.name" = select.fearure.name)
