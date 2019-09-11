setwd("D:/Ziyi/School/PMO/12.Multi_omics/7.MSPL-code/1-Simulation")

##----------------------------------
##  Code Version 1.0
##  This simulation code for multimodal self-paced learning for multimodal data analysis.
##  MSPL
##  Created by Zi-Yi Yang
##  Modified by Zi-Yi Yang on July 23, 2019
##  Concact: yangziyi091100@163.com
##----------------------------------

##---------------------
## load libraries
##---------------------
library(Matrix)
library(tseries)
library(glmnet)
library(caret)
library(ncvreg)
library(pROC)
library(ROCR)
library(ggplot2)

source("MSPL_Function.R")

#--------------------------------
# 1. Generate simulation data
#--------------------------------
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

##-----------------
## Step 1 : Initialization parameters
##-----------------
View_num <- 3
iter_num <- 100
gamma <- 0.1
lambda <- c(0.06, 0.06, 0.06)
uplambda <- 0.04
num_add <- 4
num_up <- 2

valpredmatrix <- list()
evlpredmatrix <- list()
coefmatrix <- list()
nonzerocoefmatrix <- list()
coef_idx <- list()
coef_value <- list()
marker_record <- list()
selectedidx <- list()

loss <- matrix(0, nrow = length(y_train), ncol = View_num)
v_iter <- matrix(0, nrow = length(y_train), ncol = View_num)

for(iter in 1:iter_num) {
  valpredmatrix[[iter]] <- matrix(0, nrow = length(y_train), ncol = View_num)
  evlpredmatrix[[iter]] <- matrix(0, nrow = length(y_test), ncol = View_num)
  coefmatrix[[iter]] <-  list()
  nonzerocoefmatrix[[iter]] <- matrix(0, nrow = 1, ncol = View_num)
}

val_labels <- matrix(rep(y_train, each = 3),ncol = 3, byrow = T)
evl_labels <- matrix(rep(y_test, each = 3),ncol = 3, byrow = T)
valmaps <- replicate(iter_num,0)
evlmaps <- replicate(iter_num,0)


##-------------------------------------
## Step 2.1: Initialization classifier
##-------------------------------------
for(i in 1:View_num){
  cvfit<-cv.glmnet(x = data.train[[i]],
                   y = y_train,
                   alpha = 1,
                   family = "binomial",
                   type.measure = "class")
  valpredmatrix[[1]][,i] <- predict(cvfit, data.train[[i]], type = "response", s = "lambda.min")
}

##-----------------------
## Step 2.2: Optimization
##-----------------------
for (iter in 1:iter_num){

  if(length(unlist(selectedidx)) == (length(y_train)*View_num)){break}
  cat("Starting the ",iter,"-th iteration.\n", sep = "")
  
  for(j in 1:View_num){
    # update v_view
    dev_decval = valpredmatrix[[iter]]
    v_iter = mvselfpace.rank(dev_decval = dev_decval, 
                             dev_labels = y_train, 
                             v_iter = v_iter, View_id = j, 
                             lambda = lambda, gamma = gamma, 
                             View_num = View_num,num_add = num_add)
    
    for(i in 1:View_num){
      selectedidx[[i]] = which(v_iter[,i]==1)
    }
    
    # update w_view Logistic with Lasso or Elasitc net
    train.idx <- selectedidx[[j]]
    
    cvfit<-cv.glmnet(x = data.matrix(data.train[[j]][train.idx,]),
                     y = data.matrix(y_train[train.idx]),
                     alpha = 1,
                     family = "binomial",
                     type.measure = "class") 
    
    valprediction <- predict(cvfit, data.train[[j]], type = "response", s = "lambda.min")
    tsprediction <- predict(cvfit, data.test[[j]],type = "response",s = "lambda.min")
    coefprediction <- as.vector(coef(cvfit,s = "lambda.min"))
    numbernonzerocoef <- length(which(coefprediction!=0))
    
    valpredmatrix[[iter]][,j] <- as.numeric(valprediction)
    evlpredmatrix[[iter]][,j] <- tsprediction
    coefmatrix[[iter]][[j]] <- coefprediction
    nonzerocoefmatrix[[iter]][,j] <- numbernonzerocoef
  }
  
  #evaluate the training and test error
  val_loss <- sum((valpredmatrix[[iter]] - val_labels)^2)
  evl_loss <- sum((evlpredmatrix[[iter]] - evl_labels)^2)
  valmaps[iter] <- val_loss
  evlmaps[iter] <- evl_loss
  
  # update lambda and valpredmatrix for next iteriation
  lambda <- uplambda + lambda
  num_add <- num_add + num_up
  
  valpredmatrix[[iter+1]] <- valpredmatrix[[iter]]

}


##----------------------------------------------------
# Step 3: Find the run with the best valudation map
##----------------------------------------------------
## best results ##
best.iter <- which(valmaps == min(valmaps[1:length(which(valmaps!=0))]))
best_valperf <- valpredmatrix[[best.iter]]
best_evlperf <- evlpredmatrix[[best.iter]]
best_coef <- coefmatrix[[best.iter]]
best_numcoef <- nonzerocoefmatrix[[best.iter]]

MSPL.result <- list("best.valperf" = best_valperf, 
                     "best.evlperf" = best_evlperf,
                     "best.coef" = best_coef, 
                     "best.numcoef" = best_numcoef)

##-------------------------------------------
## Step 4: Evaluate the best performance
##-------------------------------------------
pref.train <- evaluate.performance(pred_label = MSPL.result$best.valperf, true_label = y_train, View_num)
pref.test <- evaluate.performance(pred_label = MSPL.result$best.evlperf, true_label = y_test, View_num)
perf.beta <- evaluate.beta.performance(pred_beta = MSPL.result$best.coef, true_beta = beta.matrix)

MSPL.performance <- list("pref.train" = pref.train, "perf.test" = pref.test, "perf.beta" = perf.beta)


