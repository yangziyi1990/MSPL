##----------------------------------
##  Code Version 1.0
##  This real data code for multimodal self-paced learning for multiomics data analysis.
##  MSPL
##  Created by Zi-Yi Yang 
##  Modified by Zi-Yi Yang on July 23, 2019
##  Concact: yangziyi091100@163.com
##----------------------------------


## Load library
library(Matrix)
library(tseries)
library(glmnet)
library(caret)
library(ncvreg)
library(pROC)
library(ROCR)
library(ggplot2)

## Set path
setwd("D:/Ziyi/School/PMO/12.Multi_omics/7.MSPL-code/2-Benchmark")
source('MVSPL_Function.R')

## Load data
load("D:/Ziyi/School/PMO/12.Multi_omics/7.MSPL-code/2-Benchmark/1-data/SNFdatasets.RDATA")

#----------------
## Colon dataset
#----------------
# Preparing data
data <- snf_data$colon
label <- snf_group$colon

# random select samples and Setting training and testing
randidx = sample(c(1:length(label)),size=length(label))
splits = vector(mode = "numeric", length = length(label))
splits_trainidx = randidx[1:(0.7*length(randidx))]
splits_testidx = randidx[(0.7*length(randidx)+1):length(randidx)]

splits[splits_trainidx] = 0
splits[splits_testidx] = 1
splits = as.matrix(splits)

trainidx = which(splits[,1]==0)
testidx = which(splits[,1]==1)

train_mrna<-data$mrna[trainidx,]
train_mirna<-data$mirna[trainidx,]
train_cpg<-data$cpg[trainidx,]

test_mrna<-data$mrna[testidx,]
test_mirna<-data$mirna[testidx,]
test_cpg<-data$cpg[testidx,]

label[which(label=="high")] <- 1
label[which(label=="low")] <- 0

train_label <- label[trainidx]
test_label <- label[testidx]

## Generate Colon data
train_data = list(mrna=train_mrna, mirna=train_mirna, cpg=train_cpg)
train_lable = data.matrix(as.numeric(train_label))
test_data = list(mrna=test_mrna, mirna=test_mirna, cpg=test_cpg)
test_lable = data.matrix(as.numeric(test_label))

##-----------------
## Step 1 : Initialization parameters
##-----------------
View_num = 3
iter_num = 30
gamma = 0.1
lambda = c(0.06, 0.06, 0.06)
lambdaup = 0.04
num_add = 4
num_up = 2

valpredmatrix = list()
evlpredmatrix = list()
coefmatrix = list()
nonzerocoefmatrix = list()
coef_idx = list()
coef_value = list()
marker_record = list()
selectedidx = list()

loss = matrix(0, nrow = length(trainidx), ncol = View_num)
v_iter = matrix(0, nrow = length(trainidx), ncol = View_num)

for(iter in 1:iter_num) {
  valpredmatrix[[iter]] = matrix(0, nrow = length(trainidx), ncol = View_num)
  evlpredmatrix[[iter]] = matrix(0, nrow = length(testidx), ncol = View_num)
  coefmatrix[[iter]] =  list()
  nonzerocoefmatrix[[iter]] = matrix(0, nrow = 1, ncol = View_num)
}

val_labels <- matrix(rep(train_lable,each = View_num),ncol = View_num, byrow = T)
evl_labels <- matrix(rep(test_lable,each = View_num),ncol = View_num, byrow = T)
valmaps <- replicate(iter_num,0)
evlmaps <- replicate(iter_num,0)


##------------------------------------
## Step 2.1: Initialization classifier
##------------------------------------
for(i in 1:View_num){
  cvfit<-cv.glmnet(x = train_data[[i]],
                   y = train_lable,
                   alpha = 1,
                   family = "binomial",
                   type.measure = "class") 
  valpredmatrix[[1]][,i] <- predict(cvfit,train_data[[i]],type="response",s="lambda.min")
}

##-----------------------
## Step 2.2: Optimization
##-----------------------
for (iter in 1:iter_num){
  
  if(length(unlist(selectedidx)) == (length(train_lable)*View_num)){break}
  cat("Starting the ",iter,"-th iteration.\n", sep = "")
  
  for(j in 1:View_num){
    # update v_view
    dev_decval = valpredmatrix[[iter]]
    v_iter  = mvselfpace.rank(dev_decval = dev_decval, 
                              dev_labels = train_lable, 
                              v_iter = v_iter, View_id = j, 
                              lambda = lambda, gamma = gamma, 
                              View_num = View_num,num_add = num_add)
    
    
    for(i in 1:View_num){
      selectedidx[[i]] = which(v_iter[,i]==1)
    }
    
    # update w_view Logistic with Lasso or Elasitc net
    train.idx = selectedidx[[j]]
    
    cvfit<-cv.glmnet(x = data.matrix(train_data[[j]][train.idx,]),
                     y = data.matrix(train_lable[train.idx]),
                     alpha = 1,
                     family = "binomial",
                     type.measure = "class") # 1-lasso,0-ridge
    
    valprediction <- predict(cvfit, train_data[[j]], type = "response", s = "lambda.min")
    tsprediction <- predict(cvfit, test_data[[j]], type = "response", s = "lambda.min")
    coefprediction <- as.vector(coef(cvfit,s="lambda.min")[-1])
    numbernonzerocoef <- length(which(coefprediction!=0))
    
    valpredmatrix[[iter]][,j] = as.numeric(valprediction)
    evlpredmatrix[[iter]][,j] = tsprediction
    coefmatrix[[iter]][[j]] = coefprediction
    nonzerocoefmatrix[[iter]][,j] = numbernonzerocoef
  }
  
  #evaluate the training and test error
  val_loss <- sum((valpredmatrix[[iter]] - val_labels)^2)
  evl_loss <- sum((evlpredmatrix[[iter]] - evl_labels)^2)
  valmaps[iter] <- val_loss
  evlmaps[iter] <- evl_loss
  
  # update lambda and valpredmatrix for next iteriation
  lambda= lambda + lambdaup
  num_add = num_add + num_up
  
  valpredmatrix[[iter+1]]=valpredmatrix[[iter]]
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

MVSPL.result <- list("best.valperf" = best_valperf, 
                     "best.evlperf" = best_evlperf,
                     "best.coef" = best_coef, 
                     "best.numcoef" = best_numcoef)

##-------------------------------------------
## Step 4: Evaluate the best performance
##-------------------------------------------
pref.train <- evaluate.performance(pred_label = best_valperf, true_label = train_lable, View_num)
pref.test <- evaluate.performance(pred_label = best_evlperf, true_label = test_lable, View_num)

MVSPL.performance <- list("pref.train" = pref.train, "perf.test" = pref.test)
  
  


