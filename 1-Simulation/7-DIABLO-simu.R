setwd("D:/Ziyi/School/PMO/12.Multi_omics/7.MSPL-code/1-Simulation")

##----------------------------------
##  Code Version 1.0
##  This simulation code for DIABLO.
##  DIABLO
##  Created by Zi-Yi Yang
##  Modified by Zi-Yi Yang on July 23, 2019
##  Concact: yangziyi091100@163.com
##----------------------------------

##---------------------
## load libraries
##---------------------
library(mixOmics)
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


##--------------------------
## 2.1 Parameter choice
##--------------------------
data.train <- list(mrna = data.train[[1]],mirna = data.train[[2]],cpg = data.train[[3]])
data.test <- list(mrna = data.test[[1]],mirna = data.test[[2]],cpg = data.test[[3]])

design = matrix(0.1, ncol = length(data.train), nrow = length(data.train), 
                dimnames = list(names(data.train), names(data.train)))
diag(design) = 0


# 2.1.1 Tuning the number of components
sgccda.res = block.splsda(X = data.train, Y = y_train, ncomp = 5, 
                          design = design)
perf.diablo = perf(sgccda.res, validation = 'Mfold', folds = 10, nrepeat = 10)
ncomp = perf.diablo$choice.ncomp$WeightedVote["Overall.BER", "centroids.dist"]

# 2.1.2 Tuning keepX
test.keepX = list (mRNA = c(5:9, seq(10, 18, 2), seq(20,30,5)),
                   miRNA = c(5:9, seq(10, 18, 2), seq(20,30,5)),
                   proteomics = c(5:9, seq(10, 18, 2), seq(20,30,5)))
tune.TCGA = tune.block.splsda(X = data.train, Y = y_train, ncomp = ncomp, 
                              test.keepX = test.keepX, design = design,
                              validation = 'Mfold', folds = 10, nrepeat = 1,
                              cpus = 2, dist = "centroids.dist")

list.keepX = tune.TCGA$choice.keepX

##------------------------------
# 3. Final model
##------------------------------
res <- block.splsda(X = data.train, Y = y_train, ncomp = ncomp,
                    keepX = list.keepX, design = design)

##-------------------------------
# 4. Prediction
##-------------------------------

predict_train = predict(res, newdata = data.train)
confusion.mat.train = get.confusion_matrix(truth = y_train,
                                           predicted = predict_train$WeightedVote$centroids.dist[,1])

predict_test = predict(res, newdata = data.test)
confusion.mat.test = get.confusion_matrix(truth = y_test, 
                                          predicted = predict_test$WeightedVote$centroids.dist[,1])

list_selectVar_mrna = selectVar(res, comp = 1, block = "mrna")
list_selectVar_mirna <- selectVar(res, comp = 1, block = "mirna")
list_selectVar_cpg <- selectVar(res, comp = 1, block = "cpg")

Var.mrna <- gsub("X","",list_selectVar_mrna$mrna$name)
Var.mirna <- gsub("X","",list_selectVar_mirna$mirna$name)
Var.cpg <- gsub("X","",list_selectVar_cpg$cpg$name)

gsub("X","",Var.mrna)
gsub("X","",Var.mirna)
gsub("X","",Var.cpg)

Coef.idx <- list("coef.mrna" = Var.mrna,
                 "coef.mirna" = Var.mirna,
                 "coef.cpg" = Var.cpg)

##-------------------------------
# 4. Performance
##-------------------------------

perf.train <- evaluate.DIABLO.performance(confusion.mat = confusion.mat.train, true_lable = y_train, 
                                          predict_label = predict_train$WeightedVote$centroids.dist[,1])

perf.test <- evaluate.DIABLO.performance(confusion.mat = confusion.mat.test, true_lable = y_test, 
                                         predict_label = predict_test$WeightedVote$centroids.dist[,1])

perf.beta <- evaluate.beta.DIABLO.performance(coef = Coef.idx, beta = beta.matrix)


DIABLO.performance <- list("perf.train" = perf.train, "perf.test" = perf.test, "perf.beta" = perf.beta)
