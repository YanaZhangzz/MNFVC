setwd("C:/Users/yn/Desktop/MNFVC/Sigma(JOE)-final code")
library(Matrix)
library(poweRlaw)
library(MultiRNG)
library(MASS)
library(mvtnorm)
library(stats4)
library(splines)
library(KernSmooth)


source("LSE_estimation.R")
source("LSE_grad_hessian.R")
source("LSE_inference.R")
source("pred_r2_functions.R")
source("MNFVC_functions.R")
source("NFVC_functions.R")
source("FVCSIM_functions.R")

N = 150
Num.Observed = 75

min.points = Num.Observed - 10
max.points = Num.Observed
num.points = 20

Nrep = 100


pre_accu = rep(list(matrix(0,1,1)),length(Nrep))
pre_accu_NFVC = rep(list(matrix(0,1,1)),length(Nrep))
pre_accu_FVCSIM = rep(list(matrix(0,1,1)),length(Nrep))

g0=matrix(0,nrow=N,ncol=4)
beta10=rep(1/sqrt(3),3)
beta20=beta10
beta0 = cbind(beta10,beta20)

D0 = Matrix(c(0.25,-0.25,0.25,0.25), 2, 2)
B0 = Matrix(0.5, 2, 2)
Sige1 = Diagonal(n = 2, x = c(0.4,0.6))
Sige1[1,2] = 0.1
Sige1[2,1] = 0.1


for (i in 1:Nrep) {
  
  repeat{
    tryCatch({
      cat(i,"\r")
      
      dat = get.func.network(n = N, min.points = Num.Observed - 10, max.points = Num.Observed, net_type = "2")
      
      pred_accu = pred(dat)
      esti_our = pred_accu$model
      pre_accu[[i]] = c(pred_accu$R2)

      pred_accu_NFVC = pred_NFVC(dat)
      pre_accu_NFVC[[i]] = c(pred_accu_NFVC$R2)
          
      # FVCSIM 相关
      z = rmvnorm(N, rep(0, 3), diag(3))
      train_z = z[(1:(N * 4 / 5)), ]
      test_z = z[-(1:(N * 4 / 5)), ]
      pred_accu_FVCSIM = pred_Zhu(dat)
      pre_accu_FVCSIM[[i]] = c(pred_accu_FVCSIM$R2)
       
      break
    }, error=function(e){
      cat("Error occurred, restarting iteration ", i, "\n", sep = "")
    })
    
  }
}


all_result <- list(pre_accu=pre_accu,
                   pre_accu_NFVC=pre_accu_NFVC,
                   pre_accu_FVCSIM=pre_accu_FVCSIM)
filename = paste("case-b-r2-02-", N, "-", Num.Observed, ".rda", sep = "")
save(all_result, file = filename)



pre_accu=all_result$pre_accu
pre_accu_NFVC=all_result$pre_accu_NFVC
pre_accu_FVCSIM=all_result$pre_accu_FVCSIM


pre = apply(do.call(rbind,pre_accu),2,median)
pre_NFVC =  apply(do.call(rbind,pre_accu_NFVC),2,median)
pre_FVCSIM = apply(do.call(rbind,pre_accu_FVCSIM),2,median)
round(pre,4)
round(pre_NFVC,4)
round(pre_FVCSIM,4)



