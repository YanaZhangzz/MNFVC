
setwd("C:/Users/yn/Desktop/MNFVC/Sigma(JOE)-final code/MNFVC_code")
library(Matrix)
library(poweRlaw)
library(MultiRNG)
library(MASS)
library(mvtnorm)
# library(VGSM)
library(stats4)
library(splines)
library(KernSmooth)
# set.seed(124)

source("LSE_estimation-0.5.R")
source("LSE_grad_hessian.R")
source("LSE_inference.R")
source("pred_sum.R")
source("MNFVC_functions.R")
source("NFVC_functions.R")
source("FVCSIM_functions.R")


N = 150
Num.Observed = 20

min.points = Num.Observed - 10
max.points = Num.Observed
num.points = 20

Nrep = 100


pre_accu = rep(list(matrix(0,1,1)),length(Nrep))
inside_pre_accu = rep(list(matrix(0,1,1)),length(Nrep))
pre_accu_NFVC = rep(list(matrix(0,1,1)),length(Nrep))
pre_accu_FVCSIM = rep(list(matrix(0,1,1)),length(Nrep))
inside_pre_accu_NFVC = rep(list(matrix(0,1,1)),length(Nrep))
inside_pre_accu_FVCSIM = rep(list(matrix(0,1,1)),length(Nrep))

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
      inside_pred_accu = inside_pred(esti_our)
      pre_accu[[i]] = c(pred_accu$IMSE)
      inside_pre_accu[[i]] = c(inside_pred_accu$IMSE)
 
      pred_accu_NFVC = pred_NFVC(dat)
      inside_pred_accu_NFVC = inside_pred_NFVC(pred_accu_NFVC$NFVC_ans1, pred_accu_NFVC$NFVC_ans2)
      pre_accu_NFVC[[i]] = c(pred_accu_NFVC$IMSE)
      inside_pre_accu_NFVC[[i]] = c(inside_pred_accu_NFVC$IMSE)
      
   
      z = rmvnorm(N, rep(0, 3), diag(3))
      train_z = z[(1:(N * 4 / 5)), ]
      test_z = z[-(1:(N * 4 / 5)), ]
      pred_accu_FVCSIM = pred_Zhu(dat)
      inside_pred_accu_FVCSIM = inside_pred_Zhu(pred_accu_FVCSIM$esti)
      pre_accu_FVCSIM[[i]] = c(pred_accu_FVCSIM$IMSE)
      inside_pre_accu_FVCSIM[[i]] = c(inside_pred_accu_FVCSIM$IMSE)
      
      break
    }, error=function(e){
      cat("Error occurred, restarting iteration ", i, "\n", sep = "")
    })
    
  }
}

all_result <- list(pre_accu=pre_accu,inside_pre_accu=inside_pre_accu,
                   pre_accu_NFVC=pre_accu_NFVC,inside_pre_accu_NFVC=inside_pre_accu_NFVC,
                   pre_accu_FVCSIM=pre_accu_FVCSIM,inside_pre_accu_FVCSIM=inside_pre_accu_FVCSIM)
filename = paste("case-b-pred-02-", N, "-", Num.Observed, ".rda", sep = "")
save(all_result, file = filename)



####result analysis#####

pre_accu=all_result$pre_accu
pre_accu_NFVC=all_result$pre_accu_NFVC
pre_accu_FVCSIM=all_result$pre_accu_FVCSIM
inside_pre_accu=all_result$inside_pre_accu
inside_pre_accu_NFVC=all_result$inside_pre_accu_NFVC
inside_pre_accu_FVCSIM=all_result$inside_pre_accu_FVCSIM
#

######prediction accuracy####
pre = apply(do.call(rbind,pre_accu),2,median)
pre_NFVC =  apply(do.call(rbind,pre_accu_NFVC),2,median)
pre_FVCSIM = apply(do.call(rbind,pre_accu_FVCSIM),2,median)
round(pre,4)
round(pre_NFVC,4)
round(pre_FVCSIM,4)

inpre = apply(do.call(rbind,inside_pre_accu),2,median)
inpre_NFVC = apply(do.call(rbind,inside_pre_accu_NFVC),2,median)
inpre_FVCSIM = apply(do.call(rbind,inside_pre_accu_FVCSIM),2,median)
round(inpre,4)
round(inpre_NFVC,4)
round(inpre_FVCSIM,4)

