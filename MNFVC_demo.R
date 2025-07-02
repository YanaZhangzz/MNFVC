setwd("/MNFVC_code")
library(Matrix)
library(poweRlaw)
library(MultiRNG)
library(MASS)
library(mvtnorm)
library(stats4)
library(splines)
library(KernSmooth)
set.seed(1234)
source("MNFVC_functions.R")
source("LSE_estimation.R")
source("LSE_grad_hessian.R")
source("LSE_inference.R")
source("NFVC_functions.R")
source("FVCSIM_functions.R")

N = 150
Num.Observed = 75

min.points = Num.Observed - 10
max.points = Num.Observed
num.points = 20

Nrep = 100

dbeta = rep(list(matrix(0, nrow = num.points, ncol = 8)), length(Nrep))
CST11 = rep(list(matrix(0, num.points, num.points), length(Nrep)))
CST22 = rep(list(matrix(0, num.points, num.points), length(Nrep)))
NFVC_esti_coef1 = rep(list(matrix(0,nrow=max.points,ncol=3)),length(Nrep))
NFVC_esti_coef2 = rep(list(matrix(0,nrow=max.points,ncol=3)),length(Nrep))
FVCSIM_b = rep(list(matrix(0,nrow=max.points,ncol=4)),length(Nrep))

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
      
      ans = est_d_beta(dat$Y,dat$v1,dat$Tpoints,dat$W,D0,B0,Sige1,num.points)
      dbeta[[i]] = ans$est_coef
      
      CST0 = est_covariance(ans$est_coef, ans$timeline, Y1=dat$Y1, Y2=dat$Y2, X=dat$v1, Tpoints=dat$Tpoints, W=dat$W)
      CST11[[i]] = CST0$cov11
      CST22[[i]] = CST0$cov22
      
      #NFVC
      
      NFVC_coef1 = est_entire_NFVC(dat$Y1N, dat$v1, dat$Tpoints, dat$W, num.points)
      NFVC_coef2 = est_entire_NFVC(dat$Y2N, dat$v1, dat$Tpoints, dat$W, num.points)
      NFVC_esti_coef1[[i]] = NFVC_coef1$est_coef
      NFVC_esti_coef2[[i]] = NFVC_coef2$est_coef
      
      # FVCSIM
      z = rmvnorm(N, rep(0, 3), diag(3))
      FVCSIM_ans = iteration(3, 1e-5, dat$Y.all, dat$v1, z, g0)
      FVCSIM_b[[i]] = FVCSIM_ans$alpha
      
      break
    }, error=function(e){
      cat("Error occurred, restarting iteration ", i, "\n", sep = "")
    })
    
  }
}

all_result = list(dbeta=dbeta,NFVC_esti_coef1=NFVC_esti_coef1,NFVC_esti_coef2=NFVC_esti_coef2,FVCSIM_b=FVCSIM_b,CST11= CST11, CST22= CST22)
filename = paste("case-b-02-", N, "-", Num.Observed, ".rda", sep = "")
save(all_result, file = filename)



dbeta=all_result$dbeta
NFVC_esti_coef1=all_result$NFVC_esti_coef1
NFVC_esti_coef2=all_result$NFVC_esti_coef2
NFVC_esti_coef=rep(list(matrix(0, nrow = 20, ncol = 6)), length(Nrep))
for(i in 1:Nrep){
  NFVC_esti_coef[[i]] = cbind(NFVC_esti_coef1[[i]],NFVC_esti_coef2[[i]])
}

FVCSIM_b=all_result$FVCSIM_b



t=seq(0,1,length=num.points)
# true_coef = cbind(0.2*sin(2*t*pi)+0.3, 0.5 * sin(pi*(t - 0.5)),  0.5 - 1/(1 + exp(-5*(t - 0.5))), 0.2*cos(2*t*pi)+0.3,
#                   (t)^2, 5*(t-0.5)^2, (1-t)^2, (t)^0.5)
# true_coef = cbind(0.2*sin(2*t*pi)+0.3, 2*(t^2-t)+0.5,  2*(t-t^2), 0.2*cos(2*t*pi)+0.3,
#                   (t)^2, 5*(t-0.5)^2, (1-t)^2, (t)^0.5)
# true_coef = cbind(0.2*sin(2*t*pi)+0.3, 0*t,  0*t, 0.2*cos(2*t*pi)+0.3,
#                   (t)^2, 5*(t-0.5)^2, (1-t)^2, (t)^0.5)
true_coef = cbind(0.2*sin(2*t*pi)+0.3, -0.3*sin(2*pi*t)-0.2, 0.3*cos(2*pi*t)+0.2, 0.2*cos(2*t*pi)+0.3,
                  (t)^2, 5*(t-0.5)^2, (1-t)^2, (t)^0.5)
NFVC_true_coef=cbind(0.2*sin(2*t*pi)+0.3,t^2, 5*(t-0.5)^2,
                     0.2*cos(2*t*pi)+0.3,(1-t)^2, t^0.5)
tf = seq(0,1,length=max.points)
FVCSIM_true_b = cbind(tf^2,5*(tf-0.5)^2,(1-tf)^2,tf^0.5)


esti_coef = apply(simplify2array(dbeta), 1:2, median)
NFVC_escoef1 = apply(simplify2array(NFVC_esti_coef1), 1:2, median)
NFVC_escoef2 = apply(simplify2array(NFVC_esti_coef2), 1:2, median)
NFVC_escoef = cbind(NFVC_escoef1,NFVC_escoef2)
FVCSIM_esb = apply(simplify2array(FVCSIM_b), 1:2, median)

###############IMSE-SUP################

IMSE_coef =colMeans((true_coef-esti_coef)^2)
print(c('d11', 'd12','d21','d22', 'b11','b21','b12', 'b22'))
r1 = IMSE_coef
round(r1, 4)
SUP_coef = apply(((true_coef-esti_coef)^2),2,max)
round(SUP_coef,4)


###NFVC
NFVC_IMSE_coef=colMeans((NFVC_true_coef-NFVC_escoef)^2)
print(c('d11','b11','b21','d22','b12','b22'))
r2=NFVC_IMSE_coef
round(r2, 4)
NFVC_SUP_coef = apply(((NFVC_true_coef-NFVC_escoef)^2),2,max)
round(NFVC_SUP_coef,4)


####FVCSIM
FVCSIM_b_IMSE=colMeans((FVCSIM_true_b-FVCSIM_esb)^2)
r3 = FVCSIM_b_IMSE
round(r3,4)
FVCSIM_b_SUP = apply(((FVCSIM_true_b-FVCSIM_esb)^2),2,max)
round(FVCSIM_b_SUP,4)



