
getCovCST = function(s)
{
  phi11.seq = matrix(sin(2*s*pi) * sqrt(2), ncol = 1)
  phi12.seq = matrix(cos(2*s*pi) * sqrt(2), ncol = 1)
  phi1.seq = rbind(phi11.seq,phi12.seq)
  phi21.seq = matrix(cos(2*s*pi) * sqrt(2), ncol = 1)
  phi22.seq = matrix(sin(2*s*pi) * sqrt(2), ncol = 1)
  phi2.seq = rbind(phi21.seq,phi22.seq)
  
  TCST11 = 1.2*phi11.seq%*%t(phi11.seq)+0.6*phi12.seq%*%t(phi12.seq)
  TCST22 = phi21.seq%*%t(phi21.seq)+0.5*phi22.seq%*%t(phi22.seq)
  
  return(list(tcov11 = TCST11,tcov22 = TCST22))
}

gamma = function(Y,X,Tpoints,W,esti_coef,t_est){
  
  N = nrow(Y)/2
  num.points = nrow(esti_coef)
  Y1s=Y[1:N,]
  Y2s=Y[-(1:N),]
  
  W.hat = kronecker(diag(2),(W))
  X_bolang=kronecker(diag(2),X)
  
  gamma = array(1:N*2*num.points,dim=c(N,2,num.points))
  ga1 = matrix(0,nrow=N,ncol=num.points)
  ga2 = matrix(0,nrow=N,ncol=num.points)
  for (m in 1:num.points){
    T_diff = Tpoints - t_est[m]
    h = rule.thumb(T_diff[1,])/3
    Ker_mat = gaussian.kernel(T_diff, h)
    
    T_diff.hat = rbind(T_diff,T_diff)
    Ker_mat.hat = gaussian.kernel(T_diff.hat, h)
    
    ns = rowSums(!is.na(Tpoints))
    ns.hat = rowSums(!is.na(rbind(Tpoints,Tpoints)))
    fs = rowSums(Ker_mat, na.rm = T)/ns
    fs.hat = rowSums(Ker_mat.hat, na.rm = T)/ns.hat
    
    Y_Ker_mat_1 = Y1s*Ker_mat/ns/fs
    Y2_Ker_mat_1 = Y1s^2*Ker_mat/ns/fs
    Y_Ker_vec_1 = rowSums(Y_Ker_mat_1, na.rm = T) 
    Y2_Ker_vec_1 = rowSums(Y2_Ker_mat_1, na.rm = T) 
    
    Y_Ker_mat_2 = Y2s*Ker_mat/ns/fs
    Y2_Ker_mat_2 = Y2s^2*Ker_mat/ns/fs
    Y_Ker_vec_2 = rowSums(Y_Ker_mat_2, na.rm = T) 
    Y2_Ker_vec_2 = rowSums(Y2_Ker_mat_2, na.rm = T) 
    
    Y_Ker_mat.hat = rbind(Y_Ker_mat_1,Y_Ker_mat_2)
    Y2_Ker_mat.hat = rbind(Y2_Ker_mat_1,Y2_Ker_mat_2)
    #Y_Ker_mat = Y_Ker_mat.hat/ns/fs
    #Y2_Ker_mat = Y2_Ker_mat.hat/ns/fs
    Y_Ker_vec = rowSums(Y_Ker_mat.hat, na.rm = T)
    Y2_Ker_vec = rowSums(Y2_Ker_mat.hat, na.rm = T)
    
    zeta = apply(Y*Ker_mat.hat, 1, function(s){return(mean(s, na.rm = T))})#????Y*Ker_mat????ÿ???ڵ?(?????ĵ?һά??)????ֵ,2N*1
    zetat = t((zeta/fs.hat)%*%t(matrix(1, nrow = 2*N, ncol = 1)))#2N*2N
    
    wzf1 = matrix(0,nrow=N,ncol=N)
    for(i in 1:N){
      wzf1[,i] = W[,i]*(zeta/fs.hat)[1:N][i]
    }
    wzfs1 = colSums(wzf1)
    
    wzf2 = matrix(0,nrow=N,ncol=N)
    for(i in 1:N){
      wzf2[,i] = W[,i]*(zeta/fs.hat)[-(1:N)][i]
    }
    wzfs2 = colSums(wzf2)
    
    ##j=1
    ga1[,m] = (zeta/fs.hat)[1:N] - esti_coef[,1][m]*wzfs1-esti_coef[,2][m]*wzfs2-esti_coef[,5][m]*X[,1]-esti_coef[,6][m]*X[,2]
    
    ga2[,m] = (zeta/fs.hat)[-(1:N)] - esti_coef[,3][m]*wzfs1-esti_coef[,4][m]*wzfs2-esti_coef[,7][m]*X[,1]-esti_coef[,8][m]*X[,2]
    
  }
  gamma[,1,]=ga1
  gamma[,2,]=ga2
  return(gamma)
}

lk = function(K,Y,X,Tpoints,W,esti_coef,t,ecov11,ecov22){
  
  e1evalue = eigen(ecov11)$value
  e1efun = eigen(ecov11)$vectors
  e2evalue = eigen(ecov22)$value
  e2efun  = eigen(ecov22)$vectors
  
  
  N = nrow(Y)/2
  num.points=nrow(esti_coef)
  
  score = array(1:N*2*max.points,dim=c(N,2,num.points))
  sci1 = matrix(0,nrow=N,ncol=num.points)
  sci2 = matrix(0,nrow=N,ncol=num.points)
  
  gammap = gamma(Y,X,Tpoints,W,esti_coef,t)
  for(i in 1:N){
    for(k in 1:K){
      sci1[i,] = e1evalue[k]*t(e1efun[,k])%*%solve(ecov11)%*%gammap[i,1,]
      sci2[i,] = e2evalue[k]*t(e2efun[,k])%*%solve(ecov22)%*%gammap[i,2,]
    }
    score[,1,] = sci1
    score[,2,] = sci2
  }
  
  
  losum = 0
  for(m in 1:num.points){
    d11m = esti_coef[,1][m]
    d12m = esti_coef[,2][m]
    d21m = esti_coef[,3][m]
    d22m = esti_coef[,4][m]
    dm = matrix(c(d11m,d12m,d21m,d22m),nrow=2,byrow=TRUE)
    S = diag(N*2)-kronecker(t(dm),W)
    
    Y1s=Y[1:N,]
    Y2s=Y[-(1:N),]
    
    T_diff = Tpoints - t[m]
    h = rule.thumb(T_diff[1,])/3
    Ker_mat = gaussian.kernel(T_diff, h)
    
    T_diff.hat = rbind(T_diff,T_diff)
    Ker_mat.hat = gaussian.kernel(T_diff.hat, h)
    
    ns = rowSums(!is.na(Tpoints))
    ns.hat = rowSums(!is.na(rbind(Tpoints,Tpoints)))
    fs = rowSums(Ker_mat, na.rm = T)/ns
    fs.hat = rowSums(Ker_mat.hat, na.rm = T)/ns.hat
    
    Y_Ker_mat_1 = Y1s*Ker_mat/ns/fs
    Y2_Ker_mat_1 = Y1s^2*Ker_mat/ns/fs
    Y_Ker_vec_1 = rowSums(Y_Ker_mat_1, na.rm = T)
    Y2_Ker_vec_1 = rowSums(Y2_Ker_mat_1, na.rm = T)
    
    Y_Ker_mat_2 = Y2s*Ker_mat/ns/fs
    Y2_Ker_mat_2 = Y2s^2*Ker_mat/ns/fs
    Y_Ker_vec_2 = rowSums(Y_Ker_mat_2, na.rm = T)
    Y2_Ker_vec_2 = rowSums(Y2_Ker_mat_2, na.rm = T)
    
    #Y_Ker_mat = (Y*Ker_mat.hat)/ns.hat/fs.hat
    #Y2_Ker_mat = (Y^2*Ker_mat.hat)/ns.hat/fs.hat
    #Y_Ker_vec = rowSums(Y_Ker_mat, na.rm = T)
    #Y2_Ker_vec = rowSums(Y2_Ker_mat, na.rm = T)
    
    Y_Ker_mat.hat = rbind(Y_Ker_mat_1,Y_Ker_mat_2)
    Y2_Ker_mat.hat = rbind(Y2_Ker_mat_1,Y2_Ker_mat_2)
    #Y_Ker_mat = Y_Ker_mat.hat/ns/fs
    #Y2_Ker_mat = Y2_Ker_mat.hat/ns/fs
    Y_Ker_vec = rowSums(Y_Ker_mat.hat, na.rm = T)
    Y2_Ker_vec = rowSums(Y2_Ker_mat.hat, na.rm = T)
    
    zeta = apply(Y*Ker_mat.hat, 1, function(s){return(mean(s, na.rm = T))})#2N*1
    zetah = zeta/fs.hat
    
    Szeta = S%*%zetah
    
    Xhat = kronecker(diag(2),X)
    b11m = esti_coef[,5][m]
    b21m = esti_coef[,6][m]
    b12m = esti_coef[,7][m]
    b22m = esti_coef[,8][m]
    bm = c(b11m,b21m,b12m,b22m)
    
    XB = Xhat%*%bm
    
    phisum = matrix(0,nrow=2*N,ncol=1)
    for (k in 1:K){
      phi1 = score[,1,k]*e1efun[k,m]
      phi2 = score[,2,k]*e2efun[k,m]
      phi = rbind(phi1,phi2)
      phim = matrix(phi,ncol=1)
      phisum = phisum+phim
    }
    
    
    # lo = -(colSums(log((Szeta - XB - phisum)^2)))/N
    lo = -N*log(t(Szeta - XB - phisum)%*%(Szeta - XB - phisum))
  }
  
  losum = losum+lo
  lok = losum/num.points
  return(lok)
}

pred = function(dat)
{
  Y1 = dat$Y1
  Y2 = dat$Y2
  X = dat$v1
  W = dat$W
  Tpoints = dat$Tpoints
  
  t_grid <- seq(0, 1, length.out = N)
  train_size = floor(4*N/5)
  train_grid <- seq(1,train_size) 
  test_grid <- seq(train_size+1,N)  
  
  train_Y1 <- Y1[train_grid,]
  test_Y1 <- Y1[test_grid,]
  
  train_Y2 <- Y2[train_grid,]
  test_Y2 <- Y2[test_grid,]
  
  train_Y = rbind(train_Y1,train_Y2)
  test_Y = rbind(test_Y1,test_Y2)
  
  train_X = X[train_grid,]
  test_X = X[test_grid,]
  
  train_W = W[train_grid,train_grid]
  test_W = W[test_grid,test_grid]
  
  train_Tpoints <- Tpoints[train_grid,]
  test_Tpoints <- Tpoints[test_grid,]
  
  model = est_d_beta(train_Y,train_X,train_Tpoints,train_W,D0,B0,Sige1,max.points)
  ecovp = est_covariance(model$est_coef,model$timeline,test_Y1,test_Y2,test_X,test_Tpoints,test_W)
  ecov11p = ecovp$cov11
  ecov22p = ecovp$cov22
  e1evaluep = eigen(ecov11p)$values
  e1efunp = eigen(ecov11p)$vectors
  e2evaluep = eigen(ecov22p)$values
  e2efunp = eigen(ecov22p)$vectors
  
  obBIC = function(K,Y,X,Tpoints,W,esti_coef,t,ecov11,ecov22){
    return(-2*lk(K,Y,X,Tpoints,W,esti_coef,t,ecov11,ecov22)+8*log(N/5))
  }
  BIC1=c()
  for(K in 1:max.points){
    BIC1[K]=obBIC(K,test_Y,test_X,test_Tpoints,test_W,model$est_coef,model$timeline,ecov11p,ecov22p)
  }
  K = which.min(BIC1)
  
  gammap = gamma(test_Y,test_X,test_Tpoints,test_W,model$est_coef,model$timeline)
  sci1 = matrix(0,nrow=N/5,ncol=K)
  sci2 = matrix(0,nrow=N/5,ncol=K)
  for(i in 1:(N/5)){
    for(k in 1:K){
      sci1[i,k] = e1evaluep[k]*t(e1efunp[,k])%*%solve(ecov11p)%*%gammap[i,1,]
    }
    for(k in 1:K){
      sci2[i,k] = e2evaluep[k]*t(e2efunp[,k])%*%solve(ecov22p)%*%gammap[i,2,]
    }
  }
  
  phisum = matrix(0,nrow=(N*2/5),ncol=max.points)
  for (m in 1:max.points){
    for (k in 1:K){
      phi1 = sci1[,k]*e1efunp[k,m]
    }
    for (k in 1:K){
      phi2 = sci2[,k]*e2efunp[k,m]
    }
    phi = rbind(phi1,phi2)
    phim = matrix(phi,ncol=1)
    phisum[,m] = phisum[,m]+phim
  }
  
  b11 = matrix(model$est_coef[,5],nrow=1)
  b12 = matrix(model$est_coef[,7],nrow=1)
  b21 = matrix(model$est_coef[,6],nrow=1)
  b22 = matrix(model$est_coef[,8],nrow=1)
  b.mat = rbind(b11,b21,b21,b22)
 
  x1 = test_X[,1]
  x2 = test_X[,2]
  v1 = cbind(x1,x2)
  v1_hat = kronecker(diag(2),v1)
  
  d11 = matrix(model$est_coef[,1],nrow=1)
  d12 = matrix(model$est_coef[,2],nrow=1)
  d21 = matrix(model$est_coef[,3],nrow=1)
  d22 = matrix(model$est_coef[,4],nrow=1)
  d.mat = rbind(d11,d12,d21,d22)
  
  Ymat = matrix(0,nrow=2*N/5,ncol=max.points)
  for(m in 1:max.points){
    Xbeta = v1_hat%*%b.mat[,m]+phisum[,m]
    D = matrix(d.mat[,m],nrow=2,byrow=TRUE)
    DkW = kronecker(t(D),test_W)
    DkW2 = DkW%*%DkW
    DkW3 = DkW2%*%DkW
    DkW4 = DkW3%*%DkW
    DkW5 = DkW4%*%DkW
    DkW6 = DkW5%*%DkW
    
    Yvec = Xbeta+DkW%*%Xbeta+DkW2%*%Xbeta+DkW3%*%Xbeta+DkW4%*%Xbeta+DkW5%*%Xbeta+DkW6%*%Xbeta
    Ymat[,m] = Yvec
  }
  
  Yhat = Ymat

  A=test_Y-Yhat
  SSE = mean(t(A)%*%A)
  mean_testY = apply(test_Y,2,mean)
  B = test_Y-mean_testY
  SST = mean(t(B)%*%B)
  R2 = 1-SSE/SST
  
  return(list(R2=R2))
}


##########NFVC_pred#####
est_cov_NFVC <- function(est_coef, timeline,Y, X, Tpoints, W)
{
  N = dim(W)[1]
  rhott = est_coef[,1]
  betatt = est_coef[,2:dim(est_coef)[2]]
  Stt = lapply(rhott, function(rho){diag(N)-rho*W})
  
  num = length(timeline)
  
  CO = matrix(0,num,num)
  ns = rowSums(!is.na(Tpoints))
  for (s in 1:num) {
    for(t in 1:num){
      T_diff = Tpoints - timeline[s]
      h = rule.thumb(T_diff[1,])/3
      Ker_mat_s = gaussian.kernel(T_diff, h)
      Y_s = rowSums(Y*Ker_mat_s, na.rm = T)
      T_diff = Tpoints - timeline[t]
      h = rule.thumb(T_diff[1,])/3
      Ker_mat_t = gaussian.kernel(T_diff, h)
      Y_t = rowSums(Y*Ker_mat_t, na.rm = T)
      Y2_st = rowSums(Y^2*Ker_mat_s*Ker_mat_t, na.rm = T)
      Vst.n = Y_s%*%t(Y_t) - diag(Y2_st)
      K_s = rowSums(Ker_mat_s, na.rm = T)
      K_t = rowSums(Ker_mat_t, na.rm = T)
      K_st = rowSums(Ker_mat_s*Ker_mat_t, na.rm = T)
      Vst.d = K_s%*%t(K_t) - diag(K_st)
      Vst = Vst.n/Vst.d
      
      suu=0
      for(i in 1:N)
      {
        suu = suu+t(Stt[[s]][i,])%*%Vst%*%(Stt[[t]][i,])
      }
      CO[s,t] = (suu - t(X%*%betatt[s,])%*%(X%*%betatt[t,]))/N
    }
  }
  return(cov = CO)
}

gamma_NFVC = function(t_est,esti_coef,Y,X,Tpoints,W){
  
  N = nrow(Y)
  num=length(t_est)
  gamma = matrix(0,nrow=N,ncol=num)
  
  for (m in 1:num){
    T_diff = Tpoints - t_est[m]
    h = rule.thumb(T_diff[1,])/3
    Ker_mat = gaussian.kernel(T_diff, h)
    ns = rowSums(!is.na(Tpoints))
    fs = rowSums(Ker_mat, na.rm = T)/ns
    Y_Ker_mat = Y*Ker_mat/ns/fs
    Y2_Ker_mat = Y^2*Ker_mat/ns/fs
    Y_Ker_vec = rowSums(Y_Ker_mat, na.rm = T) 
    Y2_Ker_vec = rowSums(Y2_Ker_mat, na.rm = T) 
    
    zeta = apply(Y*Ker_mat, 1, function(s){return(mean(s, na.rm = T))})
    zetat = t((zeta/fs)%*%t(matrix(1, nrow = N, ncol = 1)))
    
    wzf1 = matrix(0,nrow=N,ncol=N)
    for(i in 1:N){
      wzf1[,i] = W[,i]*(zeta/fs)[i]
    }
    wzfs1 = colSums(wzf1)
    
    ##j=1
    gamma[,m] = (zeta/fs)[1:N] - esti_coef[,1][m]*wzfs1-esti_coef[,2][m]*X[,1]-esti_coef[,3][m]*X[,2]
    
  }
  return(gamma)
}

gammap_NFVC = function(NFVC_ans1,NFVC_ans2,t_est,Y,X,Tpoints,W){
  ###zeta
  N = nrow(Y)/2
  
  #Y---2n*maxpoints
  Y1s=Y[1:N,]
  Y2s=Y[-(1:N),]
  
  #bianxing
  W.hat = kronecker(diag(2),(W))
  X_bolang=kronecker(diag(2),X)#2N*4
  
  np = length(t_est)
  gamma = array(1:(N*2*np),dim=c(N,2,np))
  ga1 = matrix(0,nrow=N,ncol=np)
  ga2 = matrix(0,nrow=N,ncol=np)
  for (m in 1:np){
    
    T_diff = Tpoints - t_est[m]#N*numpointsά
    h = rule.thumb(T_diff[1,])/3
    Ker_mat = gaussian.kernel(T_diff, h)
    
    T_diff.hat = rbind(T_diff,T_diff)
    Ker_mat.hat = gaussian.kernel(T_diff.hat, h)
    
    ns = rowSums(!is.na(Tpoints))
    
    ns.hat = rowSums(!is.na(rbind(Tpoints,Tpoints)))
    
    fs = rowSums(Ker_mat, na.rm = T)/ns
    fs.hat = rowSums(Ker_mat.hat, na.rm = T)/ns.hat
    
    Y_Ker_mat_1 = Y1s*Ker_mat/ns/fs
    Y2_Ker_mat_1 = Y1s^2*Ker_mat/ns/fs
    Y_Ker_vec_1 = rowSums(Y_Ker_mat_1, na.rm = T)
    Y2_Ker_vec_1 = rowSums(Y2_Ker_mat_1, na.rm = T)
    
    Y_Ker_mat_2 = Y2s*Ker_mat/ns/fs
    Y2_Ker_mat_2 = Y2s^2*Ker_mat/ns/fs
    Y_Ker_vec_2 = rowSums(Y_Ker_mat_2, na.rm = T)
    Y2_Ker_vec_2 = rowSums(Y2_Ker_mat_2, na.rm = T)
    Y_Ker_mat.hat = rbind(Y_Ker_mat_1,Y_Ker_mat_2)
    Y2_Ker_mat.hat = rbind(Y2_Ker_mat_1,Y2_Ker_mat_2)
    
    Y_Ker_vec = rowSums(Y_Ker_mat.hat, na.rm = T)
    Y2_Ker_vec = rowSums(Y2_Ker_mat.hat, na.rm = T)
    
    zeta = apply(Y*Ker_mat.hat, 1, function(s){return(mean(s, na.rm = T))})
    zetat = t((zeta/fs.hat)%*%t(matrix(1, nrow = 2*N, ncol = 1)))
    
    wzf1 = matrix(0,nrow=N,ncol=N)
    for(i in 1:N){
      wzf1[,i] = W[,i]*(zeta/fs.hat)[1:N][i]
    }
    wzfs1 = colSums(wzf1)
    
    wzf2 = matrix(0,nrow=N,ncol=N)
    for(i in 1:N){
      wzf2[,i] = W[,i]*(zeta/fs.hat)[-(1:N)][i]
    }
    wzfs2 = colSums(wzf2)
    
    d11 = NFVC_ans1$est_coef[,1]
    d22 = NFVC_ans2$est_coef[,1]
    b11 = NFVC_ans1$est_coef[,2]
    b21 = NFVC_ans1$est_coef[,3]
    b12 = NFVC_ans2$est_coef[,2]
    b22 = NFVC_ans2$est_coef[,3]
    
    ##j=1
    ga1[,m] = (zeta/fs.hat)[1:N] - d11[m]*wzfs1-b11[m]*X[,1]-b21[m]*X[,2]
    
    ga2[,m] = (zeta/fs.hat)[-(1:N)]-d22[m]*wzfs2-b12[m]*X[,1]-b22[m]*X[,2]
    
  }
  gamma[,1,]=ga1
  gamma[,2,]=ga2
  return(gamma)
}


loss_func = function(K,t_est,esti_coef,Y,X,Tpoints,W,ecov11){
  
  e1evalue = eigen(ecov11)$value
  e1efun = eigen(ecov11)$vectors
  
  gamma1 = gamma_NFVC(t_est,esti_coef,Y,X,Tpoints,W)
  n = nrow(Y)
  sci1 = matrix(0,nrow=n,ncol=K)
  for(i in 1:n){
    for(k in 1:K){
      sci1[i,k] = e1evalue[k]*t(e1efun[,k])%*%solve(ecov11)%*%gamma1[i,]
    }
  }
  
  
  losum = 0
  num = nrow(esti_coef)
  for(m in 1:num){
    d11m = esti_coef[,1][m]
    S = diag(n)-kronecker(d11m,W)
    T_diff = Tpoints - t_est[m]
    h = rule.thumb(T_diff[1,])/3
    Ker_mat = gaussian.kernel(T_diff, h)
    ns = rowSums(!is.na(Tpoints))
    fs = rowSums(Ker_mat, na.rm = T)/ns
    Y_Ker_mat = Y*Ker_mat/ns/fs
    Y2_Ker_mat = Y^2*Ker_mat/ns/fs
    Y_Ker_vec = rowSums(Y_Ker_mat, na.rm = T) 
    Y2_Ker_vec = rowSums(Y2_Ker_mat, na.rm = T) 
    
    zeta = apply(Y*Ker_mat, 1, function(s){return(mean(s, na.rm = T))})
    zetah = zeta/fs
    
    Szeta = S%*%zetah
    
    b11m = esti_coef[,2][m]
    b21m = esti_coef[,3][m]
    bm = c(b11m,b21m)
    
    XB = X%*%bm
    
    phisum1 = matrix(0,nrow=n,ncol=1)
    for (k in 1:K){
      phi1 = sci1[,k]*e1efun[k,m]
      phim = matrix(phi1,ncol=1)
      phisum1 = phisum1+phim
    }
    
    lo = -N*log(t(Szeta - XB - phisum1)%*%(Szeta - XB - phisum1))
  }
  
  losum = losum+lo
  lok = losum/num.points
  return(lok)
}


pred_NFVC=function(dat){
  
  Y=dat$Y
  Y1 = dat$Y1
  Y2 = dat$Y2
  X = dat$v1
  W = dat$W
  Tpoints = dat$Tpoints
  
  t_grid <- seq(0, 1, length.out = N)
  train_size = floor(4*N/5)
  train_grid <- seq(1,train_size)  
  test_grid <- seq(train_size+1,N)   
  
  train_Y1 <- Y1[train_grid,]
  test_Y1 <- Y1[test_grid,]
  
  train_Y2 <- Y2[train_grid,]
  test_Y2 <- Y2[test_grid,]
  
  train_Y = rbind(train_Y1,train_Y2)
  test_Y = rbind(test_Y1,test_Y2)
  
  train_X = X[train_grid,]
  test_X = X[test_grid,]
  
  train_W = W[train_grid,train_grid]
  test_W = W[test_grid,test_grid]
  
  train_Tpoints <- Tpoints[train_grid,]
  test_Tpoints <- Tpoints[test_grid,]
  
  NFVC_ans1 = est_entire_NFVC(Y=train_Y1,X=train_X,Tpoints=train_Tpoints,W=train_W,max.points)
  NFVC_ans2 = est_entire_NFVC(Y=train_Y2,X=train_X,Tpoints=train_Tpoints,W=train_W,max.points)
  
  b11 = NFVC_ans1$est_coef[,2]
  b21 = NFVC_ans1$est_coef[,3]
  beta.seq1 = rbind(b11, b21)  
  rho.mat11 = NFVC_ans1$est_coef[,1]
  
  b12 = NFVC_ans2$est_coef[,2]
  b22 = NFVC_ans2$est_coef[,3]
  beta.seq2 = rbind(b12, b22)  
  rho.mat22 = NFVC_ans2$est_coef[,1]
  v1=test_X
  
  ecovp11 = est_cov_NFVC(NFVC_ans1$est_coef,NFVC_ans1$timeline,test_Y1,test_X,test_Tpoints,test_W)
  ecovp22 = est_cov_NFVC(NFVC_ans2$est_coef,NFVC_ans2$timeline,test_Y2,test_X,test_Tpoints,test_W)
  ecov11p = ecovp11
  ecov22p = ecovp22
  e1evaluep = eigen(ecov11p)$values
  e1efunp = eigen(ecov11p)$vectors
  e2evaluep = eigen(ecov22p)$values
  e2efunp = eigen(ecov22p)$vectors
  
  obBIC = function(K,t_est,esti_coef,Y,X,Tpoints,W,ecov11){
    return(-2*loss_func(K,t_est,esti_coef,Y,X,Tpoints,W,ecov11)+K*((N/5)^0.8))
  }
  BIC1=c()
  for(K in 1:num.points){
    BIC1[K]=obBIC(K,NFVC_ans1$timeline,NFVC_ans1$est_coef,test_Y1,test_X,test_Tpoints,test_W,ecov11p)
  }
  K1 = which.min(BIC1)
  
  BIC2=c()
  for(K in 1:num.points){
    BIC2[K]=obBIC(K,NFVC_ans2$timeline,NFVC_ans2$est_coef,test_Y2,test_X,test_Tpoints,test_W,ecov22p)
  }
  K2 = which.min(BIC2)
  
  gammap = gammap_NFVC(NFVC_ans1,NFVC_ans2,NFVC_ans1$timeline,test_Y,test_X,test_Tpoints,test_W)
  
  sci1 = matrix(0,nrow=N/5,ncol=K1)
  sci2 = matrix(0,nrow=N/5,ncol=K2)
  for(i in 1:(N/5)){
    for(k in 1:K1){
      sci1[i,k] = e1evaluep[k]*t(e1efunp[,k])%*%solve(ecov11p)%*%gammap[i,1,]
    }
    for(k in 1:K2){
      sci2[i,k] = e2evaluep[k]*t(e2efunp[,k])%*%solve(ecov22p)%*%gammap[i,2,]
    }
  }
  
  phisum = matrix(0,nrow=(N*2/5),ncol=max.points)
  for (m in 1:max.points){
    for (k in 1:K1){
      phi1 = sci1[,k]*e1efunp[k,m]
    }
    for (k in 1:K2){
      phi2 = sci2[,k]*e2efunp[k,m]
    }
    phi = rbind(phi1,phi2)
    phim = matrix(phi,ncol=1)
    phisum[,m] = phisum[,m]+phim
  }
  
  
  right.mat1 = v1%*%beta.seq1+phisum[(1:(N/5)),]
  right.mat2 = v1%*%beta.seq2+phisum[-(1:(N/5)),]
  Y1 = right.mat1+rho.mat11*(test_W%*%right.mat1)+(rho.mat11^2)*
    (test_W%*%test_W%*%right.mat1)+(rho.mat11^3)*(test_W%*%test_W%*%test_W%*%right.mat1)
  Y2 = right.mat2+rho.mat22*(test_W%*%right.mat2)+(rho.mat22^2)*
    (test_W%*%test_W%*%right.mat2)+(rho.mat22^3)*(test_W%*%test_W%*%test_W%*%right.mat2)
  Y.all=rbind(Y1,Y2)
  

  A=test_Y-Y.all
  SSE = mean(t(A)%*%A)
  mean_testY = apply(test_Y,2,mean)
  B = test_Y-mean_testY
  SST = mean(t(B)%*%B)
  R2 = 1-SSE/SST
  
  return(list(R2=R2))
}


###########FVCSIM_pred########

rule.thumb_FVCSIM = function(x){
  0.79*min(quantile(x,0.75,na.rm=T)-quantile(x,0.25,na.rm=T),sd(x,na.rm=T))*length(x)^(-1/5)
}

est_alphat=function(Y,X,z,s,bns,beta){
  n = nrow(Y)/2
  Y1 = Y[1:n,]
  Y2 = Y[-(1:n),]
  g1 = bns[,1]
  g2 = bns[,2]
  beta1=beta[,1]
  beta2=beta[,2]
  
  diff1 = Y1[,1]-z%*%beta1
  diff2 = Y2[,1]-z%*%beta2
  
  h11=rule.thumb_FVCSIM(diff1)
  h21=rule.thumb_FVCSIM(diff2)
 
  max.points=ncol(Y)
  sm=seq(0,1,length=max.points)
  Wh1 = matrix(0,nrow=4,ncol=4)
  for(i in 1:n){
    for(m in 1:max.points){
      K1 = dnorm((sm[m]-s)/h11)/h11
      
      D_im1 = c(X[i,],X[i,]*(sm[m]-s)/h11)
      Wh1 = Wh1+ K1*D_im1%*%t(D_im1)
    }
  }
  
  Wh2 = matrix(0,nrow=4,ncol=4)
  for(i in 1:n){
    for(m in 1:max.points){
      K2 = dnorm((sm[m]-s)/h21)/h21
      D_im2 = c(X[i,],X[i,]*(sm[m]-s)/h21)
      Wh2 = Wh2+ K2*D_im2%*%t(D_im2)
    }
  }
  
  absum1 = matrix(0,nrow=4,ncol=1)
  for(i in 1:n){
    for(m in 1:max.points){
      K1 = dnorm((sm[m]-s)/h11)/h11
      D_im1 = c(X[i,],X[i,]*(sm[m]-s)/h11)
      absum1 = absum1 +K1*D_im1*(Y1[i,m]-g1[i])
    }
  }
  ab1 = solve(Wh1)%*%absum1
  
  absum2 = matrix(0,nrow=4,ncol=1)
  for(i in 1:n){
    for(m in 1:max.points){
      K2 = dnorm((sm[m]-s)/h21)/h21
      D_im2 = c(X[i,],X[i,]*(sm[m]-s)/h21)
      absum2 = absum2 +K2*D_im2*(Y2[i,m]-g2[i])
    }
  }
  ab2 = solve(Wh2)%*%absum2
  
  a = cbind(ab1[1:2,],ab2[1:2,])
  
  return(a)
  
}

est_alpha<-function(Y,X,z,bns,beta){
  max.points=ncol(Y)
  s = seq(from=0,to=1,length.out=max.points)
  alpha=matrix(0,nrow=max.points,ncol=4)
  for(t in 1:max.points){
    alpha[t,] = est_alphat(Y,X,z,s[t],bns,beta)
  }
  return(alpha)
}

est_ab<-function(Y,X,z,z_est,ans,beta,verbose=F){
  
  n1 = nrow(Y)/2
  Y1=Y[1:n1,]
  Y2=Y[-(1:n1),]
  max.points=ncol(Y)
  
  alpha1s.seq = ans[,1:2]
  alpha2s.seq = ans[,3:4]
  Y1t = Y1- X%*%t(alpha1s.seq)
  Y2t = Y2- X%*%t(alpha2s.seq)
  Yt=rbind(Y1t,Y2t)
  
  beta1 = beta[,1]
  beta2 = beta[,2]
  
  diff1 = Y1t[,1]-z%*%beta1
  diff2 = Y2t[,1]-z%*%beta2
  
  h11=rule.thumb_FVCSIM(diff1)
  h21=rule.thumb_FVCSIM(diff2)

  mx1_z=cbind(rep(1,n1),(z-z_est)%*%beta1/h11)
  mx2_z=cbind(rep(1,n1),(z-z_est)%*%beta2/h21)
  k1_z=dnorm((z-z_est)%*%beta1/h11)/h11
  k2_z=dnorm((z-z_est)%*%beta2/h21)/h21
  
  ##eati a1,b1
  
  sig_sum1=matrix(0,ncol=2,nrow=2)
  for(i in 1:n1)
  {
    sig_sum1=sig_sum1+max.points*(k1_z[i])*(mx1_z[i,])%*%t(mx1_z[i,])
  }
  sig_inv1=solve(sig_sum1)
  
  re1=c(0,0)
  for(i in 1:n1){
    for(j in 1:max.points){
      re1=re1+k1_z[i]*(mx1_z[i,])*Y1t[i,j]
    }
  }
  res1=sig_inv1%*%re1
  
  ##esti a2,b2
  
  sig_sum2=matrix(0,ncol=2,nrow=2)
  for(i in 1:n1){
    sig_sum2=sig_sum2+max.points*(k2_z[i])*(mx2_z[i,])%*%t(mx2_z[i,])
  }
  sig_inv2 <- solve(sig_sum2)
  
  re2=c(0,0)
  for(i in 1:n1){
    for(j in 1:max.points){
      re2=re2+k2_z[i]*(mx2_z[i,])*Y2t[i,j]
    }
  }
  res2=sig_inv2%*%re2
  
  rest=cbind(res1,res2)
  
  return(c(res1[1,],res2[2,]))
}

est_abt<-function(Y,X,z,ans,beta){
  n=nrow(Y)/2
  hh=matrix(0,nrow=n,ncol=4)
  for(i in 1:n){
    hh[i,]=est_ab(Y,X,z,z[i,],ans,beta)
  }
  hh_a=hh[,1:2]
  hh_b=hh[,3:4]
  return(hh)
}

###esti beta
est_beta<-function(Y,X,z,ans,bns,beta){
  n1 = nrow(Y)/2
  Y1=Y[1:n1,]
  Y2=Y[-(1:n1),]
  max.points=ncol(Y)
  
  alpha1s.seq = ans[,1:2]
  alpha2s.seq = ans[,3:4]
  Y1t = Y1- X%*%t(alpha1s.seq)
  Y2t = Y2- X%*%t(alpha2s.seq)
  Yt=rbind(Y1t,Y2t)
  
  beta1=beta[,1]
  beta2=beta[,2]
  
  diff1 = Y1t[,1]-z%*%beta1
  diff2 = Y2t[,1]-z%*%beta2
  
  h11=rule.thumb_FVCSIM(diff1)
  h21=rule.thumb_FVCSIM(diff2)
  # h11=dpill(Y1t,z%*%beta1)
  # h21=dpill(Y2t,z%*%beta2)
  
  est_a=bns[,1:2]
  est_b=bns[,3:4]
  
  ##esti beta1
  Omega_1=matrix(0,nrow=3,ncol=3)
  
  for (it in 1:n1){
    sum_1 = 0
    for (i in 1:n1){
      delta_z = z[i,]-z[it,]
      Kt = as.vector(dnorm(delta_z%*%beta1/h11)/h11)
      sum_1 = sum_1+Kt
    }
  }
  
  for(i in 1:n1){
    for(j in 1:n1){
      delta_z=z[i,]-z[j,]
      k1=as.vector(dnorm(delta_z%*%beta1/h11)/h11)/(sum_1/n1)
      Omega_1=Omega_1+max.points*(est_b[i,1])^2*(k1)*delta_z%*%t(delta_z)/(h11^2)
    }
  }
  Omega_inv_1=solve(Omega_1)
  
  beta_s_1=c(0,0,0)
  
  for(i in 1:n1){
    for(j in 1:n1){
      for(t in 1:max.points){
        delta_z=z[i,]-z[j,]
        k1=as.vector(dnorm(delta_z%*%beta1/h11)/h11)/(sum_1/n1)
        beta_s_1=beta_s_1+est_b[i,1]*(delta_z/h11)*k1*(Y1t[i,t]-est_a[i,1])
      }
    }
  }
  beta_1=Omega_inv_1%*%beta_s_1
  
  ##esti beta2
  
  for (it in 1:n1){
    sum_2 = 0
    for (i in 1:n1){
      delta_z = z[i,]-z[it,]
      Kt = as.vector(dnorm(delta_z%*%beta2/h21)/h21)
      sum_2 = sum_2+Kt
    }
  }
  
  Omega_2=matrix(0,nrow=3,ncol=3)
  for(i in 1:n1){
    for(j in 1:n1){
      delta_z=z[i,]-z[j,]
      k2=as.vector(dnorm(delta_z%*%beta2/h21)/h21)/(sum_2/n1)
      Omega_2=Omega_2+max.points*(est_b[i,2])^2*(k2)*delta_z%*%t(delta_z)/(h21^2)
    }
  }
  Omega_inv_2=solve(Omega_2)
  
  beta_s_2=c(0,0,0)
  
  for(i in 1:n1){
    for(j in 1:n1){
      for(t in 1:max.points){
        delta_z=z[i,]-z[j,]
        k2=as.vector(dnorm(delta_z%*%beta2/h21)/h21)/(sum_2/n1)
        beta_s_2=beta_s_2+est_b[i,2]*(delta_z/h21)*k2*(Y2t[i,t]-est_a[i,2])
      }
    }
  }
  beta_2=Omega_inv_2%*%beta_s_2
  
  
  beta_11=as.vector(beta_1/as.vector(sqrt(t(beta_1)%*%(beta_1))))
  beta_21=as.vector(beta_2/as.vector(sqrt(t(beta_2)%*%(beta_2))))
  return(cbind(beta_11,beta_21))
  
}

iteration<-function(max_iter,tol,Y,X,z,g0){
  
  ans = est_alpha(Y,X,z=z,g0,beta=beta0)
  bns = est_abt(Y,X,z=z,ans=ans,beta=beta0)
  cns = est_beta(Y,X,z=z,ans=ans,bns=bns,beta=beta0)
  bnst = est_abt(Y,X,z=z,ans=ans,beta=cns)
  
  beta = cns
  iter = 0
  
  while(iter<max_iter){
    
    beta_old = beta
    anst = est_alpha(Y,X,z=z,bns=bnst,beta=beta_old)
    bnstt = est_abt(Y,X,z=z,ans=anst,beta=beta_old)
    cnstt = est_beta(Y,X,z=z,ans=anst,bns=bnstt,beta=beta_old)
    beta = cnstt
    
    if (sum((beta - beta_old)^2,na.rm=TRUE) < tol) {
      break
    }
    
    iter = iter+1
    
  }
  list(alpha=anst, g = bnstt ,beta = cnstt)
  
}

est_cov<-function(Y,X,z,t,t_est,ans,bns,beta){
  
  n=nrow(Y)/2
  Y1=Y[1:n,]
  Y2=Y[-(1:n),]
  max.points=ncol(Y)
  alpha1s.seq = ans[,1:2]
  alpha2s.seq = ans[,3:4]
  Y1t = Y1- X%*%t(alpha1s.seq)
  Y2t = Y2- X%*%t(alpha2s.seq)
  Yt=rbind(Y1t,Y2t)
  
  est_a=bns[,1:2]
  est_b=bns[,3:4]
  
  Yi_1=Y1t-est_a[,1]
  Yi_2=Y2t-est_a[,2]
  
  beta1=beta[,1]
  beta2=beta[,2]
  
  diff1 = Y1t[,1]-z%*%beta1
  diff2 = Y2t[,1]-z%*%beta2
  
  h12=rule.thumb_FVCSIM(diff1)
  h22=rule.thumb_FVCSIM(diff2)
  # h12=dpill(Y1t,z%*%beta1)
  # h22=dpill(Y2t,z%*%beta2)
  # 
  k1_t=dnorm((t-t_est)/h12)/h12
  mx1_t=cbind(rep(1,max.points),(t-t_est)/h12)
  k2_t=dnorm((t-t_est)/h22)/h22
  mx2_t=cbind(rep(1,max.points),(t-t_est)/h22)
  
  ##esti eta_1
  sum1=matrix(0,nrow=2,ncol=2)
  for(m in 1:max.points){
    sum1=sum1+k1_t[m]*mx1_t[m,]%*%t(mx1_t[m,])
  }
  sum_inv1=solve(sum1)
  
  sum2=matrix(0,nrow=n,ncol=2)
  for(m in 1:max.points){
    sum2=sum2+k1_t[m]*Yi_1[,m]%*%t(mx1_t[m,])
  }
  
  eta1=sum2%*%sum_inv1
  
  ##esti eta_2
  sum3=matrix(0,nrow=2,ncol=2)
  for(m in 1:max.points){
    sum3=sum3+k2_t[m]*mx2_t[m,]%*%t(mx2_t[m,])
  }
  sum_inv2=solve(sum3)
  
  sum4=matrix(0,nrow=n,ncol=2)
  for(m in 1:max.points){
    sum4=sum4+k2_t[m]*Yi_2[,m]%*%t(mx2_t[m,])
  }
  
  eta2=sum4%*%sum_inv2
  
  return(cbind(eta1[,1],eta2[,1]))
  
}

est_eta<-function(Y,X,z,t,ans,bns,beta){
  
  n=nrow(Y)/2
  max.points=ncol(Y)
  eta=array(0,dim=c(n,2,max.points))
  for(s in 1:max.points){
    ts=seq(0.01,0.99,length=max.points)
    eta[,,s]=est_cov(Y,X,z,t,ts[s],ans,bns,beta)
  }
  return(eta)
}

est_covt<-function(Y,X,z,t,ans,bns,beta){
  
  n = nrow(Y)/2
  num = nrow(ans)
  eta=array(0,dim=c(n,num,2))
  for(s in 1:num){
    ts=seq(0.01,0.99,length=num)
    eta[,s,]=est_cov(Y,X,z,t,ts[s],ans,bns,beta)
  }
  
  c11=matrix(0,num,num)
  c12=matrix(0,num,num)
  c21=matrix(0,num,num)
  c22=matrix(0,num,num)
  
  for(s in 1:num){
    for(t in 1:num){
      suu = matrix(0,2,2)
      for(i in 1:n){
        suu=suu+eta[i,s,]%*%t(eta[i,t,])
      }
      Rhat=suu/n
      c11[s,t]=Rhat[1,1]
      c12[s,t]=Rhat[1,2]
      c21[s,t]=Rhat[2,1]
      c22[s,t]=Rhat[2,2]
    }
  }
  cov1122=list(cov11=c11,cov22=c22)
  return(cov1122)
}

pred_Zhu = function(dat){
  
  Y1 = dat$Y1
  Y2 = dat$Y2
  X = dat$v1
  N=nrow(X)
  
  t_grid <- seq(0, 1, length= N)
  train_size = floor(4*N/5)
  train_grid <- seq(1,train_size)  
  test_grid <- seq(train_size+1,N)
  
  train_Y1 <- Y1[train_grid,]
  test_Y1 <- Y1[test_grid,]
  
  train_Y2 <- Y2[train_grid,]
  test_Y2 <- Y2[test_grid,]
  
  train_Y = rbind(train_Y1,train_Y2)
  test_Y = rbind(test_Y1,test_Y2)
  
  train_X = X[train_grid,]
  test_X = X[test_grid,]

  t=seq(0,1,length=max.points)
  g00 = g0[(1:N*4/5),]
  esti = iteration(3,1e-5,train_Y,train_X,train_z,g00)
  ans=esti$alpha
  beta=esti$beta
  bns=est_abt(test_Y,test_X,test_z,ans,beta)
  eta = est_eta(test_Y,test_X,test_z,t,ans,bns,beta)
  
  ecovp = est_covt(test_Y,test_X,test_z,t,ans,bns,beta)
  ecov11p = ecovp$cov11
  ecov22p = ecovp$cov22
  e1evaluep = eigen(ecov11p)$values
  e1efunp = eigen(ecov11p)$vectors
  e2evaluep = eigen(ecov22p)$values
  e2efunp = eigen(ecov22p)$vectors
  
 
  e1evalue_po = e1evaluep[e1evaluep>=0]
  total_var1 = sum(e1evalue_po)
  var1_contri= e1evalue_po/total_var1
  cumu_var1 = cumsum(var1_contri)
  e2evalue_po = e2evaluep[e2evaluep>=0]
  total_var2 = sum(e2evalue_po)
  var2_contri= e2evalue_po/total_var2
  cumu_var2 = cumsum(var2_contri)
  threshold1=0.7
  threshold2=1
  K1 = min(which(cumu_var1>=threshold1 & cumu_var1<threshold2))
  K2 = min(which(cumu_var2>=threshold1 & cumu_var2<threshold2))

  sci1 =matrix(0,nrow=N/5,ncol=max.points)
  for(m in 1:max.points){
    sci1 = sci1+eta[,1,m]%*%t(e1efunp[,m])*(1/max.points)
  }
  sci2 =matrix(0,nrow=N/5,ncol=max.points)
  for(m in 1:max.points){
    sci2 = sci2+eta[,2,m]%*%t(e2efunp[,m])*(1/max.points)
  }
  
  phisum = matrix(0,nrow=(N*2/5),ncol=max.points)
  for (m in 1:max.points){
    for (k in 1:K1){
      phi1 = sci1[,k]*e1efunp[k,m]
    }
    for (k in 1:K2){
      phi2 = sci2[,k]*e2efunp[k,m]
    }
    phi = rbind(phi1,phi2)
    phim = matrix(phi,ncol=1)
    phisum[,m] = phisum[,m]+phim
  }
  
  
  Y1 = test_X%*%t(ans[,1:2])+bns[,1]+phisum[1:(N/5),]
  Y2 = test_X%*%t(ans[,3:4])+bns[,2]+phisum[-(1:(N/5)),]
  
  Y.all=rbind(Y1,Y2)
  
  SSE=sum((test_Y-Y.all)^2)
  SST=sum((test_Y-mean(test_Y))^2)
  R2 = 1-SSE/SST
  return(list(R2=R2))
}


