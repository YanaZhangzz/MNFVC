####three network structure
getDyadW<-function(N, N1, delta, normalize = T)                                           
{
  A = matrix(0, nrow = N, ncol = N)                                                                   
  
  ind = which(upper.tri(A), arr.ind = T)                                                               
  indM = ind[sample(1:nrow(ind), N*N1),]                                                               
  A[indM] = 1                                                                                          
  A[indM[,2:1]] = 1                                                                                    
  
  ind1 = which(A==0&upper.tri(A), arr.ind = T)                                                         
  indS = ind1[sample(1:nrow(ind1), N^delta),]                                                          
  tmp = sample(1:nrow(indS), floor(N^delta/2))                                                         
  indS[tmp,] = indS[tmp, 2:1]                                                                          
  A[indS] = 1                                                                                           
  diag(A) = 0                                                                                          
  if (!normalize)
    return(A)
  W = A/rowSums(A)                                                                                     
  # W = as(W, "dgCMatrix")
  return(W)
}


getPowerLawW<-function(N, alpha, normalize = T)           
{
  Nfollowers = rpldis(N, 1, alpha)
  
  
  A = sapply(Nfollowers, function(n) {             
    vec = rep(0, N) 
    vec[sample(1:N, min(n,N))] = 1  
    return(vec)
  })
  diag(A) = 0 
  ind = which(rowSums(A)==0)                        
  for (i in ind)
  {
    A[i, sample(setdiff(1:N,i), 3)] = 1             
  }
  if (!normalize)
    return(A)
  W = A/rowSums(A)
  return(W)
}


getBlockW<-function(N, Nblock, normalize = T)                                                 
{
  if (N%%Nblock==0){                                                                
    isDiagList = rep(list(matrix(1, nrow = N/Nblock, ncol = N/Nblock)), Nblock)     
    mList = rep(list(matrix(rbinom((N/Nblock)^2, size = 1, prob = 0.9*N^{-1}),   
                            nrow = N/Nblock, ncol = N/Nblock)), Nblock)
  }
  else
  {
    
    isDiagList = rep(list(matrix(1, nrow = floor(N/Nblock), ncol = floor(N/Nblock))), Nblock-1) 
    isDiagList[[length(Nblock)]] = matrix(1, nrow = N%%Nblock, ncol = N%%Nblock)
    
    
    mList = rep(list(matrix(rbinom(floor(N/Nblock)^2, size = 1, prob = 0.9*N^{-1}),     
                            nrow = floor(N/Nblock), ncol = floor(N/Nblock))), Nblock-1)
    mList[[Nblock]] = matrix(rbinom(floor(N/Nblock)^2, size = 1, prob = 0.9*N^{-1}),   
                             nrow = floor(N/Nblock), ncol = floor(N/Nblock))
  }
  isDiag = bdiag(isDiagList)  
  
  offDiag = which(isDiag == 0, arr.ind = T)                                         
  mList = lapply(mList, function(M){
    ind = which(rowSums(M)==0)
    if (length(ind)>0)
      M[cbind(ind, sample(1:nrow(M), length(ind)))] = 1
    return(M)
  })
  bA = bdiag(mList)
  
  bA[offDiag] = rbinom(nrow(offDiag), size = 1, prob = 0.3/N)                        
  bA = as.matrix(bA)
  upperInd = which(upper.tri(bA), arr.ind = T)
  
  bA[upperInd[,2:1]] = bA[upper.tri(bA)]
  diag(bA) = 0
  
  
  ind = which(rowSums(bA)==0)                                 
  for (i in ind)
  {
    bA[i, sample(setdiff(1:N,i), 3)] = 1                      
  }
  
  if (!normalize)
    return(bA)
  W = bA/rowSums(bA)                                         
  return(W)
}

####generate data
get.func.network = function(n, min.points, max.points, net_type, Y_out_point = NULL)
{
  
  if (net_type=="1")
  {W = getDyadW(n, N1 = 10, delta = 1.2, normalize = T)}
  
  if (net_type=="2")
  {W = getBlockW(n,10)}
  if (net_type=="3")
  {W = getPowerLawW(n, alpha = 2.5, normalize = T) }
  
  s = seq(from = 0, to = 1, length.out = max.points) 
  
  Sige1 = Diagonal(n = 2, x = c(0.4,0.6))
  Sige1[1,2] = 0.1
  Sige1[2,1] = 0.1
  eps = matrix(rnorm(n*2*max.points),nrow=n*2)
  ee = eigen(Sige1)
  In = Diagonal(n = N, x = 1)
  eps = kronecker(t(sqrt(ee$values)*t(ee$vectors)), In)%*%eps
  b11 = matrix(s^2, nrow = 1)
  b12= matrix((1-s)^2, nrow = 1)
  b21= matrix(5*(s-0.5)^2, nrow = 1)
  b22 = matrix(s^0.5, nrow = 1)
  b.mat = rbind(b11,b21,b12,b22)
  
  ZSigma = 0.5^abs(outer(1:2, 1:2,"-")) 
  v1 = mvrnorm(n = N, mu = rep(0,nrow(ZSigma)), Sigma = ZSigma)
  v1_hat = kronecker(diag(2),v1)
  
  xi11 = matrix(rnorm(n, 0, 1.2), n, 1)
  xi12 = matrix(rnorm(n, 0, 0.6), n, 1)
  xi1.seq = cbind(xi11, xi12)
  xi21 = matrix(rnorm(n, 0, 1), n, 1)
  xi22 = matrix(rnorm(n, 0, 0.5), n, 1)
  xi2.seq = cbind(xi21, xi22)
  
  phi11.seq = matrix(sin(2*s*pi) * sqrt(2), nrow = 1)
  phi12.seq = matrix(cos(2*s*pi) * sqrt(2), nrow = 1)
  phi1.seq = rbind(phi11.seq,phi12.seq)
  phi21.seq = matrix(cos(2*s*pi) * sqrt(2), nrow = 1)
  phi22.seq = matrix(sin(2*s*pi) * sqrt(2), nrow = 1)
  phi2.seq = rbind(phi21.seq,phi22.seq)
  
  eta1.seq = xi1.seq%*%phi1.seq
  eta2.seq = xi2.seq%*%phi2.seq
  eta.seq = rbind(eta1.seq,eta2.seq)
  
  
  d11 = matrix(0.2*sin(2*s*pi)+0.3,nrow = 1)
  ##case_a
  # d12 = matrix(2*(s^2-s)+0.5,nrow=1)
  # d21 = matrix(2*(s-s^2),nrow=1)
  ##case_b
  d12 = matrix(-0.3*sin(2*pi*s)-0.2,nrow=1)
  d21 = matrix(0.3*cos(2*pi*s)+0.2,nrow=1)
  ##case_c
  # d12 = matrix(0.5 * sin(pi * (s - 0.5)),nrow=1)
  # d21 = matrix(0.5 - 1 / (1 + exp(-5 * (s - 0.5))),nrow=1)
  ##case_d
  # d12 = matrix(0*s,nrow=1)
  # d21 = matrix(0*s,nrow=1)
  d22 = matrix(0.2*cos(2*s*pi)+0.3,nrow=1)
  d.mat = rbind(d11,d12,d21,d22)
  
  Ymat = matrix(0,nrow=2*N,ncol=max.points)
  for(m in 1:max.points){
    Xbeta = v1_hat%*%b.mat[,m]
    eps1 = eps[,m]+Xbeta
    D = matrix(d.mat[,m],nrow=2,byrow=TRUE)
    DkW = kronecker(t(D),(W))
    DkW2 = DkW%*%DkW
    DkW3 = DkW2%*%DkW
    DkW4 = DkW3%*%DkW
    DkW5 = DkW4%*%DkW
    DkW6 = DkW5%*%DkW
    
    Yvec = eps1+DkW%*%eps1+DkW2%*%eps1+DkW3%*%eps1+DkW4%*%eps1+DkW5%*%eps1+DkW6%*%eps1
    Ymat[,m]=Yvec
  }
  
  
  Y.all=Ymat
  Y1 = Ymat[(1:N),]
  Y2 = Ymat[-(1:N),]
  
  if (!is.null(Y_out_point)){
    Y_true = get.Ytrue(xi11, xi12, xi21, xi22, t = Y_out_point, W, X, z)
  }
  else{
    Y_true = NULL
  }
  
  Y = matrix(NA, 2*n, max.points)
  
  Tpoints = matrix(NA, n, max.points)
  
  for (i in 1:n) {
    total.points = sample(min.points:max.points, 1, prob = rep(1:(max.points-min.points+1)/(max.points-min.points+1)))
    
    ind = sort(sample(1:max.points, size = total.points, replace = F))
    
    Tpoints[i, 1:total.points] = s[ind]
    
    Y[i, 1:total.points] = Y.all[i,ind]
    Y[n+i, 1:total.points] = Y.all[n+i,ind]
  }
  
  Y1NA = Y[(1:N),]
  Y2NA = Y[-(1:N),]
  list(Y = Y, Y1 = Y1, Y2 = Y2, Y1N = Y1NA,Y2N = Y2NA,Tpoints = Tpoints, v1 = v1, W = W, Y.all = Y.all, Y_true = Y_true, time_seq = s)
}

####MNFVC
est_d_beta <- function(Y,X,Tpoints,W,D,B,Sige1,num.points)
{
  # dimensions
  N = nrow(Y)/2
  
  #Y---2n*maxpoints
  Y1s=Y[1:N,]
  Y2s=Y[-(1:N),]
  
  #bianxing
  W.hat = kronecker(diag(2),(W))
  X_bolang=kronecker(diag(2),X)
  
  tlsub = seq(0,1,length=num.points)
  re = matrix(0,nrow=num.points,ncol=8)
  for(ii in 1:num.points){
    
    t_est = tlsub[ii]
    T_diff = Tpoints - t_est
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
    
    zeta = apply(Y*Ker_mat.hat, 1, function(s){return(mean(s, na.rm = T))})
    zetat = t((zeta/fs.hat)%*%t(matrix(1, nrow = 2*N, ncol = 1)))
    
    Ys = zeta/fs.hat
    
    Ym = matrix(Ys,ncol=2)
    
    result = MSAR.Lse(Ym,X,W,D,B,Sige1)
    dvec =as.vector(result$theta[(1:4)]) 
    betavec = as.vector(result$theta[-(1:4)])
    
    ddvec = c(dvec[1],dvec[3],dvec[2],dvec[4])
    re[ii,] = c(ddvec,betavec)
  }
  
  list(est_coef = re, timeline = tlsub)
  
  
}

obtain.kernel = function(Y, Tpoints, t)
{
  n = dim(Tpoints)[1]
  max.points = dim(Tpoints)[2]
  Ys = Tpoints[1,]
  h = 1*density(Ys[Ys<999])$bw
  
  K = matrix(999, n, max.points)
  for (i in 1:n) {
    ti = Tpoints[i,]
    total.points = sum(ti<999)
    K[i, 1:total.points] = gaussian.kernel.mat(t-ti[1:total.points], h)
    
  }
  list(K = K, h = h)
}

est_covariance <- function(est_coef, timeline, Y1 , Y2 , X, Tpoints, W)
{
  N = dim(W)[1]
  
  d11tt = est_coef[,1]
  d12tt = est_coef[,2]
  d21tt = est_coef[,3]
  d22tt = est_coef[,4]
  dtt = est_coef[,1:4]
  
  b11tt = est_coef[,5]
  b21tt = est_coef[,6]
  b12tt = est_coef[,7]
  b22tt = est_coef[,8]
  betatt = est_coef[,5:dim(est_coef)[2]]
  
  
  num.points = length(timeline)
  
  
  Stt=rep(list(matrix(NA,nrow=2*N,ncol=2*N)),num.points)
  for (s in 1:num.points) {Stt[[s]]=diag(2*N)-kronecker(W,cbind(rbind(d11tt[s],d12tt[s]),rbind(d21tt[s],d22tt[s])))}
  
  C110 = matrix(0, num.points, num.points)
  C120 = matrix(0, num.points, num.points)
  C210 = matrix(0, num.points, num.points)
  C220 = matrix(0, num.points, num.points)
  
  ns = rowSums(!is.na(Tpoints))
  
  for (s in 1:num.points) {
    for(t in 1:num.points){
      
      T_diff = Tpoints - timeline[s]
      h = rule.thumb(T_diff[1,])/3
      
      
      Ker_mat_s = gaussian.kernel(T_diff, h)
      Y1_s = rowSums(Y1*Ker_mat_s, na.rm = T)
      Y2_s = rowSums(Y2*Ker_mat_s, na.rm = T)
      
      
      T_diff = Tpoints - timeline[t]
      h = rule.thumb(T_diff[1,])/3
      
      Ker_mat_t = gaussian.kernel(T_diff, h)
      Y1_t = rowSums(Y1*Ker_mat_t, na.rm = T)
      Y2_t = rowSums(Y2*Ker_mat_t, na.rm = T)
      
      Y1_2_st = rowSums(Y1^2*Ker_mat_s*Ker_mat_t, na.rm = T)
      Y12_st = rowSums(Y1*Y2*Ker_mat_s*Ker_mat_t, na.rm = T)
      Y2_2_st = rowSums(Y2^2*Ker_mat_s*Ker_mat_t, na.rm = T)
      
      Vst11.n = Y1_s%*%t(Y1_t) - diag(Y1_2_st)
      Vst12.n = Y1_s%*%t(Y2_t) - diag(Y12_st)
      Vst21.n = Y2_s%*%t(Y1_t) - diag(Y12_st)
      Vst22.n = Y2_s%*%t(Y2_t) - diag(Y2_2_st)
      
      K_s = rowSums(Ker_mat_s, na.rm = T)
      K_t = rowSums(Ker_mat_t, na.rm = T)
      K_st = rowSums(Ker_mat_s*Ker_mat_t, na.rm = T)
      Vst.d = K_s%*%t(K_t) - diag(K_st)
      
      Vst_11 = Vst11.n/Vst.d
      Vst_12 = Vst12.n/Vst.d
      Vst_21 = Vst21.n/Vst.d
      Vst_22 = Vst22.n/Vst.d
      Vst=rbind(cbind(Vst_11,Vst_12),cbind(Vst_21,Vst_22))
      
      
      suu=matrix(0,2,2)
      for(i in 1:N)
      {suu = suu+(Stt[[s]][(2*i-1):(2*i),])%*%Vst%*%t(Stt[[t]][(2*i-1):(2*i),])}
      
      B_s=rbind(cbind(b11tt[s],b12tt[s]),cbind(b21tt[s],b22tt[s]))
      B_t=rbind(cbind(b11tt[t],b12tt[t]),cbind(b21tt[t],b22tt[t]))
      
      C0_st = (suu - t(X%*%B_s)%*%(X%*%B_t))/N    
      
      C110[s,t]=C0_st[1,1]
      C120[s,t]=C0_st[1,2]
      C210[s,t]=C0_st[2,1]
      C220[s,t]=C0_st[2,2]
      
    }
  }
  cov1122 = list(cov11 = C110,cov22 = C220)
  
  return(cov1122)
}

get.Ytrue <-function(xi11, xi12, xi21, xi22, s, W, v1)
{
  max.points = length(s)
  n = nrow(W)
  
  b11 = matrix(s^2, nrow = 1)
  b12= matrix((1-s)^2, nrow = 1)
  b21= matrix(5*(s-0.5)^2, nrow = 1)
  b22 = matrix(s^0.5, nrow = 1)
  b.mat = rbind(b11,b21,b12,b22)
  
  ep1 = matrix(rnorm(n*max.points, sd = 0.2), n) 
  ep2 = matrix(rnorm(n*max.points, sd = 0.15), n) 
  ep=rbind(ep1,ep2)
  
  phi11.seq = matrix(sin(2*s*pi) * sqrt(2), nrow = 1)
  phi12.seq = matrix(cos(2*s*pi) * sqrt(2), nrow = 1)
  phi1.seq = rbind(phi11.seq,phi12.seq)
  phi21.seq = matrix(cos(2*s*pi) * sqrt(2), nrow = 1)
  phi22.seq = matrix(sin(2*s*pi) * sqrt(2), nrow = 1)
  phi2.seq = rbind(phi21.seq,phi22.seq)
  eta1.seq = xi1.seq%*%phi1.seq
  eta2.seq = xi2.seq%*%phi2.seq
  eta.seq = rbind(eta1.seq,eta2.seq)
  
  d11 = matrix(0.2*sin(2*s*pi)+0.3,nrow = 1)
  ##case_a
  # d12 = matrix(2*(s^2-s)+0.5,nrow=1)
  # d21 = matrix(2*(s-s^2),nrow=1)
  ##case_b
  d12 = matrix(-0.3*sin(2*pi*s)-0.2,nrow=1)
  d21 = matrix(0.3*cos(2*pi*s)+0.2,nrow=1)
  ##case_c
  # d12 = matrix(0.5 * sin(pi * (s - 0.5)),nrow=1)
  # d21 = matrix(0.5 - 1 / (1 + exp(-5 * (s - 0.5))),nrow=1)
  ##case_d
  # d12 = matrix(0*s,nrow=1)
  # d21 = matrix(0*s,nrow=1)
  d22 = matrix(0.2*cos(2*s*pi)+0.3,nrow=1)
  d.mat = rbind(d11,d12,d21,d22)
  
  Ymat = matrix(0,nrow=2*N,ncol=max.points)
  for(m in 1:max.points){
    Xbeta = v1_hat%*%b.mat[,m]
    eps1 = ep[,m]+Xbeta+eta.seq[,m]
    D = matrix(d.mat[,m],nrow=2,byrow=TRUE)
    DkW = kronecker(t(D),(W))
    DkW2 = DkW%*%DkW
    DkW3 = DkW2%*%DkW
    DkW4 = DkW3%*%DkW
    DkW5 = DkW4%*%DkW
    DkW6 = DkW5%*%DkW
    
    Yvec = eps1+DkW%*%eps1+DkW2%*%eps1+DkW3%*%eps1+DkW4%*%eps1+DkW5%*%eps1+DkW6%*%eps1
    Ymat[,m]=Yvec
  }
  
  Ytrue=Ymat
  
  return(Ytrue)
}

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

gaussian.kernel.mat = function(u, h)
{
  
  u = u/h
  return(exp(-0.5*u^2)/sqrt(2*pi))
}

rule.thumb = function(x)
{
  1.06*sd(x, na.rm = T)*length(x)^{-1/5}
}

gaussian.kernel = function(u, h)
{
  
  u = u/h
  return(exp(-0.5*u^2)/sqrt(2*pi))
}
