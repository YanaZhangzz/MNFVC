#####NFVC
est_rho_beta <- function(t_est, Y, X, Tpoints, W, Ytrue = NULL, verbose = F)
{
  N = nrow(Y)
  T_diff = Tpoints - t_est
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
  nu_mat = Y_Ker_vec%*%t(Y_Ker_vec)
  
  I_N = Diagonal(n = N, x = 1)
  tWW = t(W)%*%W
  d_tWW = diag(tWW)
  tW_W = t(W) + W
  
  rho = 0.5
  del = 1
  iter = 0
  while (max(abs(del))>10^{-4} &  iter < 20)
  {
    St = diag(N)-rho*W
    Mt = t(St)%*%St
    Dt = diag(diag(Mt)^{-1})                                                            
    Ut = Dt%*%Mt ### Ui is the i-th column of Ut
    A = Dt%*%t(St)%*%X
    beta = solve(t(A)%*%A)%*%(t(A)%*%rowSums(Ut*zetat))   
    
    if (verbose)
      cat(del, "\n")
    S = I_N - rho*W
    tSS = t(S)%*%S
    tSS_rho = -tW_W + 2*rho*tWW 
    tSS_rho2 = 2*tWW 
    D = 1/diag(tSS)
    D1 = -2*rho*D^2*diag(tWW)
    D2 = -2*D^2*d_tWW + 8*rho^2*D^3*d_tWW^2
    U_mat = D*tSS
    U1_mat = D1*tSS - D*tW_W + 2*rho*D*tWW
    u_mat = t(U1_mat)%*%U_mat
    U2_mat = D2*tSS + 2*D1*tSS_rho+D*tSS_rho2
    u1_mat = t(U2_mat)%*%U_mat + t(U1_mat)%*%U1_mat
    
    Xbeta = X%*%beta
    xi = D*t(S)%*%Xbeta
    xi_rho = (D1*t(S) - D*t(W))%*%Xbeta
    xi_rho2 = (D2*t(S) - 2*D1*t(W))%*%Xbeta
    xi_beta = D*t(S)%*%X
    v_vec_rho = t(U1_mat)%*%xi + t(U_mat)%*%xi_rho
    v1_vec_rho = t(U2_mat)%*%xi + 2*t(U1_mat)%*%xi_rho + t(U_mat)%*%xi_rho2
    
    G1_rho = 2/N*sum(nu_mat*u_mat)
    G2_rho = -2/N*sum(Y_Ker_vec*v_vec_rho)
    G3_rho = 2/N*sum(xi_rho*xi)
    G_rho = G1_rho + G2_rho + G3_rho
    
    H1_rho = 2/N*sum(nu_mat*u1_mat)
    H2_rho = -2/N*sum(Y_Ker_vec*v1_vec_rho)
    H3_rho = 2/N*sum(xi_rho*xi_rho + xi_rho2*xi)
    H_rho = H1_rho + H2_rho + H3_rho
    
    del = H_rho^{-1}*G_rho
    rho = rho - del
    iter = iter + 1
    if(rho<0 | rho>1) {rho = runif(1)}
  }
  
  return(c(rho, beta))
  
}

est_entire_NFVC <- function(Y, X, Tpoints, W, num.points)
{
  tlsub = seq(from = 0, to = 1, length.out = num.points)
  re <- matrix(0, nrow = num.points, ncol = 3)
  for (ii in 1:num.points) {
    re[ii,] = est_rho_beta(tlsub[ii], Y, X, Tpoints, W)
  }
  list(est_coef = re, timeline = tlsub)
}

getCovCST_NFVC1 = function(s)
{
  phi11.seq = matrix(sin(2*s*pi) * sqrt(2), ncol = 1)
  phi12.seq = matrix(cos(2*s*pi) * sqrt(2), ncol = 1)
  
  TCST11 = 1.2*phi11.seq%*%t(phi11.seq)+0.6*phi12.seq%*%t(phi12.seq)
  
  return(TCST11)
}

getCovCST_NFVC2 = function(s)
{
  phi21.seq = matrix(cos(2*s*pi) * sqrt(2), ncol = 1)
  phi22.seq = matrix(sin(2*s*pi) * sqrt(2), ncol = 1)
  
  TCST22 = phi21.seq%*%t(phi21.seq)+0.5*phi22.seq%*%t(phi22.seq)
  
  return(TCST22)
}

gaussian.kernel.mat = function(u, h)
{
  u = u/h
  return(exp(-0.5*u^2)/sqrt(2*pi))
}

est_cov_NFVC <- function(est_coef, timeline,Y, X, Tpoints, W)
{
  N = dim(W)[1]
  rhott = est_coef[,1]
  betatt = est_coef[,2:dim(est_coef)[2]]
  Stt = lapply(rhott, function(rho){diag(N)-rho*W})
  CO = matrix(0, max.points, max.points)
  ns = rowSums(!is.na(Tpoints))
  for (s in 1:max.points) {
    for(t in 1:max.points){
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
