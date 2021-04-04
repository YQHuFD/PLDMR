rm(list=ls())
library(parallel)
library(mvtnorm)

#=========================================================Functions======================================================

gene_snp=function(n,p,rho)
{
  Sigma=toeplitz(rho^(0:(p-1)))
  z=rmvnorm(n,rep(0,p),Sigma,method='chol')
  u=pnorm(z)
  maf=runif(p,0.2,0.4)
  snp=do.call(cbind,lapply(1:ncol(u),function(i)qbinom(u[,i],2,maf[i])))
  colnames(snp)=paste('g',1:p,sep='')
  return(snp)
}

getsumdata_x=function(iv,exposure)
{
  n=nrow(iv)
  p=ncol(iv)
  gamma_tilde=rep(0,p)
  se_gt=rep(0,p)
  t_gt=rep(0,p)
  for (j in 1:p) {
    ols.fit=summary(lm(exposure~iv[,j]))
    gamma_tilde[j]=ols.fit$coef[2,1]
    se_gt[j]=ols.fit$coef[2,2]
    t_gt[j]=ols.fit$coef[2,3]
  }
  cbind(gamma_tilde,se_gt,t_gt)
}

getsumdata_y=function(iv,outcome)
{
  n=nrow(iv)
  p=ncol(iv)
  Gamma_tilde=rep(0,p)
  se_Gt=rep(0,p)
  t_Gt=rep(0,p)
  for (j in 1:p) {
    ols.fit=summary(lm(outcome~iv[,j]))
    Gamma_tilde[j]=ols.fit$coef[2,1]
    se_Gt[j]=ols.fit$coef[2,2]
    t_Gt[j]=ols.fit$coef[2,3]
  }
  cbind(Gamma_tilde,se_Gt,t_Gt)
}

getsumdata_x_stan=function(iv,exposure)
{
  n=nrow(iv)
  p=ncol(iv)
  gamma_tilde=rep(0,p)
  se_gt=rep(0,p)
  t_gt=rep(0,p)
  
  iv=scale(iv)
  for (j in 1:p) {
    ols.fit=summary(lm(exposure~iv[,j]))
    gamma_tilde[j]=ols.fit$coef[2,1]
    se_gt[j]=ols.fit$coef[2,2]
    t_gt[j]=ols.fit$coef[2,3]
  }
  cbind(gamma_tilde,se_gt,t_gt)
}

getsumdata_y_stan=function(iv,outcome)
{
  n=nrow(iv)
  p=ncol(iv)
  Gamma_tilde=rep(0,p)
  se_Gt=rep(0,p)
  t_Gt=rep(0,p)
  
  iv=scale(iv)
  for (j in 1:p) {
    ols.fit=summary(lm(outcome~iv[,j]))
    Gamma_tilde[j]=ols.fit$coef[2,1]
    se_Gt[j]=ols.fit$coef[2,2]
    t_Gt[j]=ols.fit$coef[2,3]
  }
  cbind(Gamma_tilde,se_Gt,t_Gt)
}

PLDMR=function(g,x,y,init,maxreplic=20){
  g=scale(g,scale=F)
  x=scale(x,scale=F)
  y=scale(y,scale=F)
  gtg=crossprod(g)
  decom=eigen(gtg,symmetric = TRUE)
  sqrt_gtg=decom$vectors%*%diag(sqrt(decom$values))%*%t(decom$vectors)
  gamma_hat=solve(gtg,crossprod(g,x))
  Gamma_hat=solve(gtg,crossprod(g,y))
  wG=sqrt_gtg%*%Gamma_hat
  wg=sqrt_gtg%*%gamma_hat
  wg=cbind(sqrt_gtg%*%rep(1,length(wg)),wg)
  fit.ldmr=summary(lm(wG~-1+wg))$coef
  fit.pldmra=summary(lm(Gamma_hat~gamma_hat))$coef
  fn=function(para){
    a0h=para[1];bth=para[2];r2=para[3];sy2=para[4] 
    -sum(log(dnorm(qt%*%Gamma_hat,qt%*%cbind(rep(1,p),gamma_hat)%*%c(a0h,bth),sqrt(sy2*(r2+lam)))+1e-9))
  }
  qt=t(decom$vectors)
  lam=1/decom$values
  mle=spg(c(init[1],init[2],init[3],init[4]),fn,lower=c(-Inf,-Inf,1e-7,1e-7),upper=c(Inf,Inf,Inf,Inf),control=list(maxit=5000),quiet=T)
  replic=0
  while (mle$convergence!=0&replic<maxreplic) {
    mle=spg(mle$par,fn,lower=c(-Inf,-Inf,1e-7,1e-7),upper=c(Inf,Inf,Inf,Inf),control=list(maxit=5000),quiet=T)
    replic=replic+1
  }
  r2=mle$par[3]
  S=solve(gtg)+diag(1,p,p)*r2
  decom_S=eigen(solve(S),symmetric = TRUE)
  w=decom_S$vectors%*%diag(sqrt(decom_S$values))%*%t(decom_S$vectors)
  wG=w%*%Gamma_hat
  wg=w%*%cbind(rep(1,p),gamma_hat)
  fit=summary(lm(wG~-1+wg))
  fit.pldmr=fit$coef
  
  sa2_true=init[3]*init[4]
  S=solve(gtg)*sigma_y^2+sa2_true*diag(1,p,p)
  decom_S=eigen(solve(S),symmetric = TRUE)
  w=decom_S$vectors%*%diag(sqrt(decom_S$values))%*%t(decom_S$vectors)
  wG=w%*%Gamma_hat
  wg=w%*%cbind(rep(1,p),gamma_hat)
  fit.pldmrt=summary(lm(wG~-1+wg))$coef
  return(list(fit.ldmr,fit.pldmra,fit.pldmr,fit.pldmrt))
}

LDA.MREgger<-function(X,Y,W){
  bX<-cbind(1,X)
  bread<-solve(crossprod(bX,W)%*%bX)
  theEsts<-bread%*%crossprod(bX,W%*%Y)
  theresid<-c(Y-theEsts[1]-X*theEsts[2])
  Sig.Est<-c(crossprod(theresid,W%*%theresid))/(length(X)-2)
  finresults<- cbind(theEsts,diag(bread)*Sig.Est)
  TestStat<-theEsts/sqrt(finresults[,2])
  Pvals<-2*pt(abs(TestStat),df = nrow(bX)-2,lower.tail = F)
  return(cbind(finresults,TestStat,Pvals))
}

TWMR<-function(beta,gamma,C,nx,ny){
  S<-t(beta)%*%solve(C)%*%beta
  H<-(1-1/sqrt(3781))*S+(1/sqrt(3781))*diag(length(S[,1]))
  alpha<-solve(H)%*%(t(beta)%*%solve(C)%*%gamma)
  
  alpha<-as.vector(alpha)
  
  C_inv <- solve(C)
  GCG_inv <- t(beta) %*% solve(C) %*% beta
  GCG_inv<-(1-1/sqrt(5000))*GCG_inv+(1/sqrt(5000))*diag(length(GCG_inv[,1]))
  GCG_inv<-solve(GCG_inv)
  
  
  df_dg <- GCG_inv %*% t(beta) %*% C_inv
  df_dG <- (GCG_inv %x% (t(gamma) %*% C_inv %*% ((beta %*% GCG_inv %*% t(beta)) %*% C_inv + diag(nrow(beta))))) + ((-t(gamma) %*% C_inv %*% beta %*% GCG_inv) %x% (GCG_inv %*% t(beta) %*% C_inv))
  J <- cbind(df_dG, df_dg)
  
  SEs<-c(rep(1/sqrt(nx),length(beta[1,])*length(beta[,1])),rep(1/sqrt(ny),length(gamma[,1])))
  R<-diag(length(beta[1,])+1)
  Sigma <- (SEs %*% t(SEs)) * (C %x% R)   
  V <- J %*% Sigma %*% t(J)
  se<- sqrt(V[1,1])
  
  N=length(beta[,1])
  Ngene=length(beta[1,])
  Z<-alpha/se
  pval<-2*pnorm(abs(Z),lower.tail=FALSE)
  return(c(alpha,se,pval))
}
#==================================================Simulation===================================================

set.seed(2020)

n=c(10000)
p=25

beta=0
gamma=runif(p+1,0.5,4)
rho_g=c(0.6)
type=c(0.1)
sa=c(0.2)
sigma_x=2
sigma_y=2

r=10000

for (i in 1:length(n)) {
  for (j in 1:length(rho_g)) {
    for (k in 1:length(type)) {
      for (l in 1:length(sa)) {
        g.all=gene_snp(n[i]+5000,p,rho_g[j])
        g=g.all[1:n[i],]
        modified.g=cbind(rep(1,n[i]),g)
        R=cor(g)
        g1=g[1:(n[i]/2),]
        g2=g[(n[i]/2+1):n[i],]
        gs=g.all[(n[i]+1):(n[i]+5000),]
        R_ts=cor(gs)
        
        cl=makeCluster(35)
        clusterExport(cl,c('PLDMR','LDA.MREgger','TWMR','getsumdata_x','getsumdata_x_stan','getsumdata_y','getsumdata_y_stan','n','p','beta','gamma','sigma_x','sigma_y','g','g1','g2','type','sa','modified.g','R','R_ts','i','k','l'))
        clusterEvalQ(cl,{
          library(mvtnorm)
          library(MendelianRandomization)
          library(BB)
          library(MR.LDP)
          library(mr.raps)})
        
        res=do.call(rbind,parLapply(cl,1:r,function(m){
          alpha=rnorm(p+1,type[k],sa[l])
          e=rmvnorm(n[i],c(0,0),matrix(c(sigma_x^2,0.5*sigma_x*sigma_y,0.5*sigma_x*sigma_y,sigma_y^2),2,2))
          x=modified.g%*%gamma+e[,1]
          y=modified.g%*%alpha+x*beta+e[,2]
          
          x1=x[1:(n[i]/2),]
          y2=y[(n[i]/2+1):n[i],]
          summary_x_ts=getsumdata_x(g1,x1)
          summary_y_ts=getsumdata_y(g2,y2)
          summary_x_ts.stan=getsumdata_x_stan(g1,x1)
          summary_y_ts.stan=getsumdata_y_stan(g2,y2)
          
          mri=mr_input(bx=summary_x_ts[,1],bxse=summary_x_ts[,2],by=summary_y_ts[,1],byse=summary_y_ts[,2],correlation=R_ts)
          res_MRE=mr_egger(mri,correl=T)
          res_PLDMR=PLDMR(g,x,y,c(type[k],beta,sa[l]^2/sigma_y^2,sigma_y^2))
          res_LDAMR=LDA.MREgger(summary_x_ts.stan[,1],summary_y_ts.stan[,1],solve(R_ts)%*%(R_ts*(summary_y_ts.stan[,2]%*%t(summary_y_ts.stan[,2])))%*%solve(R_ts))
          res_TWMR=TWMR(as.matrix(summary_x_ts[,1]),as.matrix(summary_y_ts[,2]),R_ts,n[i]/2,n[i]/2)
          res_RAPS=mr.raps.overdispersed.robust(summary_x_ts[,1],summary_y_ts[,1],summary_x_ts[,2],summary_y_ts[,2],'huber')
          
          gamma0 = rep(0.01, p)
          alpha0 = rep(0.01, p)
          sgga2 = 0.01
          sgal2 = 0.01
          beta0 = 0
          maxIter = 1000
          lam = 0.1
          epsStopLogLik = 1e-7
          SimMRLDP_Hb = MRLDP_SimPXvb(summary_x_ts[,1], summary_y_ts[,1], summary_x_ts[,2], summary_y_ts[,2], gamma0, alpha0, beta0, sgga2, sgal2, R_ts, 0, epsStopLogLik, maxIter, model = 2)
          SimMRLDP_H0 = MRLDP_SimPXvb(summary_x_ts[,1], summary_y_ts[,1], summary_x_ts[,2], summary_y_ts[,2], gamma0, alpha0, beta0, sgga2, sgal2, R_ts, 1, epsStopLogLik, maxIter, model = 2)
          beta_mrldp = SimMRLDP_Hb$beta0
          tstat = 2*(SimMRLDP_Hb$tstat - SimMRLDP_H0$tstat)
          pval_mrldp=pchisq(tstat,1,lower.tail = F)
          
          return(c(beta_mre=res_MRE@Estimate,beta_ldamr=res_LDAMR[2,1],beta_twmr=res_TWMR[1],beta_raps=res_RAPS$beta.hat,beta_mrldp=beta_mrldp,
                   beta_ldmr=res_PLDMR[[1]][2,1],beta_pldmra=res_PLDMR[[2]][2,1],beta_pldmr=res_PLDMR[[3]][2,1],beta_pldmrt=res_PLDMR[[4]][2,1],
                   p_mre=res_MRE@Causal.pval,p_ldamr=res_LDAMR[2,4],p_twmr=res_TWMR[3],p_raps=res_RAPS$beta.p.value,p_mrldp=pval_mrldp,
                   p_ldmr=res_PLDMR[[1]][2,4],p_pldmra=res_PLDMR[[2]][2,4],p_pldmr=res_PLDMR[[3]][2,4],p_pldmrt=res_PLDMR[[4]][2,4]))
        }))
        write.table(res,paste(c('simu_beta0a',type[k],'sa',sa[l],'LD',rho_g[j],'n',n[i],'.txt'),collapse=''),sep='\t',row.names=F,quote=F)
        stopCluster(cl)
      }
    }
  }
}

