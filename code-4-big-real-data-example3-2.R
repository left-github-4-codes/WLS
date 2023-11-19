#code for the big real data set by Y. Zuo on 01/19/22
#rm(list=ls())
# adding LMS in the revision on 04/27/23 by Y. Zuo
library(mvtnorm)
library(robustbase)
library(compiler)
enableJIT(1)

library(MASS)

data(Boston, package="MASS")
#attach(Boston)
options(warn=-1)

###########################################################################
##### just entire data set and calulate the total time consumed and sample variance
###########################################################################

z=Boston
z=data.matrix(z)
n=dim(z)[1]; p=dim(z)[2]

R=RepN=1000
#n=100 #tuning this n for your need
#p=10
alpha=3
N_ls=1
delta=1
Cc=6
k=6
cutoff=0.0001
#bet0=c(rep(1,5), rep(-1,p-5))

beta_aa1=beta_aa2=beta_aa3=beta_aa4=matrix(0, nrow=R, ncol=p)  #aa1-lst, aa2-wls, aa3-LTS, aa3-ls
t_aa1=t_aa2=t_aa3=t_aa4=0

for (i in 1:RepN)
{
  t1=Sys.time()  
  beta_aa1[i,]= AA1_main_new_lst_v2(z,alpha, delta, N_ls)
  t2=Sys.time()-t1
  t_aa1=t_aa1+t2
  # 
  t1=Sys.time()  
  beta_aa2[i,]=wls(z, Cc, k, cutoff)
  t2=Sys.time()-t1
  t_aa2=t_aa2+t2
  
  t1=Sys.time()  
  fit1<-ltsReg(z[,1:(p-1)], z[,p])
  beta_aa3[i,]=as.numeric(fit1$coefficients) 
  t2=Sys.time()-t1
  t_aa3=t_aa3+t2
  
  t1=Sys.time()  
  fit2<-lm(z[,p]~z[,1:(p-1)])
  beta_aa4[i,]=as.numeric(fit2$coefficients) 
  t2=Sys.time()-t1
  t_aa4=t_aa4+t2
  
  print(i)
}

beta_aa1_mean=colMeans(beta_aa1)
deviat_beta_aa1=beta_aa1-matrix(beta_aa1_mean, byrow=T, nrow=RepN, ncol=p)
beta_aa2_mean=colMeans(beta_aa2)
deviat_beta_aa2=beta_aa2-matrix(beta_aa2_mean, byrow=T,nrow=RepN, ncol=p)
beta_aa3_mean=colMeans(beta_aa3)
deviat_beta_aa3=beta_aa3-matrix(beta_aa3_mean, byrow=T,nrow=RepN, ncol=p)
beta_aa4_mean=colMeans(beta_aa4)
deviat_beta_aa4=beta_aa4-matrix(beta_aa4_mean, byrow=T,nrow=RepN, ncol=p)


EMSE_aa1=sum(deviat_beta_aa1*deviat_beta_aa1)/RepN
EMSE_aa2=sum(deviat_beta_aa2*deviat_beta_aa2)/RepN
EMSE_aa3=sum(deviat_beta_aa3*deviat_beta_aa3)/RepN
EMSE_aa4=sum(deviat_beta_aa4*deviat_beta_aa4)/RepN

print(c(R,n,p,N_ls,alpha, Cc, k, cutoff, Del, Kap) )
print( c(EMSE_aa1, EMSE_aa2, EMSE_aa3, EMSE_aa4))
print(c(t_aa1,t_aa2,t_aa2,t_aa4))
print(EMSE_aa4/c(EMSE_aa1, EMSE_aa2, EMSE_aa3,  EMSE_aa4))


############################################################################
############## sampling n1 points from the entire data set
############################################################################

#z=Boston
#Bos_minus=z[,c(-3,-7)] # delete indus and age two predictor since they are not significant with p-value <0.001
#z=Bos_minus
z=Boston
z=zz=data.matrix(z)
n=dim(z)[1]; p=dim(z)[2]

RepN=1000; R=RepN
n1=250 # tune this to 50, 100, or 200
alpha=1 # one or three?
c=0
cut_off=10^{-3}
#epsilon=cut_off
N_ls=1
delta=1
#gamma=10^2
beta_aa1=beta_aa2=beta_aa3=beta_aa4=matrix(0, nrow=R, ncol=p)  #aa1-lst, aa2-wls, aa3-LTS, aa3-ls
t_aa1=t_aa2=t_aa3=t_aa4=0

i=1
for (i in 1:RepN)
{
  index=sample(1:n, n1)
  z=Boston[index,]

  t1=Sys.time()  
  beta_aa1[i,]= AA1_main_new_lst_v2(z,alpha, delta, N_ls)
  t2=Sys.time()-t1
  t_aa1=t_aa1+t2
  # 
  t1=Sys.time()  
  beta_aa2[i,]=wls(z, Cc, k, cutoff)
  t2=Sys.time()-t1
  t_aa2=t_aa2+t2
  
  t1=Sys.time()  
  fit1<-ltsReg(z[,1:(p-1)], z[,p])
  beta_aa3[i,]=as.numeric(fit1$coefficients) 
  t2=Sys.time()-t1
  t_aa3=t_aa3+t2
  
  t1=Sys.time()  
  fit2<-lm(zz[,p]~zz[,1:(p-1)])
  beta_aa4[i,]=as.numeric(fit2$coefficients) 
  t2=Sys.time()-t1
  t_aa4=t_aa4+t2
  
  print(i)
}

beta_aa1_mean=colMeans(beta_aa1)
deviat_beta_aa1=beta_aa1-matrix(beta_aa1_mean, byrow=T, nrow=RepN, ncol=p)
beta_aa2_mean=colMeans(beta_aa2)
deviat_beta_aa2=beta_aa2-matrix(beta_aa2_mean, byrow=T,nrow=RepN, ncol=p)
beta_aa3_mean=colMeans(beta_aa3)
deviat_beta_aa3=beta_aa3-matrix(beta_aa3_mean, byrow=T,nrow=RepN, ncol=p)
beta_aa4_mean=colMeans(beta_aa4)
deviat_beta_aa4=beta_aa4-matrix(beta_aa4_mean, byrow=T,nrow=RepN, ncol=p)


EMSE_aa1=sum(deviat_beta_aa1*deviat_beta_aa1)/RepN
EMSE_aa2=sum(deviat_beta_aa2*deviat_beta_aa2)/RepN
EMSE_aa3=sum(deviat_beta_aa3*deviat_beta_aa3)/RepN
EMSE_aa4=sum(deviat_beta_aa4*deviat_beta_aa4)/RepN

print(c(R,n,p,N_ls,alpha, Cc, k, cutoff, Del, n1) )
print( c(EMSE_aa1, EMSE_aa2, EMSE_aa3, EMSE_aa4))
print(c(t_aa1,t_aa2,t_aa2,t_aa4))
print(EMSE_aa4/c(EMSE_aa1, EMSE_aa2, EMSE_aa3,  EMSE_aa4))

########################################################################################
#### Second part (a) dropping age and indus and (b) one part for estimation and the other part for prediction
########################################################################################
z=Boston
z1=z=data.matrix(z)
#Bos_minus=zz[,c(-3,-7)] # delete indus and age two predictor since they are not significant with p-value <0.001
#z=Bos_minus
n=dim(z)[1]; p=dim(z)[2]
RepN=1000; R=RepN
n1=m=250  # tune this 200, 300, 400
alpha=3 # one or three?
c=0
cut_off=10^{-3}
N=100
#gamma=10^2
beta_aa1=beta_aa2=beta_aa3=beta_aa4=matrix(0, nrow=R, ncol=p)  #aa1-lst, aa2-wls, aa3-LTS, aa3-ls
t_aa1=t_aa2=t_aa3=t_aa4=0

#i=1 #aa3 is actually replace by LTS
rss_aa1=rss_aa2=rss_aa3=rss_aa4=0
res_aa1=res_aa2=res_aa3=res_aa4=matrix(0, nrow=n-m, ncol=1)
t_aa1=t_aa2=t_aa3=t_aa4=t1=t2=0
i=1
for (i in 1:RepN)
{
  index=sample(1:n, m)
  z=z3=Boston[index,]
  z3=data.matrix(z3)
  zz=Boston[-index,] # could drop this and use all data points

  t1=Sys.time()  
 # beta_aa1[i,]=AA1_main(z, alpha, c, N, cut_off )
  beta_aa1[i,]= AA1_main_new_lst_v2(z,alpha, delta, N_ls)
  res_aa1=zz[,p]-as.matrix(cbind(matrix(1, nrow=n-m), zz[, 1:(p-1)]))%*%matrix(beta_aa1[i,], nrow=p, ncol=1)
  
  rss_temp=sum(res_aa1*res_aa1)
  rss_aa1=rss_aa1+rss_temp
  t2=Sys.time()-t1
  t_aa1=t_aa1+t2
  
  t1=Sys.time()  
  beta_aa2[i,]=wls(z, Cc, k, cutoff)
  res_aa2=zz[,p]-as.matrix(cbind(matrix(1, nrow=n-m), zz[, 1:(p-1)]))%*%matrix(beta_aa2[i,], nrow=p, ncol=1)
 
  rss_temp=sum(res_aa2*res_aa2)
  rss_aa2=rss_aa2+rss_temp
  t2=Sys.time()-t1
  t_aa2=t_aa2+t2
  
  t1=Sys.time()  
  fit1<-ltsReg(z[,1:(p-1)], z[,p])
  beta_aa3[i,]=as.numeric(fit1$coefficients) 
  res_aa3=zz[,p]-as.matrix(cbind(matrix(1, nrow=n-m), zz[, 1:(p-1)]))%*%matrix(beta_aa3[i,], nrow=p, ncol=1)
  #res_aa3=Bos[,p]-as.matrix(cbind(matrix(1, nrow=n), Bos[, 1:(p-1)]))%*%matrix(beta_aa3[i,], nrow=p, ncol=1)
 
  rss_temp=sum(res_aa3*res_aa3)
  rss_aa3=rss_aa3+rss_temp
  t2=Sys.time()-t1
  t_aa3=t_aa3+t2
  
  
  t1=Sys.time()  
  fit2<-lm(z3[,p]~z3[,1:(p-1)])
  beta_aa4[i,]=as.numeric(fit2$coefficients) 
  res_aa4=zz[,p]-as.matrix(cbind(matrix(1, nrow=n-m), zz[, 1:(p-1)]))%*%matrix(beta_aa4[i,], nrow=p, ncol=1)
  
  rss_temp=sum(res_aa4*res_aa4)
  rss_aa4=rss_aa4+rss_temp
  t2=Sys.time()-t1
  t_aa4=t_aa4+t2
  
  print(i)
}  

beta_aa1_mean=colMeans(beta_aa1)
deviat_beta_aa1=beta_aa1-matrix(beta_aa1_mean, byrow=T, nrow=RepN, ncol=p)
beta_aa2_mean=colMeans(beta_aa2)
deviat_beta_aa2=beta_aa2-matrix(beta_aa2_mean, byrow=T,nrow=RepN, ncol=p)
beta_aa3_mean=colMeans(beta_aa3)
deviat_beta_aa3=beta_aa3-matrix(beta_aa3_mean, byrow=T,nrow=RepN, ncol=p)
beta_aa4_mean=colMeans(beta_aa4)
deviat_beta_aa4=beta_aa4-matrix(beta_aa4_mean, byrow=T,nrow=RepN, ncol=p)


EMSE_aa1=sum(deviat_beta_aa1*deviat_beta_aa1)/RepN
EMSE_aa2=sum(deviat_beta_aa2*deviat_beta_aa2)/RepN
EMSE_aa3=sum(deviat_beta_aa3*deviat_beta_aa3)/RepN
EMSE_aa4=sum(deviat_beta_aa4*deviat_beta_aa4)/RepN



print(c(R, n,p, m,c, N,alpha,cut_off, n1) )
print( c(EMSE_aa1, EMSE_aa2, EMSE_aa3, EMSE_aa4))
print(c(t_aa1, t_aa2, t_aa3, t_aa4))
print(EMSE_aa4/c(EMSE_aa1, EMSE_aa2, EMSE_aa3,  EMSE_aa4))
print(c(rss_aa1, rss_aa2, rss_aa3, rss_aa4)/RepN)
