# code for wls by Y. Zuo on 10/04/23

library(mvtnorm)
library(robustbase)
library(compiler)
enableJIT(1)

#get an initial beta, could be lst of ZZ23 or LTS, or LS,
wls<-function(Zn, Cc, k, cutoff)
{# Zn is the dataset of (x', y) with the 1st p-1 columns for x and pth column for y, cutoff=10^{-12} 
  
  #INITIAL beta[0], alpha[0], v[0]
  n=dim(Zn)[1]; p=dim(Zn)[2]
  #fit1<-ltsReg(Zn[,1:(p-1)], Zn[,p]) # could be LS or LST
  #beta0=as.numeric(fit1$coefficients)

  #beta0=get_beta_from_p_points(Zn) #; print(beta0)
  #beta0[is.nan(beta0)]<-1 #; print(beta0)
  #fit2<-lm(Zn[,p]~Zn[,1:(p-1)])
  #beta0=as.numeric(fit2$coefficients)
  
  #lst_beta=AA1_main(Zn,1,0,500, 0.001) #; print(lst_beta)
  #beta0=lst_beta
  beta0=AA1_main_new_lst_v2(Zn, 3, 1, 1)
  #beta0=lst_beta #; print(beta0)
  return(beta0)
  #beta0=rep(1,p)
  
  v0=-delta_beta(beta0, Zn, Cc,k) #; print("v0"); print(v0)
  norm_delta=sqrt(sum(v0^2)) #; print("norm_delta"); print(norm_delta)
  if(norm_delta<=cutoff){return(beta0)}
  alpha0=backtracking_search(beta0, Zn, Cc, k, v0, -v0)
  #alpha0=-t(delta_beta(beta0, Zn, Cc, k))%*%v0/(t(v0)%*% temp%*%v0)

  # beta[,1] and v[,1]
  beta1=beta0+alpha0*v0 #; print(c(beta1))
  OV_new=O_beta(beta1,Zn, Cc, k); OV_old=O_beta(beta0, Zn, Cc,k)
  diff_OV=abs(OV_new-OV_old) #; print("diff_OV");print(diff_OV)
  if(diff_OV<=cutoff){return(beta0)}
  
  norm_delta=sqrt(sum(delta_beta(beta1, Zn, Cc,k)[,1]^2)) 
  if(norm_delta<=cutoff){return(beta0)}
  beta=matrix(0, nrow=p, ncol=n); alpha=rep(0, n)
  v=matrix(0, ncol=n, nrow=p)
  beta[,1]=matrix(beta1); deltaOb1=delta_beta(beta[,1], Zn, Cc, k)
  deltaOb0=delta_beta(beta0, Zn, Cc, k)
  v[,1]=-deltaOb1+ (t(deltaOb1)%*%deltaOb1/t(deltaOb0)%*%deltaOb0)[1,1]*v0
  k=1
  while(diff_OV>=cutoff | norm_delta>=cutoff)
  {
   if ((k)>=5){break}
   else{
   DObk=delta_beta(beta[,k], Zn, Cc, k)
   alpha[k]=backtracking_search(beta[,k], Zn, Cc, k, v[,k], DObk)
  # alpha[k]=-t(DObk)%*% v[,k]/t(v[,k])%*%delta_sq_beta(beta[,k],Zn,Cc,k)%*%v[,k]
 
   beta[,k+1]=beta[,k]+alpha[k]*v[,k]
   DOBkp1=delta_beta(beta[,k+1], Zn, Cc, k)
   v[,k+1]=-DOBkp1+(t(DOBkp1)%*%DOBkp1/t(DObk)%*%DObk)[1,1]*v[,k]
 
   OV_new=O_beta(beta[,k+1], Zn, Cc,k); OV_old=O_beta(beta[,k], Zn, Cc,k)
 
   diff_OV=abs(OV_new-OV_old); if(diff_OV<=cutoff){return(beta[,k+1])}
   norm_delta=sqrt(sum(delta_beta(beta[,k+1], Zn, Cc, k)[,1]^2))
   if(norm_delta<=cutoff){return(beta[,k])}
      }
   k=k+1  
  }
  return(beta[,k])
}
####################################################################################################
w<-function(r, Cc, k) #weight function
{
  w<-1*(abs(r)<=Cc)+1*(abs(r)>Cc)*(exp(-K*(1-Cc/abs(r))^2)-exp(-K))/(1-exp(-K))
}
#############################
w_deriv<-function(r, Cc, k) #first-order derivative of the weight function
{
  w1std=0
  if (abs(r)>Cc){w1std=
   1*(abs(r)>Cc)*((-2*k*Cc)/r^2*exp(-K*(1-Cc/abs(r))^2))*(1-Cc/abs(r))*sign(r)/(1-exp(-k))}
  return(w1std)
}
#########################################
w_2nd_deriv<-function(x, Cc, k) #second-order derivative of the weight function
{
  w2ndd=0
  if(abs(x)>Cc){w2ndd=
  1*(abs(x)>Cc)*(-2*k*Cc)/(1-exp(-k))/((x)^3)*exp(-K*(1-Cc/abs(x))^2)*((3*Cc/abs(x)-2)-2*k*Cc*(1-Cc/abs(x))^2/abs(x))}
  return(w2ndd)
}
###############################
O_beta1<-function(beta, Zn, Cc, k) #Capital C small k objective function
{# beta is a p by 1 vector, Zn is a data matrix with first (p-1) columns as X and the last column as y

n=dim(Zn)[1]; p=dim(Zn)[2]
x=Zn[,1:(p-1)]; y=Zn[,p] 
rn=O_beta=temp=w_vector=rep(0,n)

xmatrix=cbind(matrix(rep(1,n)), x); fitted=xmatrix %*% matrix(beta) #; print(dim(xmatrix)); print(matrix(beta))
rn = matrix(y)-fitted #;{print("rn"); print(fitted[1:10,]); print(beta)}; # print("rn"); print(c(matrix(y)[1:10]))
#if(fitted[1,]==NaN) {print(fitted[1:10,])}
#print(matrix(y)[1:5]) #; print(cbind(1,  Zn[,1:(p-1)]) %*% matrix(beta)[1:5])
#print(y[1:5]); print(cbind(1,  Zn[,1:(p-1)]) %*% matrix(beta)[1:5])
cstar=median(y*y)
#print("rn^2+cstar"); print(c(rn^2)); print(cstar)
v1=rn^2/cstar #; print(v1[1:5])
#w_vector=as.numeric(lapply(v1, w, Cc=Cc, k=k))
for (i in 1:n) {w_vector[i]=w(v1[i], Cc, k)} #; if(i<=5){print(w_vector[i])}}
#print(w_vector[1:5])
#for (i in 1:n) {temp[i]=w(rn[i]^2/cstar, Cc, k); O_beta[i]=temp[i]*rn[i]^2}

O_beta=w_vector*rn^2 #; print(c(O_beta))
O_beta=sum(O_beta)
#print("O_beta"); print(O_beta)
return (O_beta)    
}
O_beta=cmpfun(O_beta1)
###################################
delta_beta1<-function(beta, Zn, Cc, k) #first-order derivative of objective function w.r.t. beta, gradient
{
  n=dim(Zn)[1]; p=dim(Zn)[2]
  x=Zn[,1:(p-1)]; y=Zn[,p] 
  rn=rep(0, n); delta_beta=matrix(0, nrow=p, ncol=1)
  xmatrix=cbind(matrix(rep(1,n)), x); fitted=xmatrix %*% matrix(beta) 
  rn = matrix(y)-fitted
  cstar=median(y*y)
  
  for (i in 1:n){
    wi=c(1, Zn[i,1:(p-1)])
    wdi=w_deriv(rn[i]^2/cstar, Cc,k); wwi=w(rn[i]^2/cstar, Cc, k) #; print("wdi-wwi"); print(c(wdi,wwi))
    delta_beta=delta_beta+((-2)/cstar*rn[i]*(wdi*rn[i]^2+cstar*wwi))*matrix(wi) #; print(delta_beta)
  } #; print(delta_beta)
return(delta_beta) 
} 
delta_beta=cmpfun(delta_beta1)
###################################
delta_sq_beta1<-function(beta, Zn, Cc, k) #second-order derivative w.r.t. beta,  Hessian matrix
 {
  n=dim(Zn)[1]; p=dim(Zn)[2]
  x=Zn[,1:(p-1)]; y=Zn[,p] 
  rn=rep(0, n); ww=rep(0, n)
  gamma=rep(0, n); W=rep(0, nclon=n, nrow=n)
  xmatrix=cbind(matrix(rep(1,n)), x); fitted=xmatrix %*% matrix(beta) 
  rn = matrix(y)-fitted
  cstar=median(y*y)
  
  for (i in 1:n) {yi=rn[i]^2/cstar; w0=w(yi, Cc, k);w1=w_deriv(yi,Cc,k); w2=w_2nd_deriv(yi,Cc,k)
  tempi=-2*(5*yi*w1+w0+2*yi^2*w2)/cstar; ww[i]=tempi}

  W=diag(ww); x=cbind(rep(1, n), x)
  delta_sq_deriv_beta=t(x)%*%W%*%x
 }
delta_sq_beta=cmpfun(delta_sq_beta1)
#######################################
backtracking_search1<-function(beta, Zn, Cc, k, v, delta_beta) #this one is replaced by exact formula above for alpha[k]
{
#v is v_k in the iteration and delta_beta is the direction of gradient at beta_k
 a=1/4; b=1/2
 t=1; DB=delta_beta
 beta_new=matrix(beta)+t*matrix(v) #; print("beta_new+v"); print(beta_new); print(v)
 Q1=O_beta(beta_new, Zn, Cc, k)
 Q2=O_beta(beta, Zn, Cc, k)+a*t*t(DB)%*%v; Q2=Q2[1,1]
 if(Q1-Q2<=0.0001) {return(t)}
 #print("Q1,Q2");print(Q1); print(Q2)
 k=0
 while(Q1>Q2)
  { k=k+1
  if(k>10){break}
  t=b*t
  beta_new=beta+t*v
  Q1=O_beta(beta_new, Zn, Cc, k)
  Q2=O_beta(beta, Zn, Cc, k)+a*t*t(DB)%*%v; Q2=Q2[1,1]
  if(Q1-Q2<=0.0001) {return(t)}
 }
 return(t)
}
backtracking_search=cmpfun(backtracking_search1)
################################################################

get_beta_from_p_points1=function(z)
{ # input: 
  #    z: is given data matrix with n rows and p columns, 
  #       the last column represents the y coordinates. combine 1 with (p-1) vector to get p vector x
  # output:
  #    beta [determined by p points ((x_i)', y_i): y_i=(x_i)'beta]
  
  n = dim(z)[1]
  p = dim(z)[2]
  
  beta = rep(0, p); z_temp = matrix(0, nrow=p, ncol=p)
  x_temp = matrix(0, nrow=p, ncol=p)
  
  z_temp = z[sample(1:n,p),]          #sampling p rows from z
  x_temp = cbind(1, z_temp[, 1:(p-1)])#combine 1 with x_i to form p vector and take care of intercept term
  
  
  while(det(x_temp==0))               #keep samping until an invetible matrix obtained
  {
    z_temp = z[sample(1:n,p),]          #sampling p rows from z
    x_temp = cbind(1, z_temp[, 1:(p-1)])#combine 1 with x_i to form p vector,taking care of intercept term
  }
  
  y_temp = z_temp[,p]                   #corresponding p y_i in the equation y_i=(x_i)'beta
  
  beta = c(solve(x_temp)%*%y_temp)      #convert to a row vector
  
  return(beta)
}
get_beta_from_p_points=cmpfun(get_beta_from_p_points1)   #speed up when repeatedly call this function

#############################################################################


#testing

 x=c(5, 5.5, 4, 3.5, 3, 2.5, -2)
 y=c(-.5, -.5, 6, 4, 2.4, 2, .5)
 m1=Zn=cbind(x,y)
 dev.off()
 par(mfrow=c(1,1))
 plot(m1, pch=16, xlim=c(-3,8), ylim=c(-2, 7))
# abline(0,0, lty=1, col=1, lwd=1)
# abline(0, 1, lty=2, col=2, lwd=1)
 p=dim(m1)[2]
 fit1<-ltsReg(m1[,1:(p-1)], m1[,p])
 beta_ltsReg=as.numeric(fit1$coefficients) 
 abline(fit1$coefficients[[1]], fit1$coefficients[[2]], col=1, lty=1, lwd=1)
 lst_beta=AA1_main(m1,1,0,500, 0.001)
# lst_beta_1=AA2_main(m1, 1, 0,200, 0.001)
 abline(lst_beta,lty=2, col=2, lwd=1)
# abline(lst_beta_1,lty=4, col=5, lwd=1)
 fit2<-lm(m1[,2]~m1[,1]) #p is assume to be 2
 abline(fit2, col=3, lty=3, lwd=1)
 wls_beta=wls(Zn, 5, 6, 0.00001)
 abline(wls_beta, lty=4, col=4, lwd=1)
 text(x, y-0.5, as.character(seq(1:7)))
 text(7.1, 1.7, expression ('LTS'))
 text(6.7, 5.0, expression('LST'))
 text(7.1, 2.5, expression('LS'))
 text(5.7,6.9, expression('WLS'))
 
 legend(x = "topleft", legend=c("",""), cex=1, text.font=3, bg='lightblue',title="artificial seven point data set")

 legend(x = "topleft", legend=c("LTS", "LST", "LS", "WLS"),    
   col=1:4, lty=1:4, cex=1,
   title="Line types", text.font=3, bg='lightblue')