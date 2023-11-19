#code for example 1 in wls by Y. Zuo 11/06/23

#part (a)###################################
x=c(5, 5.5, 4, 3.5, 3, 2.5, -2)
y=c(-.5, -.5, 6, 4, 2.4, 2, .5)
m1=Zn=cbind(x,y)
dev.off()
par(mfrow=c(1,2))
plot(m1, pch=16, xlim=c(-3,8), ylim=c(-2, 7))
 abline(0,0, lty=1, col=1, lwd=1)
 abline(0, 1, lty=2, col=2, lwd=1)
 text(x, y-0.5, as.character(seq(1:7)))
legend(x = "topleft", legend=c("y=0", "y=x"), col=1:2, lty=1:2, cex=1, text.font=1, bg='lightblue',title="reference lines")

plot(m1, pch=16, xlim=c(-3,8), ylim=c(-2, 7))
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


legend(x = "topleft", legend=c("LTS", "LST", "LS", "WLS"),    
       col=1:4, lty=1:4, cex=1,
       title="Line types", text.font=1, bg='lightblue')

########part(b)############################################
#alpha=1     #cut-off value alpha
#N=500               # the total iteration number allowed
#c=0                #c=0,initial beta is rep(1, p);c=1,it is by LS; c=2, it is by LTS.

n=80; p=2#
epsilon=0.3; n1=floor(epsilon*n)
m1=rmvnorm(n, mean=(rep(0, p)),sigma=matrix(c(1, 0.88, 0.88, 1), ncol=2, byrow=T))
m111=m1
m1=m111

par(mfrow=c(1,2))  
plot(m1[,1], m1[,2], xlim=c(-3,8), ylim=c(-2, 7), xlab="x-axis", ylab="y-axis", pch=20)
fit0<-ltsReg(m1[,1:(p-1)], m1[,p])
abline(fit0$coefficients[[1]], fit0$coefficients[[2]], col=1, lty=1, lwd=1)
abline(AA1_main(m1,1, 0, 500, 0.001), col=2, lty=2, lwd=1)
fit2<-lm(m1[,2]~m1[,1]) #p is assume to be 2
abline(fit2, col=3, lty=3, lwd=1)
wls_beta=wls(m1, 5, 6, 0.00001)
abline(wls_beta, lty=4, col=4, lwd=1)
legend(x = "topleft", legend=c("LTS", "LST", "LS", "WLS"),    
       col=1:4, lty=1:4, cex=1,
       title="Line types", text.font=1, bg='lightblue')

if(n1>=1)
{
  m2=rmvnorm(n1, mean=c(7,-2),sigma=diag(rep(0.1, p)))
  m1[sample(1:n,n1),]<-m2 
}

m222=m1
m1=m222

plot(m1[,1], m1[,2], xlim=c(-3,8), ylim=c(-2, 7), xlab="x-axis", ylab="y-axis", pch=20)
fit1<-ltsReg(m1[,1:(p-1)], m1[,p])
abline(fit1$coefficients[[1]], fit1$coefficients[[2]], col=1, lty=1, lwd=1)
abline(AA1_main(m1,1, 0, 500, 0.001), col=2, lty=2, lwd=1)
fit2<-lm(m1[,2]~m1[,1]) #p is assume to be 2
abline(fit2, col=3, lty=3, lwd=1)
wls_beta=wls(m1, 5, 6, 0.00001)
abline(wls_beta, lty=4, col=4, lwd=1)

#text(4., -1.6, expression ('LTS'))
#text(6, 5.7, expression('LST'))
#text(4, -1, expression('LS'))

legend(x = "topleft", legend=c("LTS", "LST", "LS", "WLS"),    
       col=1:4, lty=1:4, cex=1,
       title="Line types", text.font=1, bg='lightblue')

