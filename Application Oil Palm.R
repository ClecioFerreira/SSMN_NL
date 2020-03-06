rm(list=ls(all=TRUE))
   
source('functions_ssmn_nl.R')

Oil=read.table(file="Oil.txt",header=T)

##Define non linear function
nlf<-function(x,betas){
resp<- betas[1]/(1 +betas[2]*exp(-betas[3]*x))
return(resp)
}

##Set the response y and covariate x
Y <- Oil$y
x <- Oil$x
n=length(Y)
p=3

# Derivates of nlf(x,beta) in beta

der1 <- deriv3(~beta1/(1+beta2*exp(-beta3*x)),c("beta1","beta2","beta3"),c("beta1","beta2","beta3","x"))

der_eta_beta <-function(x,beta){
# nlf = eta(beta,x)=beta1/(1+beta2*exp(-beta3*x))
n=length(x)
derbetam=matrix(0,length(beta),n)
der2betam=list()
for (i in 1:n){
   derbeta1= as.matrix(1/(1+beta[2]*exp(-beta[3]*x[i])))
   derbeta2=-as.matrix(beta[1]*exp(-beta[3]*x[i])/((1+beta[2]*exp(-beta[3]*x[i]))^2))
   derbeta3=as.matrix(beta[1]*beta[2]*x[i]*exp(-beta[3]*x[i])/((1+beta[2]*exp(-beta[3]*x[i]))^2))
   derbetam[,i]=rbind(derbeta1,derbeta2,derbeta3)# coluna
   #
   der2b1b1=0
   der2b1b2=-exp(-beta[3]*x[i])/((1+beta[2]*exp(-beta[3]*x[i]))^2)
   der2b1b3=beta[2]*x[i]*exp(-beta[3]*x[i])/((1+beta[2]*exp(-beta[3]*x[i]))^2)
   der2b2b2=2*beta[1]*exp(-2*beta[3]*x[i])/((1+beta[2]*exp(-beta[3]*x[i]))^3)
   der2b2b3=beta[1]*x[i]*exp(-beta[3]*x[i])/((1+beta[2]*exp(-beta[3]*x[i]))^3)*(1-beta[2]*exp(-beta[3]*x[i]))
   der2b3b3=beta[1]*beta[2]*(x[i]^2)*exp(-beta[3]*x[i])/((1+beta[2]*exp(-beta[3]*x[i]))^3)*(beta[2]*exp(-beta[3]*x[i])-1)
   der2beta=matrix(0,3,3)
   der2beta[1,2]=der2b1b2
   der2beta[2,1]=der2b1b2
   der2beta[1,3]=der2b1b3
   der2beta[3,1]=der2b1b3
   der2beta[2,2]=der2b2b2
   der2beta[2,3]=der2b2b3
   der2beta[3,2]=der2b2b3
   der2beta[3,3]=der2b3b3
   der2betam[[i]]=der2beta
   }
#
der_eta <-list(derbeta=derbetam,der2beta=der2betam)
return(der_eta)
}


# Initial values
beta0=c(37,44,0.7)
n=length(Y)
nls0 <- nls(Y~beta1/(1+beta2*exp(-beta3*x)), start=c(beta1=beta0[1], beta2=beta0[2],beta3=beta0[3]))
t0=summary(nls0)
beta0=as.matrix(t0$coefficients[,1])
sigma20=as.numeric(sum((Y-nlf(x,beta0))^2)/n)
resid0=t0$residuals #= Y - nlf(x,beta0)
lambda0=skewness(resid0)
theta0=c(beta0,sigma20^2,0)
p=length(beta0)

# NL - Normal
thetaN=normal.nl(Y, x, nlf, beta0)
betaN=thetaN[1:p]
EYN= nlf(x,betaN)
logN=n.logL(thetaN,p,Y,x,nlf)

# SSMN
# SN
b=sqrt(2/pi)
thetaSN <- ssmn.nl(y=Y, x=x, betas=beta0, sigma2=sigma20,shape = lambda0, nlf=nlf,family = "Skew.normal")
thetaSN=thetaSN$theta
ESN=b*sqrt(thetaSN[p+1])*thetaSN[p+2]/sqrt(1+thetaSN[p+2]^2)
betaSN=thetaSN[1:p]
sigma2SN=thetaSN[p+1]
lambdaSN=thetaSN[p+2]
viesSN=b*sqrt(sigma2SN)*lambdaSN/sqrt(1+lambdaSN^2)
EYSN=nlf(x,betaSN)+viesSN
varSN=sigma2SN-viesSN^2
logSN=snn.logL(thetaSN,p,Y,x,nlf)

# STN
thetaST <- ssmn.nl.simul(y=Y, x=x, betas=beta0, sigma2=sigma20,shape = lambda0, nu=20,nlf=nlf,family = "Skew.t")
thetaST=thetaST$theta
betaST=thetaST[1:p]
sigma2ST=thetaST[p+1]
lambdaST=thetaST[p+2]
nuST=thetaST[p+3]
logST=stn.logL(thetaST,p,Y,x,nlf)
# Esp e Var
replic=10000
V = rgamma(replic,(nuST-1)/2,nuST/2)
viesST= b*sqrt(sigma2ST)*lambdaST*sqrt(nuST/2)*gamma((nuST-1)/2)/gamma(nuST/2)*mean((V+lambdaST^2)^(-0.5)) 
EYST=nlf(x,betaST)+viesST
varST = sigma2ST*nuST/(nuST-2)-viesST^2


# Skew Slash
thetaSSL <- ssmn.nl.simul(y=Y, x=x, betas=beta0, sigma2=sigma20,shape = lambda0, nu=20,nlf=nlf,family = "Skew.slash")
thetaSSL=thetaSSL$theta
logSSL=ssl.logL(thetaSSL,p,Y,x,nlf)
betaSSL=thetaSSL[1:p]
sigma2SSL=thetaSSL[p+1]
lambdaSSL=thetaSSL[p+2]
nuSSL=thetaSSL[p+3]
replic=10000
V=rbeta(replic,1,nuSSL-1/2)
viesSSL= b*sqrt(sigma2SSL)*lambdaSSL*nuSSL/(nuSSL-1/2)*mean((V+lambdaSSL^2)^(-0.5)) 
EYSSL=nlf(x,betaSSL)+viesSSL
varSSL = sigma2SSL*nuSSL/(nuSSL-1)-viesSSL^2

# Skew Contaminated Normal
nu0=cbind(c(0.5,0.5))
thetaSCN <- ssmn.nl(y=Y, x=x, betas=beta0, sigma2=sigma20,shape = lambda0, nu=nu0,nlf=nlf,family = "Skew.cn")
thetaSCN=thetaSCN$theta
logSCN=scn.logL(thetaSCN,p,Y,x,nlf) 
betaSCN=thetaSCN[1:p]
sigma2SCN=thetaSCN[p+1]
lambdaSCN=thetaSCN[p+2]
nuSCN=thetaSCN[p+3]
gamaSCN=thetaSCN[p+4]
viesSCN= b*sqrt(sigma2SCN)*lambdaSCN*(nuSCN/sqrt(gamaSCN*(gamaSCN+lambdaSCN^2))+(1-gamaSCN)/sqrt(1+lambdaSCN^2))
EYSCN=nlf(x,betaSCN)+viesSCN
varSCN = sigma2SCN*(nuSCN/gamaSCN+1-nuSCN)-viesSCN^2

# Skew Power-exponential
thetaSPE <- ssmn.nl(y=Y, x=x, betas=beta0, sigma2=sigma20,shape = lambda0, nu=0.6,nlf=nlf,family = "Skew.pe")
thetaSPE=thetaSPE$theta
logSPE=spe.logL(thetaSPE,p,Y,x,nlf) 
betaSPE=thetaSPE[1:p]
sigma2SPE=thetaSPE[p+1]
lambdaSPE=thetaSPE[p+2]
nuSPE=thetaSPE[p+3]


################################################################################
# Standard error estimation

#SN
der1beta=matrix(0,p,n)
der2beta=list()
for (i in 1:n){
derbeta=der1(betaSN[1],betaSN[2],betaSN[3],x[i])
der1beta[,i]=t(attr(derbeta,"gradient"))
aux=attr(derbeta,"hessian")
der2beta[[i]]=matrix(aux,p,p)
}           
miSN=ssmn.nl.im3(Y,x,thetaST,p,der1beta,der2beta,family = "Skew.normal")
sqrt(diag(solve(miSN)))


#STN
der1beta=matrix(0,p,n)
der2beta=list()
for (i in 1:n){
derbeta=der1(betaST[1],betaST[2],betaST[3],x[i])
der1beta[,i]=t(attr(derbeta,"gradient"))
aux=attr(derbeta,"hessian")
der2beta[[i]]=matrix(aux,p,p)
}
miST=ssmn.nl.im3(Y,x,thetaST,p,der1beta,der2beta,family = "Skew.t")
sqrt(diag(solve(miST)))

ssmn.nl.score(Y,x,thetaST,p,der1beta,family = "Skew.t")  # OK

#SSL
der1beta=matrix(0,p,n)
der2beta=list()
for (i in 1:n){
derbeta=der1(betaSSL[1],betaSSL[2],betaSSL[3],x[i])
der1beta[,i]=t(attr(derbeta,"gradient"))
aux=attr(derbeta,"hessian")
der2beta[[i]]=matrix(aux,p,p)
}
miSSL=ssmn.nl.im3(Y,x,thetaSSL,p,der1beta,der2beta,family = "Skew.slash")
sqrt(diag(solve(miSSL)))

ssmn.nl.score(Y,x,thetaSSL,p,der1beta,family = "Skew.slash")  # OK

#SCN
der1beta=matrix(0,p,n)
der2beta=list()
for (i in 1:n){
derbeta=der1(betaSCN[1],betaSCN[2],betaSCN[3],x[i])
der1beta[,i]=t(attr(derbeta,"gradient"))
aux=attr(derbeta,"hessian")
der2beta[[i]]=matrix(aux,p,p)
}
miSCN=ssmn.nl.im3(Y,x,thetaSCN,p,der1beta,der2beta,family = "Skew.cn") 
sqrt(diag(solve(miSCN)))

ssmn.nl.score(Y,x,thetaSCN,p,der1beta,family = "Skew.cn")

#SPE
der1beta=matrix(0,p,n)
der2beta=list()
for (i in 1:n){
derbeta=der1(betaSPE[1],betaSPE[2],betaSPE[3],x[i])
der1beta[,i]=t(attr(derbeta,"gradient"))
aux=attr(derbeta,"hessian")
der2beta[[i]]=matrix(aux,p,p)
}
miSPE=ssmn.nl.im3(Y,x,thetaSPE,p,der1beta,der2beta,family = "Skew.pe") 
sqrt(diag(solve(miSPE)))

ssmn.nl.score(Y,x,thetaSPE,p,der1beta,family = "Skew.pe")  # OK



################################################################################
# Envelopes with Pearson's residuals

resSN=(Y-EYSN)/sqrt(varSN)
fit.model = lm(resSN~1)
source("envel_norm.txt") 
           
resST=(Y-EYST)/sqrt(varST)
fit.model = lm(resST~1)
source("envel_norm.txt")

resSCN=(Y-EYSCN)/sqrt(varSCN)
fit.model = lm(resSCN~1)
source("envel_norm.txt")

resSSL=(Y-EYSSL)/sqrt(varSSL)
fit.model = lm(resSSL~1)
source("envel_norm.txt")



# Envelope using quadratic forms
resSN=(Y-nlf(x,betaSN))/sqrt(sigma2SN)
envelSN(resSN,0,1)

resST=(Y-nlf(x,betaST))/sqrt(sigma2ST)
envelStN(Y,nlf(x,betaST),sigma2ST,nuST)
envelStN(resST,0,1,nuST)

resSSL=(Y-nlf(x,betaSSL))/sqrt(sigma2SSL)
envelSSL(resSSL,0,1,nuSSL)

resSCN=(Y-nlf(x,betaSCN))/sqrt(sigma2SCN)
envelSCN(resSCN,0,1,nuSCN,gamaSCN)

resSPE=(Y-nlf(x,betaSPE))/sqrt(sigma2SPE)
envelSPE(resSPE,0,1,nuSPE)

################################################################################
# Mahalanobis distance

resSN=(Y-EYSN)/sqrt(varSN)
resSN=(Y-nlf(x,betaSN))/sqrt(sigma2SN)
corte=qchisq(0.95, 1)
Ds=resSN^2
plot(Ds,ylab="Mahalanobis Distance", xlab="Index")
abline(corte,0,lty="dashed",lwd=2)
identify(Ds,n=2)

resST=(Y-EYST)/sqrt(varST)
resST=(Y-nlf(x,betaST))/sqrt(sigma2ST)
corte=qf(0.95, 1, nuST)
Ds=resST^2
plot(Ds,ylab="Mahalanobis Distance", xlab="Index",ylim=c(0,max(max(Ds),corte)))
abline(corte,0,lty="dashed",lwd=2)
identify(Ds,n=1)          
        
# SSL
fqssl=function(r,pq,nu) abs(pq-pchisq(r, 1)+2^nu*gamma(nu+0.5)/(r^nu*sqrt(pi))*pchisq(r,2*nu+1))
qssl=optimize(fqssl,c(0.001,50),0.95,nuSSL)$minimum
qssl    

resSSL=(Y-EYSSL)/sqrt(varSSL)  
resSSL=(Y-nlf(x,betaSSL))/sqrt(sigma2SSL)
corte=qssl
Ds=resSSL^2
plot(Ds,ylab="Mahalanobis Distance", xlab="Index")
abline(corte,0,lty="dashed",lwd=2)
identify(Ds,n=2)


fqsnc=function(r,pq,nu,gama) abs(nu*pchisq(gama*r, 1)+(1-nu)*pchisq(r,1)-pq)
qsnc=optimize(fqsnc,c(0.001,500),0.95,nuSCN,gamaSCN)$minimum
qsnc    


resSCN=(Y-EYSCN)/sqrt(varSCN)
resSCN=(Y-nlf(x,betaSCN))/sqrt(sigma2SCN)
corte=qsnc
Ds=resSCN^2
plot(Ds,ylab="Mahalanobis Distance", xlab="Index")
abline(corte,0,lty="dashed",lwd=2)
identify(Ds,n=2)

# SPE
fqspe=function(r,pq,nu) abs(sqrt(r)*pgamma(1/(2*nuSPE),r^nuSPE/2)/(2^(2*nuSPE))-pq)
qspe=optimize(fqspe,c(0.001,500),0.95,nuSPE)$minimum
qspe    

resSPE=(Y-nlf(x,betaSPE))/sqrt(sigma2SPE)
corte=qspe
Ds=resSPE^2
plot(Ds,ylab="Mahalanobis Distance", xlab="Index")
abline(corte,0,lty="dashed",lwd=2)
identify(Ds,n=4)

