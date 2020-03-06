source('functions_ssmn_nl.R')

##Define non linear function
nlf<-function(x,betas){
resp<- betas[1]/(1 +betas[2]*exp(-betas[3]*x))
return(resp)
}


# Derivates of nlf(x,beta) in beta
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


n=200
p=3
beta1=37
beta2=43
beta3=0.6
sigma2=0.5
lambda=-3  
nu=3
x = as.matrix(seq(0.1,15,length=n))
munl=beta1/(1+beta2*exp(-beta3*x))  
erro = as.matrix(rstn(n, 0, sigma2, lambda,nu))
Y=munl+erro


# Least squares
nls0 <- nls(Y~beta1/(1+beta2*exp(-beta3*x)), start=c(beta1=37, beta2=43,beta3=0.6))
#nls0 <- nls(Y~beta1/(1+beta2*exp(-beta3*x)), trace=TRUE)
t0=summary(nls0)
beta0=as.matrix(t0$coefficients[,1])
sigma0=as.numeric(t0$sigma)
sigma20=sigma0^2
resid0=t0$residuals #= Y - nlf(x,beta0)
lambda0=skewness(resid0)

# Fit SSMN
# SN
thetaSN <- ssmn.nl(y=Y, x=x, betas=beta0, sigma2=sigma20,shape = lambda0, nlf=nlf,family = "Skew.normal")
# STN
thetaSTN <- ssmn.nl(y=Y, x=x, betas=beta0, sigma2=sigma20,shape = lambda0, nu=20,nlf=nlf,family = "Skew.t")
# SSlash
thetaSSL <- ssmn.nl(y=Y, x=x, betas=beta0, sigma2=sigma20,shape = lambda0, nu=20,nlf=nlf,family = "Skew.slash")
# Skew Contaminated Normal
thetaSCN <- ssmn.nl(y=Y, x=x, betas=beta0, sigma2=sigma20,shape = lambda0, nu=nu0,nlf=nlf,family = "Skew.cn")
# SPE
thetaSPE <- ssmn.nl(y=Y, x=x, betas=beta0, sigma2=sigma20,shape = lambda0, nu=1,nlf=nlf,family = "Skew.pe")

