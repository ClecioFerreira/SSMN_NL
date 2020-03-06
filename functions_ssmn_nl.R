
library(moments)
library(truncdist)

# Random samples

rsnn <- function(n, location=0, scale=1, shape=0, dp=NULL) 
{
  if(!is.null(dp)) {
     if(!missing(shape)) 
        stop("You cannot set both component parameters and dp")
     location <- dp[1]
     scale <- dp[2]
     shape <- dp[3]
    }
  if(scale<0) {
     stop("Parameter scale must be positive")
   } 
  #k(u)=1: não precisa
  #
  delta<-shape/sqrt(1+shape^2)
  u1<-rnorm(n)
  u2<-rnorm(n)
  y<-location+sqrt(scale)*(delta*abs(u1)+sqrt(1-delta^2)*u2)
  return(y)
}

rstn <- function(n, location=0, scale=1, shape=0, nu=30, dp=NULL) 
{
  if(!is.null(dp)) {
     if(!missing(shape)) 
        stop("You cannot set both component parameters and dp")
     location <- dp[1]
     scale <- dp[2]
     shape <- dp[3]
     nu <- dp[4]
    }
   if (nu<0) {
    stop("Parameter nu must be positive") 
   }
   if(scale<0) {
    stop("Parameter scale must be positive")
   }
  if (nu<1) warning('Nu < 1 can generate values tending to infinite',call. = FALSE) 
  # k(u)=1/u, U ~ G(nu/2,gama/2)
  u<-rgamma(n,nu/2,nu/2)
  ku<-1/u
  shape1<-shape*sqrt(ku)
  scale1<-scale*ku
  delta<-shape1/sqrt(1+shape1^2)
  u1 <- rnorm(n)
  u2 <- rnorm(n)
  y <- location+sqrt(scale1)*(delta*abs(u1)+sqrt(1-delta^2)*u2)
  return(y)
}

rssl <- function(n, location=0, scale=1, shape=0, nu=30, dp=NULL) 
{
 if(!is.null(dp)) {
     if(!missing(shape)) 
        stop("You cannot set both component parameters and dp")
     location <- dp[1]
     scale <- dp[2]
     shape <- dp[3]                                                                 
     nu <- dp[4]
    }
   if (nu<0) {
    stop("Parameter nu must be positive") 
   }
 if(scale<0) {
    stop("Parameter scale must be positive")
  }
 v<-runif(n,0,1)
 u<-v^(1/nu)
 ku<-1/u
 shape1<-shape*sqrt(ku)
 scale1<-scale*ku
 delta<-shape1/sqrt(1+shape1^2)
 u1 <- rnorm(n)
 u2 <- rnorm(n)
 y <- location+sqrt(scale1)*(delta*abs(u1)+sqrt(1-delta^2)*u2)
 return(y)
}

rscn <- function(n, location=0, scale=1, shape=9, nu=1, gama=1, dp=NULL)
{
  if(!is.null(dp)) {
     if(!missing(shape)) 
        stop("You cannot set both component parameters and dp")
     location <- dp[1]
     scale <- dp[2]
     shape <- dp[3]
     nu <- dp[4]
     gama <- dp[5]
    }
   if (nu<0|nu>1) {
    stop("Parameter nu must be between 0.0 and 1.0") 
   }
   if(gama<0|gama>1) {
    stop("Parameter gama must be between 0.0 and 1.0")
   }
   if(scale<0) {
    stop("Parameter scale must be positive")
   }
  # k(u)=1/u, U ~ nu*I(u=gamanc)+(1-nu)*I(u=1)
  uu<-rbinom(n,1,nu)
  u<-gama*uu+1-uu
  ku<-1/u
  shape1<-shape*sqrt(ku)
  scale1<-scale*ku
  delta<-shape1/sqrt(1+shape1^2)
  u1<-rnorm(n)
  u2<-rnorm(n)
  y<-location+sqrt(scale1)*(delta*abs(u1)+sqrt(1-delta^2)*u2)
  return(y)
}

rspe <- function (n, location=0, scale=1, shape=0, nu=1, dp=NULL) 
{
   if(!is.null(dp)) {
    if(!missing(shape)) 
        stop("You cannot set both component parameters and dp")
    location <- dp[1]
    scale <- dp[2]
    shape <- dp[3]
    nu<- dp[4]
    }
   if (nu<=0.5|nu>1) {                                                        
    stop("Parameter nu must be 0.5 < nu <= 1.0") 
   }
   if(scale<=0) {
    stop("Parameter scale must be positive")
   }
 t1 <- rgamma(n,nu/2,1/2,2)
 R <- t1^(1/(2*nu))
 u <- sign(rnorm(n))
 yep <- u*R
 w <- rnorm(n)
 y <- location+sqrt(scale)*yep*sign(shape*yep-w)
 return(y)
}

# Probability density functions

dsnn <- function (x, location=0, scale=1, shape=0, dp=NULL)
{
  if(!is.null(dp)) {
     if(!missing(shape)) 
        stop("You cannot set both component parameters and dp")
     location <- dp[1]
     scale <- dp[2]
     shape <- dp[3]
    }
  if(scale<0) {
     stop("Parameter scale must be positive")
   }
   if (any(is.na(x))) x<- x[-which(is.na(x))] 
  z<-(x-location)/sqrt(scale)
  y<-2*dnorm(z)*pnorm(shape*z)/sqrt(scale)
  return(y)
}

dstn <- function(x, location=0, scale=1, shape=0, nu=30, dp=NULL)
{
   if(!is.null(dp)) {
     if(!missing(shape)) 
        stop("You cannot set both component parameters and dp")
     location <- dp[1]
     scale <- dp[2]
     shape <- dp[3]
     nu <- dp[4]
    }
   if (nu<0) {
    stop("Parameter nu must be positive") 
   }
   if(scale<0) {
    stop("Parameter scale must be positive")
   }
  if (any(is.na(x))) x<- x[-which(is.na(x))]
  z<-(x-location)/sqrt(scale)
  cte<-gamma((nu+1)/2)/gamma(nu/2)/sqrt(nu*pi)
  yt<-cte*(1+z^2/nu)^(-(nu+1)/2)
  y<-2*yt*pnorm(shape*z)/sqrt(scale)
  return(y)
}

dssl <- function(x, location=0, scale=1, shape=0, nu=30, dp=NULL)
{
 if(!is.null(dp)) {
     if(!missing(shape)) 
        stop("You cannot set both component parameters and dp")
     location <- dp[1]
     scale <- dp[2]
     shape <- dp[3]
     nu <- dp[4]
    }
   if (nu<0) {
    stop("Parameter nu must be positive") 
   }
   if(scale<0) {
    stop("Parameter scale must be positive")
   }
 if (any(is.na(x))) x<- x[-which(is.na(x))] 
 z<-(x-location)/sqrt(scale)
 n<-length(x)
 cte=nu/((2*pi*scale)^(1/2))
 d=z^2
 fdps=cte*gamma(nu+1/2)*pgamma(1,nu+1/2,scale=2/d)/((d/2)^(nu+1/2))
 fsl=2*fdps*pnorm(shape*z)
 return(fsl)
}

dscn <- function(x, location=0, scale=1, shape=0, nu=1,gama=1,dp=NULL)
{ 
    if(!is.null(dp)) {
     if(!missing(shape)) 
        stop("You cannot set both component parameters and dp")
     location <- dp[1]
     scale <- dp[2]
     shape <- dp[3]
     nu <- dp[4]
     gama <- dp[5]
    }
   if (nu<0|nu>1) {
    stop("Parameter nu must be between 0.5 and 1.0") 
   }
   if(gama<0|gama>1) {
    stop("Parameter gama must be between 0.5 and 1.0")
   }
   if(scale<0) {
    stop("Parameter scale must be positive")
   }
   if (any(is.na(x))) x<- x[-which(is.na(x))]
  z<-(x-location)/sqrt(scale)
  z2<-z*sqrt(gama)
  y<-2*(nu*dnorm(x,mean=location,sd=sqrt(scale/gama))+(1-nu)*dnorm(x,mean=location,sd=sqrt(scale)))*pnorm(shape*z)
  return(y)
}

dspe <- function(x, location=0, scale=1, shape=0, nu=1, dp=NULL)
 {
   if(!is.null(dp)) {
     if(!missing(shape)) 
        stop("You cannot set both component parameters and dp")
     location <- dp[1]
     scale <- dp[2]
     shape <- dp[3]
     nu<- dp[4]
    }
   if (nu<0.5|nu>1) {
    stop("Parameter nu must be between 0.5 and 1.0") 
   }
   if(scale<0) {
    stop("Parameter scale must be positive")
   }
    if (any(is.na(x))) x<- x[-which(is.na(x))] 
   nu1<-1/(2*nu)
   z<-(x-location)/sqrt(scale)
   d<-z^2
   cte<-nu/(2^nu1)/gamma(nu1)/sqrt(scale)
   yt<-cte*exp(-0.5*d^(nu))
   y<-2*yt*pnorm(shape*z)
   return(y)
 }
 
# Log-lilkelihood of nu (ECME-step)
 
ft<-function(nu,d2) {
  n=length(d2)
  aux=1+d2/nu
  dist=n*log(gamma((nu+1)/2))-n/2*log(nu)-n*log(gamma(nu/2))-(nu+1)/2*sum(log(aux))
  return(-dist)
}

fsl<-function(nu,sigma2,d2){
  cte=nu/(sqrt(sigma2)*sqrt(2*pi))
  v1=gamma(nu+0.5)*pgamma(1,nu+0.5,scale=2/d2)
  v1=pmax(v1,1e-300)
  f1=cte*v1/((d2/2)^(nu+0.5))
  f1=pmax(f1,1e-300)
  g=-sum(log(f1))
  #print(g)
  return(g)
}                           

fcn<-function(nugama,y,mu,sigma2,lambda)
{
  nu=nugama[1]
  gama=nugama[2]
  fy=nu*dnorm(y,mu,sqrt(sigma2/gama))+(1-nu)*dnorm(y,mu,sqrt(sigma2))
  f=-sum(log(fy))
  return(f)
}

fpe<-function(nu,d,sigma2){
  cnu<-nu/(2^(1/(2*nu))*gamma(1/(2*nu)))
  fy<-cnu*sqrt(sigma2)*exp(-0.5*d^nu)
  g=-sum(log(fy))
  return(g)
}                           

fpenu<-function(nu,d){
  cnu<-nu/(2^(1/(2*nu))*gamma(1/(2*nu)))
  fy<-cnu*exp(-0.5*d^nu)
  g=-sum(log(fy))
  return(g)
}                           
 

# Log-likelihood

n.logL <-function(theta,p,y,x,nlf) 
{
  betap<-signif(theta[1:p],digits=7)
  mu=nlf(x,betap)
  sigma2<-signif(theta[1+p],digits=7)
  sigma<-signif(sqrt(sigma2),digits=7)
  vero<-signif(dnorm(y,mu,sigma),digits=7)
  logvero<-signif(log(vero),digits=7)
  return(sum(logvero))
}



snn.logL <-function(theta,p,y,x,nlf) 
{
  beta<-signif(theta[1:p],digits=7)
  mu=nlf(x,beta)
  sigma2<-signif(theta[1+p],digits=7)
  sigma<-signif(sqrt(sigma2),digits=7)
  lambda<-signif(theta[p+2],digits=7)
  vero<-signif(2*dnorm(y,mu,sigma)*pnorm(lambda*(y-mu)/sigma),digits=7)
  logvero<-signif(log(vero),digits=7)
  return(sum(logvero))
}

stn.logL <- function(theta,p,y,x,nlf)
{
 #% Theta=[beta,sigma2,lambda,nu] %#
 beta<-signif(theta[1:p],digits=7)
 mu=nlf(x,beta)
 sigma2<-signif(theta[p+1],digits=7)
 lambda<-signif(theta[p+2],digits=7) 
 nu<-signif(theta[p+3],digits=7)
 d<-signif((y-mu)^2/sigma2,digits=7)
 aux<-signif(lambda*(y-mu)/sqrt(sigma2),digits=7)
 aux1<-signif(pmax(aux,-37),digits=7)
 cnugama<-signif(2*gamma((nu+1)/2)/(gamma(nu/2)*sqrt(nu*pi)),digits=7)
 fy<-signif(cnugama/sqrt(sigma2)*(1+d/nu)^(-(nu+1)/2)*pnorm(aux1),digits=7)
 return(sum(log(fy)))
}


ssl.logL <- function (theta,p,y,x,nlf) 
{
 #% =Theta=[beta,sigma^2,lambda,nu] #
 n<-length(y)
 beta<-theta[1:p]
 mu=nlf(x,beta)
 sigma2<-theta[p+1]
 lambda<-theta[p+2]
 nu<-theta[p+3]
 z<-(y-mu)/sqrt(sigma2)
 cte=nu/((2*pi*sigma2)^(1/2))
 d=z^2
 fdps=cte*gamma(nu+1/2)*pgamma(1,nu+1/2,scale=2/d)/((d/2)^(nu+1/2))
 fsl=2*fdps*pnorm(lambda*z)
 verot<-signif(sum(log(fsl)),digits=7)
 return(verot)
}

scn.logL <- function (theta,p,y,x,nlf) 
{
   #% Theta=[beta,sigma2,lambda,nu,gama]#
   beta<-theta[1:p]
   mu=nlf(x,beta)
   sigma2<-theta[p+1]
   lambda<-theta[p+2]
   nu<-theta[p+3]
   gama<-theta[p+4]
   d<-(y-mu)^2/sigma2
   aux<-lambda*(y-mu)/sqrt(sigma2)
   aux1<-signif(pmax(aux,-37),digits=7)
   fy<-signif(2*(nu*dnorm(y,mu,sqrt(sigma2/gama))+(1-nu)*dnorm(y,mu,sqrt(sigma2)))*pnorm(aux1),digits=7)
   return(sum(signif(log(fy),digits=7)))
}

spe.logL<- function(theta,p,y,x,nlf) 
{
  #% =Theta=[beta,sigma^2,lambda,nu] #
 n<-length(y)
 beta<-theta[1:p]
 mu=nlf(x,beta)
 sigma2<-theta[p+1]
 lambda<-theta[p+2]
 nu<-theta[p+3]
 d<-(y-mu)^2/sigma2 
 aux<-lambda*(y-mu)/sqrt(sigma2)
 aux1<-pmax(aux,-37)
 cnu<-2*nu/(2^(1/(2*nu))*sqrt(sigma2)*gamma(1/(2*nu)))
 fy<-cnu*exp(-0.5*d^nu)*pnorm(aux1)
 logfy<-signif(log(fy),digits=7)
 return(sum(logfy))
}

# Maximum likelihood estimation using EM algorithm

normal.nl <- function(y, x, nlf, betap0)
{
n=length(y)
betap<-nlm(Qbeta.n,betap0,y,x,nlf)$estimate
betap=as.matrix(betap)
res = y - nlf(x,betap)  
sigma2=as.numeric(t(res%*%res)/n)
theta=rbind(betap,sigma2) 
return(theta)
}

ssmn.nl <- function(y, x = NULL, betas = NULL, sigma2 = NULL,
    shape = NULL, nu = NULL, nlf,family = "Skew.normal",
    error = 1e-05, iter.max = 50000,lower=NULL,upper=NULL)
{
    if (ncol(as.matrix(y)) > 1)
        stop("Only univariate non linear regression supported!")
    if (length(y) != nrow(as.matrix(x)))
        stop("X variable does not have the same number of lines than y")
    if ((length(x) == 0) | (length(betas) == 0) | (length(sigma2) ==
        0) | (length(shape) == 0))
        stop("All parameters must be provided.")
    if (!is.function(nlf))
        stop("nfl parameter must be a function!")
    if ((family != "Skew.t") && (family !=
        "Skew.cn") && (family != "Skew.normal") && (family != "Skew.slash")&& (family != "Skew.pe"))
        stop("Distribution family not supported. Check documentation!")
    if (family == "Skew.normal") {
        if (length(nu) != 0)
            stop("For the Skew-Normal distribution, the Nu parameters must not be provided.")
    }
    if ((family == "Skew.t") | (family == "Skew.slash") | (family == "Skew.pe")) {
        if (length(nu) > 1)
            stop("Nu parameters must have only one parameter")
        if (length(nu) == 0)
            stop("Nu parameters must be provided.")
        if (nu <= 0)
            stop("Nu parameters must be positive.")
    }
    if (family == "Skew.cn") {
        if (length(nu) != 2)
            stop("For the Skew.cn nu must have 2 parameters")
        if (nu[1] <= 0 || nu[1] >= 1)
            stop("nu[1] must be in (0,1)")
        if (nu[2] <= 0 || nu[2] >= 1)
            stop("nu[2] must be in (0,1)")
    }
    out <- ssmn.nl.intern(y, x, betas, sigma2, shape, nu, nlf,
            family, error, iter.max,lower,upper)
    out
}

ssmn.nl.intern <- function(y, x, beta, sigma2, lambda, nu,
            nlf,family, error, iter.max,lower,upper)
{
n=length(y)
criterio=1
cont=0
theta0=rbind(as.matrix(beta),sigma2,lambda,nu)
L1<-rbind(1e-3,1e-3)
L2<-rbind(0.999,0.999)
while ((criterio > error)&&(cont<iter.max)) {
    cont=cont+1
    mu=nlf(x,beta)
    res=y-mu
    sigma=sqrt(sigma2)
    eta=lambda*res
    aux=eta/sigma
    aux1=pmax(aux,-37)
    Wphi=dnorm(aux1)/pnorm(aux1)
    t1=eta+sigma*Wphi
    t2=eta^2+sigma2+sigma*eta*Wphi
    d=res^2/sigma2
    if (family == "Skew.normal") {
       ki=as.vector(rep(1,n))
    }
    if (family == "Skew.t") {
       nu=optimize(ft,c(2,50),d)$minimum     # valores de nu>2 para estabilidade e cálculos de E(Y) e Var(Y) e tbm da MI
       ki=as.vector((nu+1)/(nu+d))
    }
    if (family == "Skew.slash") {
       nu<-optimize(fsl,c(2,30),sigma2,d)$minimum  # valores de nu>1 para estabilidade e cálculos de E(Y) e Var(Y)  e tbm da MI
       ki=as.vector((2*nu+1)/d*(pgamma(1,nu+1.5,scale=2/d)/pgamma(1,nu+0.5,scale=2/d)))
    }
    if (family == "Skew.pe") {
       nu<-optimize(fpenu,c(0.51,1),d)$minimum 
       ki=as.vector(nu*d^(nu-1))
    }
    if (family == "Skew.cn") {
       nugama0<-as.matrix(signif(theta0[c(length(theta0)-1,length(theta0))],digits=7))
       nu<-optim(nugama0,fcn,method='L-BFGS-B',lower=L1,upper=L2,y=y,mu=mu,sigma2=sigma2,lambda=lambda)$par 
       nu1<-nu[1]
       gama<-nu[2]
       aux1=pmin(1-nu1+nu1*gama^(1.5)*exp((1-gama)*d/2),1e+308)
       aux2=pmin(1-nu1+nu1*gama^(0.5)*exp((1-gama)*d/2),1e+308)
       ki<-as.vector(aux1/aux2)
     }
     if (is.null(lower)==TRUE) {
    beta<-optim(beta,Qbeta,gr = NULL,y,x,lambda,sigma2,t1,ki,nlf,method='BFGS')$par}
    else beta<-optim(beta,Qbeta,gr = NULL,y,x,lambda,sigma2,t1,ki,nlf,method='L-BFGS-B',lower=lower,upper=upper)$par
    beta=as.matrix(beta)
    sigma2=as.numeric((t(res)%*%diag(ki+lambda^2)%*%res-2*lambda*t(t1)%*%res+sum(t2))/(2*n))
    lambda=as.numeric(t(t1)%*%res/(t(res)%*%res))
    theta=rbind(beta,sigma2,lambda,nu)
    dif=theta-theta0
    criterio=sqrt(sum(dif^2))
    theta0=theta
 }
 return(list(theta=theta,cont=cont))
}

Qbeta <- function(beta,y,x,lambda,t1,ki,miv,nlf){
mu=nlf(x,beta)
res=y-mu-lambda*t1/(lambda^2+ki)
Qb=sum((lambda^2+ki)/miv*(res^2))
return(Qb)
}

# Information matrix

normal.nl.im<- function(y,x,theta,p,nlf) 
{
 # THETA: (beta,sigma2)
 n=length(y)
 beta<-theta[1:p]
 sigma2<-as.numeric(theta[p+1])
 mu=nlf(x,beta)
 res<-as.matrix(y-mu)
 derbeta12=der_eta_beta(x,beta)
 derbeta1=derbeta12$derbeta
 derbeta2=derbeta12$der2beta
 ###################### Parte Simétrica: (beta,sigma2)
    I1<-matrix(0,p+1,p+1)
    I1betaaux=matrix(0,p,p)
    I1betasigmaaux=matrix(0,p,1)
    for (i in 1:n){
        der1beta=as.matrix(derbeta1[,i])
        der2beta=derbeta2[[i]]
        I1betaaux=I1betaaux+res[i]*der2beta-der1beta%*%t(der1beta)
        I1betasigmaaux=I1betasigmaaux+res[i]*der1beta
    }
    I1[1:p,1:p]<-1/(sigma2)*I1betaaux
    I1[1:p,p+1]<--1/(sigma2^2)*I1betasigmaaux
    I1[p+1,p+1]<- n/(2*(sigma2^2))-1/(sigma2^3)*t(res)%*%res
    I1[p+1,1:p]<-t(I1[1:p,p+1])  
 ##########################################################
 IM=-I1
 return(IM)
}


ssmn.nl.im3<- function(y,x,theta,p,derbeta1,derbeta2,family = "Skew.normal") 
{
 # THETA: (beta,sigma2,lambda,nu,gama) ; nu para StN e SSL; nu e gama para SCN
 # usando derivadas em beta do R
 n=length(y)
 beta<-theta[1:p]
 sigma2<-as.numeric(theta[p+1])
 lambda<-as.numeric(theta[p+2])
 mu=nlf(x,beta)
 res<-as.matrix(y-mu)
 sigma<-sqrt(sigma2)
 aux<-lambda*res/sigma
 aux1<-pmax(aux,-37)
 Wphi<-dnorm(aux1)/pnorm(aux1)
 d<-res^2/sigma2
 ## I2(beta,sigma2,lambda): Assimétrica   (-der2 l(theta)/dtheta dthetaT)
 Wphi1<- -Wphi*(aux1+Wphi)
 dl2sigma2=-lambda/(2*(sigma^3))*res
 dl2lambda=1/sigma*res
 I2S<-matrix(0,p+2,p+2)
 for (i in 1:n){
    der2i=matrix(0,p+2,p+2)
    der1beta=as.matrix(derbeta1[,i])
    dl2beta=-lambda/sigma*der1beta
    der2beta=derbeta2[[i]]
    der2i[1:p,1:p]=-lambda/sigma*der2beta
    I2betasigma=lambda/(2*(sigma^3))*der1beta
    der2i[(p+1),1:p]=t(I2betasigma)
    der2i[1:p,(p+1)]=I2betasigma 
    I2betalambda=-1/sigma*der1beta
    der2i[(p+2),1:p]=t(I2betalambda)
    der2i[1:p,(p+2)]=I2betalambda
    der2i[(p+1),(p+1)]=3*lambda/(4*(sigma^5))*res[i]
    I2sigmalambda=-1/(2*(sigma^3))*res[i]
    der2i[(p+1),(p+2)]= I2sigmalambda
    der2i[(p+2),(p+1)]= I2sigmalambda
    #I2lambdalambda=0
    dl2i=c(dl2beta,dl2sigma2[i],dl2lambda[i])
    I2i=Wphi[i]*der2i+Wphi1[i]*dl2i%*%t(dl2i)
    I2S=I2S+I2i 
 }
 ## Parte Simétrica: (beta,sigma2, , nu, gama) (- der2 L1(theta)/dtheta dthetaT)
 if (family=="Skew.normal"){
    I2=I2S
    I1<-matrix(0,p+2,p+2)
    I1betaaux=matrix(0,p,p)
    I1betasigmaaux=matrix(0,p,1)
    for (i in 1:n){
        der1beta=as.matrix(derbeta1[,i])
        der2beta=derbeta2[[i]]
        I1betaaux=I1betaaux+res[i]*der2beta-der1beta%*%t(der1beta)
        I1betasigmaaux=I1betasigmaaux+res[i]*der1beta
    }
    I1[1:p,1:p]<-1/(sigma2)*I1betaaux
    I1[1:p,p+1]<--1/(sigma2^2)*I1betasigmaaux
    I1[p+1,p+1]<- n/(2*(sigma2^2))-1/(sigma2^3)*t(res)%*%res
    I1[p+1,1:p]<-t(I1[1:p,p+1])
 }
 if (family=="Skew.t"){
    nu<-as.numeric(theta[p+3])
    # Parte Assimétrica
    I2t=matrix(0,p+3,p+3)
    I2t[1:(p+2),1:(p+2)]=I2S
    I2=I2t
    # Parte Simétrica I1
    V<-(1+d/nu)^(-1)
    hnu<-psigamma((nu+1)/2,1)-psigamma(nu/2,1)+2/(nu^2)
    Qvbeta<- t(res)%*%diag(as.vector(V))%*%res
    I1<-matrix(0,p+3,p+3)
    for (i in 1:n){
        der1beta=as.matrix(derbeta1[,i])
        der2beta=derbeta2[[i]]
        der1i=matrix(0,p+3,p+3)
        der1i[1:p,1:p]=(nu+1)/(nu*sigma2)*V[i]*(res[i]*der2beta+(2/nu*V[i]*d[i]-1)*der1beta%*%t(der1beta))
        I1betasigma=(nu+1)/(nu*(sigma2^2))*res[i]*V[i]*(V[i]*d[i]/nu-1)*der1beta
        der1i[(p+1),1:p]=t(I1betasigma)
        der1i[1:p,(p+1)]=I1betasigma
        I1betanu=1/(nu^2*sigma2)*res[i]*V[i]*((nu+1)/nu*V[i]*d[i]-1)*der1beta
        der1i[(p+3),1:p]=t(I1betanu)
        der1i[1:p,(p+3)]=I1betanu
        I1=I1+der1i  
    }
     #I1[p+1,p+1]<- n/(2*(sigma2^2))-(nu+1)/(nu*(sigma2^3))*Qvbeta+(nu+1)/(2*(nu^2)*(sigma2^4))*t(res)%*%(diag(as.vector(res^3)))%*%diag(as.vector(V))%*%V
     I1[p+1,p+1]<- n/(2*(sigma2^2))-(nu+1)/(nu*(sigma2^2))*sum(V*d-(V*d)^2/(2*nu)) #OK
     #I1[p+1,p+3]<- 0.5/((nu*sigma2)^2)*t(res)%*%diag(as.vector(V))%*%((nu+1)/nu*diag(as.vector(V))%*%diag(as.vector(d))-diag(n))%*%res
     I1[p+1,p+3]<- 0.5/(nu^2*sigma2)*sum((nu+1)/nu*(V*d)^2-V*d) # OK
     #I1[p+3,p+3]<- n/4*hnu+0.5/(nu^2*sigma2)*Qvbeta+0.5/((nu^4)*sigma2)*t(res)%*%diag(as.vector(V))%*%((nu+1)*diag(as.vector(V))%*%diag(as.vector(d))-nu*(nu+2)*diag(n))%*%res
     I1[p+3,p+3]<- n/4*hnu+0.5/(nu^2*sigma2)*Qvbeta+0.5/(nu^4)*sum((nu+1)*(V*d)^2-nu*(nu+2)*V*d) # OK
     I1[p+3,1:p]<-I1[1:p,p+3]
     I1[p+1,1:p]<-I1[1:p,p+1]
     I1[p+3,p+1]<-I1[p+1,p+3]    
 }
 if (family=="Skew.slash"){
    nu<-as.numeric(theta[p+3])
    # Parte Assimétrica
    I2t=matrix(0,p+3,p+3)
    I2t[1:(p+2),1:(p+2)]=I2S
    I2=I2t
    # Parte Simétrica I1
    Qbeta<-t(res)%*%res
    Qwbeta<-t(res)%*%diag(as.vector(Wphi1))%*%res
    I1<-matrix(0,p+3,p+3)
    P11<-pgamma(1,nu+0.5,scale=2/d)
    P13<-pgamma(1,nu+1.5,scale=2/d)
    P15<-pgamma(1,nu+2.5,scale=2/d)
    IG1<-gamma(nu+0.5)*P11/((d/2)^(nu+0.5))
    IG3<-gamma(nu+1.5)*P13/((d/2)^(nu+1.5))
    IG5<-gamma(nu+2.5)*P15/((d/2)^(nu+2.5))
    rp<-5000
    u1<-matrix(0,n,rp)
    u2<-u1
    for (j in 1:n) {
       u1[j,]<-rtrunc(rp,"gamma",a=0,b=1,shape=nu+1/2,rate=d[j]/2)
       u2[j,]<-rtrunc(rp,"gamma",a=0,b=1,shape=nu+3/2,rate=d[j]/2)
    }
    #E11=gamma(nu+0.5)*((d/2)^(-(nu+0.5)))*P11*rowMeans(log(u1)) 
    #E13=gamma(nu+1.5)*((d/2)^(-(nu+1.5)))*P13*rowMeans(log(u2)) 
    #E21=gamma(nu+2.5)*((d/2)^(-(nu+2.5)))*P11*rowMeans(log(u1)^2) 
    E11=P11*rowMeans(log(u1)) 
    E13=P13*rowMeans(log(u2)) 
    E21=P11*rowMeans(log(u1)^2)
    I1betabeta=matrix(0,p,p)
    for (i in 1:n){
        der1beta=as.matrix(derbeta1[,i])
        Iaux=1/IG1[i]*(res[i]*IG3[i]*derbeta2[[i]]+(d[i]/IG1[i]*(IG5[i]*IG1[i]-IG3[i]^2)-IG3[i])*der1beta%*%t(der1beta))
        I1betabeta=I1betabeta+Iaux
    }
    I1[1:p,1:p]<- I1betabeta/sigma2 
    I1[1:p,p+1]<- 1/(2*sigma2^2)*derbeta1%*%diag(as.vector(res/IG1))%*%(-2*IG3+d*(IG5-IG3^2/IG1))
    I1[1:p,p+3]<- 1/(sigma2)*derbeta1%*%diag(as.vector(res/(IG1^2)))%*%(IG1*E13-IG3*E11) 
    I1[p+1,p+1]<- 1/(2*sigma2^2)*t(n+(d/IG1))%*%(-2*IG3+0.5*d*(IG5-IG3^2/IG1))
    I1[p+1,p+3]<- 1/(2*sigma2)*t(d/(IG1^2))%*%(IG1*E13-IG3*E11)
    I1[p+3,p+3]<- -n/(nu^2)+t(IG1*E21-E11^2)%*%(IG1^(-2))
    I1[p+1,1:p]<- t(I1[1:p,p+1])
    I1[p+3,1:p]<- t(I1[1:p,p+3])
    I1[p+3,p+1]<- I1[p+1,p+3]
 }
 if (family=="Skew.cn"){
     nu<-as.numeric(theta[p+3])
     gama<-as.numeric(theta[p+4])
     # Parte Assimétrica I2
     I2t=matrix(0,p+4,p+4)
     I2t[1:(p+2),1:(p+2)]=I2S
     I2=I2t
     # Parte Simétrica I1
     I1<-matrix(0,p+4,p+4)
     fi<-dnorm(y,mu,sigma)
     figama<-dnorm(y,mu,sqrt(sigma2/gama))
     fs<-nu*figama+(1-nu)*fi
     dfssigma<-1/(2*sigma2)*(nu*figama*(gama*d-1)+(1-nu)*fi*(d-1))
     dfsnu<-figama-fi
     dfsgama<-nu/(2*gama)*figama*(1-gama*d)
 #I2 = soma(1/fs[i]*der2(fs[i](theta)/der(teta_j,teta_k))-1/(fs[i]^2)*der(fs[i](theta)/der(theta_j))*der(fs[i](theta)/der(theta_j))^T)
 d2lbetabeta<-matrix(0,p,p)
 for (i in 1:n) {
    dibeta=as.matrix(derbeta1[,i])
    auxn<-nu*gama*figama[i]+(1-nu)*fi[i]
    derfsbeta=1/sigma2*auxn*res[i]*dibeta
    der2fsbb<-1/sigma2*(auxn*res[i]*derbeta2[[i]]+(-auxn+d[i]*(nu*(gama^2)*figama[i]+(1-nu)*fi[i]))*dibeta%*%t(dibeta))
    aux=1/fs[i]*der2fsbb-1/(fs[i]^2)*derfsbeta%*%t(derfsbeta)
    d2lbetabeta<-d2lbetabeta+aux
    }
 d2lsigmabeta<-matrix(0,p,1)
 for (i in 1:n) {
    dibeta=as.matrix(derbeta1[,i])
    auxn<-nu*gama*figama[i]+(1-nu)*fi[i]
    derfsbeta=1/sigma2*auxn*res[i]*dibeta
    der2fssb<-1/(2*(sigma2^2))*res[i]*(nu*gama*figama[i]*(gama*d[i]-3)+(1-nu)*fi[i]*(d[i]-3))*dibeta
    aux=1/fs[i]*der2fssb-1/(fs[i]^2)*dfssigma[i]*derfsbeta
    d2lsigmabeta<-d2lsigmabeta+aux
    }
 d2lnubeta<-matrix(0,p,1)
 for (i in 1:n) {
    dibeta=as.matrix(derbeta1[,i])
    auxn<-nu*gama*figama[i]+(1-nu)*fi[i]
    derfsbeta=1/sigma2*auxn*res[i]*dibeta
    der2fsnub<-1/sigma2*res[i]*(gama*figama[i]-fi[i])*dibeta
    aux=1/fs[i]*der2fsnub-1/(fs[i]^2)*dfsnu[i]*derfsbeta
    d2lnubeta<-d2lnubeta+aux
    }
 d2lgamabeta<-matrix(0,p,1)
 for (i in 1:n) {
    dibeta=as.matrix(derbeta1[,i])
    auxn<-nu*gama*figama[i]+(1-nu)*fi[i]
    derfsbeta=1/sigma2*auxn*res[i]*dibeta
    der2fsgamab<-nu/(2*sigma2)*res[i]*figama[i]*(3-gama*d[i])*dibeta
    aux=1/fs[i]*der2fsgamab-1/(fs[i]^2)*dfsgama[i]*derfsbeta
    d2lgamabeta<-d2lgamabeta+aux;
    }
 d2fssigmasigma=1/(2*sigma2^2)*(-2*nu*gama*figama*d-2*(1-nu)*fi*d+fs+nu/2*figama*(gama*d-1)^2+(1-nu)/2*fi*(d-1)^2)
#### 
 Isigmasigma<-1/(fs^2)*(d2fssigmasigma*fs-dfssigma^2)
 d2fsnusigma<-1/(2*sigma2)*(d*(gama*figama-fi)-figama+fi)
 Inusigma<-1/(fs^2)*(d2fsnusigma*fs-dfssigma*dfsnu)
 d2fsgamasigma<-nu/(4*sigma2)*figama*(4*d-gama*(d^2)-1/gama)
 Igamasigma<-1/(fs^2)*(d2fsgamasigma*fs-dfssigma*dfsgama)
 d2fsnunu<-0
 Inunu<- -(dfsnu/fs)^2
 d2fsnugama<-1/2*figama*(1/gama-d)
 Inugama<-1/(fs^2)*(d2fsnugama*fs-dfsnu*dfsgama)
 d2fsgamagama<-nu/4*figama*((1/gama-d)^2-2/(gama^2))
 Igamagama<-1/(fs^2)*(d2fsgamagama*fs-dfsgama^2)
 ##########################################################################
 I1[1:p,1:p]<-d2lbetabeta
 I1[1:p,p+1]<-d2lsigmabeta
 I1[1:p,p+3]<-d2lnubeta
 I1[1:p,p+4]<-d2lgamabeta
 I1[p+1,p+1]<-sum(Isigmasigma)
 I1[p+1,p+3]<-sum(Inusigma)
 I1[p+1,p+4]<-sum(Igamasigma)
 I1[p+3,p+3]<-sum(Inunu)
 I1[p+3,p+4]<-sum(Inugama)
 I1[p+4,p+4]<-sum(Igamagama)
 I1[p+1,1:p]<- t(I1[1:p,p+1])
 I1[p+3,1:p]<- t(I1[1:p,p+3])
 I1[p+4,1:p]<- t(I1[1:p,p+4])
 I1[p+3,p+1]<- I1[p+1,p+3]
 I1[p+4,p+1]<- I1[p+1,p+4]
 I1[p+4,p+3]<- I1[p+3,p+4] 
 }
 if (family=="Skew.pe"){
    nu<-as.numeric(theta[p+3])
    # Parte Assimétrica
    I2t=matrix(0,p+3,p+3)
    I2t[1:(p+2),1:(p+2)]=I2S
    I2=I2t
    # Parte Simétrica I1
    I1<-matrix(0,p+3,p+3)
    dnu<- d^(nu-1)
    d2lbetabeta<-matrix(0,p,p)
 for (i in 1:n) {
    dibeta=as.matrix(derbeta1[,i])
    aux=-nu*dnu[i]/sigma2*((2*nu-1)*dibeta%*%t(dibeta)-res[i]*derbeta2[[i]])
    d2lbetabeta<-d2lbetabeta+aux
    }
 d2lsigmabeta<-matrix(0,p,1)
 for (i in 1:n) {
    aux=-nu^2*res[i]/(sigma2^2)*dnu[i]*as.matrix(derbeta1[,i])
    d2lsigmabeta<-d2lsigmabeta+aux
    }
 d2lnubeta<-matrix(0,p,1)
 for (i in 1:n) {
    aux=res[i]/sigma2*dnu[i]*(1+nu*log(d[i]))*as.matrix(derbeta1[,i])
    d2lnubeta<-d2lnubeta+aux
    }
    I1[1:p,1:p]<- d2lbetabeta
    I1[1:p,p+1]<- d2lsigmabeta
    I1[1:p,p+3]<- d2lnubeta
    I1[p+1,p+1]<- 1/(2*(sigma2^2))*(n-nu*(nu+1)*sum(d^nu))
    I1[p+1,p+3]<- 1/(2*sigma2)*t(d^nu)%*%(matrix(1,n,1)+nu*log(d))
    nu1<- 1/(2*nu)
    I1[p+3,p+3]<- -n/(nu^2)*(1+(log(2)+psigamma(nu1))/nu+(nu1^2)*psigamma(nu1,1))-1/2*t(d^nu)%*%((log(d))^2)
    I1[p+1,1:p]<- t(I1[1:p,p+1])
    I1[p+3,1:p]<- t(I1[1:p,p+3])
    I1[p+3,p+1]<- I1[p+1,p+3]    
 }
 ##########################################################
 IM<--(I1+I2)
 return(IM)
}

# Score functions

ssmn.nl.score<- function(y,x,theta,p,derbeta1,family = "Skew.normal") 
{
 # THETA: (beta,sigma2,lambda,nu,gama) ; nu para StN e SSL; nu e gama para SCN
 # usando derivadas em beta do R
 # testa se U(theta est) = 0
 der1beta<-rbind(derbeta1)
 n=length(y)
 beta<-theta[1:p]
 sigma2<-as.numeric(theta[p+1])
 lambda<-as.numeric(theta[p+2])
 mu=nlf(x,beta)
 res<-as.matrix(y-mu)
 sigma<-sqrt(sigma2)
 aux<-lambda*res/sigma
 aux1<-pmax(aux,-37)
 Wphi<-dnorm(aux1)/pnorm(aux1)
 d<-res^2/sigma2
 ## Assimétrica 
 dl2sigma2=-lambda/(2*(sigma^3))*res
 dl2lambda=1/sigma*res
 dl2beta=-lambda/sigma*der1beta
 #dl2theta=matrix(c(dl2beta,dl2sigma2,dl2lambda),p+2,n)
 dl2theta=rbind(dl2beta,t(dl2sigma2),t(dl2lambda))
 Uasim=dl2theta%*%diag(as.vector(Wphi)) # vetor de score, para todas as observações
 ## Parte Simétrica: (beta,sigma2, , nu, gama) 
 if (family=="Skew.normal"){
    dl1beta=der1beta%*%diag(as.vector(res))/sigma2
    dl1sigma2=-1/(2*sigma2)+1/(2*(sigma2^2))*(res^2)
    dl1lambda=rep(0,n)
    dl1theta=rbind(dl1beta,t(dl1sigma2),t(dl1lambda))
 }
 if (family=="Skew.t"){
    nu<-as.numeric(theta[p+3])
    # Parte Assimétrica
    Utasim=matrix(0,p+3,n)
    Utasim[1:(p+2),]=Uasim
    Uasim=Utasim
    # Parte Simétrica
    V<-(1+d/nu)^(-1)
    h1nu<-1/2*(digamma((nu+1)/2)-digamma(nu/2)-1/nu)
    dl1beta=(nu+1)/(nu*sigma2)*der1beta%*%diag(as.vector(res*V))
    dl1sigma2=-1/(2*sigma2)+(nu+1)/(2*nu*sigma2)*V*d
    dl1lambda=rep(0,n)
    dl1nu=h1nu+1/2*(log(V)+(nu+1)/(nu^2)*V*d)
    dl1theta=rbind(dl1beta,t(dl1sigma2),t(dl1lambda),t(dl1nu))
 }
 if (family=="Skew.slash"){
    nu<-as.numeric(theta[p+3])
    # Parte Assimétrica
    Utasim=matrix(0,p+3,n)
    Utasim[1:(p+2),]=Uasim
    Uasim=Utasim
    # Parte Simétrica
    P11<-pgamma(1,nu+0.5,scale=2/d)
    P13<-pgamma(1,nu+1.5,scale=2/d)
    P15<-pgamma(1,nu+2.5,scale=2/d)
    IG1<-gamma(nu+0.5)*P11/((d/2)^(nu+0.5))
    IG3<-gamma(nu+1.5)*P13/((d/2)^(nu+1.5))
    rp<-5000
    u1<-matrix(0,n,rp)
    for (j in 1:n) {
       u1[j,]<-rtrunc(rp,"gamma",a=0,b=1,shape=nu+1/2,rate=d[j]/2)
    }
    E11=P11*rowMeans(log(u1))  #vetor n x 1
    dl1beta= 1/sigma2*der1beta%*%diag(as.vector(res*IG3/IG1))
    dl1sigma2= 1/(2*sigma2)*(-1+d*IG3/IG1)
    dl1lambda=rep(0,n)
    dl1nu= 1/nu+E11/IG1
    dl1theta=rbind(dl1beta,t(dl1sigma2),t(dl1lambda),t(dl1nu))
 }
 if (family=="Skew.cn"){
     nu<-as.numeric(theta[p+3])
     gama<-as.numeric(theta[p+4])
     # Parte Assimétrica
     Utasim=matrix(0,p+4,n)
     Utasim[1:(p+2),]=Uasim
     Uasim=Utasim
     # Parte Simétrica
     fi<-dnorm(y,mu,sigma)
     figama<-dnorm(y,mu,sqrt(sigma2/gama))
     fs<-nu*figama+(1-nu)*fi
     dfssigma<-1/(2*sigma2)*(nu*figama*(gama*d-1)+(1-nu)*fi*(d-1))
     dl1sigma2=dfssigma/fs
     dfsnu<-figama-fi
     dl1nu=dfsnu/fs
     dfsgama<-nu/(2*gama)*figama*(1-gama*d)
     dl1gama=dfsgama/fs
     dl1beta=matrix(0,p,n)
     for (i in 1:n) {
        dibeta=as.matrix(derbeta1[,i])
        auxn<-nu*gama*figama[i]+(1-nu)*fi[i]
        dl1beta[,i]=1/sigma2*auxn*res[i]*dibeta/fs[i]
    }
  #### 
    dl1lambda=rep(0,n)
    dl1theta=rbind(dl1beta,t(dl1sigma2),t(dl1lambda),t(dl1nu),t(dl1gama))
     
 }
 if (family=="Skew.pe"){
    nu<-as.numeric(theta[p+3])
    # Parte Assimétrica
    Utasim=matrix(0,p+3,n)
    Utasim[1:(p+2),]=Uasim
    Uasim=Utasim
    # Parte Simétrica
    dl1beta=matrix(0,p,n)
    dnu<- d^(nu-1)
    for (i in 1:n) {
    dibeta=as.matrix(derbeta1[,i])
    dl1beta[,i]=nu*dnu[i]/sigma2*dibeta*res[i]
    }
    dl1sigma2=1/(2*sigma2)*(-1+nu*(d^nu))
    dl1nu=1/nu+1/(nu^2*2)*(digamma(1/(2*nu))+log(2))-1/2*(d^nu)*log(d)
    #### 
    dl1lambda=rep(0,n)
    dl1theta=rbind(dl1beta,t(dl1sigma2),t(dl1lambda),t(dl1nu))
  }
 ##########################################################
 Score<-apply(dl1theta+Uasim,1,sum)  # aplicar soma de U%*%U', U_ =dl1(theta)+Wphi_1*dl2(theta)
 return(Score)
}

# Simulated envelopes

envelSN <-function(Y,mu,Sigma){
n=length(Y)
p=1
d2=(Y-mu)^2/Sigma
d2s=sort(d2)
d2s=t(d2s)
xq2 <- qchisq(ppoints(n), p)
Xsim<-matrix(0,100,n)
for(i in 1:100){Xsim[i,]<-rchisq(n, p)}
Xsim2<-apply(Xsim,1,sort)# aplica nas linhas, neste caso, coloca cada coluna em ordem crescente
d21<-rep(0,n)
d22<-rep(0,n)
for(i in 1:n){
d21[i]  <- quantile(Xsim2[i,],0.025)
d22[i]  <- quantile(Xsim2[i,],0.975)}
d2med <-apply(Xsim2,1,mean) # tira média de cada coluna; tira médias dos vetores menores...até os maiores
fy <- range(d2s,d21,d22)
plot(xq2,d2s,xlab = expression(paste("Theoretical ",chi[1]^2, " quantiles")),
     ylab="Sample values and simulated envelope",pch=20,ylim=fy)
par(new=T)
plot(xq2,d21,type="l",ylim=fy,xlab="",ylab="")
par(new=T)
plot(xq2,d2med,type="l",ylim=fy,xlab="",ylab="",lty="dashed")
par(new=T)
plot(xq2,d22,type="l",ylim=fy,xlab="",ylab="")
}

envelStN <-function(Y,mu,Sigma,nu){
n=length(Y)
p=1
d2=(Y-mu)^2/Sigma
d2s=sort(d2)
xq2 <- qchisq(ppoints(n), p)
Xsim<-matrix(0,200,n)
for(i in 1:200){Xsim[i,]<-rf(n, p,nu)}
Xsim2<-apply(Xsim,1,sort)
d21<-rep(0,n)
d22<-rep(0,n)
for(i in 1:n){
d21[i]  <- quantile(Xsim2[i,],0.025)
d22[i]  <- quantile(Xsim2[i,],0.975)}
d2med <-apply(Xsim2,1,mean)
fy <- range(d2s,d21,d22)
plot(xq2,d2s,xlab = expression(paste("Theoretical ",F(1,nu), " quantiles")),
     ylab="Sample values and simulated envelope",pch=20,ylim=fy)
par(new=T)
plot(xq2,d21,type="l",ylim=fy,xlab="",ylab="")
par(new=T)
plot(xq2,d2med,type="l",ylim=fy,xlab="",ylab="",lty="dashed")
par(new=T)
plot(xq2,d22,type="l",ylim=fy,xlab="",ylab="")
}

envelSSL <-function(Y,mu,Sigma,nu){
n=length(Y)
p=1
d2=(Y-mu)^2/Sigma
d2s1=sort(d2)
d2s=t(d2s1)
xq2 <- qchisq(ppoints(n), p)
Xsim<-matrix(0,200,n)
for(i in 1:200){
       aux=GeraSlash(n,nu)
       Xsim[i,]<-aux^2}
Xsim2<-apply(Xsim,1,sort)
d21<-rep(0,n)
d22<-rep(0,n)
for(i in 1:n){
d21[i]  <- quantile(Xsim2[i,],0.025)
d22[i]  <- quantile(Xsim2[i,],0.975)}
d2med <-apply(Xsim2,1,mean)

fy <- range(d2s,d21,d22)
plot(xq2,d2s,xlab = expression(paste("Theoretical ",chi[1]^2, " quantiles")),
     ylab="Sample values and simulated envelope",pch=20,ylim=fy)
par(new=T)
plot(xq2,d21,type="l",ylim=fy,xlab="",ylab="")
par(new=T)
plot(xq2,d2med,type="l",ylim=fy,xlab="",ylab="",lty="dashed")
par(new=T)
plot(xq2,d22,type="l",ylim=fy,xlab="",ylab="")
}

envelSCN <- function(Y,mu,Sigma,nu,gama){
n=length(Y)
p=1
d2=(Y-mu)^2/Sigma
d2s=sort(d2)
xq2 <- qchisq(ppoints(n), p)
Xsim<-matrix(0,200,n)
for(i in 1:200){
      u1<-rchisq(n, p)/gama
      u2<-rchisq(n, p)
      u3<-runif(n)
      id<-(u3<nu)
      u2[id]<-u1[id]
      Xsim[i,]=u2}
Xsim2<-apply(Xsim,1,sort)
d21<-rep(0,n)
d22<-rep(0,n)
for(i in 1:n){
d21[i]  <- quantile(Xsim2[i,],0.025)
d22[i]  <- quantile(Xsim2[i,],0.975)}
d2med <-apply(Xsim2,1,mean)

fy <- range(d2s,d21,d22)
plot(xq2,d2s,xlab = expression(paste("Theoretical ",chi[1]^2, " quantiles")),
     ylab="Sample values and simulated envelope",pch=20,ylim=fy)
par(new=T)
plot(xq2,d21,type="l",ylim=fy,xlab="",ylab="")
par(new=T)
plot(xq2,d2med,type="l",ylim=fy,xlab="",ylab="",lty="dashed")
par(new=T)
plot(xq2,d22,type="l",ylim=fy,xlab="",ylab="")
}

envelSPE <-function(Y,mu,Sigma,nu){
n=length(Y)
p=1
d2=(Y-mu)^2/Sigma
d2s=sort(d2)
n <-length(d2s)
xq2 <- qchisq(ppoints(n), p)
Xsim<-matrix(0,2000,n)
for(i in 1:2000){
       aux=GeraEP(n,nu)
       Xsim[i,]<-aux^2}
Xsim2<-apply(Xsim,1,sort)
d21<-rep(0,n)
d22<-rep(0,n)
for(i in 1:n){
d21[i]  <- quantile(Xsim2[i,],0.025)
d22[i]  <- quantile(Xsim2[i,],0.975)}
d2med <-apply(Xsim2,1,mean)

fy <- range(d2s,d21,d22)
plot(xq2,d2s,xlab = expression(paste("Theoretical ",chi[1]^2, " quantiles")),
     ylab="Sample values and simulated envelope",pch=20,ylim=fy)
par(new=T)
plot(xq2,d21,type="l",ylim=fy,xlab="",ylab="")
par(new=T)
plot(xq2,d2med,type="l",ylim=fy,xlab="",ylab="",lty="dashed")
par(new=T)
plot(xq2,d22,type="l",ylim=fy,xlab="",ylab="")
}

envel_norm_Pearson <-function(res){
# 
par(mfrow=c(1,1))
n <- length(res)
p <- 1
X=matrix(1,n,p)
H <- X%*%solve(t(X)%*%X)%*%t(X)
h <- diag(H)
si <- 1
r <- res
tsi <- r/(si*sqrt(1-h))
#
ident <- diag(n)
epsilon <- matrix(0,n,100)
e <- matrix(0,n,100)
e1 <- numeric(n)
e2 <- numeric(n)
#
for(i in 1:100){
     epsilon[,i] <- rnorm(n,0,1)
     e[,i] <- (ident - H)%*%epsilon[,i]
     u <- diag(ident - H)
     e[,i] <- e[,i]/sqrt(u)
     e[,i] <- sort(e[,i]) }
#
for(i in 1:n){
     eo <- sort(e[i,])
     e1[i] <- eo[5]#(eo[2]+eo[3])/2
     e2[i] <- eo[95]#(eo[97]+eo[98])/2 
     }
#
med <- apply(e,1,mean)
faixa <- range(tsi,e1,e2)
#
par(pty="s")
pp=qqnorm(tsi,xlab="Percentiles of N(0,1)",
ylab="Pearson Residuals", ylim=faixa, pch=16)
par(new=T)
qqnorm(e1,axes=F,xlab="",ylab="",type="l",ylim=faixa,lty=1)
par(new=T)
qqnorm(e2,axes=F,xlab="",ylab="", type="l",ylim=faixa,lty=1)
par(new=T)
qqnorm(med,axes=F,xlab="",ylab="",type="l",ylim=faixa,lty=2)
#--------------------------------------------------------------#
}

envel_chi2_Pearson <-function(res){

par(mfrow=c(1,1))
n <- length(res)
p <- 1
X=matrix(1,n,p)
H <- X%*%solve(t(X)%*%X)%*%t(X)
h <- diag(H)
si <- 1
r <- res
tsi <- r/(si*sqrt(1-h))
#
ident <- diag(n)
epsilon <- matrix(0,n,100)
e <- matrix(0,n,100)
e1 <- numeric(n)
e2 <- numeric(n)
#
for(i in 1:100){
     epsilon[,i] <- rnorm(n,0,1)
     e[,i] <- (ident - H)%*%epsilon[,i]
     u <- diag(ident - H)
     e[,i] <- e[,i]/sqrt(u)
     e[,i] <- sort(e[,i]) }
#
for(i in 1:n){
     eo <- sort(e[i,])
     e1[i] <- (eo[2]+eo[3])/2
     e2[i] <- (eo[97]+eo[98])/2 }
#
med <- apply(e,1,mean)
faixa <- range(tsi,e1,e2)
#
par(pty="s")
pp=qqnorm(tsi,xlab="Percentiles of N(0,1)",
ylab="Pearson Residuals", ylim=faixa, pch=16)
par(new=T)
qqnorm(e1,axes=F,xlab="",ylab="",type="l",ylim=faixa,lty=1)
par(new=T)
qqnorm(e2,axes=F,xlab="",ylab="", type="l",ylim=faixa,lty=1)
par(new=T)
qqnorm(med,axes=F,xlab="",ylab="",type="l",ylim=faixa,lty=2)
#--------------------------------------------------------------#
}

# Auxiliar functions

chute.inicial <- function(beta1,Y,X,nlf){
mu=nlf(X,beta1)
fmin=sum(abs(Y-mu))
return(fmin)
}

betamixt<- function (a,b,N)
{
  # gera de uma distribuiçao de misturas de betas
  wb=matrix(1,N,1)
  wt=wb
  for (j in 2:N)
  {
      wb[j]=wb[j-1]*b/(a+j-1)
  }
  for (j in 2:N)
  {
      wt[j]=wt[j-1]+wb[j]
  }
  wt=rbind(0,wt/wt[N])
  u=runif(1,0,1)
  dif=wt-u
  k=length(which(dif<0))
  Y=rbeta(1,a,k)
  return(Y)
}


GeraSlash<- function(n,nu){
  u1 <- runif(n)
  u2 <- u1^(1/(nu))
  u3<-rep(0,n)
  for(j in 1:n){u3[j] <-  rnorm(1, 0, sd=1/sqrt(u2[j]))}
  return(u3)
  }
 
GeraEP<- function(n, nu)
{ 
t=rgamma(n,0.5*nu,0.5)
R=t^(1/(2*nu))
u=sign(rnorm(n))
y=u*R
return(y)
}

GeraSEP<- function(n, lambda, nu)
{ 
# Gerando Y ~ EP(0,1,nu)
t=rgamma(n,0.5*nu,1/2)
R=t^(1/(2*nu))
u=sign(rnorm(n))
yep=u*R
# Aplicando forma de Wang, Boyer and Genton
w=rnorm(n)
yv=yep*sign(lambda*yep-w)
return(yv)
}

gamtruncrnd<- function(a,b,t,n,m)
{
 # Referência: ANNE PHILIPPE. Simulation of right and left truncated gamma
 # distributions by mixtures. Statistics and Computing (1997) 7, 173-181

 # a: parametro de forma/locaçao
 # b: parametro de taxa
 # t: truncamento a direita, f_t varia de 0 a t
 # [n,m]: tamanho da amostra
 # Se X ~ TG(a,b,t) então X/t ~ TG(a,bt,1)
 # N: nº de replicaçoes para garantir que a probabilidade de aceitação do
 # algoritmo de aceitação-rejeição P(N) > p.
  tp=qnorm(0.99,0,1)
  N=ceiling((tp+sqrt(tp^2+4*b))^2/4) # ceiling(X): inteiro maior que X
  N=min(N,171)
  K=matrix(1:N)
  M=1/sum(b^(K-1)/gamma(K))
  Y=matrix(0,n,m)
  for (j in 1:n)
  {
      for (k in 1:m)
      {
          u=runif(1,0,1)
          x=betamixt(a,b,N)
          aux=((1-x)*b)^(K-1)/gamma(K)
          rhox=1/(exp(b*x)*sum(aux))
          Y[j,k]=x
          while (u*M > rhox)
          {
              u=runif(1,0,1)
              x=betamixt(a,b,N)
              aux=((1-x)*b)^(K-1)/gamma(K)
              rhox=1/(exp(b*x)*sum(aux))
              Y[j,k]=x
          }
      }
  }
  return(Y)
}
  
  

















