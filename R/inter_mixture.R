library(survival)

########################################################################
# R Codes for implementing the analysis of the smoking cessation data #
# using mixture cure models with/without random effects. #
# SIM paper: Mixture cure model with random effects for clustered #
# interval-censored survival data #
# By Liming Xiang, Xiangmei Ma and Kelvin K.W.Yau, Nov 2010 #
########################################################################
# Use survival package
# The following R codes include two parts:
# Part 1 :Mixture cure model for smoking cessation data
# Part 2 :Mixture cure model with random effects for smoking cessation data
# The data are extracted from http://www.tibs.org/biometric/index.html
# The data set reports 223 subjects who are recruited from 51 zip codes
# in the southeastern corner of Minnesota.
# Description of variables
# SexF: 1 = Female, 0 = Male
# Duration: Duration as smoker in years
# SI/UC: 1 = Special Intervention(SI), 0 = Usual Care(UC)
# F10Cigs: Average number of cigarrettes smoked per day over
# the last 10 years (rounded)
# Noted smoking cessation data are interval censored
# left <- Timept1
# right <- Timept2
# X <- as.matrix(cbind(SexF, SI/UC, F10Cigs, Duration))
# patient <- Zip
# Description of functions in the following R codes
# loglik: Construct the negative complete-data log-likelihood function
# gradlik: Construct the first derivative of negative log-likelihood
# mixture: Use the proposed mixture cure model to fit the data,
# wreml.logit: Estimate the parameters and their asymptotic variances
# for logistic part
# loglik.effect: Construct the negative complete-data BLUP type
# log-likelihood function only for beta part
# gradlik.effect: Construct the first derivative of negative BLUP type
# log-likelihood function only for beta part
# mixture.effect: Use the proposed mixture cure model with random effects
# to fit the smoking cessation data, so as to obtain the
# estimates of parameters and their asymptotic variances
# 1.mixture(left, right, X, itmax=100)
# Fit data to the mixture cure model with random effects ignored
# 2.mixture.effect(left, right, X, patient, itmax=100)
# Fit data to the mixture cure model with random effects
##########################################################
# Part 1 :Mixture cure model with random effects ignored #
##########################################################
loglik <- function(omega, alpha, x, y) {
  #cat("LOG COME?OU!!!\n")
  #theta_time<-proc.time()
  p <- ncol(x)
  n <- nrow(x)
  q <- ncol(alpha)
  beta0 <- omega[1:p]
  gamma0 <- omega[-(1:p)]
  xb <- as.vector(x%*%beta0)
  one <- matrix(1,n,q-1)
  linb <- exp(xb*one+t(gamma0*t(one)))
  s <- exp(-t(apply(linb,1,cumsum)))
  #cat("Tempo loglik antes do g\n")
  #print(proc.time()-theta_time)
  g <- rowSums(alpha*(cbind(1,s)-cbind(s,0)))
  sum(-y*log(g))
}#end of function
########################################################
gradlik <-function(omega, alpha,x,y){
  #cat("COME?OU O GRAD EIN\n")
  #grad_time<-proc.time()
  p <- ncol(x)
  n <- nrow(x)
  q <- ncol(alpha)
  beta0 <- omega[1:p]
  gamma0 <- omega[-(1:p)]
  xb <- as.vector(x%*%beta0)
  one <- matrix(1,n,q-1)
  linb <- exp(xb*one+t(gamma0*t(one)))
  s <- exp(-t(apply(linb,1,cumsum)))
  f <- s*log(s)
  ff <- cbind(0,f)-cbind(f,0)
  g <- rowSums(alpha*(cbind(1,s)-cbind(s,0)))
  c1 <- (alpha[,-q]-alpha[,-1])*s
  c2 <- t(apply(c1,1,cumsum))
  #cat("Tempo gradlik antes dos l\n")
  #print(proc.time()-grad_time)
  c <- rowSums(c1)-cbind(0,c2[,-ncol(c2)])
  l.beta <- t(x)%*%(y*rowSums(alpha*ff)/g)
  l.gamma <- as.vector(colSums(y*(linb*c)/g))
  c(-l.beta,-l.gamma)
}#end of function
#########################################################
mixture <- function(l,r,cov_theta,cov_beta,dataset,itmax){
  X_theta<-as.matrix(dataset[,cov_theta])
  X_beta<-as.matrix(dataset[,cov_beta])
  p_0<-ncol(X_theta)
  p_1<-ncol(X_beta)
  n <- nrow(dataset) # the no. of covariates
  w <- cbind(1,X_theta)
  sdelta <- ifelse(r==Inf,0,1)
  tk <- max(r[sdelta==1])
  tau <- sort(unique(c(0,l,r))) #ordered distinct all interval censored time
  interv <- function(x,l,r) ifelse(x[1]>=l & x[2]<=r,1,0)
  tau12 <- cbind(tau[-length(tau)],tau[-1])
  alpha <- apply(tau12,1,interv,l=l,r=r)
  q <- length(tau)-1
  one <- matrix(1,n,q-1)
  b0 <- rep(0,(p_1+1))
  beta0 <- rep(0,p_0)
  #print(beta0)
  ekm <- survfit(Surv(tau[1:q],rep(1,q))~1) #library(survival)
  G<- ekm$surv
  lg <- log(G[-length(G)])
  ga <- c(0,lg[-length(lg)])-lg
  gamma <- log(ga)
  print("Tamanho gamma: ")
  print(length(gamma))
  print("GAMMA \n")
  print(gamma)
  delta0 <- as.vector(c(beta0,gamma)) #initial vatues of beta,gamma
  eta<-as.vector(X_beta%*%beta0)
  ksi <- as.vector(w%*%b0)
  pai<- 1-1/(1+exp(ksi)) # the uncured probability
  linb <- exp(eta*one+t(gamma*t(one)))
  s <- exp(-t(apply(linb,1,cumsum)))#survival function for uncured individuals
  g <- rowSums(alpha*(cbind(1,s)-cbind(s,0)))
  y <- pai*g/(1-pai+pai*g)
  y[sdelta==1] <- 1 #the latent random variable
  #cat("CHEGUEI AQUI\n")
  for(iter in 1: itmax) {
    print("ITER \n")
    print(iter)
    print("\n")
    y2 <- ifelse(l>tk,0,y)
    lb <- glm(y2~X_theta,family = binomial(link=logit))
    b <- coef(lb)
    se.b <- as.vector(coef(summary(lb))[,2])
    varB<-vcov(lb)
    ml <- optim(delta0, loglik, gradlik, method="BFGS", hessian=TRUE,
                alpha=alpha, x=X_beta, y=y)
    #maximize the complete-data log-likelihood function to estimate beta & gamma
    delta <- ml$par
    beta <- delta[1:p_1]
    gamma <- delta[-(1:p_1)]
    H <- ml$hessian
    se.beta <- sqrt(diag(solve(H[1:p_1,1:p_1])))
    eta <- as.vector(X_beta%*%beta)
    ksi <- as.vector(w%*%b)
    pai <- 1/(1+exp(-ksi))
    linb <- exp(eta*one+t(gamma*t(one)))
    s <- exp(-t(apply(linb,1,cumsum)))
    g <- rowSums(alpha*(cbind(1,s)-cbind(s,0)))
    #cat("CHEGUEI AQUI 3\n")
    y0 <- pai*g/(1-pai+pai*g)
    y0[sdelta==1] <- 1
    dy <- max(abs(y-y0))
    if(dy<0.0001)
    { cat("iter1.max=",iter, "\n","beta=",beta, "\n", "b=",b,"\n")
      break }
    else {
      delta0 <- delta
      y <- y0 }
    #cat("Iterations= ",iter, "dy= ",dy,"\n")
    #cat("beta= ",beta,"\n", "b= ",b, "\n")
    #cat("se.beta= ",se.beta,"\n", "se.b= ",se.b,"\n")
  }
  stdb1<-se.b
  stdb2<-se.beta
  bb1<-cbind(b,se.b,2*(1-pnorm(abs(b/se.b))))
  names(b)<-c("Intercept",cov_theta)
  dimnames(bb1)<-list(names(b),c("Estimate","SE","p-value"))
  bb2<-cbind(beta,se.beta,2*(1-pnorm(abs(beta/se.beta))))
  names(beta)<-cov_beta
  dimnames(bb2)<-list(names(beta),c("Estimate","SE","p-value"))
  cat("\n", "Fit of mixture cure model with random effects ignored ","\n")
  cat("\n", "Parameters in the logistic regression: ","\n")
  print(round(bb1,4))
  cat("\n", "Parameters in the survival model: ","\n")
  print(round(bb2,4))
  cat("\n","gamma=",round(gamma,4), "\n")
  cat("\n")
  cPar<-as.numeric(iter==itmax)
  return(list(beta=beta, b=b, gamma=gamma, se.beta=se.beta, se.b=se.b, varB=varB,StopC=cPar))
}#end of function
#############################################################
# Part 2 :Mixture cure model with random effects #
#############################################################
wreml.logit<-function(y,x,z,b1,v1,sig2,family="logistic",cov_names="") {
  M <- ncol(z)
  n <- length(y)
  X <- cbind(rep(1,length(y)),x)
  p1 <- ncol(X)
  zero1 <- matrix(0,ncol=p1,nrow=M)
  X1 <- rbind(X,zero1)
  Z <- rbind(z,diag(M))
  XX <- cbind(X1,Z)
  itmax <- 1000
  b0 <- c(b1,v1)
  b <- b1
  v <- v1
  flag <- 0
  for(iter in 1:itmax) {
    theta <- as.vector(X%*%b+z%*%v)
    w1 <- exp(theta)/(1+exp(theta))^2
    w2 <- c(w1,rep(1/sig2,M))
    mu <- exp(theta)/(1+exp(theta))
    ply <- c(as.vector(t(X)%*%(y-mu)),as.vector(t(z)%*%(y-mu)-v/sig2))
    w2 <- t(matrix(rep(w2,(p1+M)),ncol=(p1+M)))
    A1 <- (t(XX)*w2)%*%XX
    A <- solve(A1)
    bb <- b0+A%*%ply
    b <- bb[1:p1]
    v <- bb[(p1+1):(p1+M)]
    std.b <- sqrt(diag(A)[1:p1])
    sig2 <- as.vector(t(v)%*%v+sum(diag(A)[(p1+1):(p1+M)]))/M
    r2 <- sum(diag(A)[(p1+1):(p1+M)])/sig2
    se.sig2 <- sqrt(2*sig2^2/(M-2*r2+sum(diag(A)[(p1+1):(p1+M)]^2)/sig2^2))
    if (max(abs(bb-b0))<0.001) {flag<-1;break}
    else b0 <- bb }
  if(flag) {
    reslt <- list(b=b,v=v,sig2=sig2,std.b=std.b,se.sig2=se.sig2) }
  else stop("error:not reach the convergence")
}#end of function
###########################################################
loglik.effect <- function(omega,alpha,xl,yl,cluster,th1) {
  p <- ncol(xl)
  n <- nrow(xl)
  q <- ncol(alpha)
  m <- length(unique(cluster))
  z <- matrix(0,ncol=m,nrow=n)
  for(j in 1:m) {
    z[,j] <- ifelse(cluster==unique(cluster)[j],1,0) }
  beta0 <- omega[1:p]
  u <- omega[(p+1):(p+m)]
  gamma0 <- omega[-(1:(p+m))]
  eta <- as.vector(xl%*%beta0+z%*%u)
  one <- matrix(1,n,q-1)
  linb <- exp(eta*one+t(gamma0*t(one)))
  s <- exp(-t(apply(linb,1,cumsum)))
  g <- rowSums(alpha*(cbind(1,s)-cbind(s,0)))
  sum(-yl*log(g))+1/2*(m*log(2*pi*th1)+1/th1*sum(u*u))
}#end of function
###########################################################
gradlik.effect <- function(omega, alpha, xl, yl, cluster,th1){
  p <- ncol(xl)
  n <- nrow(xl)
  q <- ncol(alpha)
  m <- length(unique(cluster))
  z <- matrix(0,ncol=m,nrow=n)
  for(j in 1:m) {
    z[,j] <- ifelse(cluster==unique(cluster)[j],1,0) }
  beta0 <- omega[1:p]
  u <- omega[(p+1):(p+m)]
  gamma0 <- omega[-(1:(p+m))]
  eta <- as.vector(xl%*%beta0+z%*%u)
  one <- matrix(1,n,q-1)
  linb <- exp(eta*one+t(gamma0*t(one)))
  s <- exp(-t(apply(linb,1,cumsum)))
  f <- s*log(s)
  ff <- cbind(0,f)-cbind(f,0)
  g <- rowSums(alpha*(cbind(1,s)-cbind(s,0)))
  c1 <- (alpha[,-q]-alpha[,-1])*s
  c2 <- t(apply(c1,1,cumsum))
  c <- rowSums(c1)-cbind(0,c2[,-ncol(c2)])
  l.beta <- t(xl)%*%(yl*rowSums(alpha*ff)/g)
  l.u <- as.vector(t(z)%*%(yl*rowSums(alpha*ff)/g)-u/th1)
  l.gamma <- as.vector(colSums(yl*(linb*c)/g))
  c(-l.beta,-l.u,-l.gamma)
}#end of function
###############################################################
mixture.effect <- function(l,r,cov_theta,cov_beta,cluster,dataset,itmax){
  X_theta<-as.matrix(dataset[,cov_theta])
  X_beta<-as.matrix(dataset[,cov_beta])
  p_0<-ncol(X_theta)
  p_1<-ncol(X_beta)
  n <- nrow(dataset)
  w <- cbind(1,X_theta)
  sdelta <- ifelse(r==Inf,0,1)
  tk <- max(r[sdelta==1])
  tau <- sort(unique(c(0,l,r)))
  interv <- function(x,l,r) ifelse(x[1]>=l & x[2]<=r,1,0)
  tau12 <- cbind(tau[-length(tau)],tau[-1])
  alpha <- apply(tau12,1,interv,l=l,r=r)
  q <- length(tau)-1
  m <- length(unique(cluster))
  z <- matrix(0,ncol=m,nrow=n)
  for(j in 1:m) {
    z[,j] <- ifelse(cluster==unique(cluster)[j],1,0) }
  one <- matrix(1,n,q-1)
  b0 <- rep(0,(p_0+1))
  beta0 <- rep(0,p_1)
  ekm <- survfit(Surv(tau[1:q],rep(1,q))~1)
  G <- ekm$surv
  lg <- log(G[-length(G)])
  ga <- c(0,lg[-length(lg)])-lg
  gamma <- log(ga)
  u <- as.vector(rep(0,m))
  v <- as.vector(rep(0,m))
  delta0 <- as.vector(c(beta0,u,gamma))
  th1 <- 1; th2 <- 1
  eta <- as.vector(X_beta%*%beta0+z%*%u)
  ksi <- as.vector(w%*%b0+z%*%v)
  pai<- 1/(1+exp(-ksi))
  linb <- exp(eta*one+t(gamma*t(one)))
  s <- exp(-t(apply(linb,1,cumsum)))
  g <- rowSums(alpha*(cbind(1,s)-cbind(s,0)))
  y <- pai*g/(1-pai+pai*g)
  y[sdelta==1] <- 1
  cat("CHEGUEI AQUI\n")
  for(iter in 1: itmax) {
    y2 <- ifelse(l>tk,0,y)
    glb <- wreml.logit(y2,X_theta,z,b0,v,sig2=th2)
    b <- glb$b
    v <- glb$v
    th2 <- glb$sig2
    se.b <- glb$std.b
    se.th2 <- glb$se.sig2
    mle <- optim(delta0,loglik.effect,gradlik.effect,method="BFGS",hessian=TRUE,alpha=alpha,xl=X_beta, yl=y, cluster=cluster,th1=th1)
    # maximize the BLUP type log-likelihood function to estimate beta,u & gamma
    delta0 <- mle$par
    beta <- delta0[1:p_1]
    u <- delta0[(p_1+1):(p_1+m)]
    gamma <- delta0[-(1:(p_1+m))]
    H1 <- mle$hessian
    B <- H1[1:(p_1+m),1:(p_1+m)]
    if(any(is.na(B))||any(is.nan(B))||any(is.infinite(B)))
    { flag <- 1; break }
    if(qr(B)$rank<ncol(qr(B)$qr))
    { flag <- 1; break }
    else {
      flag <- 0
      H <- solve(B) }
    se.beta <-sqrt(diag(H)[1:p_1])
    th1 <- as.vector(t(u)%*%u+sum(diag(H)[(p_1+1):(p_1+m)]))/m
    r1 <- sum(diag(H)[(p_1+1):(p_1+m)])/th1
    se.th1 <- sqrt(2*th1^2/(m-2*r1+sum(diag(H)[(p_1+1):(p_1+m)]^2)/th1^2))
    eta <- as.vector(X_beta%*%beta+z%*%u)
    ksi <- as.vector(w%*%b+z%*%v)
    pai <- 1/(1+exp(-ksi))
    linb <- exp(eta*one+t(gamma*t(one)))
    s <- exp(-t(apply(linb,1,cumsum)))
    g <- rowSums(alpha*(cbind(1,s)-cbind(s,0)))
    y0 <- pai*g/(1-pai+pai*g)
    y0[sdelta==1] <- 1
    dy <- max(abs(y-y0))
    if(max(abs(y-y0))<0.001) {
      cat("iter2.max=",iter, "\n","beta=",beta, "\n","b=",b,
          "\n","th1=",th1,"th2=",th2,"\n")
      break }
    else y<-y0
    #cat("Iterations= ",iter, "dy= ",dy,"\n")
    #cat("beta= ",beta, "b= ",b, "\n")
    #cat("se.beta= ",se.beta,"\n", "se.b= ",se.b,"\n")
    #cat("th1= ",th1, "th2= ",th2, "\n")
    #cat("se.th1= ",se.th1,"\n", "se.th2= ",se.th2,"\n")
  }
  sig1<-th1
  sig2<-th2
  stdb1<-se.b
  stdb2<-se.beta
  stds1<-se.th1
  stds2<-se.th2
  ss1<-cbind(sig1,stds1)
  dimnames(ss1)<-list(names(sig1),c("Estimate","SE"))
  ss2<-cbind(sig2,stds2)
  dimnames(ss2)<-list(names(sig2),c("Estimate","SE"))
  bb1<-cbind(b,se.b,2*(1-pnorm(abs(b/se.b))))
  names(b)<-c("Intercept",cov_theta)
  dimnames(bb1)<-list(names(b),c("Estimate","SE","p-value"))
  bb2<-cbind(beta,se.beta,2*(1-pnorm(abs(beta/se.beta))))
  names(beta)<-cov_beta
  dimnames(bb2)<-list(names(beta),c("Estimate","SE","p-value"))
  cat("\n","Fit of mixture cure model with random effects","\n")
  cat("\n", "Parameters in the logistic regression: ","\n")
  print(round(bb1,4))
  cat("\n", "Random effect U= ", round(u,4), "\n")
  cat("\n", "Variance component: sig1^2 ","\n")
  cat(round(sig1, 4),"(",round(stds1,4),")", "\n")
  cat("\n", "Parameters in the survival model: ","\n")
  print(round(bb2,4))
  cat("\n", "Random effect V= ", round(v,4), "\n")
  cat("\n", " Variance component: sig2^2 ","\n")
  cat(round(sig2, 4),"(",round(stds2,4),")", "\n")
  cat("\n","gamma=",round(gamma,4), "\n")
  cat("\n")
  cPar<-as.numeric(iter==itmax)
  return(list(flag=flag, beta=beta, b=b, u=u, v=v, gamma=gamma, th1=th1,
              th2=th2, se.beta=se.beta, se.b=se.b, se.th1=se.th1, se.th2=se.th2,StopC=cPar))
}#end of function
#############################################################
# Some of the outputs obtained from the analysis of #
# the smoking cessation data #
#############################################################
# left <- Timept1
# right <- Timept2
# X <- as.matrix(cbind(SexF, SI/UC, F10Cigs, Duration))
# patient <- Zip
#######################################################
# Part 1. mixture(left, right, X, itmax=100):
# iter1.max= 27
# Fit of mixture cure model with random effects ignored
# Parameters in the logistic regression:
#   Estimate SE p-value
# Intercept 1.0573 0.7811 0.1759
# Sex 0.2463 0.3232 0.4460
# SI/UC -0.6071 0.3630 0.0945
# Cigarettes per day 0.0836 0.0162 0.0000
# Duration as smoker -0.1072 0.0233 0.0000
# Parameters in the survival model:
#   Estimate SE p-value
# Sex 0.4206 0.2713 0.1210
# SI/UC -0.0630 0.2889 0.8274
# Cigarettes per day -0.0654 0.0139 0.0000
# Duration as smoker 0.0721 0.0151 0.0000
#######################################################
# Part 2. mixture.effect(left, right, X, patient, itmax=100):
# iter2.max= 39
# Fit of mixture cure model with random effects
# Parameters in the logistic regression:
#   Estimate SE p-value
# Intercept 1.6817 0.8156 0.0392
# Sex 0.0921 0.3286 0.7793
# SI/UC -0.5154 0.3729 0.1668
# Cigarettes per day 0.0904 0.0169 0.0000
# Duration as smoker -0.1281 0.0245 0.0000
# Variance component: sig1^2
# 0.0666 ( 0.1700 )
# Parameters in the survival model:
#   Estimate SE p-value
# Sex 0.5329 0.2713 0.0500
# SI/UC -0.1741 0.2926 0.5517
# Cigarettes per day -0.0583 0.0137 0.0000
# Duration as smoker 0.0762 0.0151 0.0000
# Variance component: sig2^2
# 0.0135 ( 0.2107 ) 