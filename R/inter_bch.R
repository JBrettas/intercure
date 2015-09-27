library(ICE)
library(MASS)
library(interval)
library(Matrix)
library(compiler)
enableJIT(3)

# Create global variables s, r, m and sm1 for main function
create_sr <- function(L,R){
  eq_int <- inmost(data.frame(L,R),0)
  s <- eq_int$q
  r <- eq_int$p
  if(!(max(L[R==Inf])>max(r[r!=Inf]))) {
    cat("Warning: There is no event-free case with L_i > r_m")
  }
  m=length(r[r!=Inf])
  if(as.logical(sum(!s%in%L) | sum(!r%in%R))) cat("Alerta: algum membro do conjunto s ou r nao pertence a L ou R respectivamente (rever banco de dados)")
    s<<-s
    r<<-r
    m<<-m
    s0<-(-Inf)
    sm1<-c(s0,s)
    sm1<<-sm1[-length(sm1)]
}

# Alternative function using Aintmap below
#
# Create global variables s, r, m and sm1 for main function
# create_sr<-function(L,R){
#   eq_int<-Aintmap(L,R,Lin=T,Rin=T)$intmap
#   s<-eq_int[1,]
#   r<-eq_int[2,]
#   if(!(max(L[R==Inf])>max(r[r!=Inf]))) {
#     cat("Warning: Dataset not compatible with the model: there is no case with L_i > r_m")
#     stop
#   }
#   m=length(r[r!=Inf])
#   if(as.logical(sum(!s%in%L) | sum(!r%in%R))) cat("Alerta: algum membro do conjunto s ou r nao pertence a L ou R respectivamente (rever banco de dados)")
#   s<<-s
#   r<<-r
#   m<<-m
#   s0<-(-Inf)
#   sm1<-c(s0,s)
#   sm1<<-sm1[-length(sm1)]
# }

# Computes the X matrix for the k-th iteraction, proposed by Shen and Hao
compute_Xk<-function(theta,vector_prob,dados,L,R,covariates,F_hat){

  # Matrix with Intercepts
  Z=(data.frame(1,dados[,covariates]))
  colnames(Z)=c("intercept",covariates)

  # Evaluating gamma defined by Shen
  Xk=exp(-exp(theta%*%t(Z))%x%sapply(s[1:m],F_hat))
  Xk=Xk*(1-exp(-exp(theta%*%t(Z))%x%vector_prob))
  Xk<-rbind(Xk,exp(-exp(theta%*%t(Z))))

  # Called X_k on Shen (a_ij on Turnbull's article)
  for(i in 1:nrow(dados)){
    Xk[,i]<-as.numeric(L[i]<=s & r <=R[i])*Xk[,i]
  }
  # Defining  Xk
  Xk<-Matrix((1/colSums(Xk))*t(Xk))
  return(Xk)
}

# Computes auxiliar matrix C given Xk
compute_C<-function(Xk){
  aux_C<-diag(m)*NA
  for(j in 1:m){
    aux_C[,j]<-(s[j]<=sm1[-length(sm1)])
  }
  C<-Matrix((Xk[,1:m])%*%aux_C)
  return(C)
}

# Expected log-likelihood
log_lik<-function(theta=(1:(1+length(covariates)))*0,vector_prob=rep.int(1/m,m),dados,Xk,C,covariates,F_hat){
  Z=data.frame(1,dados[,covariates])
  colnames(Z)=c("intercept",covariates)
  l<-as.numeric(exp(theta%*%t(Z)))
  soma_1<-as.numeric(-t(l)%*%(C[][,1:m]%*%vector_prob[1:m]))
  soma_2<-sum((diag(log(1-exp(-vector_prob%*%t(l)))%*%(Xk[,1:m]))))
  soma_3<-as.numeric(-Xk[,(m+1)]%*%l)
  soma_lvero=as.numeric(soma_1+soma_2+soma_3)
  return(soma_lvero)
}


#####################################################################################################################
#  WARNING: IF THETA BIG ENOUGH -> NANs ARE OBTAINED WITH THE R PRECISION IN Xk[,1:m]/(exp(t(b_kp1%x%p))-1))        #
#####################################################################################################################

# First order derivative for p
d_fp<-function(p,theta,Xk,C,Z){
  b_kp1<-exp(theta%*%t(Z))
  dfp<-b_kp1%*%((Xk[][,1:m]/(exp(t(b_kp1%x%p))-1))-C[][,1:m])
  return(dfp)
}

# Second order derivative for p
dd_fp<-function(p,theta,Xk,C,Z){
  b_kp1<-exp(theta%*%t(Z))
  ddfp=p*0
  for(i in 1:nrow(Xk)){
    ddfp=ddfp-(Xk[][i,1:m]*exp(2*theta%*%t(Z[i,]))*exp(b_kp1[i]*p))/((exp(b_kp1[i]*p)-1)*(exp(b_kp1[i]*p)-1))
  }
  return(ddfp)
}

# Updates v
v_plus<-function(p,tau=(1/p),dfp,ddfp,Sigma){
  return(sum((1/(Sigma/mean(tau*p))+p*dfp)/(tau-p*ddfp))/sum((p/(tau-p*ddfp))))
}

# Delta p
delta_p<-function(p,tau=(1/p),v_plus,dfp,ddfp,Sigma){
  return((p*(dfp-v_plus)+(1/(Sigma/mean(tau*p))))/(tau-p*ddfp))
}

# Euclidean Distance
norm_vec <- function(x) sqrt(sum(x*x))

# Proposed function to evaluate p's convergence
Fn<-function(p,tau,v,theta,Xk,C,eta){
  dfp<-d_fp(p,theta,Xk[],C[],Z)
  v1=dfp+tau-v
  v2=(diag(tau)%*%p)-(1/eta)
  v3=1-sum(p)
  return(c(v1,v2,v3))
}

# Backtracking for p maximizing
backtracking<-function(p,delta_p,tau,delta_tau,v,delta_v,theta,Xk,C,Sigma){
  # Defining eta
  eta<-Sigma/mean(tau*p)
  pi=0.01
  rho<-0.5
  psi_max<-min(1,min(-(tau/delta_tau)[delta_tau<0]))
  psi=psi_max*0.99
  while(min(p+psi*delta_p)<=0) psi=psi*rho
  while((norm_vec(Fn(as.vector(p+psi*delta_p),as.vector(tau+psi*delta_tau),as.vector(v+psi*delta_v),theta,as.matrix(Xk[]),as.matrix(C[]),eta)))>((1-pi*psi)*norm_vec(Fn(p,tau,v,theta,as.matrix(Xk[]),as.matrix(C[]),eta)))) {
    psi=psi*rho
  }
  return(psi)
}

# Estimates p maximizing the expected log-likelihood
max_p<-function(p,theta,Xk,eps2=0.001,MAXITER=500,Sigma){
  # Stop criteria
  C<-compute_C(Xk)
  CRIT1=FALSE
  CRIT2=FALSE
  it=1

  # First and Second order derivatives
  dfp<-d_fp(p,theta,Xk,C,Z)
  ddfp<-dd_fp(p,theta,Xk,C,Z)

  # Defining initial tau
  ini_tau<-1/p
  tau<-ini_tau

  # Initial v
  v<-0.1

  while(( !CRIT1 | !CRIT2 )&it<=MAXITER){
    # Defining eta
    eta<-Sigma/mean(tau*p)

    # New v (used on delta_v)
    new_v<-v_plus(p,tau,dfp,ddfp,Sigma)

    # Delta_p
    delp<-delta_p(p,tau,new_v,dfp,ddfp,Sigma)

    # tau_p
    tau_plus<-new_v-ddfp*delp-dfp

    # Delta_tau, Delta_v
    deltau<-tau_plus-tau
    delv<-new_v-v

    # Backtracking
    psi<-backtracking(p,delp,tau,deltau,v,delv,theta,Xk,C,Sigma)
    p=as.vector(p+psi*delp)
    tau=as.vector(tau+psi*deltau)
    v=as.vector(v+psi*delv)
    dfp<-d_fp(p,theta,Xk,C,Z)
    ddfp<-dd_fp(p,theta,Xk,C,Z)

    # Updating convergence parameters
    CRIT1<-(sum(p*tau)<eps2)
    CRIT2<-(sqrt(sum((dfp+tau-v)*(dfp+tau-v)))<eps2)
    #cat("CRIT2:",sqrt(sum((dfp+tau-v)*(dfp+tau-v))))  studying convergence
    it=it+1
  }
  return(p)
}

# Simulates BCH with F= 1-exp(-t) given theta=c(alpha,beta)
Sim_BCH<-function(N,alpha,beta){
  z<-runif(N,-1,1)
  pi_Z<-rbinom(N,1,exp(-exp(alpha+beta*z)))
  U<-runif(N)
  V1=(1-exp(-exp(alpha+beta*z)))*(1-U)+exp(-exp(alpha+beta*z))
  V2=-log(1+(log(V1)/exp(alpha+beta*z)))
  Times=ifelse(pi_Z==1,Inf,V2)
  L=R=c(1:N)*0
  for(i in 1:N){
    if(pi_Z[i]==1){
      L[i]=3
      R[i]=Inf
    }
    else{
      R[i]=rexp(1,5)
      while(Times[i]<L[i] | Times[i]>R[i]){
        L[i]=R[i]
        R[i]=R[i]+rexp(1,5)
      }
      if(R[i]>=3){
        R[i]=Inf
      }
      if(L[i]>=3){
        L[i]=3
      }

    }
  }
  dados<-data.frame(Times,L,R,z)
  colnames(dados)<-c("Times","L","R","Cov")
  return(dados)
}

# Shen/Hao estimator
sh_reg<-function(dados,L,R,covariates,Sigma=20,crit_theta=0.001,crit_p=0.005,n_int=100,OUTPUT_FILE="Shen_Estimates",OUTPUT_VAR_FILE="Shen_Variances"){
  # Initializing output files
  est_file_name<-paste(OUTPUT_FILE,".txt",     sep="")
  var_file_name<-paste(OUTPUT_VAR_FILE,".txt", sep="")
  fileConn<-file(est_file_name,"w")

  # Specifying the covariates
  dados<-as.data.frame(dados)
  Z=data.frame(1,dados[,covariates])
  colnames(Z)=c("intercept",covariates)
  Z<<-Z

  # Defining global variables s, r, m and sm1
  create_sr(L,R)

  # Initial vector of probabilities (Turnbull "jumps")
  prob<-1/m #Uniform on first iteraction
  p<-rep.int(prob,m)

  # F initial estimator
  F_hat_aux<-data.frame(r[1:m],cumsum(p))
  colnames(F_hat_aux)<-c("time","cum")
  F_hat<-stepfun(F_hat_aux$time,c(0,F_hat_aux$cum))

  # Initial theta
  theta_k=c(1:ncol(Z))*0

  # Convergence flags
  CONV=FALSE
  CONV2=FALSE

  # Iteration counter
  it=1

  # ECM Loop
  while(!CONV | !CONV2){
    # Timing loop
    tic <- proc.time()

    # Computing matrix X and C
    cat("\nIT #",it)
    cat("\nComputing Xk")
    Xk<-compute_Xk(theta_k,p,dados,L,R,covariates,F_hat)
    cat("\nComputing C")
    C<-compute_C(Xk)

    # Obtaining theta and it's variance with MLE
    llk<-function(theta) log_lik(theta,p,dados,Xk,C,covariates,F_hat)
    fit_theta<-optim(theta_k,llk,method = "BFGS",control=list(fnscale=-1),hessian=T)
    theta_knew<-fit_theta$par
    theta_var<-solve(-fit_theta$hessian)
    write.table(theta_var, file=var_file_name, row.names=FALSE, col.names=FALSE)

    # Checking theta convergence
    CONV=(max(abs(theta_knew-theta_k))<crit_theta)
    cat("\nTheta Max Difference:",max(abs(theta_knew-theta_k)))
    cat("\nTHETA:", theta_knew)

    # Updating Xk for new theta
    Xk<-compute_Xk(theta_knew,p,dados,L,R,covariates,F_hat)

    # Computing p vector for new Xk and theta
    novo_p<-max_p(p,theta_knew,Xk,Sigma=Sigma)

    # Checking p convergence
    CONV2=(max(abs(novo_p-p))<crit_p)
    cat("\np Max Difference:",max(abs(novo_p-p)))
    cat("\n")

    # Updating parameters
    p=novo_p
    F_hat_aux<-data.frame(r[1:m],cumsum(p))
    colnames(F_hat_aux)<-c("time","cum")
    F_hat<-stepfun(F_hat_aux$time,c(0,F_hat_aux$cum))
    theta_k=theta_knew

    # Writing new estimates on file
    write(as.vector(theta_k),file=fileConn,append=T,sep=" ")

    # Shows iteration time
    cat("\n it time:",(proc.time()-tic))
    it=it+1

    # Check if reached iteration limit defined by user
    if(it>=n_int) {
      write("WARNING: CONVERGENCE NOT REACHED",file=fileConn,append=T,sep=" ")
      cat("\nWARNING: CONVERGENCE NOT REACHED FOR #IT=",n_int)
      break
    }
  }
  close(fileConn)

  # Check if reached the it limit
  cPar<-as.numeric(it>=n_int)

  # Remove global aux variables
  rm(Z,m,r,s,sm1, pos = ".GlobalEnv")

  # Outputs an list with useful metrics
  return(list("par"=theta_k,"p"=p,"mcov"=theta_var,"StopC"=cPar))
}


