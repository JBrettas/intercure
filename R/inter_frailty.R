library(MASS)
library(Matrix)
library(survival)
library(dplyr)
library(compiler)
library(Rcpp)
enableJIT(3)
#sourceCpp("survCpp.cpp")


# Simulate data using fragility model (Lam et al 2007)
sim_dados_LAM<-function(N){
  theta<-c(-1,1,0); beta<-c(0,0.5)
  u<-runif(N); a<-runif(N)
  C<-cbind(5,a*15)
  C<-C[,1]*(C[,1]<=C[,2])+C[,2]*(C[,1]>C[,2])
  intercept<-1
  xi1<-rbinom(N,1,0.5)
  xi2<-rnorm(N)
  cov_theta<-data.frame(intercept,xi1,xi2)
  cov_beta<-data.frame(xi1,xi2)
  eta<-exp(as.vector(theta%*%t(cov_theta)))
  K_vector<-rpois(N,eta/2)
  U_vector<-K_vector*NA
  for(i in 1:length(K_vector)){
    if(K_vector[i]==0) U_vector[i]=0
    else{
      U_vector[i]=0
      for(j in 1:K_vector[i]) U_vector[i]=U_vector[i]+rchisq(1,2,ncp=0)
    }
  }
  #U_vector<-rchisq(N,0,ncp=eta)
  beta_x<-as.vector(beta%*%t(cov_beta))
  exp_pred_beta<-exp(beta_x)
  num=-2*log(1-u)
  den=U_vector*exp_pred_beta
  tempos=ifelse(U_vector!=0,sqrt(num/den),Inf)
  Z=ifelse(tempos<C,tempos,C)
  delta=ifelse(tempos<C,1,0)
  L=R=Z*NA
  for(i in 1:N){
    if(delta[i]==0){
      L[i]=Z[i]
      R[i]=Inf
    }
    else{
      L[i]=0
      add<-runif(1,0.1,0.5)
      R[i]=add
      check=(L[i]<=Z[i] & Z[i]<R[i])
      while(!check){
        L[i]=L[i]+add
        add<-runif(1,0.1,0.5)
        R[i]=R[i]+add
        check=(L[i]<=Z[i] & Z[i]<R[i])
      }
    }  
  }
  dados=data.frame(Z,L,R,delta,xi1,xi2)
  return(dados)
}

# Fast search for interval containing x
findInterval2 <- function(x,v) {
  n = length(v)
  if (x<v[1])
    return (0)
  if (x>=v[n])
    return (n)
  i=1
  k=n
  while({j = (k-i) %/% 2 + i; !(v[j] <= x && x < v[j+1])}) {
    if (x < v[j])
      k = j
    else
      i = j+1
  }
  return (j)
}

#Creates a table with time and respective Nelson-Aalen estimates for Cum. Hazard Function
#including latent variables u_i
Nelson_Aalen_Table <- function(dataSet,eventTimes,delta,beta,covariates,u){
  factors <- u*exp(beta%*%t(dataSet[,covariates]))
  eventTimes_relev <- unique(sort(eventTimes[delta==1]))
  parc_f <- function(t) (sum(eventTimes[delta==1]==t)/sum(factors[eventTimes>=t]))
  parcels <- sapply(eventTimes_relev  , parc_f  )
  cum_haz <- cumsum(parcels)
  c_haz_f <- data.frame(eventTimes_relev,cum_haz)
  colnames(c_haz_f) <- c("time","hazard")
  return(c_haz_f)
}

# Array with piecewise Nelson-Aalen estimator
Aux_NAalen<-function(Tempos,NelsonAalen_Original){
  z<-unique(sort(Tempos))
  NAalen_z<-sapply(z,NelsonAalen_Original)
  f_NAalen=data.frame(z,NAalen_z)
  colnames(f_NAalen)=c("time","hazard")
  return(f_NAalen)
}

# Survival Function given x(0) , x(1), Theta, Beta
# Note: Precision problems may occur on the linear predictors 
S<-function(t,xi_0,xi_1,theta,beta,NelsonAalen){
  exp(-(exp(theta%*%(xi_0))/2)*((1-(1/(1+2*NelsonAalen(t)*exp(beta%*%(xi_1)))))))
}

# Inverse F function conditioned on Li and Ri (Method 2 of Lam et al 2007)
inv_F_cond_LR2<-function(w,Li,Ri,xi_0,xi_1,theta,beta,NAalen_Original){
  S_L<-S(Li,xi_0,xi_1,theta,beta,NAalen_Original)
  S_R<-S(Ri,xi_0,xi_1,theta,beta,NAalen_Original)
  if(S_L==S_R) {
    cat("Warning! Precision problem on linear predictor. (S_L=S_R) ")
    stop
  }
  k=(1-w)*S_L + w*S_R
  k_linha=-2*log(k)/((exp(theta%*%(xi_0))+2*log(k))*2*exp(beta%*%(xi_1)))
  if(k_linha==0){
    cat("Warning. Division by 0 due to precision problems for w=", w)
    stop
  }
  nAalen_Li <- NAalen_Original(Li)
  nAalen_Ri <- NAalen_Original(Ri)
  a=(nAalen_Ri*Li-nAalen_Li*Ri)/(nAalen_Ri-nAalen_Li)
  b=(Ri-Li)/(nAalen_Ri-nAalen_Li)
  a+k_linha*b
}

# Generates n observations y using the inverse transformation
gera_yh<-function(BASE,L,R,delta,cov_theta,cov_beta,theta,beta,NelsonAalen_Original){
  tam <- length(L)
  new_y=rep(NA,tam)
  intercepto=1
  tic<-proc.time()
  X_theta <- tbl_df(BASE[,cov_theta])
  X_beta <- tbl_df(BASE[,cov_beta])
  Uh <- runif(tam)
  for(i in 1:tam){
    if((delta[i]==0)|(L[i]==R[i])){new_y[i]=L[i]}
      else{
        xi_0=as.numeric(cbind(intercepto,X_theta[i,]))
        xi_1=as.numeric(X_beta[i,])
        new_y[i]=inv_F_cond_LR2(Uh[i],L[i],R[i],xi_0,xi_1,theta,beta,NelsonAalen_Original)
      }
    }
  print(proc.time()-tic)
  return(new_y)
}  

# Generates latent variables k for n observations
gera_kh<-function(y_h,BASE,delta,cov_theta,cov_beta,theta,beta,NelsonAalen){
  k_h<-y_h*NA
  intercepto=1
  xi_0=cbind(intercepto,BASE[,cov_theta])
  xi_1=BASE[,cov_beta]
  num<-exp(as.vector(theta%*%t(xi_0)))
  den<-2+4*sapply(y_h,NelsonAalen)*exp(as.vector(beta%*%t(xi_1)))
  k_h<-rpois(length(num),num/den)+delta
  return(k_h)
}

# Generates latent variables u for n observations
gera_uh<-function(y_h,k_h,BASE,R,delta,cov_beta,beta,NelsonAalen){
  u_h<-y_h*NA
  r_estrela<-max(R[delta==1])
  xi_1=BASE[,cov_beta]
  alpha_gamma<-k_h+delta
  beta_gamma<-1/(0.5+sapply(y_h,NelsonAalen)*exp(as.vector(beta%*%t(xi_1))))
  COND<-(k_h==0 | y_h>r_estrela)
  u_h<-ifelse(COND,0,rgamma(length(alpha_gamma),alpha_gamma,scale=beta_gamma))
  return(u_h)
}

# Returns covariance matrix from (5) of (Lam et al 2007)
var_matrix<-function(SUM_VAR,alpha_matrix){
  M<-nrow(alpha_matrix)
  alpha_j<-colMeans(alpha_matrix)
  matriz_soma_2=diag(ncol(alpha_matrix))*0
  for(h in 1:M){
    matriz_soma_2=matriz_soma_2+((alpha_matrix[h,]-alpha_j)%*%t(alpha_matrix[h,]-alpha_j))/(M-1)
  }
  sigma_alpha=(SUM_VAR/M)+(1+1/M)*matriz_soma_2
  return(sigma_alpha)
}

#Convergence criteria returning TRUE of FALSE
convergencia<-function(alpha_new,alpha_old,tol=0.001){
  conv=FALSE
  new_val<-c(alpha_new)
  old_val<-c(alpha_old)
  max_error<-max(abs(new_val-old_val))
  if(max_error<tol) conv=TRUE
  return(conv)
}

#ANDA (Asymptotic Normal Data Augmentation)
ANDA<-function(BASE,L,R,delta,cov_theta,cov_beta,M,b=0.001,N_INT_MAX=100,NAME_DIF="",ncores=1,burn_in=30){
  # File outputs
  est_file_name<-paste("LAM_Estimativas",NAME_DIF,".txt",     sep="")
  var_file_name<-paste("LAM_Variancias",NAME_DIF,".txt", sep="")
  fileConn<-file(est_file_name,"w")
  
  # Faster dataset
  BASE<-tbl_df(BASE)
  
  # Arrays for simulated latent variables
  y_nxM=k_nxM=u_nxM=matrix(NA,nrow=M,ncol=nrow(BASE))
  for(i in 1:M) y_nxM[i,]=ifelse(delta==1,(L+R)/2,L)
  
  # Matrix to allocate M parameter vectors
  compr_theta<-1+length(cov_theta); compr_beta<-length(cov_beta)
  compr_alpha<-compr_theta+compr_beta
  a_M<-matrix(NA,nrow=M,ncol=compr_alpha)
  rotulos<-c("intercepto",cov_theta,cov_beta)
  colnames(a_M)<-rotulos
  theta_M<-matrix(a_M[,1:compr_theta],nrow=M) 
  colnames(theta_M)<-rotulos[1:compr_theta]
  beta_M<-matrix(a_M[,(compr_theta+1):compr_alpha],nrow=M) 
  colnames(beta_M)<-rotulos[(compr_theta+1):compr_alpha]
  
  # Initial set of parameters
  alpha<-c(1:compr_alpha)*0 
  sigma_alpha<-b*diag(compr_alpha); beta<-alpha[(compr_theta+1):compr_alpha]
  
  # Initial u vector = delta vector
  u<-delta 
  
  # Initial cumulative hazard function estimate
  Vetores_NAalen<-Nelson_Aalen_Table(BASE,y_nxM[1,],delta,beta,cov_beta,u)
  NAalen_MEDIA<-stepfun(Vetores_NAalen$time,c(0,Vetores_NAalen$hazard))
  
  # Initializing convergence criteria and parameters
  conv=FALSE
  n=0
  a_M_NEW<-a_M
  
  # Iterates until convergence or user defined limit for iterations is reached
  while(!conv | n<=burn_in){
    # Set initial sum of variances as 0
    SUM_VAR_THETA=SUM_VAR_BETA=0
    
    # Loop for each h replicate
    tempoITER<-system.time({
      for(h in 1:M){
        # Generating h parameter vector
        a_M[h,]<-mvrnorm(n=1, alpha, sigma_alpha)
        theta_M[h,]<-a_M[h,1:compr_theta]; beta_M[h,]<-a_M[h,(compr_theta+1):compr_alpha]
        
        # Generating h observation vector
        y_nxM[h,]=gera_yh(BASE,L,R,delta,cov_theta,cov_beta,as.numeric(theta_M[h,]),as.numeric(beta_M[h,]),NAalen_MEDIA)
        
        # Generating h-th k vector
        k_nxM[h,]=gera_kh(y_nxM[h,],BASE,delta,cov_theta,cov_beta,theta_M[h,],beta_M[h,],NAalen_MEDIA)
        
        # Generating h-th u vector
        u_nxM[h,]=gera_uh(y_nxM[h,],k_nxM[h,],BASE,R,delta,cov_beta,beta_M[h,],NAalen_MEDIA)
        
        # Defining offset to construct paper's likelihood
        o_set<-k_nxM[h,]*0-log(2)
        
        # GLM model for cure rate
        expression_theta=paste("BASE$", cov_theta[1:length(cov_theta)] , sep = "",collapse="+")
        eq_theta<-paste("fit_theta<-glm(k_nxM[h,]~",expression_theta,"+offset(o_set),family=poisson)")
        eval(parse(text=eq_theta))
        a_M_NEW[h,1:compr_theta]=fit_theta$coef
        SUM_VAR_THETA=SUM_VAR_THETA+vcov(fit_theta)
      
        # Cox proportional hazards model for time
        expression_beta=paste("BASE$", cov_beta[1:length(cov_beta)] , sep = "",collapse="+")
        eq_beta<-paste("fit_beta<-coxph(Surv(y_nxM[h,],delta)~",expression_beta,"+offset(ifelse(log(u_nxM[h,])==-Inf,-200,log(u_nxM[h,]))),method='breslow')")
        eval(parse(text=eq_beta))
        a_M_NEW[h,(compr_theta+1):compr_alpha]=fit_beta$coef
        SUM_VAR_BETA=SUM_VAR_BETA+vcov(fit_beta)
      }
      
      # Obtaining Nelson-Aalen estimate for each h replica
      lhs2  <- paste("V",1:M,"_NAalen",     sep="")
      rhs2  <- paste("Nelson_Aalen_Table(BASE,y_nxM[",1:M,",],delta,beta_M[",1:M,",],cov_beta,u_nxM[",1:M,",])", sep="")
      eq2   <- paste(paste(lhs2, rhs2, sep="<-"), collapse=";")
      eval(parse(text=eq2))
      lhs3  <- paste("NAalen",    1:M,     sep="")
      rhs3  <- paste("stepfun(V",1:M,"_NAalen$time, c(0,V",1:M,"_NAalen$hazard))", sep="")
      eq3   <- paste(paste(lhs3, rhs3, sep="<-"), collapse=";")
      eval(parse(text=eq3))
      
      # Obtaining mean of functions defined above
      expression=paste("NAalen", 1:M,"(x)" , sep = "",collapse="+")
      eq4<-paste("NAalen_MEDIA<-function(x) (",expression,")/",M)
      eval(parse(text=eq4))
      V_NAalen<-Aux_NAalen(sort(y_nxM[,delta==1]),NAalen_MEDIA)
      NAalen_MEDIAnew<-stepfun(V_NAalen$time,c(0,V_NAalen$hazard))
      
      # Updating covariance matrix
      SUM_VAR<-as.matrix(bdiag(list(SUM_VAR_THETA,SUM_VAR_BETA)))
      VAR<-var_matrix(SUM_VAR,a_M_NEW)
      sigma_alpha<-VAR
      
      # New vector of parameters
      alpha_new=colMeans(a_M_NEW)
    })

    # Evaluating convergence
    conv=convergencia(alpha_new,alpha)
    
    # Updating vector of parameters
    alpha<-alpha_new
    
    # Inserting new estimates into outputs
    write(as.vector(alpha),file=fileConn,append=T,sep=" ")
    write.table(VAR, file=var_file_name, row.names=FALSE, col.names=FALSE)
    
    # Updating actual cumulative hazard estimate
    NAalen_MEDIA<-NAalen_MEDIAnew
    
    # Show and update iterate
    cat("Finished it ",n)
    n=n+1
    
    # Checking if iteration counter reached N_INT_MAX
    if(n==(N_INT_MAX+burn_in)) {
      write("\n Warning: Iteration Number achieved but convergence criteria not met.",file=fileConn,append=T,sep=" ")
      close(fileConn)
      cat("\n Warning: Convergence criteria not met. Estimates given for N_INT_MAX=",N_INT_MAX)
      break
    }
  }
  
  # Flag indicating iteration limit stop
  cPar<-as.numeric(n==(N_INT_MAX+burn_in))
  
  # List of outputs
  alphaList<-list(par=alpha,mcov=VAR,mcov.cura=VAR[1:(1+length(cov_theta)),1:(1+length(cov_theta))],StopC=cPar)
  print(alphaList)
  
  return(alphaList)
}