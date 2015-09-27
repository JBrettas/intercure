# devtools::use_package("foreach")
# devtools::use_package("snow")
# devtools::use_package("survival")
# devtools::use_package("MASS")
# devtools::use_package("stats4")
# devtools::use_package("Matrix")


# Generates a vector of n latent variables ksi
gera_ksih_effect <- function(y_h, BASE, delta, cov_theta, cov_beta, grupo, theta, beta, gamma, NelsonAalen){
  m <- length(unique(grupo))
  A_h = B_h = rep(NA, m)
  for(i in 1:m) {
    A_h[i] <- exp(gamma) + sum(delta[(grupo == unique(grupo)[i])])
    B_h[i] <- exp(gamma) + sum((exp(theta %*% t(cbind(1,BASE[grupo == unique(grupo)[i], cov_theta]))) / 2) *
                                 (1 - 1 / (2 * exp(beta %*% t(BASE[grupo == unique(grupo)[i], cov_beta])) *
                                             (sapply(y_h[grupo == unique(grupo)[i]],NelsonAalen)) + 1)))
  }
  rgamma(m, A_h, B_h)
}

# Generates a vector of n latent variables k
gera_kh_effect <- function(y_h, ksi_h, BASE, delta, cov_theta, cov_beta, grupo, theta, beta, NelsonAalen){
  corresp <- data.frame(clu = unique(grupo), effect = ksi_h)
  ksi_ij <- plyr::join(data.frame(clu=grupo),corresp,by="clu")[,2]
  m <- length(unique(grupo))
  n <- nrow(BASE)
  intercepto <- 1
  xi_0 <- cbind(intercepto, BASE[,cov_theta])
  xi_1 <- BASE[, cov_beta]
  num <- exp(as.vector(theta %*%t (xi_0))) * ksi_ij
  den <- 2 + 4 * sapply(y_h, NelsonAalen) * exp(as.vector(beta %*% t(xi_1)))
  k_h <- rpois(length(num), num / den) + delta
  return(list(k = k_h, ksi = ksi_ij))
}

# Log-likelihood of gamma
log_vero_gamma <- function(gamma_par, ksi_h){
  logvero <- sum((dgamma(ksi_h, shape=exp(gamma_par), rate=exp(gamma_par), log=TRUE)))
  return(-logvero)
}


# Log-likelihood of theta given y_h, k_h, covariates and theta
log_vero_theta <- function(theta, cov_theta, y_h, k_h, BASE){
  total <- 0
  intercepto <- 1
  cov <- cbind(intercepto, BASE[,cov_theta])
  total <- sum((-exp(theta %*% t(cov)) / 2) + k_h * (theta %*% t(cov)))
  return(total)
}

# Log-likelihood of beta given y_h, u_h, covariates and beta
log_vero_beta <- function(beta, cov_beta, y_h, u_h, BASE){
  BASE_completa <- BASE
  BASE_completa$u <- u_h
  delta <- BASE_completa$delta
  factors <- BASE_completa$u * exp(beta %*% t(BASE_completa[,cov_beta]))
  total <- 0
  for(i in 1:length(y_h)) {
    xi_1 <- BASE_completa[i,cov_beta]
    soma_risco_i <- sum(factors[y_h >= y_h[i]])
    if (delta[i] == 1) parcela <- delta[i] * (sum(beta * xi_1) - log(soma_risco_i))
    else parcela <- 0
    total <- total + parcela
  }
  return(total)
}


# w may need more precision (rmpfr)

S_cl_i <- function(tempo, i, j, theta, beta, gamma, cov_theta, cov_beta, grupo, NelsonAalen, dataset) {
  w <- exp(gamma)
  base_cl <- dataset[grupo == i,]
  etas_i <- as.numeric(exp(theta %*% t(cbind(1, base_cl[,cov_theta]))))
  mus_i <- as.numeric(exp(beta %*% t(base_cl[,cov_beta])))
  if (j == 1) {
    S_1_cl_i <- (1 + (etas_i[1] / 2 * w) * (1 - 1 / (2 * mus_i[1] * NelsonAalen(tempo) + 1))) ^ (-w)
    return(S_1_cl_i)
  }
  num <- (etas_i[j] / 2) * (1 - (1 / (2 * mus_i[j] * NelsonAalen(tempo) + 1)))
  den <- w + sum((etas_i[1:(j-1)] / 2) * (1 - (1 / (2 * mus_i[1:(j-1)] * NelsonAalen(tempo) + 1))))
  S_j_cl_i <- (1 + num / den)^(-w - sum(base_cl$delta[1:(j - 1)]))
  return(S_j_cl_i)
}

F_y_cond_effect <- function(tempo, l, r, i, j, theta, beta, gamma, cov_theta, cov_beta, NelsonAalen, dataset){
  S_cl <- function(t_gen) S_cl_i(t_gen, i, j, theta, beta, gamma, cov_theta, cov_beta, grupo, NelsonAalen, dataset)
  NAalen_mod <- function(t_gen) {
    (NelsonAalen(l) * (r - t_gen) + NelsonAalen(r) * (t_gen - l)) / (r - l)
  }
  S_cl_mod <- function(t_gen) S_cl_i(t_gen, i, j, theta, beta, gamma, cov_theta, cov_beta, grupo, NAalen_mod, dataset)
  num <- S_cl_mod(tempo) - S_cl(r)
  den <- S_cl(l) - S_cl(r)
  return(as.numeric(1 - (num / den)))
}

inverse = function (f, lower = 0.1, upper = 100) {
  function (y) uniroot((function (x) f(x) - y), lower = lower, upper = upper)[1]
}

#Generates a vector of n observations using the previous function
gera_yh_effect <- function(L, R, delta, cov_theta, cov_beta, grupo, theta, beta, gamma, NelsonAalen, dataset){
  new_yh <- NA
  for(i in unique(grupo)) {
    times_cl_i <- rep(NA,nrow(dataset[grupo == i,]))
    delta_cl_i <- delta[grupo == i]
    l_cl_i <- L[grupo == i]
    r_cl_i <- R[grupo == i]
    for(j in c(1:length(times_cl_i))) {
      if(delta_cl_i[j] == 0) {
        times_cl_i[j] = l_cl_i[j]
      } else {
         if (l_cl_i[j] == r_cl_i[j]) {
          times_cl_i[j] <- l_cl_i[j]
        } else{
          #print(i)
          #print(j)
          U <- runif(1)
          F_y_cond <- function(x) F_y_cond_effect(x, l_cl_i[j], r_cl_i[j], i, j, theta, beta, gamma, cov_theta, cov_beta, NelsonAalen, dataset)
          inverse_F_y_cond_effect <- inverse(function(x) F_y_cond(x), l_cl_i[j], r_cl_i[j])
          times_cl_i[j] <- inverse_F_y_cond_effect(U)
          }
        }
    }
    new_yh <- c(new_yh, times_cl_i)
  }
  new_yh <- new_yh[-1]
  return(new_yh)
}


# dataset <- read.delim("~/intercure_pkg/smoking_data.txt")
# cov_theta=cov_beta=c("SexF","Duration","SI.UC","F10Cigs")
# dataset=arrange(dataset,Zip)
# L <- dataset$Timept1
# R <- dataset$Timept2
# R <- ifelse(is.na(R),Inf,R)
# grupo=as.character(dataset$Zip)
# delta<-ifelse(R==Inf,0,1)
# M=10
# b=0.001
# N_INT_MAX=100
# BASE=dataset
# NAME_DIF=""
# ncores=1
# inter_frailty_cl(dataset,L,R,delta,cov_theta,cov_beta,grupo,10,0.0001,par_cl = NULL, outputFiles = TRUE)


#ANDA (Asymptotic Normal Data Augmentation)


#' Cure frailty model for interval censored clustered data.
#'
#' \code{inter_frailty_cl} returns a list with the estimated parameters
#' \code{par} and their covariance matrix \code{mcov}. The list also contains
#' the cure rate covariance estimates \code{mcov.cura} for cure rate part only
#' and a dummy variable \code{StopC} assuming 0 if algorithm converged and 1 if
#' a stop criteria ended the process.
#'
#' @param data_set Dataset used to fit the model.
#' @param L Vector containing the last check times before event.
#' @param R Vector containing the first check times after event.
#' @param delta Flag vector indicating failure inside interval.
#' @param cov_theta String vector containing the column names to be used on the
#'   cure rate predictor.
#' @param cov_beta String vector containing the column names to be used on the
#'   predictor associated with the hazard function.
#' @param grp Vector containing cluster identifier (numeric or string).
#' @param M Number of replicates generated by each iteration on the ANDA
#'   (Asymptotic Normal Data Augmentation) algorithm.
#' @param b Numeric for tolerance of convergence.
#' @param N_INT_MAX Maximum number of algorithm's iterations without the burn
#'   in.
#' @param par_cl Registered SOCK cluster for parallel process. If NULL (default)
#'   the program loops are executed sequentially.
#' @param burn_in Number of burn in iterations.
#' @return The \code{inter_frailty} function returns an list containing the
#'   following outputs.
#' @return \code{par} estimates of theta and beta parameters.
#' @return \code{mcov} estimates for the covariance matrix of theta and beta
#'   parameters.
#' @return \code{mcov.cura} estimates for the covariance matrix associated with
#'   the cure rate part.
#' @return \code{StopC} stop criteria indicator assuming 1 when process is
#'   stopped for a non-convergence criteria. Assumes 0 when convergence is
#'   reached.
#' @examples
#' sample_set <- sim_frailty_data(100)
#' inter_frailty(sample_set, sample_set$L, sample_set$R, sample_set$delta, c("xi1","xi2"), c("xi1","xi2"), M = 50)
#' inter_frailty(sample_set, sample_set$L, sample_set$R, sample_set$delta, c("xi1"), c("xi2"), M = 10)
#' @export
#' @import foreach
inter_frailty_cl <- function(data_set, L, R, delta, cov_theta, cov_beta, grp, M, b=0.001, N_INT_MAX=100, par_cl = NULL, burn_in=50, outputFiles = FALSE) {
  data_set <- as.data.frame(data_set)

  if (outputFiles) {
    est_file_name <- paste("LAM_Estimates.txt", sep="")
    var_file_name <- paste("LAM_Variances.txt", sep="")
    fileConn <- file(est_file_name, "w")
    write(paste("PAR:\t",paste0(c("Intercept",cov_theta,cov_beta, "log(gamma)"), collapse="\t")), file=fileConn, append=T, sep="")
  }

  #Initial values for y
  y_nxM = k_nxM = u_nxM = ksi_nxM = matrix(NA, nrow = M, ncol = nrow(data_set))
  for(i in 1:M) y_nxM[i,] <- ifelse(delta == 1, (L + R) / 2, L)

  #Auxiliar matrix for different clusters
  m <- length(unique(grp))
  z <- matrix(0, ncol = m, nrow = nrow(data_set))
  for(j in 1:m) {
    z[,j] <- ifelse(grp == unique(grp)[j], 1, 0)
  }

  #Initial values for the parameters
  compr_theta <- 1 + length(cov_theta); compr_beta <- length(cov_beta)
  compr_alpha <- compr_theta + compr_beta + 1
  a_M <- matrix(NA, nrow = M, ncol = compr_alpha)
  rotulos <- c("intercepto", cov_theta, cov_beta, "log_w")
  colnames(a_M) <- rotulos
  theta_M <- matrix(a_M[,1:compr_theta], nrow=M)
  colnames(theta_M) <- rotulos[1:compr_theta]
  beta_M <- matrix(a_M[,(compr_theta + 1):compr_alpha], nrow=M)
  colnames(beta_M) <- rotulos[(compr_theta + 1):compr_alpha]
  alpha <- c(1:compr_alpha) * 0
  sigma_alpha <- b * diag(compr_alpha)
  beta <- alpha[(compr_theta + 1):(compr_alpha - 1)]
  gamma <- alpha[compr_alpha]

  #Initial values for latent vector u
  u <- delta

  #Initial Nelson-Aalen estimator
  Vetores_NAalen <- nelson_aalen_table(data_set, y_nxM[1,], delta, beta, cov_beta, u)
  NAalen_MEDIA <- stepfun(Vetores_NAalen$time, c(0,Vetores_NAalen$hazard))

  #Initializing convergence criteria and parameters
  conv <- FALSE
  n <- 0
  a_M_NEW <- a_M

  #Iterative process (with parallel computing)
  while(!conv | n <= burn_in) {
    if(!is.null(par_cl)){
      list_reg <- foreach(iterators::icount(M), .packages=c("MASS","MLEcens","Matrix","survival","stats4","plyr"), .export=c("surv_lam","inverse_lam_f","S_cl_i","inverse","F_y_cond_effect", "gera_yh_effect","gera_ksih_effect", "gera_kh_effect","log_vero_gamma", "gera_uh"), .inorder=F) %dopar% {
        a_M <- mvrnorm(n=1, alpha, sigma_alpha)
        theta_M <- a_M[1:compr_theta]
        beta_M <- a_M[(compr_theta + 1):(compr_alpha - 1)]
        gamma_M <- a_M[compr_alpha]
        y <- as.numeric(gera_yh_effect(L, R, delta, cov_theta, cov_beta, grp, theta_M, beta_M, gamma_M, NAalen_MEDIA, data_set))
        ksi <- gera_ksih_effect(y, data_set, delta, cov_theta, cov_beta, grp, theta_M, beta_M, gamma_M, NAalen_MEDIA)
        k_aux <- gera_kh_effect(y, ksi, data_set, delta, cov_theta, cov_beta, grp, theta_M, beta_M, NAalen_MEDIA)
        k <- k_aux$k
        ksi_geral <- k_aux$ksi
        u <- gera_uh(y, k, data_set, R, delta, cov_beta, beta_M, NAalen_MEDIA)

        #Gamma Regression for w
        fit_gamma <- stats4::mle(log_vero_gamma, start = list(gamma_par = as.numeric(gamma_M)), fixed = list(ksi_h = ksi))

        #Poisson Regression for Theta
        o_set <- k * 0 - log(2) + log(ksi_geral)
        expression_theta <- paste("data_set$", cov_theta[1:length(cov_theta)] ,sep = "", collapse="+")
        eq_theta <- paste("fit_theta <- stats::glm(k~",expression_theta," + offset(o_set),family = poisson)")
        eval(parse(text = eq_theta))

        #Cox Regression for Beta
        expression_beta <- paste("data_set$", cov_beta[1:length(cov_beta)] ,sep = "", collapse="+")
        eq_beta <- paste("fit_beta <- survival::coxph(Surv(y,delta)~", expression_beta, " + offset(ifelse(log(u) == -Inf, -200, log(u))),method = 'breslow')")
        eval(parse(text = eq_beta))

        #Outputs of Parallel Computing
        out <- list(fit_theta$coef, vcov(fit_theta), fit_beta$coef, vcov(fit_beta), coef(fit_gamma)[1], vcov(fit_gamma), y, u, k, ksi_geral)
        out
      }
    } else {
      list_reg <- foreach(iterators::icount(M), .packages=c("MASS","MLEcens","Matrix","survival","stats4","plyr"), .export=c("surv_lam","inverse_lam_f","S_cl_i","inverse","F_y_cond_effect", "gera_yh_effect","gera_ksih_effect", "gera_kh_effect","log_vero_gamma", "gera_uh"), .inorder=F) %do% {
        a_M <- mvrnorm(n=1, alpha, sigma_alpha)
        theta_M <- a_M[1:compr_theta]
        beta_M <- a_M[(compr_theta + 1):(compr_alpha - 1)]
        gamma_M <- a_M[compr_alpha]
        y <- as.numeric(gera_yh_effect(L, R, delta, cov_theta, cov_beta, grp, theta_M, beta_M, gamma_M, NAalen_MEDIA, data_set))
        ksi <- gera_ksih_effect(y, data_set, delta, cov_theta, cov_beta, grp, theta_M, beta_M, gamma_M, NAalen_MEDIA)
        k_aux <- gera_kh_effect(y, ksi, data_set, delta, cov_theta, cov_beta, grp, theta_M, beta_M, NAalen_MEDIA)
        k <- k_aux$k
        ksi_geral <- k_aux$ksi
        u <- gera_uh(y, k, data_set, R, delta, cov_beta, beta_M, NAalen_MEDIA)

        #Gamma Regression for w
        fit_gamma <- stats4::mle(log_vero_gamma, start = list(gamma_par = as.numeric(gamma_M)), fixed = list(ksi_h = ksi))

        #Poisson Regression for Theta
        o_set <- k * 0 - log(2) + log(ksi_geral)
        expression_theta <- paste("data_set$", cov_theta[1:length(cov_theta)] ,sep = "", collapse="+")
        eq_theta <- paste("fit_theta <- stats::glm(k~",expression_theta," + offset(o_set),family = poisson)")
        eval(parse(text = eq_theta))


        #Cox Regression for Beta
        expression_beta <- paste("data_set$", cov_beta[1:length(cov_beta)] ,sep = "", collapse="+")
        eq_beta <- paste("fit_beta <- survival::coxph(Surv(y,delta)~", expression_beta, " + offset(ifelse(log(u) == -Inf, -200, log(u))),method = 'breslow')")
        eval(parse(text = eq_beta))

        #Outputs of Parallel Computing
        out <- list(fit_theta$coef, vcov(fit_theta), fit_beta$coef, vcov(fit_beta), coef(fit_gamma)[1], vcov(fit_gamma), y, u, k, ksi_geral)
        out
      }
    }

    #Allocating M new and auxiliary parameter vectors
    SUM_VAR_GAMMA = SUM_VAR_THETA = SUM_VAR_BETA = 0
    for(h in 1:M) {
      SUM_VAR_THETA <- SUM_VAR_THETA + list_reg[[h]][[2]]; SUM_VAR_BETA =SUM_VAR_BETA + list_reg[[h]][[4]]
      a_M_NEW[h,1:compr_theta] <- list_reg[[h]][[1]]
      a_M_NEW[h,(compr_theta + 1):(compr_alpha - 1)] <- list_reg[[h]][[3]]
      a_M_NEW[h,compr_alpha] <- list_reg[[h]][[5]]
      SUM_VAR_GAMMA <- SUM_VAR_GAMMA + list_reg[[h]][[6]]
      y_nxM[h,] <- list_reg[[h]][[7]]
      u_nxM[h,] <- list_reg[[h]][[8]]
      k_nxM[h,] <- list_reg[[h]][[9]]
      ksi_nxM[h,] <- list_reg[[h]][[10]]
    }

    #Matrix of the M beta vectors
    beta_M <- matrix(a_M_NEW[,(compr_theta + 1):(compr_alpha - 1)], nrow=M)

    #Obtaining new Nelson-Aalen estimator for Cum. Hazard function
    if(!is.null(par_cl)){
      step_list <- foreach::foreach(h=1:M, .export="nelson_aalen_table", .inorder=F) %dopar% {
        V_NAalen <- nelson_aalen_table(data_set, y_nxM[h,], delta, beta_M[h,], cov_beta, u_nxM[h,])
        step_list <- stepfun(V_NAalen$time, c(0, V_NAalen$hazard))
        step_list
      }
    } else {
      step_list <- foreach::foreach(h=1:M, .export="nelson_aalen_table", .inorder=F) %do% {
        V_NAalen <- nelson_aalen_table(data_set, y_nxM[h,], delta, beta_M[h,], cov_beta, u_nxM[h,])
        step_list <- stepfun(V_NAalen$time, c(0, V_NAalen$hazard))
        step_list
      }
    }

    expression <- paste("step_list[[", 1:M,"]](x)", sep = "", collapse="+")
    eq4 <- paste("NAalen_MEDIA<-function(x) (",expression,") / ", M)
    eval(parse(text = eq4))

    #Creating new times/survival table and a more efficient estimator
    #V_NAalen<-aux_naalen(sort(y_nxM[,delta==1]),NAalen_MEDIA,cl)
    V_NAalen <- aux_naalen(sort(y_nxM[,delta==1]), NAalen_MEDIA, par_cl)
    NAalen_MEDIAnew <- stepfun(V_NAalen$time, c(0, V_NAalen$hazard))

    #Calculating the new covariance matrix
    SUM_VAR <- as.matrix(bdiag(list(SUM_VAR_THETA, SUM_VAR_BETA, SUM_VAR_GAMMA)))
    VAR <- var_matrix(SUM_VAR, a_M_NEW)
    sigma_alpha <- VAR

    #New vector of estimates
    alpha_new <- Matrix::colMeans(a_M_NEW)
    cat("")
    cat("Alpha NEW:", alpha_new)


    #Checking convergence
    conv <- convergence_lam(alpha_new, alpha)

    #Setting new alpha as old one for iteractive process
    alpha <- alpha_new

    #Writing alpha values
    if (outputFiles) {
      write(paste("IT",n + 1,":\t",paste0(alpha, collapse="\t")), file=fileConn, append=T, sep="")
      write.table(VAR, file=var_file_name, row.names=FALSE, col.names=FALSE)
    }

    #Setting new baseline cum. hazard estimator as old one for iteractive process
    NAalen_MEDIA <- NAalen_MEDIAnew

    #Updating the iteration counter
    n <- n + 1

    #Checking if iteration counter reached N_INT_MAX
    if(n == (N_INT_MAX + burn_in)) {
      if (outputFiles) {
        write("\nWarning: Iteration Number achieved but convergence criteria not met.",file=fileConn,append=T, sep=" ")
        close(fileConn)
      }
      cat("\nWarning: Convergence criteria not met. Estimates given for N_INT_MAX=", N_INT_MAX)
      break
    }
  }
  cPar <- as.numeric(n == (N_INT_MAX + burn_in))
  alphaList <- list(par = alpha, mcov = VAR, mcov.cura = VAR[1:(1 + length(cov_theta)), 1:(1 + length(cov_theta))], StopC = cPar)
  if (outputFiles) close(fileConn)
  return(alphaList)
}
