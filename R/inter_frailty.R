# devtools::use_package("foreach")
# devtools::use_package("snow")
# devtools::use_package("survival")
# devtools::use_package("MASS")
# devtools::use_package("Matrix")
# devtools::use_testthat()

# Fast search for interval containing x
findInterval2 <- function(x,v) {
  n = length(v)
  if (x < v[1])
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
nelson_aalen_table <- function(dataSet, eventTimes, delta, beta, covariates, u) {
  factors <- u * exp(beta %*% t(dataSet[,covariates]))
  event_times_relev <- unique(sort(eventTimes[delta == 1]))
  parc_f <- function(t) (sum(eventTimes[delta == 1] == t) / sum(factors[eventTimes >= t]))
  parcels <- sapply(event_times_relev  , parc_f )
  cum_haz <- cumsum(parcels)
  c_haz_f <- data.frame(event_times_relev, cum_haz)
  colnames(c_haz_f) <- c("time", "hazard")
  return(c_haz_f)
}

#Auxiliar function to Nelson-Aalen estimates table
aux_naalen <- function(tempos, naalen_f, par_cl = NULL) {
  z <- unique(sort(tempos))
  if (!is.null(par_cl)) {
    vector_naalen <- snow::parSapply(par_cl, z, naalen_f)
  } else {
    vector_naalen <- sapply(z, naalen_f)
  }
  pieceAalen <- data.frame(z, vector_naalen)
  colnames(pieceAalen) <- c("time", "hazard")
  return(pieceAalen)
}


# Survival Function given x(0) , x(1), Theta, Beta
# Note: Precision problems may occur on the linear predictors
surv_lam <- function(t, xi_0, xi_1, theta, beta, nelson_aalen_function){
  exp(-(exp(theta %*% (xi_0)) / 2) * ((1 - (1 / (1 + 2 * nelson_aalen_function(t)*exp(beta %*% (xi_1)))))))
}

# Inverse F function conditioned on Li and Ri (Method 2 of Lam et al 2007)
inverse_lam_f <- function(w, Li, Ri, xi_0, xi_1, theta, beta, naalen_original){
  S_L <- surv_lam(Li, xi_0, xi_1, theta, beta, naalen_original)
  S_R <- surv_lam(Ri, xi_0, xi_1, theta, beta, naalen_original)
  if (S_L == S_R) {
    cat(" Warning! Precision problem on linear predictor. (S_L = S_R) ")
    stop
  }
  k <- (1 - w) * S_L + w * S_R
  k_linha <- -2 * log(k) / ((exp(theta %*% (xi_0)) + 2 * log(k)) * 2 * exp(beta %*% (xi_1)))
  if (k_linha == 0) {
    cat(" Warning. Division by 0 due to precision problems for w= ", w)
    stop
  }
  nAalen_Li <- naalen_original(Li)
  nAalen_Ri <- naalen_original(Ri)
  a <- (nAalen_Ri * Li - nAalen_Li * Ri) / (nAalen_Ri - nAalen_Li)
  b <- (Ri - Li) / (nAalen_Ri - nAalen_Li)
  a + k_linha * b
}

# Generates n observations y using the inverse transformation
gera_yh <- function(data_set, L, R, delta, cov_theta, cov_beta, theta, beta, naalen_original) {
  tam <- length(L)
  new_y <- rep(NA, tam)
  intercepto <- 1
  x_theta <- data.frame(data_set[,cov_theta])
  x_beta <- data.frame(data_set[,cov_beta])
  Uh <- runif(tam)
  for (i in 1:tam) {
    if ((delta[i] == 0) | (L[i] == R[i])) {new_y[i] <- L[i]}
    else {
      xi_0 <- as.numeric(cbind(intercepto, x_theta[i,]))
      xi_1 <- as.numeric(x_beta[i,])
      new_y[i] <- inverse_lam_f(Uh[i], L[i], R[i], xi_0, xi_1, theta, beta, naalen_original)
    }
  }
  return(new_y)
}

# Generates latent variables k for n observations
gera_kh <- function(y_h, data_set, delta, cov_theta, cov_beta, theta, beta, nelson_aalen_f) {
  k_h <- y_h * NA
  intercepto <- 1
  xi_0 <- cbind(intercepto, data_set[,cov_theta])
  xi_1 <- data_set[,cov_beta]
  num <- exp(as.vector(theta %*% t(xi_0)))
  den <- 2 + 4 * sapply(y_h, nelson_aalen_f) * exp(as.vector(beta %*% t(xi_1)))
  k_h <- rpois(length(num), num / den) + delta
  return(k_h)
}

# Generates latent variables u for n observations
gera_uh <- function(y_h, k_h, data_set, R, delta, cov_beta, beta, nelson_aalen_f) {
  u_h <- y_h * NA
  r_estrela <- max(R[delta == 1])
  xi_1 <- data_set[,cov_beta]
  alpha_gamma <- k_h + delta
  beta_gamma <- 1 / (0.5 + sapply(y_h, nelson_aalen_f) * exp(as.vector(beta %*% t(xi_1))))
  cond_u <- (k_h == 0 | y_h > r_estrela)
  u_h <- ifelse(cond_u, 0, rgamma(length(alpha_gamma), alpha_gamma, scale=beta_gamma))
  return(u_h)
}


# Returns covariance matrix from (5) of (Lam et al 2007)
var_matrix <- function(sum_var, alpha_matrix) {
  M <- nrow(alpha_matrix)
  alpha_j <- Matrix::colMeans(alpha_matrix)
  matriz_soma_2 <- diag(ncol(alpha_matrix)) * 0
  for (h in 1:M) {
    matriz_soma_2 <- matriz_soma_2 + ((alpha_matrix[h,] - alpha_j) %*% t(alpha_matrix[h,] - alpha_j)) / (M - 1)
  }
  sigma_alpha <- (sum_var / M) + (1 + 1 / M) * matriz_soma_2
  return(sigma_alpha)
}

# Convergence criteria returning TRUE of FALSE
convergence_lam <- function(alpha_new, alpha_old, tol = 0.001) {
  conv <- FALSE
  new_val <- c(alpha_new)
  old_val <- c(alpha_old)
  max_error <- max(abs(new_val - old_val))
  if (max_error < tol) conv <- TRUE
  return(conv)
}


#' Fits cure rate frailty model for interval censored data.
#'
#' \code{inter_frailty} returns a list with the estimated parameters \code{par}
#' and their covariance matrix \code{mcov}. The list also contains the cure rate
#' covariance estimates \code{mcov.cura} for cure rate part only and a dummy
#' variable \code{stop_c} assuming 0 if algorithm converged and 1 if a stop
#' criteria ended the process.
#'
#' @param data_set Dataset used to fit the model.
#' @param L Vector containing the last check times before event.
#' @param R Vector containing the first check times after event.
#' @param delta Flag vector indicating failure inside interval.
#' @param cov_theta String vector containing the column names to be used on the
#'   cure rate predictor.
#' @param cov_beta String vector containing the column names to be used on the
#'   predictor associated with the hazard function.
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
#' @return \code{stop_c} stop criteria indicator assuming 1 when process is
#'   stopped for a non-convergence criteria. Assumes 0 when convergence is
#'   reached.
#' @examples
#' sample_set <- sim_frailty_data(100)
#' inter_frailty(sample_set, sample_set$L, sample_set$R, sample_set$delta, c("xi1","xi2"), c("xi1","xi2"), M = 50)
#' inter_frailty(sample_set, sample_set$L, sample_set$R, sample_set$delta, c("xi1"), c("xi2"), M = 10)
#' @export
#' @import foreach
inter_frailty <- function(data_set, L, R, delta, cov_theta, cov_beta, M, b=0.001, N_INT_MAX=100, par_cl = NULL, burn_in=30, outputFiles = FALSE) {
# Defining output variables
  if (outputFiles) {
    est_file_name <- paste("LAM_Estimates.txt", sep="")
    var_file_name <- paste("LAM_Variances.txt", sep="")
    fileConn <- file(est_file_name, "w")
    write(paste("PAR:\t",paste0(c("Intercept",cov_theta,cov_beta), collapse="\t")), file=fileConn, append=T, sep="")
  }


  # Initial values for y
  y_nxM = k_nxM = u_nxM = matrix(NA, nrow=M, ncol = nrow(data_set))
  for(i in 1:M) y_nxM[i,] <- ifelse(delta == 1, (L + R) / 2, L)

  # Initial values for the parameters
  compr_theta <- 1 + length(cov_theta); compr_beta <- length(cov_beta)
  compr_alpha <- compr_theta + compr_beta
  a_M <- matrix(NA, nrow=M, ncol=compr_alpha)
  rotulos <- c("intercepto", cov_theta, cov_beta)
  colnames(a_M) <- rotulos
  theta_M <- matrix(a_M[,1:compr_theta], nrow=M)
  colnames(theta_M) <- rotulos[1:compr_theta]
  beta_M <- matrix(a_M[,(compr_theta + 1):compr_alpha], nrow=M)
  colnames(beta_M) <- rotulos[(compr_theta + 1):compr_alpha]
  alpha <- c(1:compr_alpha) * 0
  sigma_alpha <- b * diag(compr_alpha); beta <- alpha[(compr_theta + 1):compr_alpha]

  # Initial values for latent vector u
  u <- delta

  # Initial Nelson-Aalen estimator
  Vetores_NAalen <- nelson_aalen_table(data_set, y_nxM[1,], delta, beta, cov_beta, u)
  NAalen_MEDIA <- stepfun(Vetores_NAalen$time, c(0,Vetores_NAalen$hazard))

  # Initializing convergence criteria and parameters
  conv <- FALSE; n <- 0
  a_M_NEW <- a_M

  # Iterative process (with parallel computing)
  while(!conv | n <= burn_in) {
    cat("ITER#", (n + 1))
    tempoITER <- system.time({
    if(!is.null(par_cl)){
      list_reg <- foreach(iterators::icount(M), .packages=c("MASS","Matrix","survival"), .export=c("surv_lam","inverse_lam_f", "gera_yh", "gera_kh", "gera_uh"), .inorder=F) %dopar% {
        a_M <- mvrnorm(n=1, alpha, sigma_alpha)
        theta_M <- a_M[1:compr_theta]; beta_M <- a_M[(compr_theta + 1):compr_alpha]
        y <- gera_yh(data_set, L, R, delta, cov_theta, cov_beta, as.numeric(theta_M), as.numeric(beta_M), NAalen_MEDIA)
        k <- gera_kh(y, data_set, delta, cov_theta, cov_beta, theta_M, beta_M, NAalen_MEDIA)
        u <- gera_uh(y, k, data_set, R, delta, cov_beta, beta_M, NAalen_MEDIA)

        # Poisson Regression for Theta
        o_set <- k * 0 - log(2)
        expression_theta <- paste("data_set$", cov_theta[1:length(cov_theta)], sep = "", collapse="+")
        eq_theta <- paste("fit_theta <- glm(k~", expression_theta, "+offset(o_set),family=poisson)")
        eval(parse(text = eq_theta))

        # Cox Regression for Beta
        expression_beta <- paste("data_set$", cov_beta[1:length(cov_beta)] ,sep = "", collapse="+")
        eq_beta <- paste("fit_beta <- coxph(Surv(y,delta)~", expression_beta," + offset(ifelse(log(u)==-Inf,-200,log(u))),method='breslow')")
        eval(parse(text = eq_beta))

        # Outputs of Parallel Computing
        out <- list(fit_theta$coef, vcov(fit_theta), fit_beta$coef, vcov(fit_beta), y, u)
        out
      }
    } else {
      list_reg <- foreach(iterators::icount(M), .packages=c("MASS","Matrix","survival"), .export=c("surv_lam","inverse_lam_f", "gera_yh", "gera_kh", "gera_uh"), .inorder=F) %do% {
        a_M <- mvrnorm(n=1, alpha, sigma_alpha)
        theta_M <- a_M[1:compr_theta]; beta_M <- a_M[(compr_theta + 1):compr_alpha]
        y <- gera_yh(data_set, L, R, delta, cov_theta, cov_beta, as.numeric(theta_M), as.numeric(beta_M), NAalen_MEDIA)
        k <- gera_kh(y, data_set, delta, cov_theta, cov_beta, theta_M, beta_M, NAalen_MEDIA)
        u <- gera_uh(y, k, data_set, R, delta, cov_beta, beta_M, NAalen_MEDIA)

        # Poisson Regression for Theta
        o_set <- k * 0 - log(2)
        expression_theta <- paste("data_set$", cov_theta[1:length(cov_theta)], sep = "", collapse="+")
        eq_theta <- paste("fit_theta <- stats::glm(k~", expression_theta, "+offset(o_set),family=poisson)")
        eval(parse(text = eq_theta))

        # Cox Regression for Beta
        expression_beta <- paste("data_set$", cov_beta[1:length(cov_beta)] ,sep = "", collapse="+")
        eq_beta <- paste("fit_beta <- survival::coxph(Surv(y,delta)~", expression_beta," + offset(ifelse(log(u)==-Inf,-200,log(u))),method='breslow')")
        eval(parse(text = eq_beta))

        # Outputs of Parallel Computing
        out <- list(fit_theta$coef, vcov(fit_theta), fit_beta$coef, vcov(fit_beta), y, u)
        out
      }
    }



    # Allocating M new and auxiliary parameter vectors
    SUM_VAR_THETA = SUM_VAR_BETA = 0
    for (h in 1:M) {
      SUM_VAR_THETA = SUM_VAR_THETA + list_reg[[h]][[2]]; SUM_VAR_BETA = SUM_VAR_BETA + list_reg[[h]][[4]]
      a_M_NEW[h,1:compr_theta] = list_reg[[h]][[1]]
      a_M_NEW[h,(compr_theta + 1):compr_alpha] = list_reg[[h]][[3]]
      y_nxM[h,] = list_reg[[h]][[5]]
      u_nxM[h,] = list_reg[[h]][[6]]
    }

    # Matrix of the M beta vectors
    beta_M <- matrix(a_M_NEW[,(compr_theta + 1):compr_alpha], nrow=M)

    # Obtaining new Nelson-Aalen estimator for Cum. Hazard function
    if (!is.null(par_cl)) {
      step_list <- foreach(h=1:M, .export="nelson_aalen_table", .inorder=F) %dopar% {
        V_NAalen <- nelson_aalen_table(data_set, y_nxM[h,], delta, beta_M[h,], cov_beta, u_nxM[h,])
        step_list <- stepfun(V_NAalen$time, c(0, V_NAalen$hazard))
        step_list
      }
    } else {
      step_list <- foreach(h=1:M, .export="nelson_aalen_table", .inorder=F) %do% {
        V_NAalen <- nelson_aalen_table(data_set, y_nxM[h,], delta, beta_M[h,], cov_beta, u_nxM[h,])
        step_list <- stepfun(V_NAalen$time, c(0, V_NAalen$hazard))
        step_list
      }
    }

    expression <- paste("step_list[[", 1:M,"]](x)", sep = "", collapse = "+")
    eq4 <- paste("NAalen_MEDIA <- function(x) (",expression,")/", M)
    eval(parse(text = eq4))

    # Creating new times/survival table and a more efficient estimator
    V_NAalen <- aux_naalen(sort(y_nxM[,delta==1]), NAalen_MEDIA, par_cl)
    NAalen_MEDIAnew <- stepfun(V_NAalen$time, c(0, V_NAalen$hazard))

    #Calculating the new covariance matrix
    SUM_VAR <- as.matrix(Matrix::bdiag(list(SUM_VAR_THETA, SUM_VAR_BETA)))
    VAR <- var_matrix(SUM_VAR, a_M_NEW)
    sigma_alpha <- VAR

    #New vector of estimates
    alpha_new <- Matrix::colMeans(a_M_NEW)
    })
    print(tempoITER)

    #Checking convergence
    conv <- convergence_lam(alpha_new ,alpha)

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
    if (n == (N_INT_MAX + burn_in)) {
      if (outputFiles) {
        write("Warning: Iteration Number achieved but convergence criteria not met.", file=fileConn, append=T, sep=" ")
        close(fileConn)
      }
      cat("Warning: Convergence criteria not met. Estimates given for N_INT_MAX=", N_INT_MAX)
      cat("\n")
      break
    }
  }
  #Kills the parallel proccess
  cPar <- as.numeric(n == (N_INT_MAX + burn_in))
  alphaList <- list(par = alpha, mcov = VAR, mcov.cura = VAR[1:(1 + length(cov_theta)), 1:(1 + length(cov_theta))], stop_c = cPar)
  if (outputFiles) close(fileConn)
  return(alphaList)
}

