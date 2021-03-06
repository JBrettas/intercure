#' Generates a interval censored clustered dataset using the frailty cure rate
#' model
#'
#' \code{sim_frailty_cl} returns a dataset generated by the frailty cure rate
#' model for clustered data.
#' @param N Size of the sample to be generated.
#' @param theta Three parameters associated with the cure linear predictor.
#' @param beta Two parameters associated with the hazard function.
#' @param A A positive number representing a fixed right censoring.
#' @param B A positive number which multiplies an exponential random variable
#'   with mean 1, defining another right censoring case.
#' @param prob Probability that individual presents treatment T1 (baseline is
#'   T0).
#' @param nclus Number of clusters to generate with balanced sizes.
#' @param w Shape and rate parameters value for the Gamma distribution with mean
#'   1.
#' @return A generated dataset with columns: \code{Z}, the actual event time;
#'   \code{L}, the leftmost limit of the censored interval; \code{R}, the
#'   rightmost limit of the censored interval; \code{delta}, the failure
#'   indicator; \code{xi1}, the treatment covariate assuming 1 with probability
#'   \code{prob} and 0 otherwise; \code{xi2}, second variable generated by a
#'   standard normal distribution; \code{U}, the frailty for each individual,
#'   with value 0 for cured individuals. \code{clus}, representing the cluster
#'   id for each observation.
#' @examples
#' sim_frailty_cl(50)
#' @export
sim_frailty_cl <- function(N, theta = c(-1,1,0),
                        beta = c(0,0.5), A = 5, B = 15,
                        prob = 0.5,
                        nclus = 2,
                        w = exp(-0.5)) {
  u <- stats::runif(N)
  a <- stats::rexp(N)
  C <- cbind(A, a * B)
  C <- C[,1] * (C[,1] <= C[,2]) + C[,2] * (C[,1] > C[,2])
  intercept <- 1
  xi1 <- stats::rbinom(N,1,prob)
  xi2 <- stats::rnorm(N)
  cov_theta <- data.frame(intercept, xi1, xi2)
  cov_beta <- data.frame(xi1, xi2)
  clus <- sample(1:nclus, N, replace=TRUE,
                 prob = rep(1/nclus, nclus))
  clus <- as.character(clus)
  Ksis <- stats::rgamma(nclus, shape = w, rate = w)
  maptab <- data.frame(clus = 1:nclus, value = Ksis)
  maptab$clus <- as.character(maptab$clus)
  clus_effect <- maptab$value[match(clus, maptab$clus)]
  eta <- exp(as.vector(theta %*% t(cov_theta)))
  K_vector <- stats::rpois(N, (eta*clus_effect) / 2)
  U_vector <- K_vector * NA
  for(i in 1:length(K_vector)) {
    if(K_vector[i] == 0) U_vector[i] <- 0
    else{
      U_vector[i] <- 0
      for (j in 1:K_vector[i])
        U_vector[i] <- U_vector[i] + stats::rchisq(1, 2, ncp = 0)
    }
  }
  beta_x <- as.vector(beta %*% t(cov_beta))
  exp_pred_beta <- exp(beta_x)
  num <- -1 * log(1 - u)
  den <- U_vector * exp_pred_beta
  tempos <- ifelse(U_vector != 0, num / den, Inf)
  Z <- ifelse(tempos < C, tempos, C)
  delta <- ifelse(tempos < C, 1, 0)
  L <- R <- Z * NA
  for (i in 1:N) {
    if(delta[i] == 0) {
      L[i] <- Z[i]
      R[i] <- Inf
    }
    else {
      L[i] <- 0
      add <- stats::runif(1, 0.1, 0.5)
      R[i] <- add
      check <- (L[i] <= Z[i] & Z[i] < R[i])
      while(!check) {
        L[i] <- L[i] + add
        add <- stats::runif(1, 0.1, 0.5)
        R[i] <- R[i] + add
        check <- (L[i] <= Z[i] & Z[i] < R[i])
      }
    }
  }
  dados <- data.frame(Z, L, R, delta, xi1, xi2, U = U_vector, clus)
  return(dados)
}
