# devtools::use_package("doSNOW")

library(intercure)
set.seed(30)
data_test <- sim_frailty(30)


# one covariate
fit <- inter_frailty(data_test, data_test$L, data_test$R, data_test$delta, c("xi1"), c("xi1"), M = 10, max_n = 5, burn_in = 0)
fitx2 <- inter_frailty(data_test, data_test$L, data_test$R, data_test$delta, c("xi2"), c("xi2"), M = 10, max_n = 5, burn_in = 0)
test_that("just one covariate works on inter_frailty", {
  expect_is(fit$par,"numeric")
  expect_is(fit$stop_c,"numeric")
  expect_is(fitx2$par,"numeric")
  expect_is(fitx2$stop_c,"numeric")
})

# different covariate for each predictor
fit <- inter_frailty(data_test, data_test$L, data_test$R, data_test$delta, c("xi1"), c("xi2"), M = 10, max_n = 5, burn_in = 0)
fitx2 <- inter_frailty(data_test, data_test$L, data_test$R, data_test$delta, c("xi2"), c("xi1"), M = 10, max_n = 5, burn_in = 0)
test_that("program works for different covariates for each predictor", {
  expect_is(fit$par,"numeric")
  expect_is(fit$stop_c,"numeric")
  expect_is(fitx2$par,"numeric")
  expect_is(fitx2$stop_c,"numeric")
})

# parallelism
cl <- snow::makeCluster(2,type="SOCK")
doSNOW::registerDoSNOW(cl)
fit2 <- inter_frailty(data_test, data_test$L, data_test$R, data_test$delta, c("xi1"), c("xi1"),
                     M = 10, max_n = 4, burn_in = 1, par_cl = cl)

test_that("parallel mode is working", {
  expect_is(fit2$par,"numeric")
  expect_is(fit2$stop_c,"numeric")
})

snow::stopCluster.default(cl)


# M = 1 gives error
test_that("M = 1 gives error", {
  expect_error(inter_frailty(data_test, data_test$L, data_test$R, data_test$delta, c("xi1"), c("xi1"), M = 1, max_n = 5, burn_in = 0))
})


# M = 250 works
test_that("works for M = 250", {
  expect_is(inter_frailty(data_test, data_test$L, data_test$R, data_test$delta, c("xi1"), c("xi1"), M = 250, max_n = 5, burn_in = 0)$par, "numeric")
})


# error if delta contains different of 1 or 0
test_that("delta contains different from 1 or 0", {
  expect_is(inter_frailty(data_test, data_test$L, data_test$R, (1+data_test$delta), c("xi1"), c("xi1"), M = 250, max_n = 5, burn_in = 1)$par, "numeric")
})


