test_that("cpt_dist is zero for identical changepoint configurations", {
  tau <- c(20, 35, 70)
  expect_equal(cpt_dist(tau1 = tau, tau2 = tau, N = 100), 0)
})

test_that("cpt_dist handles empty changepoint configurations", {
  tau <- c(25, 50, 75)
  
  expect_equal(cpt_dist(tau1 = NULL, tau2 = tau, N = 100), 3)
  expect_equal(cpt_dist(tau1 = tau, tau2 = NULL, N = 100), 3)
  expect_warning(
    expect_equal(cpt_dist(tau1 = NULL, tau2 = NULL, N = 100), 0),
    "Both configurations are NULL."
  )
})

test_that("cpt_dist is symmetric", {
  tau1 <- c(20, 35, 70, 80, 90)
  tau2 <- c(25, 50, 75)
  
  expect_equal(
    cpt_dist(tau1 = tau1, tau2 = tau2, N = 100),
    cpt_dist(tau1 = tau2, tau2 = tau1, N = 100)
  )
})

test_that("cpt_dist returns a finite numeric scalar", {
  out <- cpt_dist(tau1 = c(20, 40), tau2 = c(25, 50), N = 100)
  
  expect_true(is.numeric(out))
  expect_length(out, 1)
  expect_true(is.finite(out))
})