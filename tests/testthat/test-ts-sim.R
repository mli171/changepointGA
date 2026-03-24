test_that("ts_sim is reproducible with the same seed", {
  Ts <- 100
  XMatT <- matrix(1, nrow = Ts, ncol = 1)
  colnames(XMatT) <- "intercept"
  
  sim1 <- ts_sim(
    Ts = Ts,
    beta = 0.5,
    XMat = XMatT,
    sigma = 1,
    Delta = c(2, -2),
    CpLoc = c(25, 75),
    seed = 1234
  )
  
  sim2 <- ts_sim(
    Ts = Ts,
    beta = 0.5,
    XMat = XMatT,
    sigma = 1,
    Delta = c(2, -2),
    CpLoc = c(25, 75),
    seed = 1234
  )
  
  expect_equal(sim1, sim2)
})

test_that("ts_sim returns numeric output with expected length", {
  Ts <- 100
  XMatT <- matrix(1, nrow = Ts, ncol = 1)
  colnames(XMatT) <- "intercept"
  
  sim <- ts_sim(
    Ts = Ts,
    beta = 0.5,
    XMat = XMatT,
    sigma = 1,
    Delta = c(2, -2),
    CpLoc = c(25, 75),
    seed = 1234
  )
  
  expect_true(is.numeric(sim))
  expect_length(sim, Ts)
})

test_that("ts_sim works with ARMA terms", {
  Ts <- 100
  XMatT <- matrix(1, nrow = Ts, ncol = 1)
  colnames(XMatT) <- "intercept"
  
  sim <- ts_sim(
    Ts = Ts,
    beta = 0.5,
    XMat = XMatT,
    sigma = 1,
    phi = c(0.5, -0.5),
    theta = 0.8,
    d = 0,
    Delta = c(2, -2),
    CpLoc = c(25, 75),
    seed = 1234
  )
  
  expect_true(is.numeric(sim))
  expect_length(sim, Ts)
  expect_equal(attr(sim, "diffEff"), 0)
})

test_that("ts_sim works with ARIMA(1,1,1) terms", {
  Ts <- 100
  XMatT <- matrix(1, nrow = Ts, ncol = 1)
  colnames(XMatT) <- "intercept"
  
  sim <- ts_sim(
    Ts = Ts,
    beta = 0.5,
    XMat = XMatT,
    sigma = 1,
    phi = 0.5,
    theta = -0.3,
    d = 1,
    Delta = c(2, -2),
    CpLoc = c(25, 75),
    seed = 1234
  )
  
  expect_true(is.numeric(sim))
  expect_length(sim, Ts)
  expect_equal(attr(sim, "arEff"), 0.5)
  expect_equal(attr(sim, "maEff"), -0.3)
  expect_equal(attr(sim, "diffEff"), 1)
})