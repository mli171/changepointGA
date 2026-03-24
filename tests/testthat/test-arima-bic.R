test_that("arima_bic returns a finite numeric scalar", {
  Ts <- 200
  XMatT <- matrix(1, nrow = Ts, ncol = 1)
  colnames(XMatT) <- "intercept"
  
  Xt <- ts_sim(
    Ts = Ts,
    beta = 0.5,
    XMat = XMatT,
    sigma = 1,
    phi = 0.5,
    theta = NULL,
    Delta = c(2, -2),
    CpLoc = c(50, 150),
    seed = 1234
  )
  
  chromosome <- c(2, 50, 150, Ts + 1)
  
  out <- arima_bic(chromosome, XMat = XMatT, Xt = Xt)
  
  expect_true(is.numeric(out))
  expect_length(out, 1)
  expect_true(is.finite(out))
})

test_that("arima_bic is deterministic for fixed input", {
  Ts <- 200
  XMatT <- matrix(1, nrow = Ts, ncol = 1)
  colnames(XMatT) <- "intercept"
  
  Xt <- ts_sim(
    Ts = Ts,
    beta = 0.5,
    XMat = XMatT,
    sigma = 1,
    phi = 0.5,
    theta = NULL,
    Delta = c(2, -2),
    CpLoc = c(50, 150),
    seed = 1234
  )
  
  chromosome <- c(2, 50, 150, Ts + 1)
  
  out1 <- arima_bic(chromosome, XMat = XMatT, Xt = Xt)
  out2 <- arima_bic(chromosome, XMat = XMatT, Xt = Xt)
  
  expect_equal(out1, out2)
})