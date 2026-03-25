test_that("arima_bic_order_pdq returns a numeric value", {
  Ts <- 200
  XMatT <- matrix(1, nrow = Ts, ncol = 1)
  colnames(XMatT) <- "intercept"
  
  Xt <- ts_sim(
    Ts = Ts,
    beta = 0.5,
    XMat = XMatT,
    sigma = 1,
    phi = 0.5,
    theta = 0.3,
    d = 1,
    Delta = c(2, -2),
    CpLoc = c(50, 150),
    seed = 1234
  )
  
  chromosome <- c(2, 1, 1, 1, 50, 150)
  
  bic_obj <- arima_bic_order_pdq(
    chromosome = chromosome,
    plen = 3,
    XMat = XMatT,
    Xt = Xt
  )
  
  expect_true(is.numeric(bic_obj))
  expect_length(bic_obj, 1)
  expect_true(is.finite(bic_obj))
})

test_that("arima_bic_order_pdq works when there is no changepoint", {
  Ts <- 200
  XMatT <- matrix(1, nrow = Ts, ncol = 1)
  colnames(XMatT) <- "intercept"
  
  Xt <- suppressWarnings(
    ts_sim(
      Ts = Ts,
      beta = 0.5,
      XMat = XMatT,
      sigma = 1,
      phi = 0.5,
      theta = 0.3,
      d = 1,
      seed = 1234
    )
  )
  
  chromosome <- c(0, 1, 1, 1)
  
  bic_obj <- arima_bic_order_pdq(
    chromosome = chromosome,
    plen = 3,
    XMat = XMatT,
    Xt = Xt
  )
  
  expect_true(is.numeric(bic_obj))
  expect_length(bic_obj, 1)
  expect_true(is.finite(bic_obj))
})

test_that("arima_bic_order_pdq rounds model orders to integers", {
  Ts <- 200
  XMatT <- matrix(1, nrow = Ts, ncol = 1)
  colnames(XMatT) <- "intercept"
  
  Xt <- ts_sim(
    Ts = Ts,
    beta = 0.5,
    XMat = XMatT,
    sigma = 1,
    phi = 0.5,
    theta = 0.3,
    d = 1,
    Delta = c(2, -2),
    CpLoc = c(50, 150),
    seed = 1234
  )
  
  bic1 <- arima_bic_order_pdq(
    chromosome = c(2, 1, 1, 1, 50, 150),
    plen = 3,
    XMat = XMatT,
    Xt = Xt
  )
  
  bic2 <- arima_bic_order_pdq(
    chromosome = c(2, 1.2, 0.8, 1.1, 50, 150),
    plen = 3,
    XMat = XMatT,
    Xt = Xt
  )
  
  expect_equal(bic1, bic2)
})

test_that("arima_bic_order_pdq ignores invalid changepoint locations", {
  Ts <- 200
  XMatT <- matrix(1, nrow = Ts, ncol = 1)
  colnames(XMatT) <- "intercept"
  
  Xt <- ts_sim(
    Ts = Ts,
    beta = 0.5,
    XMat = XMatT,
    sigma = 1,
    phi = 0.5,
    theta = 0.3,
    d = 1,
    Delta = c(2, -2),
    CpLoc = c(50, 150),
    seed = 1234
  )
  
  bic_obj <- arima_bic_order_pdq(
    chromosome = c(2, 1, 1, 1, -5, 250),
    plen = 3,
    XMat = XMatT,
    Xt = Xt
  )
  
  expect_true(is.numeric(bic_obj))
  expect_length(bic_obj, 1)
})

test_that("arima_bic_order_pdq returns penalty when fitting fails", {
  Ts <- 50
  XMatT <- matrix(1, nrow = Ts, ncol = 1)
  colnames(XMatT) <- "intercept"
  
  Xt <- rep(1, Ts)
  
  bic_obj <- arima_bic_order_pdq(
    chromosome = c(1, 5, 2, 5, 25),
    plen = 3,
    XMat = XMatT,
    Xt = Xt
  )
  
  expect_equal(bic_obj, 1e10)
})