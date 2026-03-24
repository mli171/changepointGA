test_that("arima_bic_order_pdq returns numeric BIC without changepoints", {
  Ts <- 200
  XMatT <- matrix(1, nrow = Ts, ncol = 1)
  colnames(XMatT) <- "intercept"
  
  Xt <- suppressWarnings(
    ts_sim(
      Ts = Ts,
      beta = 0.5,
      XMat = XMatT,
      sigma = 1,
      phi = 0.4,
      theta = -0.2,
      d = 1,
      seed = 1234
    )
  )
  
  chromosome <- c(0, 1, 1, 1, Ts + 1)
  
  bic_obj <- arima_bic_order_pdq(
    chromosome = chromosome,
    plen = 3,
    XMat = XMatT,
    Xt = Xt
  )
  
  expect_true(is.numeric(bic_obj))
  expect_length(bic_obj, 1)
  expect_false(is.na(bic_obj))
})

test_that("arima_bic_order_pdq returns numeric BIC with changepoints", {
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
    d = 0,
    Delta = c(2, -2),
    CpLoc = c(50, 150),
    seed = 1234
  )
  
  chromosome <- c(2, 1, 0, 1, 50, 150, Ts + 1)
  
  bic_obj <- arima_bic_order_pdq(
    chromosome = chromosome,
    plen = 3,
    XMat = XMatT,
    Xt = Xt
  )
  
  expect_true(is.numeric(bic_obj))
  expect_length(bic_obj, 1)
  expect_false(is.na(bic_obj))
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
    theta = -0.3,
    d = 0,
    Delta = c(2, -2),
    CpLoc = c(50, 150),
    seed = 1234
  )
  
  chromosome1 <- c(2, 1, 0, 1, 50, 150, Ts + 1)
  chromosome2 <- c(2, 1.2, 0.2, 1.1, 50, 150, Ts + 1)
  
  bic1 <- arima_bic_order_pdq(
    chromosome = chromosome1,
    plen = 3,
    XMat = XMatT,
    Xt = Xt
  )
  
  bic2 <- arima_bic_order_pdq(
    chromosome = chromosome2,
    plen = 3,
    XMat = XMatT,
    Xt = Xt
  )
  
  expect_equal(bic1, bic2)
})

test_that("arima_bic_order_pdq handles invalid changepoint locations", {
  Ts <- 200
  XMatT <- matrix(1, nrow = Ts, ncol = 1)
  colnames(XMatT) <- "intercept"
  
  Xt <- suppressWarnings(
    ts_sim(
      Ts = Ts,
      beta = 0.5,
      XMat = XMatT,
      sigma = 1,
      phi = 0.4,
      theta = -0.2,
      d = 1,
      seed = 1234
    )
  )
  
  chromosome_bad <- c(2, 1, 1, 1, -5, 250, Ts + 1)
  
  bic_obj <- arima_bic_order_pdq(
    chromosome = chromosome_bad,
    plen = 3,
    XMat = XMatT,
    Xt = Xt
  )
  
  expect_true(is.numeric(bic_obj) || is.na(bic_obj))
  expect_length(bic_obj, 1)
})