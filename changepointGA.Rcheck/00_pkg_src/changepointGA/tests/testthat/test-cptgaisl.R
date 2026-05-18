# tests/testthat/test-cptgaisl.R

test_that("cptgaisl returns a cptgaisl object with sensible slots", {
  d <- make_demo_data()
  
  res <- cptgaisl(
    ObjFunc = arima_bic,
    N = d$N,
    XMat = d$XMatT,
    Xt = d$Xt,
    popSize = 20,
    numIslands = 2,
    maxMig = 2,
    maxgen = 3,
    maxconv = 1,
    parallel = FALSE,
    monitoring = FALSE,
    seed = 99
  )
  
  expect_s4_class(res, "cptgaisl")
  expect_equal(res@N, d$N)
  expect_equal(res@popSize, 20)
  expect_equal(res@numIslands, 2)
  expect_equal(res@Islandsize, floor(20 / 2))
  expect_true(is.numeric(res@overbestfit))
  expect_length(res@overbestfit, 1)
  expect_true(is.finite(res@overbestfit))
  expect_true(is.numeric(res@overbestchrom))
  expect_true(length(res@overbestchrom) >= 2)
})

test_that("cptgaisl is reproducible with a fixed seed", {
  d <- make_demo_data()
  
  res1 <- cptgaisl(
    ObjFunc = arima_bic,
    N = d$N,
    XMat = d$XMatT,
    Xt = d$Xt,
    popSize = 20,
    numIslands = 2,
    maxMig = 2,
    maxgen = 3,
    maxconv = 1,
    parallel = FALSE,
    monitoring = FALSE,
    seed = 123
  )
  
  res2 <- cptgaisl(
    ObjFunc = arima_bic,
    N = d$N,
    XMat = d$XMatT,
    Xt = d$Xt,
    popSize = 20,
    numIslands = 2,
    maxMig = 2,
    maxgen = 3,
    maxconv = 1,
    parallel = FALSE,
    monitoring = FALSE,
    seed = 123
  )
  
  expect_equal(res1@overbestchrom, res2@overbestchrom)
  expect_equal(res1@overbestfit, res2@overbestfit)
})

test_that("cptgaisl summary and plot methods work", {
  d <- make_demo_data()
  
  res <- cptgaisl(
    ObjFunc = arima_bic,
    N = d$N,
    XMat = d$XMatT,
    Xt = d$Xt,
    popSize = 20,
    numIslands = 2,
    maxMig = 2,
    maxgen = 3,
    maxconv = 1,
    parallel = FALSE,
    monitoring = FALSE,
    seed = 99
  )
  
  summary_out <- capture.output(summary(res))
  expect_true(length(summary_out) > 0)
  
  tf <- tempfile(fileext = ".pdf")
  grDevices::pdf(tf)
  on.exit(grDevices::dev.off(), add = TRUE)
  
  expect_no_error(plot(res, data = d$Xt))
})