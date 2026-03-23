test_that("cptga returns a cptga object with sensible slots", {
  d <- make_demo_data()
  
  res <- cptga(
    ObjFunc = arima_bic,
    N = d$N,
    XMat = d$XMatT,
    Xt = d$Xt,
    popSize = 12,
    maxgen = 4,
    maxconv = 2,
    parallel = FALSE,
    monitoring = FALSE,
    seed = 99
  )
  
  expect_s4_class(res, "cptga")
  expect_equal(res@N, d$N)
  expect_true(is.numeric(res@overbestfit))
  expect_length(res@overbestfit, 1)
  expect_true(is.finite(res@overbestfit))
  expect_true(is.numeric(res@bestfit))
})

test_that("cptga is reproducible with a fixed seed", {
  d <- make_demo_data()
  
  res1 <- cptga(
    ObjFunc = arima_bic,
    N = d$N,
    XMat = d$XMatT,
    Xt = d$Xt,
    popSize = 12,
    maxgen = 4,
    maxconv = 2,
    parallel = FALSE,
    monitoring = FALSE,
    seed = 123
  )
  
  res2 <- cptga(
    ObjFunc = arima_bic,
    N = d$N,
    XMat = d$XMatT,
    Xt = d$Xt,
    popSize = 12,
    maxgen = 4,
    maxconv = 2,
    parallel = FALSE,
    monitoring = FALSE,
    seed = 123
  )
  
  expect_equal(res1@overbestchrom, res2@overbestchrom)
  expect_equal(res1@overbestfit, res2@overbestfit)
})

test_that("cptga summary and plot methods work", {
  d <- make_demo_data()
  
  res <- cptga(
    ObjFunc = arima_bic,
    N = d$N,
    XMat = d$XMatT,
    Xt = d$Xt,
    popSize = 12,
    maxgen = 4,
    maxconv = 2,
    parallel = FALSE,
    monitoring = FALSE,
    seed = 99
  )
  
  s <- capture.output(summary(res))
  expect_true(length(s) > 0)
  expect_true(any(grepl("GA results|Optimal Fitness value", s)))
  
  tf <- tempfile(fileext = ".pdf")
  grDevices::pdf(tf)
  on.exit(grDevices::dev.off(), add = TRUE)
  
  expect_no_error(plot(res, data = d$Xt))
})