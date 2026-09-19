test_that("cptgascisl returns a cptgascisl object with sensible slots", {
  d <- make_demo_data()
  
  res <- cptgascisl(
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
  
  expect_s4_class(res, "cptgascisl")
  expect_equal(res@N, d$N)
  expect_equal(res@popSize, 20)
  expect_equal(res@numIslands, 2)
  expect_equal(res@Islandsize, floor(20 / 2))
  expect_true(is.numeric(res@overbestfit))
  expect_length(res@overbestfit, 1)
  expect_true(is.finite(res@overbestfit))
  expect_true(is.numeric(res@overbestchrom))
  expect_true(length(res@overbestchrom) >= 2)
  
  expect_true(res@consensus)
  expect_true(res@native_consensus)
  expect_equal(res@consensus_radius, 2)
  expect_equal(res@consensus_lambda, 0.5)
  expect_equal(res@consensus_candidates, 5)
  expect_equal(res@local_refine, "periodic")
  expect_equal(res@local_every, 5)
  expect_equal(res@local_radius, 5)
  expect_equal(res@local_max_passes, 1)
  expect_true(is.numeric(res@n_eval_local))
})


test_that("cptgascisl is reproducible with a fixed seed", {
  d <- make_demo_data()
  
  res1 <- cptgascisl(
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
  
  res2 <- cptgascisl(
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


test_that("cptgascisl summary and plot methods work", {
  d <- make_demo_data()
  
  res <- cptgascisl(
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
  expect_true(any(grepl("Structured and Consensus-Guided", summary_out)))
  
  tf <- tempfile(fileext = ".pdf")
  grDevices::pdf(tf)
  on.exit(grDevices::dev.off(), add = TRUE)
  
  expect_no_error(plot(res, data = d$Xt))
})


test_that("cptgascisl works without consensus and local refinement", {
  d <- make_demo_data()
  
  res <- cptgascisl(
    ObjFunc = arima_bic,
    N = d$N,
    XMat = d$XMatT,
    Xt = d$Xt,
    popSize = 20,
    numIslands = 2,
    maxMig = 2,
    maxgen = 3,
    maxconv = 1,
    consensus = FALSE,
    local_refine = "none",
    parallel = FALSE,
    monitoring = FALSE,
    seed = 99
  )
  
  expect_s4_class(res, "cptgascisl")
  expect_false(res@consensus)
  expect_equal(res@local_refine, "none")
  expect_equal(res@n_eval_local, 0)
  expect_true(is.finite(res@overbestfit))
})


test_that("cptgascisl supports simultaneous model order selection", {
  d <- make_demo_data()
  
  prange <- list(ar = c(0, 2), ma = c(0, 2))
  
  res <- suppressWarnings(cptgascisl(
    ObjFunc = arima_bic_order_pq,
    N = d$N,
    prange = prange,
    option = "both",
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
  ))
  
  expect_s4_class(res, "cptgascisl")
  expect_equal(res@option, "both")
  expect_equal(res@prange, prange)
  expect_true(is.finite(res@overbestfit))
  
  m <- res@overbestchrom[1]
  ar_order <- res@overbestchrom[2]
  ma_order <- res@overbestchrom[3]
  
  expect_true(ar_order %in% 0:2)
  expect_true(ma_order %in% 0:2)
  expect_true(m >= 0)
  
  summary_out <- capture.output(summary(res))
  expect_true(any(grepl("Model hyperparameters", summary_out)))
})