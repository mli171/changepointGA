test_that("mutation_birth_death_relocate returns a valid chromosome", {
  N <- 80
  lmax <- 42
  mmax <- 40
  
  child <- c(2, 20, 60, N + 1, rep(0, lmax - 4))
  
  set.seed(123)
  res <- mutation_birth_death_relocate(
    child = child,
    prange = NULL,
    minDist = 5,
    pchangepoint = 0.01,
    lmax = lmax,
    mmax = mmax,
    N = N
  )
  
  expect_true(is.matrix(res))
  expect_equal(dim(res), c(lmax, 1))
  
  K <- res[1, 1]
  expect_true(K %in% c(1, 2, 3))
  expect_equal(res[K + 2, 1], N + 1)
  
  if (K > 0) {
    tau <- res[2:(K + 1), 1]
    expect_true(all(tau > 1))
    expect_true(all(tau <= N))
    expect_equal(tau, sort(tau))
  }
})


test_that("mutation_birth_death_relocate preserves model order parameters", {
  N <- 80
  lmax <- 44
  mmax <- 40
  prange <- list(ar = c(0, 2), ma = c(0, 2))
  
  child <- c(2, 1, 1, 20, 60, N + 1, rep(0, lmax - 6))
  
  set.seed(123)
  res <- mutation_birth_death_relocate(
    child = child,
    prange = prange,
    minDist = 5,
    pchangepoint = 0.01,
    lmax = lmax,
    mmax = mmax,
    N = N
  )
  
  expect_equal(res[2:3, 1], c(1, 1))
  
  K <- res[1, 1]
  expect_equal(res[length(prange) + K + 2, 1], N + 1)
})


test_that("zero consensus weight reproduces unguided BDR mutation", {
  N <- 80
  lmax <- 42
  mmax <- 40
  child <- c(2, 20, 60, N + 1, rep(0, lmax - 4))
  consensus_score <- seq(0, 1, length.out = N)
  
  set.seed(123)
  res1 <- mutation_birth_death_relocate(
    child = child,
    prange = NULL,
    minDist = 5,
    pchangepoint = 0.01,
    lmax = lmax,
    mmax = mmax,
    N = N
  )
  
  set.seed(123)
  res2 <- mutation_birth_death_relocate(
    child = child,
    prange = NULL,
    minDist = 5,
    pchangepoint = 0.01,
    lmax = lmax,
    mmax = mmax,
    N = N,
    consensus_score = consensus_score,
    consensus_lambda = 0
  )
  
  expect_equal(res1, res2)
})