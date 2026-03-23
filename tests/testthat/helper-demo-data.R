make_demo_data <- function(N = 80) {
  XMatT <- matrix(1, nrow = N, ncol = 1)
  colnames(XMatT) <- "intercept"
  
  Xt <- ts_sim(
    beta = 0.5,
    XMat = XMatT,
    sigma = 1,
    phi = 0.5,
    theta = NULL,
    Delta = c(2, -2),
    CpLoc = c(20, 60),
    seed = 1234
  )
  
  list(N = N, XMatT = XMatT, Xt = Xt)
}