.default.GA_param = list(
  popsize      = 200,
  Pcrossover   = 0.95,
  Pmutation    = 0.15,
  Pchangepoint = 0.01,
  minDist      = 2,
  mmax         = 499,
  lmax         = 501,
  maxgen       = 100000,
  maxconv      = 1000,
  option       = "cp",
  monitoring   = FALSE,
  parallel     = FALSE,
  nCore        = NULL,
  tol          = 1e-5,
  seed         = NULL
)


.default.IslandGA_param = list(
  subpopsize   = 40,
  Islandsize   = 5,
  Pcrossover   = 0.95,
  Pmutation    = 0.15,
  Pchangepoint = 0.01,
  minDist      = 2,
  mmax         = 499,
  lmax         = 501,
  maxMig       = 500,
  maxgen       = 100,
  maxconv      = 100,
  option       = "cp",
  monitoring   = FALSE,
  parallel     = FALSE,
  nCore        = NULL,
  tol          = 1e-5,
  seed         = NULL
)

.default.operators = list(population = "random_population",
                          selection  = "selection_linearrank",
                          crossover  = "uniformcrossover",
                          mutation   = "mutation")
