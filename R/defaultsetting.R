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

.default.GA_operators = list(population = "random_population_cpp",
                            selection  = "selection_linearrank_cpp",
                            crossover  = "offspring_uniformcrossover_cpp",
                            mutation   = "mutation")


.default.IslandGA_param = list(
  popsize      = 40,
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

.default.IslandGA_operators = list(population = "random_population_cpp",
                                   selection  = "selection_linearrank_cpp",
                                   crossover  = "offspring_uniformcrossover_cpp",
                                   mutation   = "mutation")
