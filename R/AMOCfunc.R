AMOCpopulation = function(popsize, p.range, N, minDist, Pchangepoint, mmax, lmax){

  tauclc = floor(0.05*N):ceiling(0.95*N)

  pop = matrix(0, nrow=lmax, ncol=popsize)
  pop[1,] = rep(1, popsize)
  pop[2,] = sample(tauclc, size=popsize)
  pop[3,] = N+1

  return(pop)
}

AMOCselection = function(pop, popFit){
  return(selection_linearrank_cpp(pop, popFit))
}

AMOCcrossover = function(mom, dad, p.range, minDist, lmax, N){

  child = 1
  child[2] = round((dad[2]+mom[2])/2)
  child[3] = N + 1

  return(child)
}

AMOCmutation = function(child, p.range, minDist, Pchangepoint, lmax, mmax, N){

  tmptau = 1

  while(tmptau < floor(0.05*N) | tmptau > ceiling(0.95*N)){
    tmpsign = sample(x=c(-1,1), size = 1)
    tmptau = child[2] + tmpsign*minDist
  }
  child[2] = tmptau

  return(child)
}
