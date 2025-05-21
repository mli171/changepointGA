#-------------------------- Check Convergence
checkConv = function(a, maxconv, tol){
  # This function is used to check convergence criterion for overall GA
  # some inputs ++++++++++++++++++
  #   a= the best fitness values for maxconv consecutive migrations
  #   maxconv= if maxconv consecutive migrations, the overall best does not change, then stop
  #   tol= tolerance level for iterations
  # outputs ++++++++++++++++++
  #   decision= 1 means continue and 0 means stop
  i<-1
  repeat{
    diff = abs(a[i+1]-a[i])
    if (diff < tol){
      i = i+1
      if (i >= maxconv) {
        decision = 1
        break
      }
    }
    else{
      decision = 0
      break
    }
  }
  return(decision)
}

#' The default mutation operator in genetic algorithm
#'
#' In a certain probability, the \code{mutation} genetic operator can be applied
#' to generate a new \code{child}. By default, the new child selection can be
#' down by the similar individual selection method in population initialization,
#' \code{\link{selectTau}}.
#'
#' @param child The child chromosome resulting from the \code{crossover} genetic
#' operator.
#' @param p.range Default is \code{NULL} for only changepoint detection. If
#' \code{p.range} is specified as a list object, which contains the range of
#' each model order parameters for order selection (integers). The number of
#' order parameters must be equal to the length of \code{p.range}.
#' @param minDist The required minimum distance between two adjacent changepoints.
#' @param Pb The probability of changepoints for every time series.
#' @param lmax The user specified maximum number of changepoints, by default,
#' as \code{N/2 - 1}.
#' @param mmax The user specified maximum length of individual chromosome,
#' by default, as \code{2+N/2 + 1}.
#' @param N The sample size of the time series.
#' @details
#' A function can apply mutation to the produced child with the specified
#' probability \code{Pmutation} in \code{\link{GA_param}} and
#' \code{\link{IslandGA_param}}. If order selection is not requested
#' (\code{option = "cp"} in \code{GA_param} and \code{Island_GA}), the default
#' \code{\link{mutation}} operator function uses \code{selectTau} to select
#' a completely new individual with a new chromosome as the mutated child.
#' For details, see \code{\link{selectTau}}. If order selection is needed
#' (\code{option = "both"} in \code{GA_param} and \code{Island_GA}), we first
#' decide whether to keep the produced child's model order with a probability
#' of 0.5. If the child's model order is retained, the \code{selectTau}
#' function is used to select a completely new individual with a new chromosome
#' as the mutated child. If a new model order is selected from the candidate
#' model order set, there is a 0.5 probability to either select a completely new
#' individual with new changepoint locations or retain the original child's
#' changepoint locations for the mutated child. Note that the current model
#' orders in the child's chromosome are excluded from the set to avoid redundant
#' objective function evaluation. Finally, the function returns a vector
#' containing the modified chromosomes for mutated \code{child}.
#' @return The resulting child chromosome representation.
#' @import Rcpp
#' @import stats
#' @import graphics
#' @useDynLib changepointGA
#' @export
mutation = function(child, p.range=NULL, minDist, Pb, lmax, mmax, N){

  plen = length(p.range)

  if(plen>0){
    childMut = matrix(0, nrow=lmax, 1)
    a1 = runif(1)
    if(a1 > 0.5){
      # 1. order from child
      childMut[2:(plen+1),] = child[2:(plen+1)]
      # 1.1 cpt from new
      tmpchildMut = selectTau(N=N, prange=NULL, minDist=minDist, Pb=Pb,
                                  mmax=mmax, lmax=lmax)
      childMut[1,] = tmpchildMut[1]
      childMut[(plen+2):(plen+tmpchildMut[1]+2),] = tmpchildMut[2:(tmpchildMut[1]+2)]
    }else{
      # 2. order from new
      new.p.range = rep(NA, plen)
      for(ii in 1:plen){
        tmp.p.range = setdiff(p.range[[ii]][1]:p.range[[ii]][2], child[2+ii-1])
        new.p.range[ii] = sample(tmp.p.range, 1)
      }
      childMut[2:(plen+1),] = new.p.range
      a2 = runif(1)
      if(a2 > 0.5){
        # 2.1 cpt from new
        tmpchildMut = selectTau(N=N, prange=NULL, minDist=minDist, Pb=Pb,
                                    mmax=mmax, lmax=lmax)
        childMut[1,] = tmpchildMut[1]
        childMut[(plen+2):(plen+tmpchildMut[1]+2),] = tmpchildMut[2:(tmpchildMut[1]+2)]
      }else{
        # 2.2 cpt from child
        childMut[1,] = child[1]
        childMut[(plen+2):(plen+child[1]+2),] = child[(plen+2):(plen+child[1]+2)]
      }
    }
  }else{
    tmpchildMut = selectTau(N=N, prange=NULL, minDist=minDist, Pb=Pb,
                                mmax=mmax, lmax=lmax)
    childMut = tmpchildMut
  }

  return(childMut)
}

NewpopulationIsland = function(ObjFunc, selection, crossover, mutation, pop, fit, minDist, lmax, mmax, Pc, Pm, Pb, maxgen, N, p.range, ...){
  # This function is used to form new population
  # some inputs ++++++++++++++++++
  #   pop= population
  #   fit= fitness evaluated for population
  #   minDist= minimum distances between two adjacent changepoints
  #   lmax= max length of chromosome
  #   mmax= max number of changepoints
  #   Pc= prob of crossover
  #   Pm= prob of mutation
  #   Pb= prob of changepoints for every time series
  #   maxgen= for each subpopulation, after maxgen then apply migration
  #   N= sample size
  #   X_hour= categorical time series
  #   XMat= Design matrix including covariate other than changepoint
  #   penalty= selection criterion to choose
  # outputs ++++++++++++++++++
  #   pop= updated population
  #   fit= updated population fitness
  #   bestfit  = currnt minimum of fitness values
  #   bestchrom = chromosome representation of the individual associated with bestfit

  plen = length(p.range)

  count = 1
  repeat{
    # indicator for c("crossover", "mutation")
    #     flag[1]=1 indicating no cross-over
    #     flag[2]=1 indicating no mutation
    flag = rep(0, 2)
    ##### step 2: parents selection
    parents = selection(pop, fit)
    dad = parents$dad
    mom = parents$mom

    ##### step 3: crossover
    a1 = runif(1)
    if(a1 <= Pc){
      child = crossover(mom, dad, p.range, minDist, lmax, N)
    }else{
      child = dad
      flag[1] = 1
    }

    ## step 4-2: mutation
    a2 = runif(1)
    if(a2 <= Pm){
      child = mutation(child, p.range, minDist, Pb, lmax, mmax, N)
    }else{
      flag[2] = 1
    }

    ## step 5: form new generation
    #  steady state method:
    #     replace the least fit in the current pop with child if child is better.
    flagsum = flag[1] + flag[2]
    if (flagsum<2){
      # flagsum < 2 indicating new individual produced and fitness evaluation needed
      fitChild = do.call(ObjFunc, c(list(child[1:(child[1]+plen+2)], plen, ...)))
      # fitChild = do.call(ObjFunc, c(list(child[1:(child[1]+plen+2)], plen, XMat, Xt)))
      leastfit = max(fit) # with largest fitness value

      if (fitChild < leastfit) {
        # indicating child is better than the worst one and replace
        pp = which.max(fit)
        pop[,pp] = child
        fit[pp] = fitChild
      }
      count = count + 1
    }

    # check: after every maxgen generations, apply migration in GA.Main()
    if (count >= maxgen){break}
  }

  return(rbind(fit, pop))
}
