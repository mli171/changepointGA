// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>

using namespace Rcpp;
using namespace arma;

int m, i, j, k, id, jp;
int popsize, islandSize;
int lmax;

// [[Rcpp::export]]
IntegerVector rank_asR(NumericVector x, bool decreasing = false)
{
  IntegerVector rank = match(x, clone(x).sort());
  if(decreasing) rank = rank.length() + 1 - rank;
  return rank;
}

//' Randomly select the chromosome
//'
//' Randomly select the changepoint configuration for population initialization.
//' The selected changepoint configuration represents a changepoint chromosome.
//' The first element of the chromosome represent the number of changepoints
//' and the last non-zero element always equal to the length of time series + 1.
//'
//' @param N The length of time series.
//' @param prange A list object containing the possible range for other
//' pre-defined model parameters, i.e. AR/MA order of ARMA models.
//' @param minDist The minimum length between two adjacent changepoints.
//' @param Pb Same as \code{Pchangepoint}, the probability that a changepoint has occurred.
//' @param mmax The maximum possible number of changepoints in the data set.
//' @param lmax The maximum possible length of the chromosome representation.
//' @return A single changepoint configuration format as above.
//' @export
// [[Rcpp::export]]
arma::vec selectTau(int N, List prange, int minDist, double Pb, int mmax, int lmax){

  m = 0;
  double a;
  int plen = prange.length();
  arma::vec tau(lmax, fill::zeros);

  if(plen > 0){
    for(k=0; k<plen; k++){
      arma::vec tmp = prange[k];
      tau(k+1) = randi(distr_param(tmp[0],tmp[1]));
    }
  }

  i = 1 + minDist;
  while(i < N - minDist){
    a = runif(1)[0];
    if(a <= Pb){
      m = m + 1;
      tau(plen+m) = i;
      i = i + minDist;
    }else{
      i = i+1;
    }
  }

  tau(0) = m;
  tau(plen+m+1) = N+1;
  return(tau);
}

//' Random population initialization
//'
//' Randomly generate the individuals' chromosomes (changepoint confirgurations)
//' to construct the first generation population.
//'
//' @param popsize An integer represents the number of individual in each
//' population for GA (or subpopulation for IslandGA).
//' @param prange Default is \code{NULL} for only changepoint detection. If
//' \code{prange} is specified as a list object, which contains the range of
//' each model order parameters for order selection (integers). The number of
//' order parameters must be equal to the length of \code{prange}.
//' @param N The length of time series.
//' @param minDist The minimum length between two adjacent changepoints.
//' @param Pb Same as \code{Pchangepoint}, the probability that a changepoint has occurred.
//' @param mmax The maximum possible number of changepoints in the data set.
//' @param lmax The maximum possible length of the chromosome representation.
//' @details
//' The default population initialization uses \code{\link{selectTau}} to
//' select the chromosome for the first generation population. Each column from
//' the produced population matrix represent an chromosome of an individual.
//' The first element of every chromosome represent the number of changepoints
//' and the last non-zero element always equal to the length of time series
//' plus one (N+1).
//' @return A matrix that contains each individual's chromosome.
//' @export
// [[Rcpp::export]]
arma::mat random_population(int popsize, List prange, int N, int minDist, double Pb, int mmax, int lmax){

  arma::mat pop(lmax, popsize, fill::zeros);

  for(j=0;j<popsize;j++){
    pop.col(j) = selectTau(N, prange, minDist, Pb, mmax, lmax);
  }

  return(pop);
}

//' Uniform crossover to produce offsprings
//'
//' In uniform crossover, typically, each bit is chosen from either parent with
//' equal probability. Other mixing ratios are sometimes used, resulting in
//' offspring which inherit more genetic information from one parent than the
//' other. In a uniform crossover, we don’t divide the chromosome into segments,
//' rather we treat each gene separately. In this, we essentially flip a coin
//' for each chromosome to decide whether or not it will be included in the
//' off-spring. If model order selection is requested, each child's model order
//' has the equal probability (0.5) from \code{dad} and \code{mom}.
//' @param mom Among two selected individuals, \code{mom} represents the selected
//' chromosome representation with lower fitness function value.
//' @param dad Among two selected individuals, \code{dad} represents the selected
//' chromosome representation with larger fitness function value.
//' @param prange Default value is \code{NULL} for only changepoint detection. If
//' \code{prange} is specified as a list object, which contains the range of
//' each model order parameters for order selection (integers). The number of
//' order parameters must be equal to the length of \code{prange}.
//' @param minDist The required minimum distance between two adjacent changepoints.
//' @param lmax The user specified maximum number of changepoints, by default,
//' as \code{N/2 - 1}.
//' @param N The length of time series.
//' @return The child chromosome that produced from \code{mom} and \code{dad} for
//' next generation.
//' @export
// [[Rcpp::export]]
arma::vec uniformcrossover(arma::vec& mom, arma::vec& dad, List prange, int minDist, int lmax, int N){

  int plen = prange.length();

  arma::vec child(lmax, fill::zeros);

  // 2). obtain order from dad or mom with equal probability
  if(plen > 0){
    double ap = runif(1)[0];
    if(ap > 0.5){
      child.subvec(1,plen) = dad.subvec(1,plen);
    }else{
      child.subvec(1,plen) = mom.subvec(1,plen);
    }
  }

  int child_mmax = dad(0) + mom(0);
  if(child_mmax == 0){
    // no changepoint for both mom and dad
    child(0) = 0;
    child(1+plen) = N+1;
  }else{
    // 3). combine dad's and mom's changepoints
    int mom_m = mom(0);
    int dad_m = dad(0);
    arma::vec parents(child_mmax, fill::zeros);
    if(mom_m > 0){
      parents.subvec(0,mom_m-1) = mom.subvec(plen+1,plen+mom_m);
    }
    if(dad_m >0){
      parents.subvec(mom_m,child_mmax-1) = dad.subvec(plen+1,plen+dad_m);
    }

    // 2). remove same changepoints and sort
    parents = sort(unique(parents));
    int parents_count = parents.size();

    // 3). select subset from parents
    double a;
    int child_m = 0;
    int tmpm = 1;
    int tmptau = parents(0);
    int temp;
    do{
      a = runif(1)[0];
      if(a > 0.5){
        child_m = child_m + 1;
        child(plen+child_m) = tmptau;
        tmptau = tmptau + 2*minDist;
        if(tmptau >= N-minDist){break;} // boundary changepoints limits
      }

      if(tmpm >= parents_count-1){break;} // if the number of changepoints of child is larger than parents, break

      temp = parents(tmpm);

      if(tmptau < temp){
        tmptau = temp;
      }else{
        tmpm = tmpm + 1;
        if(tmpm >= child_mmax){break;}
        tmptau = parents(tmpm);
      }

      tmpm = tmpm + 1;
      if(tmptau >= N-minDist){break;}
    }
    while(1<2);

    child(0) = child_m;
    child(plen+child_m+1) = N+1;
  }

  return(child);
}


//' The default parents selection genetic algorithm operator
//'
//' The genetic algorithm require to select a pair of chromosomes, representing
//' \code{dad} and \code{mom}, for the \code{crossover} operator to
//' produce offspring (individual for next generation). The parents chromosomes
//' are randomly selectd from the initialized population by a linear ranking
//' method according to each individual's fittness in the input argument
//' \code{popFit}. By default, the dad has better fit/smaller fitness function
//' value/larger rank than \code{mom}.
//' @param pop A matrix contains the chromosomes for all individuals. The number of
//' rows is equal to \code{lmax} and the number of columns is equal to the
//' \code{popsize}.
//' @param popFit A vector contains the objective function value (population fit)
//' being associated to each individual chromosome from above.
//' @return A list contains the chromosomes for \code{dad} and \code{mom}.
//' @export
// [[Rcpp::export]]
List selection_linearrank(arma::mat& pop, arma::vec& popFit){

  popsize = popFit.size();
  arma::vec myorder = arma::conv_to<arma::vec>::from(arma::sort_index(arma::sort_index(popFit))+1);
  arma::vec selectrank = popsize - myorder;
  arma::vec probs0 = 2.0*selectrank/(popsize*(popsize-1));

  arma::vec popID = linspace(0, popsize-1, popsize);
  arma::vec parentI = Rcpp::RcppArmadillo::sample(popID, 2, FALSE, probs0);

  arma::vec dad = pop.col(parentI(0));
  double dadfit = popFit(parentI(0));
  arma::vec mom = pop.col(parentI(1));
  double momfit = popFit(parentI(1));

  // dad has (better fit/smaller value/larger rank) than mom
  if(dadfit > momfit){
    arma::vec tmp = dad;
    double tmpfit = dadfit;
    dad = mom;
    dadfit = momfit;
    mom = tmp;
    momfit = tmpfit;
  }

  List res;
  res["dad"] = dad;
  res["mom"] = mom;
  return(res);
}
