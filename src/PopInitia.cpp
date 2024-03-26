// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>

using namespace Rcpp;
using namespace arma;

int m, i, j, id, jp;
int popsize, islandSize;
int lmax;

//' @export
// [[Rcpp::export]]
IntegerVector rank_asR(NumericVector x, bool decreasing = false)
{
  IntegerVector rank = match(x, clone(x).sort());
  if(decreasing) rank = rank.length()+ 1 - rank;
  return rank;
}

//' @export
// [[Rcpp::export]]
arma::vec selectTau_cpp(int n, int minDist, double Pb, int mmax, int lmax){

  m = 0;
  double a;
  arma::vec tau(lmax, fill::zeros);
  i = 1 + minDist;

  while(i < n - minDist){
    a = runif(1)[0];
    if(a <= Pb){
      m = m + 1;
      tau(m) = i;
      i = i + minDist;
    }else{
      i = i+1;
    }
  }
  tau(0) = m;
  tau(m+1) = n+1;
  return(tau);
}

//' @export
// [[Rcpp::export]]
arma::mat random_population_cpp(int popsize, int n, int minDist, double Pb, int mmax, int lmax){

  arma::mat pop(lmax, popsize, fill::zeros);

  for(j=0;j<popsize;j++){
    pop.col(j) = selectTau_cpp(n, minDist, Pb, mmax, lmax);
  }

  return(pop);
}

// //' @export
// // [[Rcpp::export]]
// void InitialPopCpp(arma::cube& Island, int n, int minDist, double Pb, int mmax){
//
//   lmax = size(Island)[0];
//   popsize = size(Island)[1];
//   islandSize = size(Island)[2];
//
//   for(id=0;id<islandSize;id++){
//     for(jp=0;jp<popsize;jp++){
//       if(lmax > mmax+2){
//         // TBD for other parameters except changepoint
//         Island.slice(id).col(jp).subvec(0, mmax+1) = SelectTauCpp(n, minDist, Pb, mmax);
//       }else{
//         Island.slice(id).col(jp).subvec(0, mmax+1) = SelectTauCpp(n, minDist, Pb, mmax);
//       }
//     }
//   }
//
// }

//' @export
// [[Rcpp::export]]
arma::vec offspring_uniformcrossover_cpp(arma::vec& mom, arma::vec& dad, int minDist, int lmax, int n){
  //  In uniform crossover, typically, each bit is chosen from either parent with equal probability.
  //  Other mixing ratios are sometimes used, resulting in offspring which inherit more genetic
  // information from one parent than the other. In a uniform crossover, we don’t divide the
  // chromosome into segments, rather we treat each gene separately. In this, we essentially
  // flip a coin for each chromosome to decide whether or not it will be included in the off-spring.

  // Rcout << "mom:" << mom << std::endl;
  // Rcout << "dad:" << dad << std::endl;

  arma::vec child(lmax, fill::zeros);

  int child_mmax = dad(0) + mom(0);
  if(child_mmax == 0){
      // no changepoint for both mom and dad
      child(0) = 0;
      child(1) = n+1;
    }else{
      double a;
      // 1). combine dad's and mom's changepoints
      int mom_m = mom(0);
      int dad_m = dad(0);
      arma::vec parents(child_mmax, fill::zeros);
      if(mom_m > 0){
        parents.subvec(0,mom_m-1) = mom.subvec(1,mom_m);
      }
      if(dad_m >0){
        parents.subvec(mom_m,child_mmax-1) = dad.subvec(1,dad_m);
      }

      // 2). remove same changepoints and sort
      parents = sort(unique(parents));
      int parents_count = parents.size();

      // 3). select subset from parents
      int child_m = 0;
      int tmpm = 1;
      int tmptau = parents(0);
      int temp;
      do{
        a = runif(1)[0];
        if(a > 0.5){
          child_m = child_m + 1;
          child(child_m) = tmptau;
          tmptau = tmptau + 2*minDist;
          if(tmptau >= n-minDist){break;} // boundary changepoints limits
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
        if(tmptau >= n-minDist){break;}
      }
      while(1<2);

      child(0) = child_m;
      child(child_m+1) = n+1;
    }

  return(child);
}

//' @export
// [[Rcpp::export]]
List selection_linearrank_cpp(arma::mat& pop, arma::vec& popFit){

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
