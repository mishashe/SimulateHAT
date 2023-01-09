#include <Rcpp.h>
#include <iostream>
#include <random>
#include <chrono>
#include <math.h>
#include <stdint.h>

// [[Rcpp::depends(RcppArmadillo)]]
// #define ARMA_WARN_LEVEL 0
// #include <RcppArmadillo.h>
using namespace Rcpp;
// using namespace arma;
using namespace std;

// [[Rcpp::export]]
std::map<int,double> getMLD(int& L, int& K1, double& rho1, double& rho2, IntegerVector& times, NumericVector& minus_log_runif_L, int& nT, double& dt)
{
  std::map<int,double> rTable;
  double tmax=(double)(*max_element(times.begin(),times.end()));
  double tmin=0.1*(double)(*min_element(minus_log_runif_L.begin(),minus_log_runif_L.end()));
  times(L-1) = 10000.0;

  minus_log_runif_L(L-1) = -10000000.0;
  NumericVector t2B(nT+1); for (int i = 0; i < nT+1; ++i) t2B(i) = pow(10.0,log10(tmin)+(1.0*i)/(1.0*nT)*(log10(tmax)-log10(tmin)));
  NumericVector t2V(nT); for (int i = 0; i < nT; ++i) t2V(i) = sqrt(t2B(i)*t2B(i+1));
  NumericVector P(nT); for (int i = 0; i < nT; ++i) P(i) = (t2B(i+1)-t2B(i))*dt*rho2*exp(-t2V(i)*dt*rho2);
  double P0 = exp(-rho2*tmax*dt);


  int rprev0 = -1;
  IntegerVector rprev(nT); for (int i = 0; i < nT; ++i) rprev(i) = -1;

  int tseg = times(0);
  int tsegloc = 0;
  int imin;
  for (int i = 0; i < L; ++i)
  {
    tseg = times(i);
    // if (i-tsegloc > K1-1)
    // {
    //   imin=std::distance(times.begin(), min_element(times.begin()+i-K1+1,times.begin()+i));
    //   tseg = times(imin);
    //   tsegloc = imin;
    // }
    // else if (times(i)<tseg)
    // {
    //   tseg = times(i);
    //   tsegloc = i;
    // }
    if (minus_log_runif_L(i) < tseg)
    {
      rTable[i-rprev0-1] += P0;
      rprev0 = i;
      bool keeploop = true;
      for (int j = nT-1; j >=0 && keeploop; j--)
      {
        if (minus_log_runif_L(i) < t2V(j))
        {
          rTable[i-rprev(j)-1] += P(j);
          rprev(j) = i;
        }
        else
        {
          keeploop=false;
        }
      }
    }
  }
  return(rTable);
}

// [[Rcpp::export]]
IntegerVector SimulateOnce(long int &L, long int& T, long int& K, long int seed)
{
  std::random_device dev;
  std::mt19937 rng(seed);
  std::uniform_int_distribution<std::mt19937::result_type> distL(0,L-K);

  IntegerVector Genome(L);

  long int loc;
  for (long int t = 0; t < T; t++)
  {
    loc = distL(rng);
    long int kmax = loc + K - 1;
    for (long int k = loc; k <= kmax; k++) Genome(k) = t;
  }
  return(Genome);
}

// [[Rcpp::export]]
IntegerVector SimulateIndirect(long int &L, long int& T, long int& K, long int seed)
{
  std::random_device dev;
  std::mt19937 rng(seed);
  std::uniform_int_distribution<std::mt19937::result_type> distL(0,L-K);

  IntegerVector Genome12(L);
  IntegerVector Genome13(L);
  IntegerVector Genome23(L);

  long int loc;
  long int ran;

  for (long int t = 0; t < T; t++)
  {
    loc = distL(rng);
    long int kmax = loc + K - 1;
    ran = distL(rng) %4;
    if (ran==0) for (long int k = loc; k <= kmax; k++) {Genome13(k) = t; Genome23(k) = Genome12(k);} // 1 to 3
    else if (ran==1) for (long int k = loc; k <= kmax; k++) {Genome13(k) = t; Genome12(k) = Genome23(k);} // 3 to 1
    else if (ran==2) for (long int k = loc; k <= kmax; k++) {Genome23(k) = t; Genome13(k) = Genome12(k);} // 2 to 3
    else if (ran==3) for (long int k = loc; k <= kmax; k++) {Genome23(k) = t; Genome12(k) = Genome13(k);} // 3 to 2
  }
  return(Genome12);
}
