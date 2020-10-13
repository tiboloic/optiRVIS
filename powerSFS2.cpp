// LT 13/12/2019

// multinomial powerlaw with random effect on slope
// for detection of change of slope among genes as an intolerance metric

// version 2 using padé approximant of power law probabilities

#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  
  // DATA
  
  // observations: counts per category per gene
  DATA_MATRIX(obs);
  
  // Padé approximation of power probabilities
  DATA_MATRIX(pade);
  
  // average slope beta
  PARAMETER(b);
  
  //log-stdev of random effect on slope
  PARAMETER(lsigb);
  
  // random effects on slope
  PARAMETER_VECTOR (bdevs);
  
  
  //Type nLL=0.0;
  
  Type beta;

  parallel_accumulator<Type> nll(this);
  
  // add likelihood contribution for gene level random effects
  nll-= dnorm(bdevs, Type(0), Type(1), true).sum();
  
  // likelihood contribution per gene
  for(int i=0;i<obs.rows();i++){
    

    // build beta
    beta = b + exp(lsigb) * bdevs(i);
    
    // vector of probabilities
    vector<Type> ps(obs.cols());
    
    for (int j=0; j<obs.cols(); j++) {
      ps(j) = pade(j, 0) + pade(j, 1) * beta + pade(j, 2) * pow(beta, 2) + pade(j, 3) * pow(beta, 3);
      ps(j) = ps(j) / (1 + pade(j, 4) * beta + pade(j, 5) * pow(beta, 2) + pade(j, 6) * pow(beta, 3));
    }
    
    Type psum = ps.sum();
    ps = ps / psum;
    
    vector<Type> counts = obs.row(i);
    
    // multinomial sampling
    nll -= dmultinom(counts, ps, true);
    
  }
  return nll;
}
