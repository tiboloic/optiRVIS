// LT 13/12/2019

// multinomial powerlaw with random effect on slope
// for detection of change of slope among genes as an intolerance metric


#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  
  // DATA
  
  // observations: counts per category per gene
  DATA_MATRIX(obs);
  
  // average slope beta
  PARAMETER(b);
  
  //log-stdev of random effect on slope
  PARAMETER(lsigb);
  
  // random effects on slope
  PARAMETER_VECTOR (bdevs);
  
  
  Type nLL=0.0;
  
  Type beta;

  
  // add likelihood contribution for gene level random effects
  nLL -= dnorm(bdevs, Type(0), Type(1), true).sum();
  

  // likelihood contribution per gene
  for(int i=0;i<obs.rows();i++){
    

    // build beta
    beta = b + exp(lsigb) * bdevs(i);
    
    // vector of probabilities
    vector<Type> ps(obs.cols());
    for (int j=0; j<obs.cols(); j++) {
      ps(j) = Type(pow(j+1, -beta));
    }
    vector<Type> counts = obs.row(i);
    
    vector<Type> ps2 = ps / ps.sum();
    //REPORT(ps2);
    
    // multinomial sampling
    nLL -= dmultinom(counts, ps2, true);
    
  }
  return nLL;
}
