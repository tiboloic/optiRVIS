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
      ps(j) = Type(0.0);
    }
    vector<Type> counts = obs.row(i);
    
    // buildup vector of probabilities by aggregating power law over bins
    // should I work on the log scale to avoid underflow ?
    //for (int j=0; j<251496; j++) {
    
    // try vector operation
    //std::vector<int> is(251496);
    //std::iota(is.begin(), is.end(), 1);
    //vector<Type> powtest = Type(pow(is, -beta));
    
    for (int j=0; j<251496; j++) {
      ps(floor(log(j+1)/log(2))) +=  Type(pow(j+1, -beta));
    }
    
    //REPORT(beta)
    //REPORT(ps);
    vector<Type> ps2 = ps / ps.sum();
    //REPORT(ps2);
    
    // multinomial sampling
    nll -= dmultinom(counts, ps2, true);
    
  }
  return nll;
}
