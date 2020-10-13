// LT 19/08/2020

// version for variable allele numbers
// use asymptotic approximation 

#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  
  // DATA
  
  // observations: counts per category per gene
  DATA_MATRIX(obs);
  
  // Average allele number for each gene
  DATA_VECTOR(an);
  
  // average slope beta
  PARAMETER(b);
  
  //log-stdev of random effect on slope
  PARAMETER(lsigb);
  
  // random effects on slope
  PARAMETER_VECTOR (bdevs);
  
  // Bernouilli polynomials at 0
  const Type Bs[] = {Type(1/6), Type(-1/30), Type(1/42), Type(-1/30),
                     Type(5/66), Type(-691/2730), Type(7/6), Type(-3617/510)};
    
  Type nLL=0.0;
  
  Type beta;
  
 /* Type zeta(Type n, Type beta) {
    // sum_{i=n}^\infty 1/i^beta
    // using asymptotic approximation
    Type ans = pow(n,1-beta)/(beta  - 1) + 0.5 / pow(n, beta);
    Type term = beta/2 / pow(n,beta-1);
    ans = ans + term * Bs[1];
    
    for (int i=2;i<9;i++) {
      term = term * (beta + 2*i-3) * (beta + 2*i -2) / ((2*i -1) * 2*i * pow(n, 2));
      ans = ans + term * Bs[i];
    }
    return(ans);
  }
*/
  //parallel_accumulator<Type> nll(this);
  
  // add likelihood contribution for gene level random effects
  nll-= dnorm(bdevs, Type(0), Type(1), true).sum();
  
  // likelihood contribution per gene
  for(int i=0;i<obs.rows();i++){
    

    // build beta
    beta = b + exp(lsigb) * bdevs(i);
    
    // vector of probabilities
    vector<Type> ps(obs.cols(), Type(0.0));
    
    // initialize to 0
    //for (int j=0; j<obs.cols(); j++) {
    //  ps(j) = Type(0.0);
    //}
    
    // calculate the first terms up to 15
    for (int j=0; j<15; j++) {
      ps(floor(log(j+1)/log(2))) +=  Type(pow(j+1, -beta));
    }
    
    // use asympotic approximation from 16 to AN. 16 = 2^4
    // we iterate on the power of 2
    for (int j=4; j<floor(log(an(i))/log(2)) + 1) {
      ps(j) = zeta(pow(2,j), beta);
    }

    for (int j=4; j<floor(log(an(i))/log(2))) {
      ps(j) = ps(j) - ps(j+1);
    }
    ps(floor(log(an(i))/log(2))) = ps(floor(log(an(i))/log(2))) - zeta(an(i), beta);

    Type psum = ps.sum();
    ps = ps / psum;
    
    vector<Type> counts = obs.row(i);
    
    // multinomial sampling
    nll -= dmultinom(counts, ps, true);
    
  }
  return nll;
}
