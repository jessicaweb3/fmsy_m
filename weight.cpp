#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  
  DATA_VECTOR(age);
  DATA_VECTOR(wprop);
  
  PARAMETER(k);
  PARAMETER(b);
  PARAMETER(Winf);
  PARAMETER(logsigma);
  
  vector<Type> pred = Winf * pow((1-exp(-k * age)), b);
  Type nll = -sum(dnorm(wprop, pred, exp(logsigma), true));
  ADREPORT(pred);
  return nll;
}



