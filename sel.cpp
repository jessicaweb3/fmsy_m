#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(age);
  DATA_VECTOR(sel);
  
  PARAMETER(s50);
  PARAMETER(ss);
  PARAMETER(logsigma);
  
  vector<Type> pred = 1/(1+exp(-((age-s50)/ss)));
  Type nll = -sum(dnorm(sel, pred, exp(logsigma), true));
  ADREPORT(pred);
  return nll;
}