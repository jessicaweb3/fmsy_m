#include <TMB.hpp>
template<class Type>
  Type objective_function<Type>::operator() ()
{
 
  DATA_VECTOR(age);
  DATA_VECTOR(mprop);
  
  PARAMETER(a50);
  PARAMETER(env);
  PARAMETER(logsigma);
  
  vector<Type> pred = 1/(1+exp(-((age-a50)/env)));
  Type nll = -sum(dnorm(mprop, pred, exp(logsigma), true));
  ADREPORT(pred);
   return nll;
  }


