#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{

DATA_VECTOR(SSB);
DATA_VECTOR(Rec);
  
  PARAMETER(alpha);
  PARAMETER(beta);
  PARAMETER(logsigma);
  
  vector<Type> pred = alpha*SSB*exp(beta*SSB);
  Type nll = -sum(dnorm(Rec, pred, exp(logsigma), true));
  ADREPORT(pred);
return nll;
}
