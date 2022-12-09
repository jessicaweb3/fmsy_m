#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(ssb);
  DATA_VECTOR(logR);

  PARAMETER(loga);
  PARAMETER(logb);
  PARAMETER(logsigma);
  Type sigma=exp(logsigma);
//  Type a=exp(loga);
  vector<Type> pred = loga+log(ssb)-log(Type(1)
                      +exp(logb)*ssb);
  Type nll = -sum(dnorm(logR,pred,sigma,true));
  ADREPORT(sigma);
//  ADREPORT(a);
  return nll;
}

