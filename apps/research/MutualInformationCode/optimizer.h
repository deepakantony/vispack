/*

optimizer.h

Sarah Geneser
written: 06-25-03
last updated: 07-25-03

*/

#include <matrix.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
using namespace std;


//#define DEBUG 1
#define EPS 1e-14
#define EPSDecrease 1e-7
#define numSlowDecreases 4

class Optimizer
{
 public:
  Optimizer() {;}
};

template <class M>
class GradientDescent: public Optimizer
{
 private:
  float _lambda;
  int _iter;
  int _endIter;
  
 public:
  GradientDescent() 
  {
    _lambda = 0.01; _iter=0; _endIter = 1000;
  }
  virtual void setLambda(float lambda) {_lambda=lambda;}
  virtual void setEndIter(int endIter) {_endIter=endIter;}
  virtual int getIter() {return _iter;}
  virtual float getLambda() {return _lambda;}
  virtual void restart() {_iter=0;}
  boolean notFinished() {return _iter<_endIter;}
  virtual void print()
  {
    cout << "Gradient Descent Parameters:" << endl;
    cout << "number of iterations = " << _endIter << endl;
    cout << "lambda = " << _lambda << endl;

  }
  void optimize(M* metric);
};

template <class M>
class Marquardt: public Optimizer
{
 private:
  int _iter;
  int _endIter;
  float _lambda;
  float _lambdaMin;
  float _lambdaMax;
  int _slowDecrease;
  boolean _done;
  float _bestValue;
  VISVector _oldParams;

 public:
  Marquardt() 
  {
    _iter=0;
    _lambda=1;_lambdaMin=1e-7;_lambdaMax=1e13;
    _slowDecrease = 0;
    _bestValue = 1e16;
    _done = false;
  }
  virtual void setLambda(float lambda) {_lambda=lambda;}
  virtual void setLambdaMin(float lambda) {_lambdaMin=lambda;}
  virtual void setLambdaMax(float lambda) {_lambdaMax=lambda;}
  virtual void setEndIter(int endIter) {_endIter=endIter;}
  virtual void restart() {_done=false;_iter=0;_slowDecrease=0;_bestValue=1e16;}
  virtual int getIter() {return _iter;}
  virtual float getLambda() {return _lambda;}
  virtual float getBestValue(){return _bestValue;}
  boolean notFinished() {return !_done;}
  virtual void print()
  {
    cout << "Marquardt Parameters:" << endl;
    cout << "lambda = " << _lambda << endl;
    cout << "number of iterations = " << _endIter << endl << endl;
  }

  void optimize(M* metric);
 
};


/////////////////////////////////////////////
// gradient descent method of maximization //
/////////////////////////////////////////////
template <class M>                           
void GradientDescent<M>::optimize(M* metric) 
{
  float old_value, new_value;
  old_value = metric->getValue();
  (metric->getTransform())->setParams((metric->getTransform())->getParams() \
				      + _lambda*(metric->getDeriv()).getParams());
  metric->getMetric();
  new_value = metric->getValue();
  _iter++;
}

////////////////////////////
// Marquardt maximization //
////////////////////////////
template <class M>                    
void Marquardt<M>::optimize(M* metric) 
{

  VISVector new_params;
  int i, j;
  
  _oldParams = (metric->getTransform())->getParams();
  int n = _oldParams.n();
  VISMatrix hessian(n, n), ident = VISIdentity(n);
  VISVector derivative(n), difference(n);
  
  float old_value, new_value;
  
  //NOTE: the metric *must* have been previously initialized

  old_value = (metric->getValue()); 
  if (_iter == 0) 
    _bestValue = old_value;
  derivative = (metric->getDeriv()).getParams(); 
  hessian = (metric->getHessian()); 

  //set the new transform in the metric
  (metric->getTransform())->setParams(new_params = _oldParams 
				      + solvLinearEqn(-hessian+_lambda*ident, derivative));

  //compute the new metric values and get the new value of the metric
  metric->getMetric();  new_value = (metric->getValue());
  
  //if the metric is greater than the previous best value
  //else if((new_value-_bestValue) > EPS)
  if(new_value > _bestValue)
  {
    _bestValue = new_value;
    _lambda =VISmax(_lambda/2.0f, _lambdaMin);
    _oldParams = new_params;
    //check to see if there has been a slow decrease
    if ( (new_value-old_value) < EPSDecrease)
      _slowDecrease++;
  }
  //if the metric increased, but not above the best value
  //or if it decreased
  else
  {
    _lambda = VISmin(_lambda*2.0f, _lambdaMax);
    //set the transform back to the old transform
    (metric->getTransform())->setParams(_oldParams);
    if (_lambda >= _lambdaMax)
    {
      _done = true;
      cout << "reached lamMax." << endl;
    }  

    if((new_value-EPS) > _bestValue)
      _bestValue = new_value;
  }

  //increment the iteration counter
  _iter++;
  //check to see if reached end iteration
  if (_iter>=_endIter)
  {
    _done = true;
    cout << "reached end iter." << endl;
  }

  //check to see if there have been too many slow decreases
  //and if so, exit
  if (_slowDecrease >= numSlowDecreases)
  {
    _done = true;
    cout << "halting due to stagnation." << endl;
  }  


}
