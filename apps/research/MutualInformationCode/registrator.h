/*
  registrator.h
  
  Sarah Geneser
  written: 06-08-03
  last updated: 07-25-03
  
  takes a metric and transform and optimizer and registrates two images
*/


#include "metric.h"
#include "optimizer.h"

#include <math.h>
#include "image/image.h"
#include "image/imagefile.h"
#include "util/mathutil.h"
#include <matrix.h>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
using namespace std;
#include <volume.h>
#include <volumefile.h>
#include <imagefile.h>

using namespace std;


#define DEBUG 1
#define BLUR .5
//#define graphs 1

template <class M, class T, class O>
class Registrator
{
 protected:
  int _EndIter;

  // weighting the parameters of the transform to restrict or encourage movement 
  // in each parameter space
  VISVector _paramWeights;

  //for multi-resolution registration
  int _numModes;

  //test and reference images
  int _refW, _refH, _refD;


  //what type of metric to use
  M* _metric;
  //the value of the metric at the current transform
  float _metricValue;
  //the current transform
  T* _transform;
  //the derivative of the metric w.r.t. the transform
  T _derivative;
  //the second derivative (Hessian) of the metric w.r.t the transform
  VISMatrix _hessian;
  // what type of optimizer to use 
  O* _optimizer;
    
  // the largest value of MI that has been found
  float _bestMetricValue;   
  // the transform that corresponds to the larges MI
  T* _bestTransform;

  // for region of interest _upperROI contains (minX, minY) and 
  // _lowerROI contains (maxX, maxY) in IMAGE CENTERED COORDINATES
  VISVector _upperROI, _lowerROI;
  

 public:
  Registrator() 
  {
    _bestMetricValue = -9999999999999;
    _upperROI = VISVector(3); 
    _upperROI = 0.0f;
    _lowerROI = VISVector(3);
    _lowerROI = 0.0f;
  }

  // the following functions can be used to set
  // specific registration routine parameters

  virtual void setMetric(M* metric) {_metric = metric;}
  virtual void setInitialTransform(T* transform)
  {
    _transform = transform;
    _bestTransform = new T();
    int n = (_transform->getParams()).n();
    _hessian = VISMatrix(n,n);
    //give the metric access to the transform
    _metric->setTransform(_transform);
  }
  
  // does some normalization preprocessing to the reference 
  // volume and sets the volume to the metric
  virtual void setReferenceVolume(VISVolume<float> ref)
  {
    //make sure pixel values range from 0 to 1
    ref -= ref.min();
    ref /= ref.max();
    _metric->setReferenceVolume(ref);
    _refW = ref.width(); _refH = ref.height(); _refD = ref.depth();

    //update the values for the region of interest
    //assume the entire image as the ROI
    _upperROI = 0.0f;
    _lowerROI.poke(0) = ref.width()-1; 
    _lowerROI.poke(1) = ref.height()-1;
    _lowerROI.poke(2) = ref.depth()-1;
        
    //transform them to image centered coords
    _upperROI = volumeToCentered(_upperROI, _refW, _refH, _refD);
    _lowerROI = volumeToCentered(_lowerROI, _refW, _refH, _refD);
  }  

  // does some normalization preprocessing to the test 
  // volume and sets the metric test volume
  virtual void setTestVolume(VISVolume<float> test)
  {
    // make sure the pixel values range from 0 to 1
    test -= test.min();
    test /= test.max();
    _metric->setTestVolume(test);
  }


  // takes the upper and lower corner of the ROI in 
  // image coordinates and transforms them to 
  // image centered coordinates
  // NOTE: ROI must be set *after* the test and ref
  // images
  virtual void setROI(VISVector x, VISVector y)
  {
    if (x.n()==3)
      _upperROI = volumeToCentered(x, _refW, _refH, _refD);
    else 
      cout << "WARNING: you have specified an unacceptable " \
	   << "x vector for the region of interest" << endl;
    if (y.n()==3)
      _lowerROI = volumeToCentered(y, _refW, _refH, _refD);
    else
      cout << "WARNING: you have specified an unacceptable " \
	   << "y vector for the region of interest" << endl;
    _metric->setROI(x, y);
  }

  virtual void setOptimizer(O* optimizer)
  {
    _optimizer = optimizer;
  }
  virtual void setIters(int iters){_optimizer->setEndIter(iters);} 
  virtual void setStepSize(float lambda)
  {
    _optimizer->setLambda(lambda);
  }
  virtual float getBestMetric(){return _bestMetricValue;}
  virtual T* getBestTransform(){return _bestTransform;}
  
  virtual void print()
  {
    cout << "Registration Paramters:" << endl;

  }

  T* registrate();
  

};


////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////
// Completes the registration process with the //
// parameters that have been set               //
/////////////////////////////////////////////////
template <class M, class T, class O>
T* Registrator<M, T, O>::registrate()
{
  //print out the parameters of the registrator
  print();
  //print out the metric parameters
  _metric->print();
  //print out the optimizer parameters
  _optimizer->print();
  
  //compute the metric values
  _metric->getMetric();
  
  //write out the information for the current iteration
  cout << "Iteration " << 0 << endl;
  cout << "Metric = " << _metric->getValue() << endl;
  cout << "initial transform : " << endl;
  _transform->print();
  cout << endl << endl;

    
#ifdef graphs
  //write metricGraph the metric information for each iteration to file
  ofstream metricGraph, lambdaGraph;
  metricGraph.open("/scratch/research/MutualInformationCode/graphGradDescent/metric.txt");
  lambdaGraph.open("/scratch/research/MutualInformationCode/graphGradDescent/lambda.txt");
  metricGraph << _metric->getValue() << endl;
  lambdaGraph << _optimizer->getLambda() << endl;
#endif

  while (_optimizer->notFinished())
  { 
    //use the optimizer to determine the next transform
    _optimizer->optimize(_metric);
    _metricValue = _metric->getValue();
    //check to see if the current transform is the best so far
    //update the best variables
    if (_metricValue > _bestMetricValue)
    {
      _bestTransform->setParams(_transform->getParams());
      _bestMetricValue = _metricValue;
    }    
        
    //write out the information for the current iteration
    cout << "Iteration " << _optimizer->getIter() << endl;
    cout << "Metric = " << _metric->getValue() << endl;
    cout << "New transform: " << endl;
    _transform->print();
    cout << endl << endl;

#ifdef graphs
    metricGraph << _metric->getValue() << endl;
    lambdaGraph << _optimizer->getLambda() << endl;
#endif

  }
#ifdef graphs
  metricGraph.close(); lambdaGraph.close();
#endif
  print();
  cout << "best transform = " << endl;
  _bestTransform->print();
  cout << "best metric = " << _bestMetricValue << endl;
  return _transform;
}
