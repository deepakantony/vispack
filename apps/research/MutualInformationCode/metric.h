/*
  metric.h

  Sarah Geneser
  written: 06-08-03
  last updated: 07-23-03
  
*/

#include <math.h>
#include "image/image.h"
#include "image/imagefile.h"
#include "util/mathutil.h"
#include <matrix.h>
#include <stdio.h>
#include <stdlib.h>
#include <volume.h>
#include <volumefile.h>
#include <imagefile.h>
#include "transform.h"


//#define WRITEDERIVS
//#define WRITENORM
//#define WRITERAND
#define DEBUG 1
//#define WARNING 1


#define BLUR .5
#define SLICE 0

#define TWOPI 6.283185307179586
#define eps .000000001

//#define moreAccurateHessian 1

class Metric
{
 protected:

  // current value of the metric
  float _metric;

  // test and reference images and derivatives in the x, y, and z
  // directions of the test image
  VISVolume<float> _ref, _dxRef, _dyRef, _dzRef;
  VISVolume<float> _test, _dxTest, _dyTest, _dzTest;
  
  // for region of interest _upperROI contains (maxX, maxY) and 
  // _lowerROI contains (minX, minY) in IMAGE CENTERED COORDINATES
  VISVector _upperROI, _lowerROI;

  // width, height, and depth of the reference image
  int _refW, _refH, _refD;

  
 public:
  Metric() : _metric( 0.0f ) {} 
  virtual ~Metric() {;}
  virtual VISVolume<float> processVolume(VISVolume<float> volume)=0;

  virtual void setReferenceVolume(VISVolume<float> ref)
  {
    _ref = processVolume(ref);

    //update the values for width, height, and depth
    _refW = ref.width(); _refH = _ref.height(); _refD = _ref.depth();
    _dxRef = -1.0f*ref.dx(); 
    _dyRef = 1.0f*ref.dy();
    _dzRef = ref.dz();
    
    //update the values for the region of interest
    //assume the entire image as the ROI
    _upperROI = 0.0f;
    _lowerROI.poke(0) = _refW-1; 
    _lowerROI.poke(1) = _refH-1; 
    _lowerROI.poke(2) = _refD-1;    
    //transform them to image centered coords
    _upperROI = volumeToCentered(_upperROI, _refW, _refH, _refD);
    _lowerROI = volumeToCentered(_lowerROI, _refW, _refH, _refD);
  }

  virtual void setTestVolume(VISVolume<float> test)
  {
    _test = processVolume(test);

    //update values for the derivatives
    _dxTest = -1.0f*test.dx(); 
    _dyTest = 1.0f*test.dy();
    _dzTest = test.dz();
  }
  
  // takes the upper and lower corner of the ROI in 
  // image coordinates and transforms them to 
  // image centered coordinates
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
  }

  virtual VISVolume<float> getReferenceVolume(){return _ref;}
  virtual VISVolume<float> getTestVolume() {return _test;}
  virtual VISVolume<float> getDxTestVolume() {return _dxTest;}  
  virtual VISVolume<float> getDyTestVolume() {return _dyTest;}
  virtual VISVolume<float> getDzTestVolume() {return _dzTest;}
  
  virtual VISVector getLowerROI() {return _lowerROI;}
  virtual VISVector getUpperROI() {return _upperROI;}

  virtual float getValue() {return _metric;}
};

//////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////
// Correlation metric computed for two images and a transform //
// that takes one image into the other image's coordinate     //
// system                                                     //
////////////////////////////////////////////////////////////////
template <class T>
class Correlation: public Metric
{
 protected:

  // the current transform
  T* _transform;
  //the number of transform parameters
  int _numParams;

  public:

  Correlation()
  {  
    _numParams = 0;
    _upperROI = VISVector(3);
    _lowerROI = VISVector(3);
   }

  virtual ~Correlation() {;}
  
  virtual void setTransform(T* transform) 
  {
    _transform = transform;
    _numParams = (_transform->getParams()).n();
  }
  
  virtual T* getTransform() {return _transform;}
  
  virtual VISVolume<float> processVolume(VISVolume<float> volume)
  {
#ifdef BLUR
    volume = volume.gauss(BLUR);
#endif
    volume = volume - volume.average();
    volume /= pow((float)volume.variance(), 0.5f);
    return volume;
  }

  // prints out the parameters of the registration process
  virtual void print()
  {
    cout << "Correlation Metric:" << endl;
  }
  virtual void printLong()
  {
    print();
    cout << "Region of interest = (" << _lowerROI.peek(0) << \
      ", " << _upperROI.peek(0) << ", " << _lowerROI.peek(1) << \
      ", " << _upperROI.peek(1) << ")" << endl;
    cout << "transform = " << endl;
    _transform->print();
  }
  void getMetric();
};

template <class T>
void Correlation<T>::getMetric()
{
  VISVolume<float> testTransformed, metric;
  testTransformed = _transform->transformVolume(_test);
  metric = _ref*testTransformed;
  _metric = metric.sum();
}


///////////////////////////////////////////////////////////////////////
// Mutual information metric computed for two images and a transform //
// that takes one image into the other image's coordinate system     //
///////////////////////////////////////////////////////////////////////
template <class T>                         
class MutualInformation: public Metric
{
 protected:
  // mutual information parameters
  int _nA, _nB;
  // covariance values
  VISMatrix _Phi_uv, _Phi_uvINV;
  float _Phi_v, _Phi_vINV, _Phi_vvINV;  
  // weighting the parameters of the transform to restrict or encourage movement 
  // in each parameter space
  VISVector _paramWeights;
  
  // the current transform
  T* _transform;
  //the number of transform parameters
  int _numParams;
  // the derivative of mutual information
  T _derivative;
  // the second derivative (Hessian) of mutual information
  VISMatrix _hessian;

  //store a matrix of random points (ref image space)
  VISMatrix _xis, _xjs;
  //store the transformation of the random points (test image space)
  VISMatrix _Txis, _Txjs;
  //store the pixel values of random points in test and ref image space
  VISVector _Vxis, _Vxjs, _Uxis, _Uxjs;
  //store the dx and dy's in the test ima2ge space
  VISVector _dxTestXis, _dxTestXjs, _dyTestXis, _dyTestXjs, _dzTestXis, _dzTestXjs;
  //determines if create new random points each iteration
  boolean _newRands;
  
  void generateRandomPoints(); 
  float Gphi_v(float x);
  float Gphi_uv(VISMatrix x);
  float W_v(float x, float y);
  float W_uv(VISMatrix x, VISMatrix y);
  
 public:

  MutualInformation()
  {  
    int _nA=50, _nB=50;
    VISMatrix _Phi_uv, _Phi_uvINV;
    VISVector _paramWeights;
    float _Phi_v, _Phi_vINV, _Phi_vvINV;  
    _numParams = 0;
    _upperROI = VISVector(3);
    _lowerROI = VISVector(3);
    _newRands = true;
    _derivative.setParams(0.0f);
    _hessian = 0.0f;
  }

  virtual ~MutualInformation() {;}

  virtual void setNA(int nA)
  {
    _nA = nA;
    VISMatrix temp(_nA,3);
    _xis=temp; _xis = 1.0f;
    _Txis = _xis;
    VISVector temp2(_nA);
    temp2 = 0.0f;
    _Vxis = temp2;
    _Uxis = temp2;
    _dxTestXis = temp2;
    _dyTestXis = temp2;
    _dzTestXis = temp2;
  }

  virtual void setNB(int nB)
  {
    _nB = nB;
    VISMatrix temp(_nB,3);
    _xjs = temp; _xjs = 1.0f;
    _Txjs = _xjs;
    VISMatrix temp2(_nB); temp2 = 0.0f;
    _Vxjs = temp2;
    _Uxjs = temp2;
    _dxTestXjs = temp2;
    _dyTestXjs = temp2;
    _dzTestXjs = temp2;
  }

  virtual void setPhi_uv(VISMatrix Phi_uv)
  {
    _Phi_uv = Phi_uv;
    _Phi_uvINV = _Phi_uv.inv();
    _Phi_vvINV = _Phi_uvINV.peek(1,1);
  }

  virtual void setPhi_v(float Phi_v)
  {
    _Phi_v = Phi_v;
    _Phi_vINV = 1/Phi_v;
  }

  virtual void setParamWeights(VISVector paramWeights)
  {
    if (paramWeights.n() == _numParams)
      _paramWeights = paramWeights;
    else 
      cout << "incorrect size of param weights" << endl;
  }

  virtual void setTransform(T* transform)
  {
    _transform = transform;
    _numParams = (_transform->getParams()).n();
    _hessian = VISMatrix(_numParams, _numParams);
  }
  
  virtual VISVolume<float> processVolume(VISVolume<float> volume)
  {
#ifdef BLUR
    volume = volume.gauss(BLUR);
#endif
    volume = volume - volume.average();
    volume /= pow((float)volume.variance(), 0.5f);
    return volume;
  }

  // the following functions return the values of the requested
  // registration routine parameters
  virtual int getNA(){return _nA;}
  virtual VISMatrix getPhi_uv() {return _Phi_uv;}
  virtual VISMatrix getPhi_uvINV() {return _Phi_uvINV;}
  virtual float getPhi_v() {return _Phi_v;}
  virtual float getPhi_vINV() {return _Phi_vINV;}

  virtual T* getTransform() {return _transform;}
  virtual T getDeriv() {return _derivative;}
  virtual VISMatrix getHessian() {return _hessian;}
  
  // prints out the parameters of the registration process
  virtual void print()
  {
    cout << "Mutual Information Metric Parameters:" << endl;
    cout << "samples for set A = " << _nA <<endl;
    cout << "samples for set B = " << _nB << endl;
    cout << "1D covariance = " << _Phi_v << endl;
    cout << "2D covariance Matrix = " << _Phi_uv ;
  }
  virtual void printLong()
  {
    print();
    cout << "Region of interest = (" << _lowerROI.peek(0) << \
      ", " << _upperROI.peek(0) << ", " << _lowerROI.peek(1) << \
      ", " << _upperROI.peek(1) << ")" << endl;
    cout << "parameter weights = " << _paramWeights << endl;
    cout << "transform = " << endl;
    _transform->print();
  }
  void getMetric();
};

///////////////////////////////////////////////////////////////////////
// Experimental New metric computed for two images and a transform   //
// that takes one image into the other image's coordinate system     //
///////////////////////////////////////////////////////////////////////
template <class T>                         
class Experimental: public Metric
{
 protected:

  // weighting the parameters of the transform to restrict or encourage movement 
  // in each parameter space
  VISVector _paramWeights;
  
  // the current transform
  T* _transform;
  //the number of transform parameters
  int _numParams;
  // the derivative of the metric
  T _derivative;
  // the second derivative (Hessian) of the metric
  VISMatrix _hessian;

 public:

  Experimental()
  {  
    VISVector _paramWeights;
    _numParams = 0;
    _upperROI = VISVector(3);
    _lowerROI = VISVector(3);
    _derivative.setParams(0.0f);
    _hessian = 0.0f;
  }

  virtual ~Experimental() {;}

  virtual void setParamWeights(VISVector paramWeights)
  {
    if (paramWeights.n() == _numParams)
      _paramWeights = paramWeights;
    else 
      cout << "incorrect size of param weights" << endl;
  }

  void setTransform(T* transform)
  {
    _transform = transform;
    _numParams = (_transform->getParams()).n();
    _hessian = VISMatrix(_numParams, _numParams);
  }

  virtual VISVolume<float> processVolume(VISVolume<float> volume)
  {
#ifdef BLUR
    volume = volume.gauss(BLUR);
#endif
    volume = volume - volume.average();
    volume /= pow((float)volume.variance(), 0.5f);
    return volume;
  }
  
  // the following functions return the values of the requested
  // registration routine parameters
  virtual T* getTransform() {return _transform;}
  virtual T getDeriv() {return _derivative;}
  virtual VISMatrix getHessian() {return _hessian;}
  
  // prints out the parameters of the registration process
  virtual void print()
  {
    cout << "Experimental Metric Parameters:" << endl;
  }
  virtual void printLong()
  {
    print();
    cout << "Region of interest = (" << _lowerROI.peek(0) << \
      ", " << _upperROI.peek(0) << ", " << _lowerROI.peek(1) << \
      ", " << _upperROI.peek(1) << ")" << endl;
    cout << "parameter weights = " << _paramWeights << endl;
    cout << "transform = " << endl;
    _transform->print();
  }
  void getMetric();
};



/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////
// gausian density with covariance Phi //
/////////////////////////////////////////
template <class T>
float MutualInformation<T>::Gphi_v(float x)
{
  float n = -.5;  //this really means that n = 1
  float rtNormPhi = pow(fabs(_Phi_v), -.5f);
  float temp = exp(-x*_Phi_vINV*x/2.0f);
#ifdef WARNING
  if (fabs(temp) < eps)
  {
    cout << "Warning: one of your weigts is zero" << endl << \
      "         consider altering the 1D covariance" << endl;
    cout << "x = " << x << endl;
  }
#endif  
  return  pow((float)TWOPI, (float)n) * rtNormPhi * temp;
}

/////////////////////////////////////////
// gausian density with covariance Phi //
/////////////////////////////////////////
template <class T>
float MutualInformation<T>::Gphi_uv(VISMatrix x)
{
  float n = -.5;   //this really means that n = 1
  float rtNormPhi = pow(_Phi_uv.det(), -.5);
  VISMatrix temp = -x.t()*_Phi_uvINV*x/2.0f;  
  float temp2 = exp(temp.peek(0,0)); 
#ifdef WARNING
  if (fabs(temp2) < eps)
  {
    cout << "Warning: one of your weigts is zero" << endl <<\
      "         consider altering the 2D covariance" << endl;
  }
#endif 
  return pow((float)TWOPI, (float)n) * rtNormPhi * temp2;
}

///////////////////////////
// 1D Weighting function //
///////////////////////////
template <class T>
float MutualInformation<T>::W_v(float x, float y)
{
  float temp = 0.0f, uk, vk;
  float ret  = Gphi_v(x-y);
  int j,a,b;
  float wk;
  for (j = 0; j<_nB; j++)
  {
    //get the pixel value corresponding to xk of the test image
    vk = _Vxjs.peek(j);
    temp += Gphi_v(x-vk);  

    /*    cout << "bottom sum = " << temp << endl;
    cout << "vk = " << vk << endl;
    cout << "x = " << x << endl; 
    */
  }
  if (fabs(temp) < eps)
  {
    cout << "ERROR: DIVIDE BY ZERO (in W_v) !!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
    cout << "x = " << x << " y = " << y << endl;
  }
  ret /=  temp;
  return ret;
}

///////////////////////////
// 2D Weighting function //
///////////////////////////
template <class T>
float MutualInformation<T>::W_uv(VISMatrix x, VISMatrix y)
{
  float temp= 0.0f, uk, vk;
  float ret = Gphi_uv(x-y);
  int j;
  VISMatrix wk(2,1);
  
  for (j=0; j<_nB; j++)
  {
    //get the pixel value of xk of the ref image
    uk = _Uxjs.peek(j);
    //get the pixel value corresponding to xk of the test image
    vk = _Vxjs.peek(j);
    wk.poke(0,0) = uk; wk.poke(1,0) = vk;
    temp += Gphi_uv(x-wk);  
    /*
    cout << "bottom sum = " << temp << endl;
    cout << "wk = " << wk << endl;
    cout << "x = " << x << endl; 
    */
  }
  if (fabs(temp) < eps)
  {
    cout << "ERROR: DIVIDE BY ZERO (in W_uv)  !!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
    cout << "x = " << x << "y = " << y << endl;
  }
  ret  /= temp;
  return ret;
}

//////////////////////////////////////////////////////////////
// function to generate a new set of random points for each //
// iteration also builds vectors of pixel values            // 
// for test and ref images                                  // 
//////////////////////////////////////////////////////////////
template <class T>
void MutualInformation<T>::generateRandomPoints()
{
  VISMatrix xi(4,1), xj(4,1);
  xi = 1.0f; xj = 1.0f;
  VISVolume<float> randRef, randTest;
  VISVolumeFile vol_file;
  VISImage<float> randRefImage, randTestImage;
  VISImageFile im_file;

  randRef = _ref; randTest = _test;
  _Uxis = 0.0f; _Uxjs = 0.0f;
  float white = _ref.max(), black = _ref.min();
  for (int i=0; i<_nA; i++)
  {  
    if (_newRands) 
    {
      //generate points within the region of interest
      _xis.poke(i,0) = (rand1())*(_upperROI.peek(0)-_lowerROI.peek(0))\
	+_lowerROI.peek(0);
      _xis.poke(i,1) = (rand1())*(_upperROI.peek(1)-_lowerROI.peek(1))\
	+_lowerROI.peek(1);
      _xis.poke(i,2) = (rand1())*(_upperROI.peek(2)-_lowerROI.peek(2))\
	+_lowerROI.peek(2);
    }
    //grab xi from the matrix of random points (volume centered coords)
    xi.poke(0,0) = _xis.peek(i,0);
    xi.poke(1,0) = _xis.peek(i,1);
    xi.poke(2,0) = _xis.peek(i,2);
    xi.poke(3,0) = 1.0f; 
    //transform to image coordinates (vispack volume coordinates)
     xi = centeredToVolume(xi, _refW, _refH, _refD);
    //update the reference pixel value vector
    if (_ref.checkBounds((float)xi.peek(0,0), (float)xi.peek(1,0), (float)xi.peek(2,0)))
    {
      _Uxis.poke(i) = _ref.interp(xi.peek(0,0), xi.peek(1,0), xi.peek(2,0));
#ifdef DEBUG
      //change pixel value of reference image to white for each random point    
      randRef.poke((int)floor(xi.peek(0,0)), (int)floor(xi.peek(1,0)), \
		   (int)floor(xi.peek(2,0))) = white;
#endif
    }
    else
      _Uxis.poke(i) = 0.0f;
    
    //grab xi from the matrix of random points (volume centered coords)
    xi.poke(0,0) = _xis.peek(i,0);
    xi.poke(1,0) = _xis.peek(i,1);
    xi.poke(2,0) = _xis.peek(i,2);
    xi.poke(3,0) = 1.0f; 
    
     
    //transform the coordinate to the test image (volume centered coords)
    xi = _transform->transformPoint(xi); xi /= xi.peek(3,0);
    
    //update the transformed points matrix (volume centered coords)
    _Txis.poke(i,0) = xi.peek(0,0);
    _Txis.poke(i,1) = xi.peek(1,0);
    _Txis.poke(i,2) = xi.peek(2,0);
    
    //transform the coordinate into image coordinates
    xi = centeredToVolume(xi, _refW, _refH, _refD);
    //check to see if point is within the image

    if(_test.checkBounds((float)xi.peek(0,0), (float)xi.peek(1,0), (float)xi.peek(2,0)))
    {
      //update the test pixel value vector
      _Vxis.poke(i) = _test.interp(xi.peek(0,0), xi.peek(1,0), xi.peek(2,0));
      //update the dx and dy vectors
      _dxTestXis.poke(i) = _dxTest.interp(xi.peek(0,0), xi.peek(1,0), xi.peek(2,0));
      _dyTestXis.poke(i) = _dyTest.interp(xi.peek(0,0), xi.peek(1,0), xi.peek(2,0));
      _dzTestXis.poke(i) = _dzTest.interp(xi.peek(0,0), xi.peek(1,0), xi.peek(2,0));
#ifdef DEBUG
      randTest.poke((int)floor(xi.peek(0,0)), (int)floor(xi.peek(1,0)), \
		    (int)floor(xi.peek(2,0))) = white;
#endif
    }
    else
    {
      _Vxis.poke(i) = 0.0f;
      _dxTestXis.poke(i) = 0.0f;
      _dyTestXis.poke(i) = 0.0f;
      _dzTestXis.poke(i) = 0.0f;
    }
  }
#ifdef WRITERAND
  vol_file.write_float(randRef, "/scratch/research/MutualInformationCode/rand/randrefXi.float");
  randRefImage = randRef.image(SLICE);
  im_file.write(randRefImage, "/scratch/research/MutualInformationCode/rand/randrefXi.fits");
  vol_file.write_float(randTest, "/scratch/research/MutualInformationCode/rand/randtestXi.float");
  randTestImage = randTest.image(SLICE);
  im_file.write(randTestImage, "/scratch/research/MutualInformationCode/rand/randtestXi.fits");
#endif  
  randRef = _ref;
  randTest = _test;
  for (int j=0; j<_nB; j++)
  {
    if (_newRands)
    {
      //image centered random points
      _xjs.poke(j,0) = (rand1())*(_upperROI.peek(0)-_lowerROI.peek(0))\
	+_lowerROI.peek(0);
      _xjs.poke(j,1) = (rand1())*(_upperROI.peek(1)-_lowerROI.peek(1))\
	+_lowerROI.peek(1);
      _xjs.poke(j,2) = (rand1())*(_upperROI.peek(2)-_lowerROI.peek(2))\
	+_lowerROI.peek(2);
    }

    //grab xj from the matrix of random points
    xj.poke(0,0) = _xjs.peek(j,0);
    xj.poke(1,0) = _xjs.peek(j,1);
    xj.poke(2,0) = _xjs.peek(j,2);
    xj.poke(3,0) = 1.0f; 
    //transform the image centered point into image coordinates
    xj = centeredToVolume(xj, _refW, _refH, _refD);
    //update the reference pixel value vector
    if (_ref.checkBounds((float)xj.peek(0,0), (float)xj.peek(1,0), (float)xj.peek(2,0)))
    {
      _Uxjs.poke(j) = _ref.interp(xj.peek(0,0), xj.peek(1,0), xj.peek(2,0));
#ifdef DEBUG
      //change grey scale value to white for random point
      randRef.poke((int)floor(xj.peek(0,0)), (int)floor(xj.peek(1,0)), 
		   (int)floor(xj.peek(2,0))) = white;
#endif
    }
    else
      _Uxjs.poke(j) = 0.0f;

    //grab xj from the matrix of random points
    xj.poke(0,0) = _xjs.peek(j,0);
    xj.poke(1,0) = _xjs.peek(j,1);
    xj.poke(2,0) = _xjs.peek(j,2);
    xj.poke(3,0) = 1.0f; 
    //transform the coordinate to the test image (volume centered coords)
    xj = _transform->transformPoint(xj); xj /= xj.peek(3,0);
    //update the transformed points matrix (image centered coords)
    _Txjs.poke(j,0) = xj.peek(0,0);
    _Txjs.poke(j,1) = xj.peek(1,0);
    _Txjs.poke(j,2) = xj.peek(2,0);
    //change from image centered to image coordinates
    xj = centeredToVolume(xj, _refW, _refH, _refD);
    //check to see if point is within the image
    if(_test.checkBounds((float)xj.peek(0,0), (float)xj.peek(1,0), (float)xj.peek(2,0)))
    {
      _Vxjs.poke(j) = _test.interp(xj.peek(0,0), xj.peek(1,0), xj.peek(2,0));
      //update the dx and dy vectors
      _dxTestXjs.poke(j) = _dxTest.interp(xj.peek(0,0), xj.peek(1,0), xj.peek(2,0));
      _dyTestXjs.poke(j) = _dyTest.interp(xj.peek(0,0), xj.peek(1,0), xj.peek(2,0));
      _dzTestXjs.poke(j) = _dzTest.interp(xj.peek(0,0), xj.peek(1,0), xj.peek(2,0));
#ifdef WRITERAND
      randTest.poke((int)floor(xj.peek(0,0)), (int)floor(xj.peek(1,0)), \
		    (int)floor(xj.peek(2,0))) = white;
#endif
    }
    else
    {
      _Vxjs.poke(j) = 0.0f;
      _dxTestXjs.poke(j) = 0.0f;
      _dyTestXjs.poke(j) = 0.0f;
      _dzTestXjs.poke(j) = 0.0f;
    }
  }
#ifdef WRITERAND
  vol_file.write_float(randRef, "/scratch/research/MutualInformationCode/rand/randrefXj.float");
  randRefImage = randRef.image(SLICE);
  im_file.write(randRefImage, "/scratch/research/MutualInformationCode/rand/randrefXj.fits");
  vol_file.write_float(randTest, "/scratch/research/MutualInformationCode/rand/randtestXj.float");
  randTestImage = randTest.image(SLICE);
  im_file.write(randTestImage, "/scratch/research/MutualInformationCode/rand/randtestXj.fits");
#endif
}


////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////
// Calculates the Mutual Information metric, gradient and hessian //
// with respect to the current transform                          //
////////////////////////////////////////////////////////////////////
template <class T>
void MutualInformation<T>::getMetric()
{
  VISVector xi(4), xj(4);
  VISMatrix wi(2,1), wj(2,1), wk(2,1), wij(2,1), wik(2,1);
  xi = 1.0f; xj = 1.0f;
  int i, j, k, m, n;
  
  VISVector delTesti(3), delTestj(3);
  delTesti = 0.0f; delTestj= 0.0f;
  
  float ui, vi, uj, vj, vk, vik, vij, tempTest, tempBoth, tempRef, H1D, H2D;
  float partA, partB, partC, partD, multA, multB, multC, multD;
  float weighting_factor, dVxijm;
  VISMatrix partOne, partTwo;
  VISVector dVxidParams, dVxjdParams, dVxijdParams;

  generateRandomPoints();
  _newRands = false;
  
  //get new random points for A and B
  generateRandomPoints();
  
  //reset MI, derivative, and hessian to zero
  _metric = 0.0f; 
  _derivative.setParams(0.0f);
  _hessian = 0.0f;
  partOne = _hessian;
  partTwo = _hessian;

  for (i=0; i<_nA; i++)
  {
    //get xi
    xi.poke(0) = _xis.peek(i,0);
    xi.poke(1) = _xis.peek(i,1);
    xi.poke(2) = _xis.peek(i,2);
    //get the pixel value at xi of the reference image
    ui = _Uxis.peek(i);
    //get the pixel value corresponding to xi of the test image
    vi = _Vxis.peek(i);
    //update the wi vector
    wi.poke(0,0) = ui; wi.poke(1,0) = vi;
    
    tempRef = 0.0f; tempTest = 0.0f; tempBoth = 0.0f;
    for (j=0; j<_nB; j++)
    {
      //get xj
      xj.poke(0) = _xjs.peek(j,0);
      xj.poke(1) = _xjs.peek(j,1);
      xj.poke(2) = _xjs.peek(j,2);
      //get the pixel value at xj of the reference image
      uj = _Uxjs.peek(j);
      //get the pixel value corresponding to xj of the test image
      vj = _Vxjs.peek(j);
      //update the wj vector
      wj.poke(0,0) = uj; wj.poke(1,0) = vj;

      //save on computations
      vij = (vi-vj); wij = (wi-wj); 

      //compute the approximation to mutual information
      /////////////////////////////////////////////////
      tempRef += Gphi_v(ui-uj);
      tempTest += Gphi_v(vij);
      tempBoth += Gphi_uv(wij);
      
      //save on computations
      weighting_factor = (W_v(vi,vj)*_Phi_vINV-W_uv(wi,wj)*_Phi_vvINV);
      dVxidParams = _transform->dVdParams((float)_dxTestXis.peek(i), \
					  (float)_dyTestXis.peek(i), \
					  (float)_dzTestXis.peek(i), xi);
      dVxjdParams = _transform->dVdParams((float)_dxTestXjs.peek(j), \
					  (float)_dyTestXjs.peek(j), \
					  (float)_dzTestXjs.peek(j), xj);
      dVxijdParams = (dVxidParams-dVxjdParams);

      // compute the approximation to the derivative (gradient)
      // of mutual information wrt the transform
      ////////////////////////////////////////////////////////////
      _derivative.setParams( _derivative.getParams() \
			     + (vij)*weighting_factor \
			     * (dVxijdParams));
            
      // compute the approximation to the second derivative (hessian)
      // of the mutual information wrt the transform
      // note that: because order of differentiation does not matter, 
      // we only compute the lower triangular part of the hessian and 
      // then fill in the rest later.
      ///////////////////////////////////////////////////////////////////

      //calculate the first part of the approx to the hessian
      ///////////////////////////////////////////////////////
      partOne += weighting_factor*VISMatrix((dVxijdParams)) \
	*(dVxijdParams).t();
      
      //calculate the rest of the approx to the hessian
      /////////////////////////////////////////////////
#ifdef moreAccurateHessian
      multA = _Phi_vINV*_Phi_vINV*Gphi_v(vij);
      multB = _Phi_vvINV*_Phi_vvINV*Gphi_uv(wij);
      
      for (m=0; m<_numParams; m++)
      {
	dVxijm = dVxijdParams.peek(m);
	for (n=0; n<=m; n++)
	{  
      	  partA = 0.0f; 
	  partB = 0.0f;
	  partC = 0.0f;
	  partD = 0.0f;
	  for (k=0;k<_nB;k++)
	  {
	    vk = _Vxjs.peek(k);
	    wk.poke(0,0) = _Uxjs.peek(k);
	    wk.poke(1,0) = vk;
	    
	    //save on computations
	    vik = (vi-vk); wik = (wi-wk);
	    partB += Gphi_v(vik);	    
	    partA += (vik)*partB*dVxijm;
	    partC += Gphi_uv(wik);
	    partD += (vik)*partC*dVxijm;
	  }  //end k for loop
	   
	  partTwo.poke(m,n) += (vij)*(dVxijdParams.peek(n))\
	    *(multA*(partA-partB*vij*(dVxijm))/(pow(partB,2))\
	      +multB*(partC*vij*(dVxijm)-partD)/pow(partC,2));

	} //end n for loop
      }  //end m for loop
      //first part is symmetric, so fill in rest
      for (m=0; m<_numParams; m++)
	for (n=0; n<m; n++)
	  partTwo.poke(n,m) = partTwo(m,n);
#endif      

    } //end for j
        
    tempTest = log(tempTest/_nB); 
    tempBoth = log(tempBoth/_nB);
    tempRef = log(tempRef/_nB);
    _metric += tempBoth - tempTest - tempRef;
  }  //end for i

  _metric /= _nA;
  _derivative.setParams(_derivative.getParams()/_nA);

  //multiply the new directions by the parameterWeights
  //_derivative.setParams(_derivative.getParams().elementMult(_paramWeights));
  //_hessian = _hessian.elementMult(VISMatrix(_paramWeights)*VISMatrix(_paramWeights.t()));
  
  //add the first part of the approximation to the second part 
  //(if the second has been calculated)
  _hessian += partOne;

#ifdef moreAccurateHessian
  _hessian += partTwo;
#endif
  _hessian /= _nA;
}

////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////


  


//////////////////////////////////////////////////////////////
// Calculates the Experimental metric, gradient and hessian //
// with respect to the current transform                    //  
//////////////////////////////////////////////////////////////
template <class T>
void Experimental<T>::getMetric()
{

  _metric = 0;
  VISVector delF(3), delG(3), x(4);
  float q;
  
  for (int k; k<= _refD; k++)
  {
    for (int j; j<= _refH; j++)
    {
      for (int i; i<=_refW; i++)
      {
	
	x.poke(0) = i;
	x.poke(1) = j;
	x.poke(2) = k;
	x.poke(3) = 1;
	x = _transform->transformPoint(x);

	//build up del f
	if (_dxTest.checkBounds((int)floor(x.peek(0)), 
				(int)floor(x.peek(1)),
				(int)floor(x.peek(2))))
	{
	  delF.poke(0) = _dxTest.interp(x.peek(0), x.peek(1), x.peek(2));
	  delF.poke(1) = _dyTest.interp(x.peek(0), x.peek(1), x.peek(2)); 
	  delF.poke(2) = _dzTest.interp(x.peek(0), x.peek(1), x.peek(2));
	  
	  //build up del g
	  delG.poke(0) = _dxRef.peek(i,j,k);
	  delG.poke(1) = _dyRef.peek(i,j,k);
	  delG.poke(2) = _dzRef.peek(i,j,k);
	  
	  q = delF.dot(delG) / delF.dot(delF);
	  _metric += (q*delF-delG).dot(q*delF-delG);
	}
      }
    }
  }

}
