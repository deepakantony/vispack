/*
  transform.h
  Sarah Geneser
  written: 05-30-03
  last updated: 08-09-03
  
  specifies several types of transformations, including:
  
  *2D affine transformations
  *2D euler rigid transformations
  *2D translations

*/


#include<math.h>
#include<matrix.h>
#include<image.h>

#define SLICE 0
 
///////////////////////////////////////////////////
// converts from vispack to centered coordinates //
///////////////////////////////////////////////////
VISVector volumeToCentered(VISVector x, int w, int h, int d)
{
  x.poke(0) = x.peek(0) - w/2.0f;
  x.poke(1) = h/2.0f - x.peek(1);
  x.poke(2) = x.peek(2) - d/2.0f;
  return x;
}

///////////////////////////////////////////////////
// converts from centered to vispack coordinates //
///////////////////////////////////////////////////
VISVector centeredToVolume(VISVector x, int w, int h, int d)
{
  x.poke(0) = x.peek(0) + w/2.0f;
  x.poke(1) = h/2.0f - x.peek(1);
  x.poke(2) = x.peek(2) + d/2.0f;
  return x;
} 

class Transform
{
public:
  virtual VISVector getParams() const {return VISVector();}
  virtual void setParams(const VISVector &params) {;}
  virtual int numParams() const {return(0);}
  virtual VISVector dVdParams(float dVdX, float dVdY, VISVector &x){return VISVector();}
  virtual VISMatrix dXdParams(const VISVector &x) const {return VISMatrix();}
  virtual VISVector transformPoint(const VISVector &x) const {return x;}
  Transform() {;}

};


//////////////////////////////////////////////////////////////////////////////
// AffineTransforms are assumed to act upon 3D image centeredoordinates   //
// in the image domain (not the normalized image coordinates)               //
//////////////////////////////////////////////////////////////////////////////
class AffineTransform: public Transform
{
protected:
  VISMatrix _T;
  VISMatrix _Tinv;
  VISVector _params;
public:
  AffineTransform() 
  {
    _T = VISMatrix(4,4);
    _T = 0.0f; _T.poke(0,0) = _T.poke(1,1) = _T.poke(2,2) = _T.poke(3,3) = 1.0f;
    _Tinv = _T.inv();
    VISVector _params(12); 
    _params = 0.0f; _params.poke(0) = _params.poke(4) = 1.0f;
  }

  AffineTransform(VISVector params)
  {
    //if there are the appropriate number of parameters
    if(params.n() == 12)
    {
      //fill in the parameter vector
      VISVector _params;
      _params = params;
      
      //fill in the transform matrix
      _T = VISMatrix(4, 4);
      //assume the parameters to be in the appropriate order
      _T.poke(0,0) = params.peek(0);
      _T.poke(0,1) = params.peek(1);
      _T.poke(0,2) = params.peek(2);
      _T.poke(0,3) = params.peek(3);
      _T.poke(1,0) = params.peek(4);
      _T.poke(1,1) = params.peek(5);
      _T.poke(1,2) = params.peek(6);
      _T.poke(1,3) = params.peek(7);
      _T.poke(2,0) = params.peek(8);
      _T.poke(2,1) = params.peek(9);
      _T.poke(2,2) = params.peek(10);
      _T.poke(2,3) = params.peek(11);

      //fill in the last row for homogenous coordinates
      _T.poke(3,0) = _T.poke(3,1) = _T.poke(3,2) = 0.0f; _T.poke(3,3) = 1.0f;
      _Tinv = _T.inv();
    }
    else 
      cout << "ERROR: specified the incorrect number of transform parameters" \
	   << endl;
  }

  virtual VISVector getParams() const {return _params;}
  
  virtual void setParams(const float param) 
  {
    _params = VISVector(12);
    _params = param;
    _T = param;
    _Tinv = _T.inv();
  }
  
  virtual void setParams(const VISMatrix T)
  {
    if (T.c()==4 && T.r()==4)
    {
      _T = T;
      _Tinv = _T.inv();
      _params = VISVector(12);
      _params.poke(0) = _T.peek(0,0);
      _params.poke(1) = _T.peek(0,1);
      _params.poke(2) = _T.peek(0,2);
      _params.poke(3) = _T.peek(0,3);
      _params.poke(4) = _T.peek(1,0);
      _params.poke(5) = _T.peek(1,1);
      _params.poke(6) = _T.peek(1,2);
      _params.poke(7) = _T.peek(1,3);
      _params.poke(8) = _T.peek(2,0);
      _params.poke(9) = _T.peek(2,1);
      _params.poke(10) = _T.peek(2,2);
      _params.poke(11) = _T.peek(2,3);
    }
    else
      cout << "ERROR: specified transform of incorrect size" << endl;
  }

  virtual void setParams(const VISVector params) 
  {

    //if there are the appropriate number of parameters
    if(params.n() == 12)
    {
      //fill in the parameter vector
      _params = params;
      
      //assume the parameters to be in the appropriate order
      _T.poke(0,0) = params.peek(0);
      _T.poke(0,1) = params.peek(1);
      _T.poke(0,2) = params.peek(2);
      _T.poke(0,3) = params.peek(3);
      _T.poke(1,0) = params.peek(4);
      _T.poke(1,1) = params.peek(5);
      _T.poke(1,2) = params.peek(6);
      _T.poke(1,3) = params.peek(7);
      _T.poke(2,0) = params.peek(8);
      _T.poke(2,1) = params.peek(9);
      _T.poke(2,2) = params.peek(10);
      _T.poke(2,3) = params.peek(11);

      //fill in the last row for homogenous coordinates
      _T.poke(3,0) = _T.poke(3,1) = _T.poke(3,2) = 0.0f; _T.poke(3,3) = 1.0f;
      _Tinv = _T.inv();
    }
    else 
      cout << "ERROR: specified the incorrect number of transform parameters" \
	   << endl;
  }

  virtual int numParams() const {return(12);}

  virtual VISMatrix dXdParams(const VISVector &x) const
  {
    VISMatrix dXdP(4,12);
    // x y z 0 0 0 0 0 0 1 0 0
    // 0 0 0 x y z 0 0 0 0 1 0
    // 0 0 0 0 0 0 x y z 0 0 1
    // 0 0 0 0 0 0 0 0 0 0 0 0 
    dXdP = 0.0f;
    dXdP.poke(0,0) = x.peek(0);
    dXdP.poke(0,1) = x.peek(1);
    dXdP.poke(0,2) = x.peek(2);
    dXdP.poke(0,9) = 1.0f;
    dXdP.poke(1,4) = x.peek(0);
    dXdP.poke(1,5) = x.peek(1);
    dXdP.poke(1,6) = x.peek(2);
    dXdP.poke(1,10) = 1.0f;
    dXdP.poke(2,6) = x.peek(0);
    dXdP.poke(2,7) = x.peek(1);
    dXdP.poke(2,8) = x.peek(2);
    dXdP.poke(2,11) = 1.0f;
    
    return dXdP;
  }
  
  virtual VISVector dVdParams(float dVdX, float dVdY, float dVdZ, VISVector &x)
  {
    VISMatrix dVdP(4,4);
    VISMatrix delV(4,1);
    VISVector ret(12);
    delV.poke(0,0) = dVdX;
    delV.poke(1,0) = dVdY;
    delV.poke(2,0) = dVdZ;
    delV.poke(3,0) = 0.0f;
    if (SLICE == 0) 
    {
      x.poke(2) = 0.0f;
      delV.poke(2,0) = 0;
    }
    
    dVdP = delV * x.t();

    //unroll the matrix into a vector
    ret.poke(0)  = dVdP.peek(0,0);
    ret.poke(1)  = dVdP.peek(0,1);
    ret.poke(2)  = dVdP.peek(0,2);
    ret.poke(3)  = dVdP.peek(0,3);
    ret.poke(4)  = dVdP.peek(1,0);
    ret.poke(5)  = dVdP.peek(1,1);
    ret.poke(6)  = dVdP.peek(1,2);
    ret.poke(7)  = dVdP.peek(1,3);
    ret.poke(8)  = dVdP.peek(2,0);
    ret.poke(9)  = dVdP.peek(2,1);
    ret.poke(10) = dVdP.peek(2,2);
    ret.poke(11) = dVdP.peek(2,3);
    return ret;
  }
  
  virtual VISVector transformPoint(const VISVector &x) const {return(_T*x);}

  virtual VISVolume<float> transformVolume(const VISVolume<float> &vol) const
  {
    VISVector x(4);
    x = 1.0f;
    int h = vol.height(), w = vol.width(), d = vol.depth(), i, j, k;
    VISVolume<float> ret(w,h,d); ret = 0.0f;
    for (j=0; j<h; j++)
    {
      for (i=0; i<w; i++)
      {
	for (k=0; k<d; k++)
	{
	  //fill in the coordinates of the vector
	  x.poke(0) = i;  x.poke(1) = j; x.poke(2) = k;
	  //transform the point into image centered coordinates
	  x = volumeToCentered(x,w,h,d);
	  //transform the point via the affine transform
	  x = _T*x;
	  //transform the point into vispack image coordinates
	  x = centeredToVolume(x,w,h,d);
	  
	  if (vol.checkBounds((float)x.peek(0),(float)x.peek(1),(float)x.peek(2)))
	    ret.poke(i,j,k) = \
	      vol.interp((float)x.peek(0),(float)x.peek(1),(float)x.peek(2));
	}
      }
    }
    return ret;
  }

  virtual void print() const {cout << "T= " << _T << endl;}
};

//////////////////////////////////////////////////////////////////////////////
// EulerRigidTransforms are assumed to act upon 3D homogenous coordinates //
//////////////////////////////////////////////////////////////////////////////
class EulerRigidTransform: public Transform
{
protected:
  VISMatrix _T, _Tinv;
  float _xTrans, _yTrans, _zTrans;
  float _xTheta, _yTheta, _zTheta;
  VISVector _params;
public:
  EulerRigidTransform() 
  {
    _xTheta = _yTheta = _zTheta = _xTrans = _yTrans = _zTrans = 0.0f;
    _params = VISVector(6);
    _T = VISMatrix(4,4);
    _T = 0.0f; 
    _T.poke(0,0) = _T.poke(1,1) = _T.poke(2,2) = _T.poke(3,3) = 1.0f;
    _Tinv = _T.inv();
  }

  EulerRigidTransform(VISVector params)
  {
    //if there are the appropriate number of parameters
    if(params.n() == 6)
    {
      //fill in the parameter vector
      VISVector _params;
      _params = params;
      
      //assume the parameters to be in the appropriate order
      _xTrans = params.peek(0);
      _yTrans = params.peek(1);
      _zTrans = params.peek(2);
      _xTheta = params.peek(3);
      _yTheta = params.peek(4);
      _zTheta = params.peek(5);
      
      //fill in the transform matrix
      _T = VISMatrix(4,4); _T = 0.0f;
      
      float cosX = cos(_xTheta), cosY = cos(_yTheta), cosZ = cos(_zTheta);
      float sinX = sin(_xTheta), sinY = sin(_yTheta), sinZ = sin(_zTheta);
           
      _T.poke(0,0) = cosY*cosZ;
      _T.poke(0,1) = sinX*sinY*cosZ + cosX*sinZ;
      _T.poke(0,2) = -cosX*sinY*cosZ + sinX*sinZ;
      _T.poke(0,3) = _xTrans;
      _T.poke(1,0) = -cosY*sinZ;
      _T.poke(1,1) = -sinX*sinY*sinZ + cosX*cosZ;
      _T.poke(1,2) = cosX*sinY*sinZ + sinX*cosZ;
      _T.poke(1,3) = _yTrans;
      _T.poke(2,0) = sinY;
      _T.poke(2,1) = -sinX*cosY;
      _T.poke(2,2) = cosX*cosY;
      _T.poke(2,3) = _zTrans;

      //fill in the last row for homogenous coordinates
      _T.poke(3,0) = _T.poke(3,1) = _T.poke(3,2) = 0.0f; _T.poke(3,3) = 1.0f;
      _Tinv = _T.inv();
    }
    else 
      cout << "ERROR: specified the incorrect number of transform parameters" \
	   << endl;
  }
  
  virtual void setParams(const VISVector &params) 
  {
    if(params.n() == 6)
    {
      //if there are the appropriate number of parameters
      //fill in the parameter vector
      _params = params;
      
      //assume the parameters to be in the appropriate order
      _xTrans = params.peek(0);
      _yTrans = params.peek(1);
      _zTrans = params.peek(2);
      _xTheta = params.peek(3);
      _yTheta = params.peek(4);
      _zTheta = params.peek(5);
      
      //fill in the transform matrix
      _T = VISMatrix(4,4); _T = 0.0f;
     
      float cosX = cos(_xTheta), cosY = cos(_yTheta), cosZ = cos(_zTheta);
      float sinX = sin(_xTheta), sinY = sin(_yTheta), sinZ = sin(_zTheta); 
                
      _T.poke(0,0) = cosY*cosZ;
      _T.poke(0,1) = sinX*sinY*cosZ + cosX*sinZ;
      _T.poke(0,2) = -cosX*sinY*cosZ + sinX*sinZ;
      _T.poke(0,3) = _xTrans;
      _T.poke(1,0) = -cosY*sinZ;
      _T.poke(1,1) = -sinX*sinY*sinZ + cosX*cosZ;
      _T.poke(1,2) = cosX*sinY*sinZ + sinX*cosZ;
      _T.poke(1,3) = _yTrans;
      _T.poke(2,0) = sinY;
      _T.poke(2,1) = -sinX*cosY;
      _T.poke(2,2) = cosX*cosY;
      _T.poke(2,3) = _zTrans;

      //fill in the last row for homogenous coordinates
      _T.poke(3,0) = _T.poke(3,1) = _T.poke(3,2) = 0.0f; _T.poke(3,3) = 1.0f;
      _Tinv = _T.inv();
    }
    else 
      cout << "ERROR: specified the incorrect number of transform parameters" \
	   << endl;
  }

  virtual void setParams(const float param) 
  {
    _params = param;
    
    //assume the parameters to be in the appropriate order
    _xTrans = _yTrans = _zTrans = param;
    _xTheta = _yTheta = _zTheta = param;
    
    //fill in the transform matrix
    _T = VISMatrix(4,4); _T = 0.0f;
        
    float cosX = cos(_xTheta), cosY = cos(_yTheta), cosZ = cos(_zTheta);
    float sinX = sin(_xTheta), sinY = sin(_yTheta), sinZ = sin(_zTheta); 
        
    _T.poke(0,0) = cosY*cosZ;
    _T.poke(0,1) = sinX*sinY*cosZ + cosX*sinZ;
    _T.poke(0,2) = -cosX*sinY*cosZ + sinX*sinZ;
    _T.poke(0,3) = _xTrans;
    _T.poke(1,0) = -cosY*sinZ;
    _T.poke(1,1) = -sinX*sinY*sinZ + cosX*cosZ;
    _T.poke(1,2) = cosX*sinY*sinZ + sinX*cosZ;
    _T.poke(1,3) = _yTrans;
    _T.poke(2,0) = sinY;
    _T.poke(2,1) = -sinX*cosY;
    _T.poke(2,2) = cosX*cosY;
    _T.poke(2,3) = _zTrans;
    
    //fill in the last row for homogenous coordinates
    _T.poke(3,0) = _T.poke(3,1) = _T.poke(3,2) = 0.0f; _T.poke(3,3) = 1.0f;
    _Tinv = _T.inv();
  }
  
  virtual VISVector getParams() const {return _params;}

  virtual int numParams() const {return(3);}

  virtual VISMatrix dXdParams(const VISVector &x) const
  {
    
    VISMatrix dXdP(4,4);
    
    //////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////
    // 1 0  -sin(theta)*x+cos(theta)*y
    // 0 1  -cos(theta)*x-sin(theta)*y
    // 0 0               0

    dXdP = 0.0f;
    
    return dXdP;
  }
  
  virtual VISVector dVdParams(float dVdX, float dVdY,  float dVdZ, VISVector &x)
  {
    VISImageFile im;
    VISVector dVdP(6);
    dVdP = 0.0f;
    
    float cosX = cos(_xTheta), cosY = cos(_yTheta), cosZ = cos(_zTheta);
    float sinX = sin(_xTheta), sinY = sin(_yTheta), sinZ = sin(_zTheta);
    
    //compute d(v(T(xi))/dT
    dVdP.poke(0) = dVdX;
    dVdP.poke(1) = dVdY;
    dVdP.poke(2) = dVdZ;
    dVdP.poke(3) = (dVdX*(x.peek(1)*(cosX*sinY*cosZ-sinX*sinZ) \
			  +x.peek(2)*(sinX*sinY*cosZ+cosX*sinZ)) \
		    + dVdY*(x.peek(1)*(-cosX*sinY*sinZ-sinX*cosZ) \
			    +x.peek(2)*(-sinX*sinY*sinZ+cosX*cosZ)) \
		    +dVdZ*(x.peek(1)*(-cosX*cosY)+x.peek(2)*(-sinX*cosY)));
    
    dVdP.poke(4) = (dVdX*(x.peek(0)*(-sinY*cosZ)+x.peek(1)*(sinX*cosY*cosZ) \
			  +x.peek(2)*(-cosX*cosY*cosZ)) \
		    + dVdY*(x.peek(0)*(sinY*sinZ)+ x.peek(1)*(-sinX*cosY*sinZ) \
			    +x.peek(2)*(cosX*cosY*sinZ)) \
		    +dVdZ*(x.peek(0)*cosY+x.peek(1)*(sinX*sinY) \
			   +x.peek(2)*(-cosX*sinY)));

    dVdP.poke(5) = (dVdX*(x.peek(0)*(-cosY*sinZ) \
			  + x.peek(1)*(-sinX*sinY*sinZ+cosX*cosZ) \
			  +x.peek(2)*(cosX*sinY*sinZ+sinX*sinZ)) \
		    + dVdY*(x.peek(0)*(-cosX*cosZ) + \
			    x.peek(1)*(-sinX*sinY*sinZ-cosX*sinZ) \
			    +x.peek(2)*(cosX*sinY*cosZ-sinX*sinZ)) \
		    +dVdZ*(x.peek(1)*(sinX*sinZ)));

    return dVdP;
  }

  virtual VISVector transformPoint(const VISVector &x) const {return(_T*x);}

  virtual VISVolume<float> transformVolume(const VISVolume<float> &vol) const
  {
    VISVector x(4);
    x = 1.0f;
    int h = vol.height(), w = vol.width(), d = vol.depth(), i, j, k;
    VISVolume<float> ret(w,h,d); ret = 0.0f;
    for (j=0; j<h; j++)
    {
      for (i=0; i<w; i++)
      {
	for (k=0; k<d; k++)
	{
	  
	  //fill in the coordinates of the vector
	  x.poke(0) = i;  x.poke(1) = j; x.poke(2) = k;
	  
	  //transform the point into image centered coordinates
	  x = volumeToCentered(x,w,h,d);
	  
	  //transform the point via the affine transform
	  x = _T*x;
	  
	  //transform the point into vispack image coordinates
	  x = centeredToVolume(x,w,h,d);
	  
	  if (vol.checkBounds((float)x.peek(0),(float)x.peek(1),(float)x.peek(2)))
	    ret.poke(i,j,k) = \
	      vol.interp((float)x.peek(0),(float)x.peek(1),(float)x.peek(2));
	}
      }
    }
    return ret;
  }
  
  virtual void print() const 
  {
    cout << "x translation = " << _xTrans << endl;
    cout << "y translation = " << _yTrans << endl;
    cout << "z translation = " << _zTrans << endl;
    cout << "x axis rotation = " << _xTheta << endl;
    cout << "y axis rotation = " << _yTheta << endl;
    cout << "z axis rotation = " << _zTheta << endl;
    cout << "transform matrix = " << _T << endl; 
  }

};

////////////////////////////////////////////////////////////////////////////
// TranslateTransform are assumed to act upon 3D homogenous coordinates //
////////////////////////////////////////////////////////////////////////////
class TranslateTransform: public Transform
{
protected:
  float _xTrans;
  float _yTrans;
  float _zTrans;
  VISVector _params;
public:
  TranslateTransform() 
  {
    _xTrans = _yTrans = _zTrans = 0.0f;
    _params = VISVector(3);
    _params = 0.0f;
  }

  TranslateTransform(VISVector params)
  {
    //if there are the appropriate number of parameters
    if(params.n() == 3)
    {
      //fill in the parameter vector
      VISVector _params;
      _params = params;
      
      //assume the parameters to be in the appropriate order
      _xTrans = params.peek(0);
      _yTrans = params.peek(1);
      _zTrans = params.peek(2);
    }
    else 
      cout << "ERROR: specified the incorrect number of transform parameters" \
	   << endl;
  }

  virtual VISVector getParams() const {return _params;}

  virtual void setParams(const float param) 
  {
    _params = param;
    _xTrans = param;
    _yTrans = param;
    _zTrans = param;
  }
  
  virtual void setParams(const VISVector &params) 
  {
    if(params.n() == 3)
    {
      //if there are the appropriate number of parameters
      //fill in the parameter vector
      _params = params;
      
      //assume the parameters to be in the appropriate order
      _xTrans = params.peek(0);
      _yTrans = params.peek(1);
      _zTrans = params.peek(2);
    }
    else 
      cout << "ERROR: specified the incorrect number of transform parameters" \
	   << endl;
  }

  virtual int numParams() const {return(3);}

  virtual VISMatrix dXdParams(const VISVector &x) const
  {
    VISMatrix dXdP(3,3);
    // 1 0 0
    // 0 1 0
    // 0 0 1

    dXdP = 0.0f;
    dXdP.poke(0,0) = 1.0f;
    dXdP.poke(1,1) = 1.0f;
    dXdP.poke(2,2) = 1.0f;
    return dXdP;
  }
  
  virtual VISVector dVdParams(float dVdX, float dVdY, float dVdZ, VISVector &x)
  {
    VISImageFile im;
    VISVector dVdP(3);
    dVdP = 0.0f;
    //compute d(v(T(xi))/dT
    dVdP.poke(0) = dVdX;
    dVdP.poke(1) = dVdY;
    dVdP.poke(2) = dVdZ;
    return dVdP;
  }

  //transforms a homogenous point
  virtual VISVector transformPoint(const VISVector &x) const 
  {
    VISVector ret(4);
    ret.poke(0) = x.peek(0) + _xTrans;
    ret.poke(1) = x.peek(1) + _yTrans;
    ret.poke(2) = x.peek(2) + _zTrans;
    ret.poke(3) = 1.0f;
    return(ret);
  }

  virtual VISVolume<float> transformVolume(const VISVolume<float> &vol) const
  {
    VISVector x(3);
    x = 1.0f;
    int h = vol.height(), w = vol.width(), d= vol.depth(), i, j, k;
    VISVolume<float> ret(w,h,d);  ret = 0.0f;
    for (j=0; j<h; j++)
    {
      for (i=0; i<w; i++)
      {
	for (k=0; k<d; k++)
	{
	  //fill in the coordinates of the vector
	  x.poke(0) = i;  x.poke(1) = j; x.poke(2) = k;
	  
	  //translate the point
	  x.poke(0) = x.peek(0) + _xTrans;
	  x.poke(1) = x.peek(1) + _yTrans;
	  x.poke(2) = x.peek(2) + _zTrans;
	  
	  if (vol.checkBounds((float)x.peek(0),(float)x.peek(1),(float)x.peek(2)))
	    ret.poke(i,j,k) = vol.interp((float)x.peek(0),(float)x.peek(1),(float)x.peek(2));
	}
      }
    }
    return ret;
  }

  virtual void print() const {cout << _params << endl;}

};
