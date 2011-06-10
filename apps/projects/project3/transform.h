
#include <matrix.h>
#include <image.h>
#include <array.h>

// This is a declaration so that pointer can be used in other classes
class Transform;


// This object takes transformations and gives the ability to
// warp/resample images
// It is also a "wrapper" or "envelope" for a coordinate
// transformation 
class Warp
{
protected:
  // This object does the coordinate transformation.  
  // It is managed through a "copy on write" protocol
  Transform *_transform;
public:
  // If the _transform object is missing, we do an identity
  // transformation
  Warp() {_transform = NULL;}

  // In case you want to transform a single point
  virtual VISVector transform(const VISVector& coord) const;

  // Resample a whole image with the transformation
  // w, h are dimensions of the new image
  // scale and offset give the ability to resize and translate the
  // output relative to the new image
  // im is the input image
  VISImage<float> resample(int w, int h, float scale, VISVector offset, 
			    const VISImage<float> &im)  const;

  // The asignment operator.  Takes care of sharing the transform and
  // copy on write if necessary
  const Warp& operator=(const Warp& other);

  // Copy constructor---allways define this for an object.
  Warp(const Warp& other) {operator=(other);}

  // Clean up the _transform, if necessary
  ~Warp();

  // Give the _transform new parameters and do copy on write if
  // necessary 
  virtual void setParams(const VISArray<float> &params);
  // Get the current set of parameters from the _transform object
  virtual VISArray<float>  getParams() const;
};

class Transform
{
protected:
  int _ref_count;
public: 
  // stuff for managing memory.  Ref count, copy-on-write, etc.
  int ref() {return(++_ref_count);}
  int refCount() {return(_ref_count);}
  int unref() {return(--_ref_count);}

  // This is the key method (virtual).  Redefine it in sub classes to
  // do what you want.  At this level it is the identity.
  virtual VISVector transform(const VISVector& coord) const {return(coord);}

  // Not much to do (no data in this class) except initialize memory
  // management 
  Transform() {_ref_count = 0;}
  // Same as above 
  Transform(const Transform &t) {_ref_count = 0;}

  // This is important for maintaining the "copy on write" protocol.
  // You must declare/define this method just like this for derived
  // classes 
  virtual Transform* newCopy() {return(new Transform(*this));}

  // Gets called by Warp (the wrapper object). 
  // Plug a new set parameters into the object
  virtual void setParams(const VISArray<float> &params){;}
  // Tell the outside world the state of your parameters.
  virtual VISArray<float> getParams() const {return VISArray<float>();}
};

// This is a special type of Warp object, that knows what type of
// transform it uses.  
// This is a convenience class---it goes ahead and creates the
// transform object and plugs it in.
// This is the object you will use most of the time (see
// testtransform.cxxx).  Either this or the CascadeWarp object below
template<class T> 
class WarpType : public Warp
{
public:
  // create and plug in the appropriate transform "T"
  WarpType() {_transform = new T(); _transform->ref(); }

  // always define this operator (it does not inherit)
  const WarpType& operator=(const WarpType& other)
  {Warp::operator=(other);}
  
  // always define this constructor (it does not inherit)
  WarpType(const WarpType& other) {operator=(other);}
};

// This is a special type of Warp object that allows you to combine
// transformations, by giving it two other Warp objects.
class CascadeWarp:  public Warp
{
  // Store the other warp objects.
protected:
  Warp _t1, _t2;
public:
  // Don't apply your own transform, but rather call two
  // transformations from the Warps you have been given. 
  virtual VISVector transform(const VISVector& coord) const
  {return(_t1.transform(_t2.transform(coord)));}

  // The constructor requires you to give it two Warp objects.  
  CascadeWarp(const Warp& T1, const Warp& T2)
  {setWarps(T1, T2);}

  void setWarps(const Warp& T1, const Warp& T2)
    {
      _t1 = T1; _t2 = T2;
      //      cout << "in set " << _t1._transform << " " << _t2._transform << endl;
    }

  // Must define this in order to get transforms from the other object
  const CascadeWarp& operator=(const CascadeWarp& other)
  {
    _t1 = other._t1;
    _t2 = other._t2;
  }

  CascadeWarp(const CascadeWarp& other) {operator=(other);}
};


// ------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------
// This is an example of how you would create a new transformation object
// ------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------

// **********************   Rotation **********************   
//
// This is a rotation transformation that allows you to choose theta (radians)
// and a center of rotation
//
// Rotation has the protocal of three params.  Theta is first, and
// then come the x and y offset
//
// New tranform objects should be derived from "Transform" but follow
// the format given here
//
class Rotation: public Transform
{
protected:
  //
  // these are data kept around that are specific to the transform
  //
  // angle of rotation in radians
  float _theta;
  // center of rotation as specified and the offset that gets added to
  // the rotated image to acheive that (review your linear algebra).
  VISVector _offset, _center;
  // keep around a rotation matrix, because it will get used a lot.
  VISMatrix _rotation;

public:
  // the object gets an array of parameters, it is up to the object to pull 
  // out the appropriate parameters and fill in the protected data above
  virtual void setParams(const VISArray<float> &params);

  // the object fills in the appropriate parameters and hands them
  // back out.
  virtual VISArray<float> getParams() const;

  // this is the key thing.  you must define this method and it should
  // perform the appropriate transformation
  virtual VISVector transform(const VISVector &coord) const
  {return(VISVector(_rotation*coord + _offset));}

  // this is necessary to make sure Cascade and things work.  
  // you must have this method, and it should look like this, except 
  // replace Rotation with the name of the class
  virtual Transform* newCopy() {cout << "new copy in rotation"; return(new Rotation(*this));}
  
};

// ********************** Perspective **********************   
//
// This is a perspective transformation
// Takes in 8 perspective projection parameters 
// + current image number + image number of the transformed space
//
class Perspective: public Transform
{
protected:
	VISMatrix _perspective;
	VISArray<float> params; // 8 p values + current image number + image number of the transformed space

public:
  // the object gets an array of parameters, it is up to the object to pull 
  // out the appropriate parameters and fill in the protected data above
  virtual void setParams(const VISArray<float> &params);

  // the object fills in the appropriate parameters and hands them
  // back out.
  virtual VISArray<float> getParams() const {return params;}

  // this is the key thing.  you must define this method and it should
  // perform the appropriate transformation
  virtual VISVector transform(const VISVector &coord) const;

  // this is necessary to make sure Cascade and things work.
  virtual Transform* newCopy() {cout << "new copy in perspective"; return(new Perspective(*this));}
  
};

// ********************** Inverse Perspective **********************   
//
// This is a inverse perspective transformation
// Takes in 8 perspective projection parameters 
// + current image number + image number of the transformed space
//
class InversePerspective: public Transform
{
protected:
	VISMatrix _INVperspective;
	VISArray<float> params; // 8 p values + current image number + image number of the transformed space

public:
  // the object gets an array of parameters, it is up to the object to pull 
  // out the appropriate parameters and fill in the protected data above
  virtual void setParams(const VISArray<float> &params);

  // the object fills in the appropriate parameters and hands them
  // back out.
  virtual VISArray<float> getParams() const {return params;}

  // this is the key thing.  you must define this method and it should
  // perform the appropriate transformation
  virtual VISVector transform(const VISVector &coord) const;

  // this is necessary to make sure Cascade and things work.
  virtual Transform* newCopy() {cout << "new copy in inverse perspective"; return(new InversePerspective(*this));}
  
};

