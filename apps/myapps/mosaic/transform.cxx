#include "transform.h"


// Produces a new image that has been warped according to the tranform
// type and the input parameters.
VISImage<float> Warp::resample(int w, int h, float scale, VISVector offset, 
					  const VISImage<float> &im)  const
{
  // creat the new image
  VISImage<float>  ret(w, h);
  // coords in source and target
  VISVector x(2), x_prime(2); 

  int i, j;

  // for every pixel in the new image
  for (i = 0; i < w; i++)
    for (j = 0; j < w; j++)
      {
	// create a vector object for this pixel location (incl. scale
	// and offset
	x = scale*VISVector((float)i, (float)j) + offset;
	// transform to other image space
	x_prime = transform(x);
	// make sure it lands in the other image and interpolate
	if (im.checkBounds((float)x_prime.peek(0), (float)x_prime.peek(1)))
	  ret.poke(i, j) = im.interp((float)x_prime.peek(0), (float)x_prime.peek(1));
	else
	  ret.poke(i, j);
      }

  return(ret);
}


// If you have a valid transform object, use it to transform the
// coordinate, otherwise, use identity.
inline VISVector Warp::transform(const VISVector& coord) const 
{
  //  cout << "in Warp::transform got transform " << _transform << endl;
  if (_transform == NULL) 
    return(coord); 
  else return(_transform->transform(coord)); 
}

//
// The assignment operator gets rid of the old transform and takes the
// new one
//
const Warp& Warp::operator=(const Warp& other)
  {
    // if you have a transform object and noone else needs it, delete it
    if ((_transform != NULL)&&(_transform->unref() < 1))
      delete _transform;
    // get the transform from the other object, and increase it's ref_count
    _transform = other._transform;
    _transform->ref();
  }

 Warp::~Warp()
  {
    // if you have a transform object and noone else needs it, delete it
    if ((_transform != NULL)&&(_transform->unref() < 1))
      delete _transform;
  }

void Warp::setParams(const VISArray<float> &params)
  {
    // if someone else is sharing your transform, deref it and create new one
    if ((_transform)&&(_transform->refCount()  > 1))
      {
	_transform->unref();
	_transform = _transform->newCopy();
	_transform->ref();
      }
    // give the new params to the transform
    _transform->setParams(params);
  }

VISArray<float>  Warp::getParams() const 
  {
    if (_transform)
      return(_transform->getParams());
    else 
      return VISArray<float>();
  }

// Rotation has the protocal of three params.  Theta is first, and
// then come the x and y offset
VISArray<float> Rotation::getParams() const
{
  VISArray<float> ret;
  ret.poke(0) = _theta;
  ret.poke(1) = _center.peek(0);
  ret.poke(2) = _center.peek(1);
  return(ret);
}

// Rotation has the protocal of three params.  Theta is first, and
// then come the x and y offset
void Rotation::setParams(const VISArray<float> &params)
{
  // must have 3 parameters
  if (params.n() < 3)
    cout << "Rotation:setParams -- set parameter warning" << endl;
  _theta = params.peek(0);
  _center = VISVector(params.peek(1), params.peek(2));

  // build a rotation matrix
  _rotation = VISMatrix(2, 2);
  _rotation.poke(0, 0) = _rotation.poke(1, 1) = cos(_theta);
  _rotation.poke(0, 1) = -sin(_theta);
  _rotation.poke(1, 0) = sin(_theta);

  // this is how you rotate about a given point (review your linear
  // algebra)   
  _offset = _center  - _rotation*_center;
}

