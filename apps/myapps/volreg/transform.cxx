#include "transform.h"


// Produces a new image that has been warped according to the tranform
// type and the input parameters.
VISVolume<float> Warp::resample(int w, int h, int d, float scale, 
				const VISVector &origin_target, 
				const VISVector &origin_source, 
				const VISVolume<float> &im)  const
{
  // creat the new volume
  VISVolume<float>  ret(w, h, d);
  // coords in source and target
  VISVector x(3), x_prime(3); 

  int i, j, k;

  // for every pixel in the new volume
  for (k = 0; k < d; k++)
  //  for (k = d/2 -1 ; k < d/2 + 1; k++)
    {
      for (i = 0; i < w; i++)
	for (j = 0; j < h; j++)
	  {
	    // create a vector object for this pixel location (incl. scale
	// and offset
	    x = scale*(VISVector((float)i, (float)j, (float)k)- origin_target);
	// transform to other volume space
	    x_prime = transform(x) + origin_source;
	    //	    if ((i == j)&&(j==k))
	    //	      {
	    //		cout << x << endl;
	    //		cout << x_prime << endl;
	    //	      }
	// make sure it lands in the other volume and interpolate
	if (im.checkBounds((float)x_prime.peek(0), 
			   (float)x_prime.peek(1), 
			   (float)x_prime.peek(2)))
	  ret.poke(i, j, k) = im.interp((float)x_prime.peek(0), 
					(float)x_prime.peek(1), 
					(float)x_prime.peek(2));
	else
	  ret.poke(i, j, k) = 0.0f;
	  }
      cout << "done slice " << k << endl;
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
    return(*this);
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
VISArray<float> Affine::getParams() const
{
  int i, j, k = 0;
  VISArray<float> ret;
  for (i = 0; i < 3; i++)
    {
      for (j = 0; j < 3; j++)
	ret.poke(k++) = _matrix.peek(i, j);
      ret.poke(k++) = _trans.peek(i);
    }
  return(ret);
}

// Rotation has the protocal of three params.  Theta is first, and
// then come the x and y offset
void Affine::setParams(const VISArray<float> &params)
{
  _matrix = VISMatrix(3,3);
  _trans = VISVector(3);
  // must have 3 parameters
  int i, j, k = 0;
  if (params.n() < 12)
    cout << "Rotation:setParams -- set parameter warning" << endl;
  else
    for (i = 0; i < 3; i++)
      {
	for (j = 0; j < 3; j++)
	  _matrix.poke(i, j) = params.peek(k++);
	_trans.poke(i) = params.peek(k++);
      }
  cout << _matrix << endl;
  cout << _trans << endl;
}

