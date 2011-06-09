#include <array.h>
#include <mathutil.h>
#include <matrix.h>

#define SQ_SIZE (1.0)
#define SQ_BOT (0.5f*SQ_SIZE)
#define SQ_TOP (-0.5f*SQ_SIZE)
#define SQ_LEFT (-0.5f*SQ_SIZE)
#define SQ_RIGHT (0.5f*SQ_SIZE)

class Model
{
public:
  virtual VISArray< VISVector > intersections(const VISMatrix &line) const
  {return VISArray< VISVector >();}
  virtual VISMatrix dXdParams(const VISVector &x) const {return VISMatrix();}
  virtual VISVector dXdParamsNormal(const VISVector &x) const {return VISVector();}
  virtual VISVector normal(const VISVector &X) const {return VISVector();}
  virtual VISVector getParams() const {return VISVector();}
  virtual void setParams(const VISVector &params) {;}
  // params are theta, x0, and y0
  virtual int numParams() const {return(0);}
};


class Square: public Model
{
public:
  VISMatrix _transform;
  VISMatrix _transform_inv;
  float _theta, _x0, _y0, _scale;
  VISMatrix _points;
  virtual VISArray< VISVector > intersections(const VISMatrix &line) const;
  virtual VISMatrix dXdParams(const VISVector &x) const;
  //virtual VISVector dXdParamsNormal(const VISVector &x) const;
  virtual VISVector normal(const VISVector &X) const;
  virtual VISVector getParams() const
    {
      VISVector ret(4);
      ret.poke(0) = _theta;
      ret.poke(1) = _x0;
      ret.poke(2) = _y0;
      ret.poke(3) = _scale;
      return(ret);
    }
  virtual void setParams(const VISVector &params);
  // params are theta, x0, and y0
  virtual int numParams() const {return(4);}
  Square();
  typedef enum {THETA=0, X0=1, Y0=2, SCALE=3} Params;
};

class Projection
{
public:
  int _size;
  virtual VISMatrix lineOfSight(int i) const {return VISMatrix();}
  VISVector project(const Square &sq) const;
  float getProjection(const VISArray< VISVector > &inters) const;
  Projection(int size) {_size = size;}
  Projection() {_size = 0;}
  int n() const {return(_size);}
};

class PerspectiveProjection: public Projection
{
public:
  float _origin;
  float _focal_length;
  // returns a line of sight in terms of cosine angles and origin
  virtual VISMatrix lineOfSight(int i) const;
  PerspectiveProjection(int size, int origin, float f):Projection(size) {
  _origin = origin; _focal_length = f;}
  PerspectiveProjection():Projection()  { 
  _origin = 0.0f; 
  _focal_length = 0.0f;}
  PerspectiveProjection(int size):Projection(size) {
  _origin = _size/2; 
  _focal_length = _size/2;
  }
};

class OrthoProjection: public Projection
{
public:
  float _origin;
  // returns a line of sight in terms of cosine angles and origin
  virtual VISMatrix lineOfSight(int i) const;
  OrthoProjection(int size, int origin):Projection(size) {
    _origin = origin;}
  OrthoProjection():Projection()  { 
  _origin = 0.0f; 
}
  OrthoProjection(int size):Projection(size) {
  _origin = _size/2; 
  }
};


class System
{
protected:
  int _current_scale;
  int _num_scales;
  VISArray< float > _errors;
public:
  virtual VISVector getParams() const {return VISVector();}
  virtual void setParams(const VISVector &p) {;}
  virtual int numPoints() const {return(0);}
  virtual VISVector dEdParams(int i) const {cout << "got parent dE" << endl; 
return VISVector();}
  virtual VISVector error() const {return VISVector();}
  virtual float errorSq() const {return(0.0f);}
  virtual float errorX(int i) const {return(0.0f);};
  virtual int numScales() const {return(0);}
  virtual void setScale(int i) {;}
};


template <class TProjection, class TModel>
class SystemFitModel: public System
{
protected:
  TProjection _projection;
  TModel _model;
  VISVector _data, _errors;
public:
  SystemFitModel(const TModel &model, 
			 const TProjection &p, 
			 const VISVector &d)
    {
      _model = model;  _projection = p; _data = d;
      _errors = _projection.project(_model) - _data;
    }
  virtual float errorSq() const;
  // right now only the square moves
  virtual VISVector getParams() const {return(_model.getParams());}
  virtual void setParams(const VISVector &p)
    {
      _model.setParams(p);
      _errors = _projection.project(_model) - _data;
    }
  virtual int numPoints() const {return(_projection.n());}
  virtual VISVector dEdParams(int i) const;
  virtual VISVector error() const {return(_errors);}
  virtual float errorX(int i) const
    {
      return(_errors.peek(i));
    }
  virtual int numScales() const {return(1);}
  virtual void setScale(int i) {_current_scale = 0;}
};


class SystemProjectionSquare: public System
{
protected:
  OrthoProjection _projection;
  Square _square;
  VISVector _data, _errors;
public:
  SystemProjectionSquare(const Square &sq, 
			 const OrthoProjection &p, 
			 const VISVector &d)
    {
      _square = sq;  _projection = p; _data = d;
      _errors = _projection.project(_square) - _data;
    }
  virtual float errorSq() const;
  // right now only the square moves
  virtual VISVector getParams() const {return(_square.getParams());}
  virtual void setParams(const VISVector &p)
    {
      _square.setParams(p);
      _errors = _projection.project(_square) - _data;
    }
  virtual int numPoints() const {return(_projection.n());}
  virtual VISVector dEdParams(int i) const;
  virtual VISVector error() const {return(_errors);}
  virtual float errorX(int i) const
    {
      return(_errors.peek(i));
    }
  virtual int numScales() const {return(1);}
  virtual void setScale(int i) {_current_scale = 0;}
};


#define LAMBDA_MIN (1.0e-6)
#define LAMBDA_MAX (1.0e10)
#define ITERATIONS_MAX (10)

class Solver
{
protected:
  float _lambda, _lambda_init, _lambda_max, _lambda_min; 
public:
  Solver()
  {
    _lambda_init = 100.0; 
    _lambda_max = LAMBDA_MAX; 
    _lambda_min = LAMBDA_MIN; 
  }
  //  Solver(VPF::ParameterFile pf);
  void reset(float lambda, int scale) {_lambda = lambda;}
  int iterate(System &sys);
  int solve(System &sys);
};

#include <projection.txx>

