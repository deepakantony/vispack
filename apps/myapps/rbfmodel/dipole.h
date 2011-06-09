#include <image.h>
#include <array.h>
#include <matrix.h>
#include <param.h>

#define DEFAULT_LENGTH (3.0f)
#define DEFAULT_WEIGHT (1.0f)

#include <imageRGBA.h>

#define DOTRADIUS 3

void placeDot(VISImageRGBA &im, const VISVector &p, VISColorChannel c);


class ScalarFunction
{
public:
//    static evaluate(const VISVector &p, const ScalarFunction *F) const
//    {
//      return(F->evaluate(p));
//    }
  virtual double evaluate(const VISVector &p) const {return(0.0f);}
};

class Dipole
{
  VISVector _position, _orientation;
  double _length, _weight;
public:
  Dipole(const VISVector &position, const VISVector &orientation, 
	 double length = DEFAULT_LENGTH)
  {
    _position = position; _orientation = orientation; _length = length;
    _weight = DEFAULT_WEIGHT;
  }

  Dipole()
  {
    _position = VISVector(0.0f, 0.0f);
    _orientation = VISVector(0.0f, 0.0f);
    _weight = _length = 0.0f;
  }
	 
  VISVector inside() const {return(_position + _length*_orientation);}
  VISVector outside() const {return(_position - _length*_orientation);}
  VISVector position() const {return(_position);}
  VISVector normal() const {return(_orientation);}

  
  
  double insideConstraint() const {return(_length*_weight);}
  double outsideConstraint() const {return(-_length*_weight);}

  // use a root traping to find the zero of the given function between the two ends of the dipole
  VISVector findZero(const ScalarFunction &F) const;

  void print() const {cout << "Diple: " << _position << _orientation << endl;}
  void print(const ScalarFunction &F) const 
  {
    cout << "Diple: " << position() << ":::inside" << F.evaluate(inside()) << " middle " 
	 << F.evaluate(position()) << " outside " 
	 << F.evaluate(outside()) << endl;
  }

};

class RBF
{
public:
  static double value(const VISVector &p);
  static VISVector gradient(const VISVector &p);
  static VISMatrix hessian(const VISVector &p);
};

class ThinPlateSplineInterp: public ScalarFunction
{
  VISMatrix  _points;
  VISVector _coefficients;
  int _num_RBFs;
public:
  ThinPlateSplineInterp() {_num_RBFs = -1;}
  int numRBFs() {return(_num_RBFs);}
  virtual double evaluate(const VISVector &p) const;
  VISVector evaluateGradient(const VISVector &p);
  VISMatrix evaluateHessian(const VISVector  &p);
  double evaluateCurvature(const VISVector  &point);
  void solve(const VISMatrix &constraints, const VISVector &values);
};

#define DEFAULT_NUM_OUTSIDE_CONSTRAINTS  (16)
#define DEFAULT_OUTSIDE_CONSTRAINTS_VALUE  (-1.0f*DEFAULT_WEIGHT)

class ImplicitModel
{

  VISMatrix _outside_constraints;

  ThinPlateSplineInterp _spline;
  VISMatrix boundingBox();
  int _num_dipoles;
  VISArray< Dipole > _dipoles;
  VISArray<VISMatrix> readParams(VPF::ParameterFile F, char* prefix);
  VISMatrix readMatrix(VPF::ParameterFile pf, const char* string) const;
  static char* _keywords[3];
  VISMatrix outsideConstraints(int num) const;
  
public:

  ImplicitModel(VPF::ParameterFile F);
  VISImage<float> getImage(int w, int h);
  VISImageRGBA getImageModel(int w, int h);
  VISImage<float> getImageCurve(int w, int h);
  void printCurvatures();
  void setDipoles(const VISMatrix &points, const VISMatrix &normals);
  void setDipoles(const VISArray<Dipole> &dipoles);
  void adjustDipoles();
  void updateDipoles();
  void rebuildSpline();
  void deform();

  typedef enum {NUMDIPOLES= 0, DIPOLENORMALS = 1, 
		DIPOLEPOINTS = 2} KeyWordNum;

};


