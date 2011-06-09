#include "dipole.h"

double RBF::value(const VISVector &p)
{ 
  double r2 = p.dot(p); 
  if (r2 >0.0) 
    return r2*log(r2); 
  else 
    return(0.0f);
}

VISVector RBF::gradient(const VISVector &p)
{
  double r2 = p.dot(p), r = sqrt(r2);
  VISVector ret(2);
  if (!(r > 0.0)) return((ret = 0.0f));
  return(p*(2.0f*(1.0f + log(r2))));
}

VISMatrix RBF::hessian(const VISVector &p)
{
  double r2 = p.dot(p), r = sqrt(r2);
  VISMatrix ret;
  if (!(r > 0)) 
    {
      ret = VISMatrix(2,2);
      return((ret = 0.0f));
    }
  ret = p.exterior(p)*(2.0f/r2) + VISIdentity(2)*(1.0f + log(r2));
  return(ret);
}

void ThinPlateSplineInterp::solve(const VISMatrix &points, const VISVector &values)
{
  if ((points.c() != values.n())&&(points.r() != 2))
    {
      cout << "Warning:  ThinPlateSplineInterp:solve - constraints not sized correctly" << endl;
      return;
    }
  _num_RBFs = points.c();
  _points = points;
  int i, j, k;
  VISMatrix A(_num_RBFs + 3, _num_RBFs + 3);
  VISVector b(_num_RBFs + 3);
  for (i = 0; i < _num_RBFs; i++)
    {
      b.poke(i) = values.peek(i);
      for (j = 0; j < _num_RBFs; j++)
	if (i != j)
	  A.poke(j, i) = RBF::value(points.vec(i) - points.vec(j));
	else 
	  A.poke(j, i) = 0.0f;
      
      A.poke(j, i) = points.peek(0, i);
      A.poke(j + 1, i) = points.peek(1, i);
      A.poke(j + 2, i) = 1.0f;
    }

  for (i = 0; i < _num_RBFs; i++)
    {
      A.poke(i, _num_RBFs) = points.peek(0, i);
      A.poke(i, _num_RBFs + 1) = points.peek(1, i);
      A.poke(i, _num_RBFs + 2) = 1.0f;
    }
  for (i = _num_RBFs; i < _num_RBFs + 3; i++)
    {
      A.poke(i, _num_RBFs) =  A.poke(i, _num_RBFs + 1) = A.poke(i, _num_RBFs + 2) = 0.0f;
    }
  //  _coefficients = solvLinearEqn(A, b);
  //  cout << "solve residuals SVD" << (A*_coefficients - b) << endl;
  _coefficients = A.inverseGJ()*b;
  cout << "solve residuals GJ" << (A*_coefficients - b) << endl;
  //  cout << _coefficients << endl;
}

double ThinPlateSplineInterp::evaluate(const VISVector &p) const
{
  double total = 0.0f;
  int i;
  for (i = 0; i < _num_RBFs; i++)
    total += _coefficients.peek(i)*RBF::value(p - _points.vec(i));
  total += _coefficients.peek(i++)*p.peek(0);
  total += _coefficients.peek(i++)*p.peek(1);
  total += _coefficients.peek(i++);
  return(total);
}

VISVector ThinPlateSplineInterp::evaluateGradient(const VISVector &p)
{
  VISVector total(2);
  total = 0.0f;
  int i;
  for (i = 0; i < _num_RBFs; i++)
    total += _coefficients.peek(i)*RBF::gradient(p - _points.vec(i));
  total.poke(0) += _coefficients.peek(i++);
  total.poke(1) += _coefficients.peek(i++);
  return(total);
}

VISMatrix ThinPlateSplineInterp::evaluateHessian(const VISVector &p)
{
  VISMatrix total(2, 2);
  total = 0.0f;
  int i;
  for (i = 0; i < _num_RBFs; i++)
    total += _coefficients.peek(i)*RBF::hessian(p - _points.vec(i));
  return(total);
}

double ThinPlateSplineInterp::evaluateCurvature(const VISVector &point)
{
  VISVector grad = evaluateGradient(point);
  VISMatrix P, hessian = evaluateHessian(point);
  P = VISIdentity(2) - grad.exterior(grad)/grad.dot(grad);
  return((P*hessian*P).trace()/grad.norm(10.0e-6));
}


char* ImplicitModel::_keywords[3] = 
{
  "NUM_DIPOLES", "DIPOLE_NORMALS", "DIPOLE_POINTS", 
};


  #define TRYKEYWORDREQ(var, pref, key) { \
sprintf(keyword, "%s%s", pref, _keywords[key]); \
      if (VPF::set(var, pf[keyword][0]) != VPF::VALID) { \
      cout << "Scan.C::myload -- missing required keyword : " << keyword << endl; \
      exit(-1); } \
}


VISMatrix ImplicitModel::readMatrix(VPF::ParameterFile pf, 
const char* string) const
{
  int r, c, j;
  VISMatrix matrix;
  double double_tmp;
  try
    {
      VPF::set(r, pf[string][0]);
      VPF::set(c, pf[string][1]);
      matrix = VISMatrix(r,c);
      for (j = 0; j < r*c; j++)
	{
	  VPF::set(double_tmp,
		   pf[string][j+2]);
	  matrix.poke(j/c, j%c) = double_tmp;
	}
    }
  catch (VPF::Exception e)
    {
      // igore matrix
      cout << "ImplicitModel::bad matrix "
	   << endl;
    }
  return(matrix);
}

VISArray< VISMatrix > ImplicitModel
::readParams(VPF::ParameterFile pf, char* prefix = "")
{
  char keyword[80];
  VISArray< VISMatrix > ret;
  VISMatrix points, normals;
  TRYKEYWORDREQ(_num_dipoles, prefix, NUMDIPOLES);
  sprintf(keyword, "%s%s", prefix, _keywords[DIPOLEPOINTS]);  
  points = (readMatrix(pf, keyword)).t();
  sprintf(keyword, "%s%s", prefix, _keywords[DIPOLENORMALS]);  
  normals = (readMatrix(pf, keyword)).t();
  if (!points.compareSize(normals))
    cout << "Error: ImplicitModel - read params points and normals don't match" << endl;
  ret.poke(0) = points;
  ret.poke(1) = normals;
  return(ret);
}

ImplicitModel::ImplicitModel(VPF::ParameterFile F)
{
  VISArray<VISMatrix> data = readParams(F);
  setDipoles(data.peek(0), data.peek(1));
}

void ImplicitModel::setDipoles(const VISMatrix &points, 
			       const VISMatrix &normals)
{
  VISArray<Dipole> dipoles;
  int num_dipoles = points.c();
  int i;
  for (i = 0; i < num_dipoles; i++)
    dipoles.poke(i) = Dipole(points.col(i), normals.col(i));
  setDipoles(dipoles);
}

void ImplicitModel::setDipoles(const VISArray<Dipole> &dipoles) 
{
  int i;
  _num_dipoles = dipoles.n();
  _dipoles = dipoles;
  rebuildSpline();
  //  adjustDipoles();
  //  rebuildSpline();
}

void ImplicitModel::rebuildSpline()
{
  int i;
  int num_outside_constraints = DEFAULT_NUM_OUTSIDE_CONSTRAINTS;
  double outside_value = DEFAULT_OUTSIDE_CONSTRAINTS_VALUE;

  VISMatrix constraint_points(2, _num_dipoles*2 + num_outside_constraints);
  VISVector constraint_values(_num_dipoles*2 + num_outside_constraints);
  //  cout << "dipoles before fitting:" << endl;
  //  if (_spline.numRBFs() > 0)
  //    for (i = 0; i < _num_dipoles; i++)
  //      (_dipoles.peek(i)).print(_spline);

  for (i = 0; i < _num_dipoles; i++)
    {
      //      (_dipoles.peek(i)).print();
      constraint_points.pokeROI(0, 2*i,(_dipoles.peek(i).inside()));
      constraint_values.poke(2*i) = (_dipoles.peek(i).insideConstraint());
      constraint_points.pokeROI(0, 2*i +1,(_dipoles.peek(i).outside()));
      constraint_values.poke(2*i+1) = (_dipoles.peek(i).outsideConstraint());
    }
  
  // set the outside constraints
  for (i = 0; i < num_outside_constraints; i++)
    constraint_values.poke(i+2*_num_dipoles) = outside_value;
  constraint_points.pokeROI(0, 2*_num_dipoles, 
			 outsideConstraints(num_outside_constraints));    

  //  cout << "constraint_points: " << constraint_points.t() << endl;
  //  cout << "constraint_values: " << constraint_values.t() << endl;
  _spline.solve(constraint_points, constraint_values);

  //  cout << "dipoles after fitting:" << endl;
  //  if (_spline.numRBFs() > 0)
  //  for (i = 0; i < _num_dipoles; i++)
  //      (_dipoles.peek(i)).print(_spline);
}


#define IMAGE_BB_PADDING_PERCENT (0.1)

VISImage<float> ImplicitModel::getImage(int w, int h)
{
  // construct a bounding box
  double min_x= FLT_MAX, min_y = FLT_MAX, 
    max_x = -FLT_MAX, max_y = -FLT_MAX;
  double width, height;
  double x, y;
  VISVector position;
  int i, j;
  
  for (i = 0; i < _num_dipoles; i++)
    {
      position = (_dipoles.peek(i)).position();
      x = position.peek(0);
      y = position.peek(1);
      max_x = VISmax(max_x, x);
      max_y = VISmax(max_y, y);
      min_x = VISmin(min_x, x);
      min_y = VISmin(min_y, y);
    }

  width = max_x - min_x;
  height = max_y - min_y;

  min_x -= width*IMAGE_BB_PADDING_PERCENT; 
  min_y -= height*IMAGE_BB_PADDING_PERCENT; 
  max_x += width*IMAGE_BB_PADDING_PERCENT; 
  max_y += height*IMAGE_BB_PADDING_PERCENT; 

  VISImage<float> ret(w, h);

  double scale = VISmin((double)w/(double)width, (double)h/(double)height)
    /(1.0f + 2.0f*IMAGE_BB_PADDING_PERCENT);
  
  double x_offset, y_offset;

  //  x_offset = -scale*min_x;
  //  y_offset = -scale*min_y;
  

  for (i = 0; i < w; i++)
    for (j = 0; j < h; j++)
      {
	ret.poke(i, j) = _spline.evaluate(VISVector(i/scale + min_x, j/scale + min_y));
      }
  return(ret);
}


VISImageRGBA ImplicitModel::getImageModel(int w, int h)
{
  // construct a bounding box
  double min_x= FLT_MAX, min_y = FLT_MAX, 
    max_x = -FLT_MAX, max_y = -FLT_MAX;
  double width, height;
  double x, y;
  VISVector position;
  int i, j;
  
  for (i = 0; i < _num_dipoles; i++)
    {
      position = (_dipoles.peek(i)).position();
      x = position.peek(0);
      y = position.peek(1);
      max_x = VISmax(max_x, x);
      max_y = VISmax(max_y, y);
      min_x = VISmin(min_x, x);
      min_y = VISmin(min_y, y);
    }

  width = max_x - min_x;
  height = max_y - min_y;

  min_x -= width*IMAGE_BB_PADDING_PERCENT; 
  min_y -= height*IMAGE_BB_PADDING_PERCENT; 
  max_x += width*IMAGE_BB_PADDING_PERCENT; 
  max_y += height*IMAGE_BB_PADDING_PERCENT; 

  VISImageRGBA ret(w, h);
  VISImage<float> im_float(w, h);

  double scale = VISmin((double)w/(double)width, (double)h/(double)height)
    /(1.0f + 2.0f*IMAGE_BB_PADDING_PERCENT);
  
  double x_offset, y_offset;

  //  x_offset = -scale*min_x;
  //  y_offset = -scale*min_y;
  

  for (i = 0; i < w; i++)
    for (j = 0; j < h; j++)
      {
	im_float.poke(i, j) = _spline.evaluate(VISVector(i/scale + min_x, j/scale + min_y));
	//	if (i == j)
	//	  cout << "eval vector" << VISVector(i/scale + min_x, j/scale + min_y) << endl;
      }

  
  im_float = im_float.zeroCrossings();
  float the_min = im_float.min(), the_max = im_float.max();
  im_float = (255.0f/(the_max - the_min))*(im_float - the_min);
  ret = VISImage<rgba>(im_float);

  for (i = 0; i < _num_dipoles; i++)
    {
      position = scale*((_dipoles.peek(i)).outside() - VISVector(min_x, min_y));
      //      cout << position << endl;
      placeDot(ret, position, R);
      position = scale*((_dipoles.peek(i)).inside() - VISVector(min_x, min_y));
      //      cout << position << endl;
      placeDot(ret, position, B);
    }
  return(ret);
}


VISImage<float> ImplicitModel::getImageCurve(int w, int h)
{
  // construct a bounding box
  double min_x= FLT_MAX, min_y = FLT_MAX, 
    max_x = -FLT_MAX, max_y = -FLT_MAX;
  double width, height;
  double x, y;
  VISVector position;
  int i, j;
  
  for (i = 0; i < _num_dipoles; i++)
    {
      position = (_dipoles.peek(i)).position();
      x = position.peek(0);
      y = position.peek(1);
      max_x = VISmax(max_x, x);
      max_y = VISmax(max_y, y);
      min_x = VISmin(min_x, x);
      min_y = VISmin(min_y, y);
    }

  width = max_x - min_x;
  height = max_y - min_y;

  min_x -= width*IMAGE_BB_PADDING_PERCENT; 
  min_y -= height*IMAGE_BB_PADDING_PERCENT; 
  max_x += width*IMAGE_BB_PADDING_PERCENT; 
  max_y += height*IMAGE_BB_PADDING_PERCENT; 

  VISImage<float> ret(w, h);

  double scale = VISmin((double)w/(double)width, (double)h/(double)height)
    /(1.0f + 2.0f*IMAGE_BB_PADDING_PERCENT);
  
  double x_offset, y_offset;

  //  x_offset = -scale*min_x;
  //  y_offset = -scale*min_y;

  VISVector grad, point;
  double grad_mag_sq;
  VISMatrix hessian, P;
  

  for (i = 0; i < w; i++)
    for (j = 0; j < h; j++)
      {
	ret.poke(i, j) = 
	  _spline.evaluateCurvature(VISVector(i/scale + min_x, j/scale + min_y));
      }
  return(ret);
}

void ImplicitModel::printCurvatures()
{
  int i;
  VISVector position, orientation;
  cout << "Implicit model curvatures at dipole centers:" << endl;
  for (i = 0; i < _num_dipoles; i++)
    {
      position = (_dipoles.peek(i)).position();
      orientation = (_dipoles.peek(i)).normal();
      cout << "pos: " << position.peek(0) << " : " << position.peek(1) << "  orient: " << 
	orientation.peek(0) << " : "  << orientation.peek(1) << " curve:  " <<
	_spline.evaluateCurvature(position) << endl;
	
    }
  
}



VISVector Dipole::findZero(const ScalarFunction &F) const
{
  double distance;
  VISVector p1 = inside(), p2 = outside(), p3;
  double v1 = F.evaluate(p1), v2 =  F.evaluate(p2), v3, range = fabs(v1-v2);

  //  cout << "in find zero p1 " << p1 << endl;
  //  cout << "in find zero v1 " << v1 << endl;
  //  cout << "in find zero p2 " << p2 << endl;
  //  cout << "in find zero v2 " << v2 << endl;

  // special input cases
  if (v1==v2)
    return VISVector();
  if (v1==0.0) return(p1);
  if (v2==0.0) return(p2);

  boolean done = false;
  int iterations = 0;
  // general case
  while(!done)
    {
      p3 = p1 + (p2-p1)*v1/(v1 - v2);
      //      cout << "in find zero p3 " << p3 << endl;
      v3 =  F.evaluate(p3);
      //      cout << "v3 and range " << v3 << " " << range << " v1 " << v1 << " v2 " << v2 << endl;
      if ((fabs(v3) < 0.01*range)||(iterations++ > 10)) 
	done = TRUE;
      else
	if (v3*v1 < 0.0f)
	  {
	    p2 = p3; 
	    v2 = v3;
	  }
	else 
	  {
	    p1 = p3; 
	    v1 = v3;
	  }
    }

  //  cout << "in find zero p3 " << p2 << endl;

  return(p3);
}

void ImplicitModel::adjustDipoles()
{
  int i;
  Dipole dipole;
  VISVector position, normal;

  for (i = 0; i < _num_dipoles; i++)
    {
      dipole = _dipoles.peek(i);
      //      dipole.print();
      position = dipole.findZero(_spline);
            normal = _spline.evaluateGradient(position);
      normal /= normal.norm();
      _dipoles.poke(i) = Dipole(position, normal);
    }
}


VISMatrix ImplicitModel::outsideConstraints(int num) const
{

  VISVector middle(0.0, 0.0), farthest;
  double farthest_dist = 0.0, this_dist;
  VISMatrix ret = VISMatrix(2, num);
  int i;


  // find the radius and position of the smallest circle that encloses all
  // (centered at mean)
  for (i = 0; i < _num_dipoles; i++)
    middle += (_dipoles.peek(i)).position();
  middle /= _num_dipoles;
  for (i = 0; i < _num_dipoles; i++)
    if ((this_dist = ((_dipoles.peek(i)).position() - middle).norm()) > farthest_dist)
      {
	farthest = (_dipoles.peek(i)).position();
	farthest_dist = this_dist;
      }

  // set the radius to be twice of what it needs to be  
  double radius = farthest_dist*3.0f;

  for (i = 0; i < num; i++)
    ret.pokeROI(0, i, radius*VISVector(cos(2.0*M_PI*(double)i/(double)num), 
				       sin(2.0*M_PI*(double)i/(double)num))
		+ middle);
  return ret;
}

void ImplicitModel::updateDipoles()
{
  int i;
  Dipole dipole;
  VISVector position, normal;
  double curve, dt = 400.0;

  for (i = 0; i < _num_dipoles; i++)
    {
      dipole = _dipoles.peek(i);
      position = dipole.position();
      //      normal = dipole.normal();
      normal = _spline.evaluateGradient(position);
      normal  /= normal.norm(10.0e-6);
      curve = _spline.evaluateCurvature(position);
      position -= (dt*curve)*normal;
      normal = _spline.evaluateGradient(position);
      normal /= normal.norm(10.0e-6);
      _dipoles.poke(i) = Dipole(position, normal);
    }
}

void ImplicitModel::deform()
{
  //    adjustDipoles();
  //    rebuildSpline();
    updateDipoles();
    rebuildSpline();
    //    adjustDipoles();
    //    rebuildSpline();
}


void placeDot(VISImageRGBA &im, const VISVector &p, VISColorChannel c)
{
  int i, j;
  int u = (int)rint(p.peek(0)), v = (int)rint(p.peek(1));
  for (i =0; i < DOTRADIUS; i++)
    for (j =0; j < DOTRADIUS; j++)
      {
	im.poke(u+i, v+j) = rgba(0);
	im.put(u+i, v+j, c, 255);
      }
}

