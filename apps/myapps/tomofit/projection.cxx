#include "projection.h"


VISVector Square::normal(const VISVector &X) const
{
  VISVector pt_tmp = _transform_inv*X;
  float x = pt_tmp.peek(0), y = pt_tmp.peek(1);
  // must fix the scaling parameter here!!!!!!!!
  if (fabs(x) >= fabs(y))
    return _transform*(VISVector((float)tsign(x), 0.0f, 0.0f))/_scale;
  else
    return _transform*(VISVector(0.0f, (float)tsign(y), 0.0f))/_scale;
}

VISMatrix Square::dXdParams(const VISVector &x) const
{
  VISMatrix ret(3, 4);
  VISMatrix dR(3, 3), R(3, 3);
  VISVector X = _transform_inv*x;

  R.poke(0,0) = R.poke(1,1) = cos(_theta);
  R.poke(1,0) = sin(_theta);
  R.poke(0,1) = -sin(_theta);
  R.poke(2, 0) = R.poke(2, 1) = 0.0f;
  // zero or one here...?
  R.poke(2, 2) = 1.0f;

  dR.poke(0, 0) = (dR.poke(1, 1) =  -_scale*sin(_theta)); 
  dR.poke(1,0) = _scale*cos(_theta);
  dR.poke(0,1) = -_scale*cos(_theta);
  dR.poke(2, 0) = dR.poke(2, 1) = 0.0f;
  // zero or one here...?
  dR.poke(2, 2) = 0.0f;
  ret.pokeROI(0, 0, dR*X);
  ret.pokeROI(0, 1, VISVector(1.0f, 0.0f , 0.0f));
  ret.pokeROI(0, 2, VISVector(0.0f, 1.0f, 0.0f));
  ret.pokeROI(0, 3, R*X);
  return(ret);
}

Square::Square()
{
  _points = VISMatrix(3, 4);
  _points.poke(0, 0) = SQ_LEFT;
  _points.poke(0, 1) = SQ_RIGHT;
  _points.poke(0, 2) = SQ_RIGHT;
  _points.poke(0, 3) = SQ_LEFT;

  _points.poke(1, 0) = SQ_TOP;
  _points.poke(1, 1) = SQ_TOP;
  _points.poke(1, 2) = SQ_BOT;
  _points.poke(1, 3) = SQ_BOT;

  _points.poke(2, 0) =_points.poke(2, 1) = _points.poke(2, 2) = 
    _points.poke(2, 3) = 1.0f;
}


VISArray< VISVector > Square::intersections(const VISMatrix &line) const
{
  // transform the line into the object space  
  VISMatrix line_tmp = _transform_inv*line;

  //  cout << "line" << line << endl;
  //  cout << "line_tmp" << line_tmp << endl;

  float lambda[2], lambda_min, lambda_max, lambda_tmp;
  int intersection = 0;
  float coord_tmp;

  //  check the four sides for intersections
  if (line_tmp.peek(0, 0) != 0)
    {
      // solve for x = SQ_LEFT;
      //      cout << (SQ_LEFT-line_tmp.peek(0, 1))/line_tmp.peek(0, 0) << endl;
      lambda_tmp = (SQ_LEFT-line_tmp.peek(0, 1))/line_tmp.peek(0, 0);
      coord_tmp = line_tmp.peek(1, 0)*lambda_tmp + line_tmp.peek(1, 1);
      //      cout << "coord_tmp x" << coord_tmp << endl;
      //      cout << "lambda tmp " << lambda_tmp << endl;
      if ((coord_tmp >= SQ_TOP)&&(coord_tmp <= SQ_BOT))
	{
	  lambda[intersection] = lambda_tmp;
	  intersection++;
	}
      lambda_tmp = (SQ_RIGHT-line_tmp.peek(0, 1))/line_tmp.peek(0, 0);
      coord_tmp = line_tmp.peek(1, 0)*lambda_tmp + line_tmp.peek(1, 1);
      if ((coord_tmp >= SQ_TOP)&&(coord_tmp <= SQ_BOT))
	{
	  lambda[intersection] = lambda_tmp;
	  intersection++;
	}
    }
  if ((intersection < 2)&&(line_tmp.peek(1, 0) != 0))
    {
      // solve for y = 0
      lambda_tmp = (SQ_TOP-line_tmp.peek(1, 1))/line_tmp.peek(1, 0);      
      coord_tmp = line_tmp.peek(0, 0)*lambda_tmp + line_tmp.peek(0, 1);
      //      cout << "coord_tmp y " << coord_tmp << endl;
      //      cout << "lambda tmp " << lambda_tmp << endl;
      if ((coord_tmp >= SQ_LEFT)&&(coord_tmp <= SQ_RIGHT))
	{
	  lambda[intersection] = lambda_tmp;
	  intersection++;
	}
      if (intersection < 2)
	{
	  lambda_tmp = (SQ_BOT - line_tmp.peek(1, 1))/line_tmp.peek(1, 0);
	  coord_tmp = line_tmp.peek(0, 0)*lambda_tmp + line_tmp.peek(0, 1);
	  if ((coord_tmp >= SQ_LEFT)&&(coord_tmp <= SQ_RIGHT))	  
	    {
	      lambda[intersection] = lambda_tmp;
	      intersection++;
	    }
	}
    }
  
  if (intersection == 0)
    return VISArray< VISVector >();
  VISArray< VISVector > ret;
  VISVector pt(3);
  lambda_min = VISmin(lambda[0], lambda[1]);
  lambda_max = VISmax(lambda[0], lambda[1]);
  //  cout << "lambdas " << lambda[0] << " " <<  lambda[1] << endl;
  ret.poke(0) = _transform*(lambda_min*line_tmp.vec(0) + line_tmp.vec(1));
  ret.poke(1) = _transform*(lambda_max*line_tmp.vec(0) + line_tmp.vec(1));
  //  cout << "intersections" << endl;
  //  cout << ret.poke(0) << endl;
  //  cout << ret.poke(1) << endl;
  return(ret);
}

void Square::setParams(const VISVector &params)
{
  if (params.n() != 4)
    {
      cout << "Warning: Square::setParams---bad number of params" << endl; 
   }
  _theta = params.peek(0);
  _x0 = params.peek(1);
  _y0 = params.peek(2);
  _scale = params.peek(3);

  _transform = VISMatrix(3, 3);
  
  _transform.poke(0,0) = _transform.poke(1,1) = _scale*cos(_theta);
  _transform.poke(1,0) = _scale*sin(_theta);
  _transform.poke(0,1) = -_scale*sin(_theta);
  
  _transform.poke(2, 0) = _transform.poke(2, 1) = 0.0f;

  _transform.poke(0, 2) = _x0;
  _transform.poke(1, 2) = _y0;
  _transform.poke(2, 2) = 1.0f;

  _transform_inv = _transform.inv();

  //  cout << _transform << endl;
  //  cout << _transform_inv << endl;

}

VISMatrix PerspectiveProjection::lineOfSight(int i) const
{
  VISMatrix ret(3, 2);
  ret = 0.0f;
  float x = (i - _origin), 
    hyp = sqrt(x*x + _focal_length*_focal_length);
  ret.poke(0, 0) = x/hyp;
  ret.poke(1, 0) = _focal_length/hyp;
  ret.poke(2, 1) = 1.0f;
  return(ret);
}

VISMatrix OrthoProjection::lineOfSight(int i) const
{
  VISMatrix ret(3, 2);
  ret.poke(0, 0) = 0.0f;
  ret.poke(1, 0) = 1.0f;
  ret.poke(2, 0) = 0.0f;
  ret.poke(0, 1) = (float)i - (float)_origin;
  ret.poke(1, 1) = 0.0f;
  ret.poke(2, 1) = 1.0f;
  return(ret);
}

float Projection::getProjection(const VISArray< VISVector > &inters) const
{
  float total = 0.0f;
  if ((inters.n()%2) != 0)
    {
      cout << "Projection::getProjection - bad number of intersections" << endl;
      return(0.0f);
    }
  int n = inters.n()/2;
  for (int i = 0; i < n; i++)
    {
      //      cout << "projection pts" << inters.peek(2*i)  << inters.peek(2*i+1) << endl;
      //      cout << "projection length" << (inters.peek(2*i) - inters.peek(2*i + 1)).norm() << endl;
      total += (inters.peek(2*i) - inters.peek(2*i + 1)).norm();
    }
  return(total);
}


VISVector Projection::project(const Square &sq) const
{
  int i;
  VISVector ret(_size);
  VISMatrix los;
  for (i = 0; i < _size; i++)
    {
      los = lineOfSight(i);
      //      cout << "about to do " << i << endl;
      //      cout << los << endl;
      ret.poke(i) = getProjection(sq.intersections(los));
    }
  return(ret);
}

int Solver::iterate(System &sys)
{
  VISVector params_old, params_new, dX;
  float lambda_old = _lambda;
  
  int i, j;

  params_old = sys.getParams();
  int n = params_old.n();
  VISMatrix hessian(n, n), ident = VISIdentity(n);
  VISVector grad(n), difference(n);
  hessian = 0.0f;
  grad = 0.0f;
  VISVector delta_params;

  int num_points = sys.numPoints();
  float old_error, new_error;

  cout << "orig errors " << sys.error() << endl;

  //  cout << "in iterate" << num_points << endl;
  // build the hessian and the grad vector;
  for (i = 0; i < num_points; i++)
    {
      dX = sys.dEdParams(i);
      // doube check this!!!!!!!!!!!!!!!!!!!!!!!
      grad -= dX*sys.errorX(i);
      hessian += dX.exterior(dX);
    }

  old_error = sys.errorSq();
  cout << "old error " << old_error << endl;
  sys.setParams(params_new = params_old 
		+ (delta_params = solvLinearEqn(hessian+_lambda*ident, grad)));
  new_error =sys.errorSq();
  cout << "grad " << grad << endl;
  cout << "lambda " << _lambda << endl;
  cout << "delta params " << delta_params << endl;
  cout << "params new " << params_new << endl; 	
  cout << "new error " << new_error << endl;

  boolean done = false;
  if (new_error > old_error)
    while (!done)
      {
	_lambda *= 2.0f;
	sys.setParams(params_new = params_old 
		      + (delta_params = solvLinearEqn(hessian+_lambda*ident, grad)));
	new_error = sys.errorSq();
	cout << "new error " << new_error << endl;
	cout << "delta params " << delta_params << endl;	
	cout << "params new " << params_new << endl; 	
	cout << "lambda " << _lambda << endl;
	cout << "errors " << sys.error() << endl;
	if ((new_error < old_error)||(_lambda >= _lambda_max))
	  done = true;
      }
  else
    _lambda = VISmax(_lambda/2.0f, _lambda_min);
//    while (!done)
//      {
//        _lamda = VISmax(_lambda/2.0f, _lambda_min);
//        sys.setParams(params_new = solveLinearEqn(hessian+_lambda*ident, grad));
//        new_error = sys.error();
//        if (new_error > old_error)
//  	{
//  	  done = true();
//  	  _lambda *= 2.0f;
//  	  sys.setParams(solveLinearEqn(hessian+_lambda*ident, grad));
//  	}
//        if (_lambda == _lambda_min)
//  	done = true();      
//      }
  if (_lambda >= _lambda_max)
    return(1);
  else
    return(0);
}

int Solver::solve(System &sys)
{
  reset(_lambda_init, 0);
  int i = 0;
  while ((!iterate(sys))&&(i++ < ITERATIONS_MAX));
}

// ******************************************************************
// **********************  SystemProjectionSquare *******************
// ******************************************************************

float SystemProjectionSquare::errorSq() const
{
  int n = numPoints();
  int i;
  float ret = 0.0f;
  for (i = 0; i < n; i++)
    ret += power(errorX(i), 2);
  return(ret);
}

VISVector SystemProjectionSquare::dEdParams(int i) const
{
  cout << "SystemProjectionSquare::dEdParams(int i)" << endl;
  VISMatrix los = _projection.lineOfSight(i);
  cout << "los" << los << endl;
  VISArray< VISVector > intersects = _square.intersections(los);
  cout << "intersects.n " << intersects.n() << endl;
  int n = intersects.n();
  VISVector ret(_square.numParams()), point, normal;
  ret = 0.0f;
  for (int i = 0; i < n; i++)
    {
      point = intersects.peek(i);
      normal = _square.normal(point);
      cout << "point" << point << endl;
      cout << "normal" << normal << endl;
      cout << "dxdParams" << _square.dXdParams(point) << endl;
      ret += ((_square.dXdParams(point)).t())*normal;
    }
  return(ret);
}








