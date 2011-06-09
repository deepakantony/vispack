#include "model.h"
#include "image/imagefile.h"
#include <fstream.h>

ModelFit::ModelFit(Model *model)
{
    _num_scales = 1;
    _force_window = DEFAULT_FORCE_WINDOW;
    _partition_factor = DEFAULT_PARTITION_FACTOR;
    _model = model;
}

ModelFit::ModelFit(Model *model, int num_scales)
{
    _num_scales = num_scales;
    _force_window = DEFAULT_FORCE_WINDOW;
    _partition_factor = DEFAULT_PARTITION_FACTOR;
    _model = model;
}

float ModelFit::calculateError(const VISMatrix &normal_weights, int scale)
{
  int num_points = _model->numPoints();
  int num_scans = _scans.n();
  //    cout << "num scans is " << num_scans << endl << flush;
  //    cout << "num points is " << num_points << endl << flush;
  VISMatrix 
    points = _model->point(), 
    normals = _model->normal();
  VISVector this_point;
  //	this_normal;
  float total_error = 0.0f;

  // partition(this_normal.dot(
  //   los(this_point, j)))

  int i, j;
  float nw_tmp, conf_total = 0.0f;
  float conf_tmp;
  float nw_total = 0;

  //    cout << "points " << points;
			
  for (i = 0; i < num_points; i++)
    {
      this_point = points.vec(i);
      //	    this_normal = points.col(i);
      
      //      cout << "normal weights size:" << normal_weights.c() << endl;

      for (j = 0; j < num_scans; j++)
	{
	  if ((nw_tmp = normal_weights(j, i)) != 0.0f)
	    {
//  	      cout << "this point " << i << " " 
//  		   << this_point ;
//  	      cout << this_point << flush;
//  	      // don't you have to transform for the scan here
//  	      cout << " err " << g_integral(this_point, j, scale, conf_tmp)
// 	       		   << endl << flush;
//  	      cout << "nw_tmp " << i << " " 
//  		   << nw_tmp;
	      total_error += nw_tmp
		*g_integral(this_point, j, 
			    scale, conf_tmp);
	      conf_total += conf_tmp*nw_tmp;
	      //conf_total += conf_tmp;

	    }
	  //		    else
	  //			{
	  //			    cout << "nw_tmp " << nw_tmp << endl << flush;
	  //			}
	}
    }

  //    return(total_error/num_points);
    
  //  cout << "conf_total is " << conf_total << endl;
    //    cout << "total error" << total_error << endl;

  if (conf_total == 0.0f)
       return(FLT_MAX);
     else
       return(total_error/conf_total);
  //return(total_error - normal_weights.sum());
  //return(total_error - conf_total);
  //  cout << "normal weight" << normal_weights.sum() << endl;
  //cout << "total conf" << conf_total << endl;
    //  if ((nw_tmp = normal_weights.sum()) == 0.0f)
  //    return(FLT_MAX);
  //  else
  //  return(total_error/nw_tmp);
  //  return(total_error);
}

VISMatrix ModelFit::normalWeights() const
{
    int num_points = _model->numPoints();
    int num_scans = _scans.n();
    VISMatrix normal_weights(num_scans, num_points);
    VISMatrix 
	points = _model->point(), 
      normals = _model->normal();
    VISVector
	this_point, 
	this_normal;

    int i, j;

//    cout << "normal weights" << endl;

    for (i = 0; i < num_points; i++)
	{
	    this_point = points.vec(i);
	    this_normal = normals.vec(i);

	    for (j = 0; j < num_scans; j++)
		{
		  //		  cout << "in NW" << endl;
		  //		  cout << this_normal << endl;
		  //		  cout << _scans.peek(j)
		  //		    ->lineOfSight(this_point) << endl;

		    normal_weights.poke(j, i) = 
			partition((float)this_normal.dot(
			    (_transforms_inv.peek(j)
			     *(_scans.peek(j))
			     ->lineOfSight(this_point))));
// 		    if (normal_weights.peek(j, i) == 0.0f)
// 			{
// 			    cout << "point" << endl << this_point;
// 			    cout << "normal" << endl << this_normal;
// 			    cout << "los" << endl 
// 				 << (_scans.peek(j))->lineOfSight(this_point);
// 			    cout << normal_weights.peek(j, i) << endl;
// 			}
		}
	}
    return(normal_weights);
}

VISMatrix ModelFit::dotWeights() const
{
  int num_points = _model->numPoints();
  int num_scans = _scans.n();
  VISMatrix normal_weights(num_scans, num_points);
  VISMatrix 
    points = _model->point(), 
    normals = _model->normal();
  VISVector
    this_point, 
    this_normal;

  int i, j;

  //    cout << "normal weights" << endl;

  for (i = 0; i < num_points; i++)
    {
      this_point = points.vec(i);
      this_normal = normals.vec(i);

      for (j = 0; j < num_scans; j++)
	{
	  normal_weights.poke(j, i) = 
	    MAX((float)this_normal
		.dot((_transforms_inv.peek(j)
		      *(_scans.peek(j))
		      ->lineOfSight(this_point))), 0.0f);
	  //		    if (normal_weights.peek(j, i) == 0.0f)
	  //			{
	  //			    cout << "point" << endl << this_point;
	  //			    cout << "normal" << endl << this_normal;
	  //			    cout << "los" << endl 
	  //				 << (_scans.peek(j))->lineOfSight(this_point);
	  //		    cout << normal_weights.peek(j, i) << endl;
	  //			}
	}
    }
  return(normal_weights);
}

// returns FALSE in the case of no problem
int ModelFit::fit()
{
  boolean r = FALSE;
  int i = 0, j;
  _lambda = LAMBDA_INIT;
  _model->setResolution(0, _model->betaRaw());
  for (j = (_num_scales - 1); (j >= 0)&&(!r); j--)
  //  j = (_num_scales - 1);
  {
    _lambda = LAMBDA_INIT;
    while (!iterate(j)&&(i++ < MAX_ITERATIONS))
      {
	//		    cout << "done iteration " << i << endl;
      }

    //    cout << "done scale " << j << endl;

    if (i >= MAX_ITERATIONS)
      r = TRUE;
  }
  // this is just the output we need to keep track of convergence
  //    cout << i << "\t";
  if (r)
    return(-1);
  else
    return(i);
}

//
// this is if the model also has associated levels of resolution
//
int ModelFit::fitScales()
{
  boolean r = FALSE;
  int i = 0, j;
  _lambda = LAMBDA_INIT;
  for (j = (_num_scales - 1); (j >= 0)&&(!r); j--)
    //    j = (_num_scales - 1);
    {
      //      if (j == (_num_scales - 1))
      //	
      //      else
      //	_model->setResolution(j);
      _model->setResolution(j, _model->betaRaw());
      _lambda = LAMBDA_INIT;
      while (!iterate(j)&&(i++ < MAX_ITERATIONS))
	{
	  //		    cout << "done iteration " << i << endl;
	  //	  cout << i << " " << transformToEuler(quaternionToTransform(beta())) << endl;
	  //	  cout << i << " " << beta() << endl;
	}

      //    cout << "done scale " << j << endl;

      if (i >= MAX_ITERATIONS)
	r = TRUE;
	    
      //	    cout << "i is " << i << endl;
    }

  if (r)
    return(-1);
  else 
    return(i);
	    
}


boolean ModelFit::iterate(int scale)
{

  _model->updatePoints();
  _model->updateNormals();
  int num_points = _model->numPoints();
  int num_params = _model->numParams();
  int num_scans = _scans.n();
  VISMatrix 
    points = _model->point(),
    normals = _model->normal(),
    ddE_ddB(num_params, num_params), 
    ddE_ddB_tmp;
  VISVector 	
    dE_dB(num_params), 
    d_N, d_g,
    this_point, 
    this_normal;

  //  cout << "points " << points.t() << endl;

  float this_g;

  VISMatrix normal_weights = normalWeights();
  ///    VISMatrix normal_weights(num_scans, num_points);

  int i, j;
  float nw_tmp;


  VISMatrix dS_N_dB = _model->dS_dot_dB(normals);
  VISMatrix grad_g_pts(_model->rangeDim() + 1, num_points); 
  VISVector this_grad_g(_model->rangeDim() + 1);
  VISVector g_pts(num_points);
  // cout << "dS_N_dB.t()" << dS_N_dB.t() << endl;
  //   VISMatrix dS_N_dB(num_params, num_points);
  //   dS_N_dB = 0.0f;
  int num_misses = 0;

  //    cout << points << flush;
  
  

  for (i = 0; i < num_points; i++)
    {
      this_point = points.vec(i);
      this_normal = normals.vec(i);

      this_g = 0.0f;
      this_grad_g = 0.0f;

//        for (j = 3; j >= 0; j--)
//  	{
//  	  cout  << "scale " << j << endl;
//  	  cout << "point " <<  this_point << endl;
//  	  cout  << "depth " << (_scans.peek(0))->depth(this_point, j) << endl;
//  	  cout  << "confidence " << (_scans.peek(0))->confidence(this_point, j) << endl;
//  	  cout  << "distance " << ((_scans.peek(0))->distance(this_point) - (_scans.peek(0))->depth(this_point, j)) << endl;
//  	}

      for (j = 0; j < num_scans; j++)
	{
	  if ((nw_tmp = normal_weights(j, i)) != 0.0f)
	    {
	      this_g +=
		nw_tmp*g(this_point, j, scale);
	      this_grad_g += 
		nw_tmp*(_transforms_inv.peek(j)
			*grad_g(this_point, j, scale));
	      //	  cout << " got one " << endl;
	      //	      cout << "nw" << nw_tmp << endl;
	      //	      cout << "this_g " << this_g << endl;
	    }
	  else
	    {
	      //               cout << "error:  this point fail" << this_point << endl;
// 	      cout << "nw " << nw_tmp << endl;
// 	      cout << "los " <<
// 		(_scans.peek(0))->lineOfSight(this_point)
// 		   << endl;
// 	      cout << "this normal " << this_normal << endl;
// 	      cout << flush;
	      //	      num_misses++;
	    }
	}
      //      cout << "this_grad_g "  << this_grad_g;
      grad_g_pts.pokeROI(0, i, this_grad_g);
      g_pts.poke(i) = this_g;
    }
    

  //    cout << "dS_N_dB" << dS_N_dB.t() << endl;
  //  cout << "grad_g_pts" << grad_g_pts.t() << endl;
  //  cout << "g_pts" << g_pts.t() << endl;

  //  cout << "num misses " << num_misses << endl << flush;


  VISMatrix dS_g_dB = _model->dS_dot_dB(grad_g_pts);
  dE_dB = 0.0f;
  ddE_ddB = 0.0f;

  for (i = 0; i < num_points; i++)
    {
      //      cout << "this ds " <<  _model->dS_N_dB(i) << endl << flush;
      d_N = dS_N_dB.vec(i);
      d_g = dS_g_dB.vec(i);
      dE_dB += g_pts.peek(i)*d_N;

      //      cout << "this d_N " <<  d_N << endl << flush;      
      //      cout << "this_g " << g_pts.peek(i) << endl;
      //      cout << "this d_g" << d_g << endl;
      //      cout << "this_d "  << d_model;

      // make sure the order is right!!!!!!!!!!!!!	    

      //ddE_ddB += power(g_pts.peek(i), 2)*exterior(d_N, d_N);
      ddE_ddB += exterior(d_N, d_g);

      // the first variation has the normals pointed inward --- but we store them outward
      //            ddE_ddB += exterior(-1.0f*d_N, d_g);
	    //ddE_ddB += exterior(d_N, d_g);
      //ddE_ddB += exterior(d_g, d_N);
      //      cout << "d_N" << d_N << endl;
      //      cout << "d_g" << d_g << endl;
      //      cout << "ext " << exterior(d_N, d_g);
    }


  boolean 
    done = false, 
    done_iterations = false;
    
  int which_direction = 0;
  VISVector new_beta, this_beta = _model->betaRaw(), delta_beta;
  float tmp_error, new_error, last_error;

  //  cout << "beta raw " << beta;

  VISMatrix I = VISIdentity(num_params);
  float lambda = _lambda;

  //  VISMatrix dot_weights = dotWeights();
  //last_error = new_error = calculateError(dot_weights, scale);
  last_error = new_error = calculateError(dotWeights(), scale);

  //  cout <<"actual beta " <<  _model->beta() << endl;
  //    cout << " about to enter iterations " << endl << flush;

  // debugging
  //  for (i = 0; i < 10; i++)
  //  {
  //      beta = _model->beta();
  //      beta.poke(0) += -0.1;
  //      _model->beta(beta);
  //      _model->updatePoints();
  //      cout << "error at x = " <<   beta.peek(0) << ":  "  
  //	   << calculateError(dotWeights(), scale) << endl;
  //    }
  //
  //  beta.poke(0) += 1.0f;

  while (!done)
    //    while (1)
    {
      ddE_ddB_tmp = ddE_ddB + lambda*I;
      //	    cout << "ddE_ddB " << ddE_ddB << endl;
      //	    cout << "ddE_ddB_tmp " << ddE_ddB_tmp << endl;
      //
      delta_beta = ddE_ddB_tmp.inverseSVD()*dE_dB;
      //	    delta_beta = (1.0f/lambda)*dE_dB;
      //
      // sign depends on how you make the normals face
      //      cout << " delta beta " << delta_beta << endl << flush;

      new_beta = this_beta + delta_beta;
      beta(new_beta);
      //      _model->beta(new_beta);
      _model->updatePoints();

      // this will slow things down
      // for DEBUGGING!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      //	    normal_weights = normalWeights();

      tmp_error = new_error;

      //      new_error = calculateError(dot_weights, scale);
      new_error = calculateError(dotWeights(), scale);
      
//              cout <<"lambda " << lambda << endl;
//              cout << " old beta " << this_beta << " new beta " 
//        	   << new_beta << " last error " << last_error 
//        	   << " new error " << new_error << endl << flush;
//                cout << " delta beta " << delta_beta << endl << flush;
//                cout << ddE_ddB << endl;
//                cout << ddE_ddB_tmp << endl;
//                cout << dE_dB << endl;
//  	      cout << new_beta;
      
      if (which_direction == 0)
	{
	  if (new_error > last_error)
	    {
	      lambda = ::VISmax(::VISmin(lambda*LAMBDA_FACTOR, 
					 (float)MAX_LAMBDA), 
				(float)MIN_LAMBDA);
	      which_direction = -1;
	    }
	  else if (new_error < last_error)
	    {
	      //			    lambda = ::VISmax(lambda/LAMBDA_FACTOR, 
	      //					 (float)MIN_LAMBDA);
	      //				    which_direction = 1;
	      lambda = lambda/LAMBDA_FACTOR;
	      done = true;
	    }
	  else if (new_error == last_error)
	    {
	      //				    cout << "got equal error" << endl;
	      //				    lambda = max(lambda/LAMBDA_FACTOR, 
	      //						 (float)MIN_LAMBDA);
	      //				    if (lambda == MIN_LAMBDA)
	      //					done = true;
	      //				    which_direction = 1;
	      done = true;
	      done_iterations = true;
	    }
	}

      else if (which_direction == -1)
	{
	  if (new_error > last_error)
	    lambda = ::VISmin(lambda*LAMBDA_FACTOR, 
			      (float)MAX_LAMBDA);
	  else
	    done = true;
	}
      else if (which_direction == 1)
	{
	  if (new_error < last_error)
	    lambda /= LAMBDA_FACTOR;
	  else
	    {
	      done = true;
	      lambda *= LAMBDA_FACTOR;
	      new_beta = this_beta + lambda*dE_dB;
	      new_error = tmp_error;
	    }
	}
      else
	{
	  cout << "got a problem on which_direction \n";
	}
      //	    cout << "old_error " << last_error << " this_error "
      //		 <<  new_error << " diff " <<
      //		new_error - last_error << endl;
      if (lambda >= (float)MAX_LAMBDA)
	{
	  //	  _model->beta(this_beta);
	  beta(this_beta);
	  done = true;
	  done_iterations = true;
	  //		    cout << "got done iterations" << endl;
	}
    }
  if (lambda < MAX_LAMBDA)
    {
      this_beta = new_beta;
      last_error = new_error;
    }
  _lambda = lambda;

  //    cout << "MATRIX TOTAL SIZE is " << MATRIX_TOTAL_SIZE << endl;
  //    VISImageMemory();
  return(done_iterations);
  ///    return(false);
}


// how do you make the relationship between a point and 
// a map
float ModelFit::g(const VISVector &point, int map_num)
{
  VISVector this_point = _transforms.peek(map_num)*point;
    float r = (_scans.peek(map_num))->depth(this_point);
    float d = (_scans.peek(map_num))->distance(this_point);
    float c = (_scans.peek(map_num))->confidence(this_point);
//    cout << "force window is " << _force_window;
    d = d - r;
    if (fabs(d) > _force_window)
      return(0);
    else
      return(c*d);
}

float ModelFit::g(VISVector point, int map_num, int scale)
{
  VISVector this_point;
  //    cout << point << flush;
  this_point = _transforms.peek(map_num)*point;
  float r = (_scans.peek(map_num))->depth(this_point, scale);
  float d = (_scans.peek(map_num))->distance(this_point);
  float c = (_scans.peek(map_num))->confidence(this_point, scale);
  //   cout << "point " << point;
  //   cout << "this_point " << this_point;
  //   cout << "r " << r << " d " << d << " c " << c << endl << flush;
   //   cout << "g " << c*(r - d)*gaussFast((r - d)/_force_window) 
   //        << endl << flush;
   //   cout << "force window is " << _force_window;
  d = d - r;
  if (fabs(d) > _force_window)
    return(0);
  else
    return(c*d);
}


VISVector ModelFit::grad_g(VISVector point, int map_num, int scale)
{
    VISVector this_point;
//    cout << point << flush;
    this_point = _transforms.peek(map_num)*point;
    float r = (_scans.peek(map_num))->depth(this_point, scale);
    float d = (_scans.peek(map_num))->distance(this_point);
    float c = (_scans.peek(map_num))->confidence(this_point, scale);
//    cout << "point " << point;
//    cout << "this_point " << this_point;
//    cout << "r " << r << " d " << d << " c " << c << endl << flush;
//    cout << "g " << c*(r - d)*gaussFast((r - d)/_force_window) 
//	 << endl << flush;
    //    d = (d - r)/_force_window;
    
//    cout << "force window is " << _force_window;
// the -1 in front counteracts the "-r" and the -1 in front of line of sight

//    cout << "los" << (_scans.peek(map_num))
//          ->lineOfSight(point) << endl;
    //    cout << "g depth" << (_scans.peek(map_num))
    //      ->grad_depth(this_point, scale) 
    //      	 << endl;
    //    cout << "grad_g" << c*(1.0f - d*d)*gaussFast(d)
    //      *((_scans.peek(map_num))->lineOfSight(this_point)
    //	+ (_scans.peek(map_num))->grad_depth(this_point, scale) 
    //	)
    //	 << endl;

    // the factor of -1 for los is to account fo the fact that it points back toward the origin

    d = d - r;

    if (fabs(d) > _force_window)
      return(0.0f*(_scans.peek(map_num))->lineOfSight(this_point));
    else
      return(c*((_scans.peek(map_num))->lineOfSight(this_point)
	     + (_scans.peek(map_num))->grad_depth(this_point, scale)));
    //    return(c*(1.0f - d*d)*gaussFast(d)
    //	   *((_scans.peek(map_num))->lineOfSight(this_point)
    //	     + (_scans.peek(map_num))->grad_depth(this_point, scale) 
    //	     ));
}


float  ModelFit::g_integral(const VISVector &point, int map_num)
{
    VISVector this_point = _transforms.peek(map_num)*point;
    float r = (_scans.peek(map_num))->depth(this_point);
    float d = (_scans.peek(map_num))->distance(this_point);
    float c = (_scans.peek(map_num))->confidence(this_point);
//    cout << "force window is " << _force_window;
    //    return(c*(1.0f - gaussFast((r - d)/_force_window)));

    d = d - r;

    if (fabs(d) > _force_window)
      return(c*_force_window*_force_window);
    else
      return(c*d*d);
}

float  ModelFit::g_integral(const VISVector &point, int map_num, int scale)
{
  VISVector this_point = _transforms.peek(map_num)*point;
  float r = (_scans.peek(map_num))->depth(this_point, scale);
  float d = (_scans.peek(map_num))->distance(this_point);
  float c = (_scans.peek(map_num))->confidence(this_point, scale);
  //      cout << "point " << point;
  //    cout << "this_point " << this_point;
  //      cout << "r " << r << " d " << d << " c " << c << endl << flush;
  //      cout << "g_int " << c*(1.0f - gaussFast((r - d)/_force_window)) << endl << flush;
  //    cout << "force window is " << _force_window;
  //  return(c*(1.0f - gaussFast((r - d)/_force_window)));

  d = d - r;

  if (fabs(d) > _force_window)
    return(c*_force_window*_force_window);
  else
    return(c*d*d);
}


float  ModelFit::g_integral(const VISVector &point,
			    int map_num, int scale, 
			    float &c)
{
    VISVector this_point = _transforms.peek(map_num)*point;
    //    float this_conf;
    float r = (_scans.peek(map_num))->depth(this_point, scale);
    float d = (_scans.peek(map_num))->distance(this_point);
    //    this_conf = (_scans.peek(map_num))->confidence(this_point, scale);
    c = (_scans.peek(map_num))->confidence(this_point, scale);
// must update the total confidence, which gets normalized.
//         cout << "point " << point;
//         cout << "this_point " << this_point;
//         cout << "r " << r << " d " << d << " c " << c << endl << flush;
//         cout << "g_int " << c*(1.0f - gaussFast((r - d)/_force_window))
//     	 << endl << flush;
// does this make more sense --- penalize things that don't hit
//    if (r == (_scans.peek(map_num))->errorCode())
//	{
//	    c = 1.0f;
//	    return(1.0f);
//	}
    //    if (this_conf == 0.0) 
    //      c = 0.0f;
    //    else
    //      c = 1.0f;
    //    if (r == (_scans.peek(map_num))->errorCode())
    //          return(1.0f);
    //    return(c*(1.0f - gaussFast((r - d)/_force_window)));

    d = d - r;

    if (fabs(d) > _force_window)
      return(c*_force_window*_force_window);
    else
      return(c*d*d);

}


Scan::Scan(const VISImage<float> &range_map)
{
  //  	cout << "have entered scan create in" << endl;
	_range_map_orig = range_map;
	_range_map_orig_dx = (range_map.dx()).setBorder(0.0f, 1);
	_range_map_orig_dy = (range_map.dy()).setBorder(0.0f, 1);
	_num_scales = 0;
	_scale_factor = NULL;
//	_dx_kernel = dx_kernel(1);
//	_dy_kernel = dy_kernel(1);
}

Scan::Scan(const VISImage<float> &range_map, float normals_scale)
{
  //	cout << "have entered scan create in" << endl;
  _range_map_orig = range_map;
  VISImage<float> dx_kernel = gauss_dx_kernel(1, normals_scale),
    dy_kernel = gauss_dy_kernel(1, normals_scale);
  //  VISImage<float> dx_kernel = ::dx_kernel(1),
  //    dy_kernel = ::dy_kernel(1);
  
  //  cout << "dx kernel" << endl;
  //    dx_kernel.printData();

  //  cout << "dx kernel" << endl;
  //  (::gauss_dx_kernel(1, 2.0f)).printData();

  //  cout << "dx kernel" << endl;
  //  (::dx_kernel(1)).printData();

  if (range_map.height() > 1)
    _range_map_orig_dx = (range_map.convolve(dx_kernel))
      .setBorder(0.0, dx_kernel.width()/2);
  else
    // this is a hack to accomodate the 2d range maps
    _range_map_orig_dx = range_map.convolve(dx_kernel);
  //	_range_map_orig.printData();
  //	_range_map_orig_dx.printData();
  //	_range_map_orig_dy = range_map.dy(normals_scale);
  //	
  //	_range_map_orig_dy = range_map
  //	    .convolve(gauss_dy_kernel(1, normals_scale));

  // this is a hack to accomodate the 2d range maps
  if (range_map.height() > 1)
    {
      _range_map_orig_dy = (range_map.convolve(dy_kernel))
	.setBorder(0.0, dy_kernel.height()/2);
    }
  else
    {
    _range_map_orig_dy = range_map;
    }

  _num_scales = 0;
  _scale_factor = NULL;
  //	_dx_kernel = dx_kernel(1);
  //	_dy_kernel = dy_kernel(1);
}


void Scan::createScales(int num_scales)
{
    int k;
    _num_scales = num_scales;
    _range_maps.clear();
    _range_maps_dx.clear();
    _range_maps_dy.clear();
    _range_maps.appendItem(_range_map_orig);
    _range_maps_dx.appendItem(_range_map_orig_dx);
    _range_maps_dy.appendItem(_range_map_orig_dy);
    _scale_factor = new float[num_scales];
    _scale_factor[0] = 1.0f;
    char filename[80];
    VISImageFile im_file;
    for (k = 1; k < num_scales; k++)
	{
	    _range_maps.appendItem(subSampleAverage((_range_maps.peek(k - 1))));
	    _range_maps_dx.appendItem(subSampleAverage((_range_maps_dx.peek(k - 1))));
	    _range_maps_dy.appendItem(subSampleAverage((_range_maps_dy.peek(k - 1))));
	    _scale_factor[k] = 2.0f*_scale_factor[k - 1];
	    sprintf(filename, "sub%d.fit", k);
	    im_file.write(_range_maps_dx.peek(k), filename);
	}
}

void SphericalScan::createScales(int num_scales)
{

  Scan::createScales(num_scales);

  char filename[80];
  VISImageFile im_file;
  
  _conf_maps.clear();
  _conf_maps.appendItem(_conf_map);
  //  im_file.write(_conf_maps.peek(0), "sub0.fit");
  for (int k = 1; k < num_scales; k++)
    {
      _conf_maps.appendItem(subSampleMin((_conf_maps.peek(k - 1))));
      //      sprintf(filename, "sub%d.fit", k);
      //     im_file.write(_conf_maps.peek(k), filename);
    }
}


VISImage<float> subSampleAverage(VISImage<float> im)
{
    int w = im.width(), h = im.height();
    int new_w, new_h;

    VISImage<float> 
	r((new_w = w/2) + w%2, 
	  (new_h = h/2) + h%2),
	tmp(new_w + w%2, h);
    
    int i, j, this_x, this_y;
    
    for (j = 0; j < h; j++)
	{
	    tmp.poke(0, j) = 0.5f*(im.peek(0, j) 
				   + im.peek(1, j));
	    
	    for (i = 1; i < new_w; i++)
		{
		    this_x = 2*i;
		    tmp.poke(i, j) 
			= 0.25f*(im.peek(this_x - 1, j) + 
				 im.peek(this_x + 1, j)) + 
			0.5f*im.peek(this_x, j);
		}

	    if (w%2)
		tmp.poke(new_w, j) 
		    = 0.5f*(im.peek(w - 1, j) + im.peek(w - 2, j));
	}

    new_w = new_w + w%2;

    if (h < 2)
	{
	    return(tmp);
	}
    
    for (i = 0; i < new_w; i++)
	{
	    r.poke(i, 0) = 0.5f*(tmp.peek(i, 0) 
				   + tmp.peek(i, 1));
	    
	    for (j = 1; j < new_h; j++)
		{
		    this_y = 2*j;
		    r.poke(i, j) 
			= 0.25f*(tmp.peek(i, this_y - 1) + 
				 tmp.peek(i, this_y + 1)) + 
			0.5f*tmp.peek(i, this_y);
		}

	    if (h%2)
		r.poke(i, new_h) 
		    = 0.5f*(tmp.peek(i, h - 1) + tmp.peek(i, h - 2));
	}
    
    return(r);
}

VISImage<float> subSampleMin(VISImage<float> im)
{
  int w = im.width(), h = im.height();
  int new_w, new_h;

  VISImage<float> 
    r((new_w = w/2) + w%2, 
      (new_h = h/2) + h%2),
    tmp(new_w + w%2, h);
    
  int i, j, this_x, this_y;
    
  for (j = 0; j < h; j++)
    {
      tmp.poke(0, j) = VISmin(im.peek(0, j), im.peek(1, j));
	    
      for (i = 1; i < new_w; i++)
	{
	  this_x = 2*i;
	  tmp.poke(i, j) 
	    = VISmin(VISmin(im.peek(this_x - 1, j), im.peek(this_x + 1, j)), 
		     im.peek(this_x, j));
	}

      if (w%2)
	tmp.poke(new_w, j) 
	  = VISmin(im.peek(w - 1, j), im.peek(w - 2, j));
    }

  new_w = new_w + w%2;

  if (h < 2)
    {
      return(tmp);
    }
    
  for (i = 0; i < new_w; i++)
    {
      r.poke(i, 0) = VISmin(tmp.peek(i, 0), tmp.peek(i, 1));
	    
      for (j = 1; j < new_h; j++)
	{
	  this_y = 2*j;
	  r.poke(i, j) 
	    = VISmin(VISmin(tmp.peek(i, this_y - 1), tmp.peek(i, this_y + 1)),  tmp.peek(i, this_y));
	}

      if (h%2)
	r.poke(i, new_h) 
	  = VISmin(tmp.peek(i, h - 1), tmp.peek(i, h - 2));
    }
    
  return(r);
}


// *****************************************************
// *****************************************************
// ****************** TwoDOrthoScane *******************
// *****************************************************
// *****************************************************


//*********************************************************************
//*********************************************************************
//********************** 3D Spherical Scan ****************************
//*********************************************************************
//*********************************************************************


SphericalScan::SphericalScan(const VISImage<float> &range_image):Scan(range_image)
{
//    cout << "about to enter sphere create in" << endl;
  setC0(range_image.width()/2.0f);
  setR0(range_image.height()/2.0f);
  setDeltaTheta(DEFAULT_FOV/range_image.width()); 
  setDeltaPhi(DEFAULT_FOV/range_image.height());
  setRangeOffset(0.0f);
  setRangeDelta(1.0f);
  VISImage<float> conf(range_image.width(), range_image.height());
  conf = 1.0f;
  setConfidence(conf);
}

SphericalScan::SphericalScan(const VISImage<float> &range_image, float normals_scale)
  :Scan(range_image, normals_scale)
{
//    cout << "about to enter sphere create in" << endl;
  setC0(range_image.width()/2.0f);
  setR0(range_image.height()/2.0f);
  setDeltaTheta(DEFAULT_FOV/range_image.width()); 
  setDeltaPhi(DEFAULT_FOV/range_image.height());
  setRangeOffset(0.0f);
  setRangeDelta(1.0f);
  VISImage<float> conf(range_image.width(), range_image.height());
  conf = 1.0f;
  setConfidence(conf);
}



VISVector TwoDOrthoScan::getSurfaceNormal(int i, int j) const
{
    VISVector r(4);
    float di = _range_map_orig_dx(i, j);
    float dj = _range_map_orig_dy(i, j);
    float range = _range_map_orig.peek(i, j);

    float norm = 1.0f/sqrt(1.0f + pow(di,2) * pow(dj,2));

    r.poke(0) = di;
    r.poke(1) = dj;
    r.poke(2) = -1.0f;
    r.poke(3) = 0.0f;

    return(r);
}

VISVector TwoDOrthoScan::getSurfaceNormal(int i, int j, int scale) const
{

    VISVector r(4);
    //    VISImage<float> range_map = _range_maps.peek(scale);
    float di = (_range_maps_dx.peek(scale)).peek(i, j);
    float dj = (_range_maps_dy.peek(scale)).peek(i, j);
    float this_scale = _scale_factor[scale];
    float range = (_range_maps.peek(scale)).interp(i/this_scale, j/this_scale);

    float norm = 1.0f/sqrt(1.0f + pow(di,2) * pow(dj,2));

    r.poke(0) = di;
    r.poke(1) = dj;
    r.poke(2) = -1.0f;
    r.poke(3) = 0.0f;

    return(r);
}



float TwoDOrthoScan::depth(const VISVector &p, int scale) const
{
  float x = p.peek(0)/_scale_factor[scale];
  VISImage<float> this_map =  _range_maps.peek(scale);
  if (!this_map.checkBounds(x, 0.0f))
	{
//	    cout << "check bounds bad " << p << "scale " 
//		 << scale << " x " << x << endl;
	    return(errorCode());
	}
    else
	return(this_map.interp(x, 0.0f));
}

VISVector TwoDOrthoScan::grad_depth(const VISVector &p, int scale) const
{
  float x = p.peek(0)/_scale_factor[scale];
  VISImage<float> this_map_dx =  _range_maps_dx.peek(scale);
  //    cout << "grad depth at " << x << " is " << 
  //	this_map_dx.interp(x, 0.0f) << endl;
  if (!this_map_dx.checkBounds(x, 0.0f))
    {
      //	    cout << "check bounds bad " << p << "scale " 
      //		 << scale << " x " << x << endl;
      return VISVector(0.0f, 0.0f, 0.0f);
    }
  else
    return VISVector(this_map_dx.interp(x, 0.0f), 0.0f, 0.0f);
}

float TwoDOrthoScan::depth(const VISVector &p) const
{
    float x = p.peek(0);
    if (!_range_map_orig.checkBounds(x, 0.0f))
	return(errorCode());
    else
	return(_range_map_orig.interp(x, 0.0f));
}

VISVector TwoDOrthoScan::lineOfSight(const VISVector &p) const
{
  return(VISVector(0.0f, -1.0f, 0.0f));
}

//float TwoDOrthoScan::distance(const VISVector &p)
//{
//    return(p.peek(0));
//}



void TwoCircles::initialize(unsigned num_pts, float radius, 
			    float separation)
{

    _radius = radius; 
    _length = separation;

    int i;

    _range_dim = 2;
    _num_params = 4;
    _num_points = num_pts;
    _beta = VISVector(_num_params);


// use homogeneous coordinates
    _points_obj = VISMatrix(_range_dim + 1, _num_points);
    _normals_obj = VISMatrix(_range_dim + 1, _num_points);

    _points = VISMatrix(_range_dim + 1, _num_points);
    _normals = VISMatrix(_range_dim + 1, _num_points);

    for (i = 0; i < num_pts/2; i++)
	{
	    _points_obj.poke(0, i) 
		= _radius*cos((float)i/((float)num_pts/2)
			      *2.0f*M_PI) - _length;
	    _normals_obj.poke(0, i) 
		= cos((float)i/((float)num_pts/2)
		      *2.0f*M_PI);

	    _points_obj.poke(1, i) 
		= _radius*sin((float)i/((float)num_pts/2)
			      *2.0f*M_PI);
	    _normals_obj.poke(1, i) = sin((float)i/((float)num_pts/2)
					  *2.0f*M_PI);

	    _points_obj.poke(2, i) = 1.0f;
	    _normals_obj.poke(2, i) = 0.0f;
	}
    for (; i < num_pts; i++)
	{
	    _points_obj.poke(0, i) = _radius*cos((float)i/((float)num_pts/2)
						 *2.0f*M_PI) +
		_length;
	    _normals_obj.poke(0, i) = cos((float)i/((float)num_pts/2)
					  *2.0f*M_PI);

	    _points_obj.poke(1, i) = _radius*sin((float)i/((float)num_pts/2)
						 *2.0f*M_PI);
	    _normals_obj.poke(1, i) = sin((float)i/((float)num_pts/2)
					  *2.0f*M_PI);

	    _points_obj.poke(2, i) = 1.0f;
	    _normals_obj.poke(2, i) = 0.0f;
	}
//    cout << _points_obj << flush;
//    cout << _normals_obj << flush;
}


// this should be done as a matrix for efficiency
VISMatrix TwoCircles::dS_dB(int i) const
{
    VISVector p = _points_obj.peekROI(0, 1, i, i);
    
    VISMatrix R(2, 2), dR(2, 2);
    VISMatrix r(_range_dim + 1, _num_params);

    float theta = _beta.peek(3), s = _beta.peek(0);

    R.poke(0, 0) = cos(theta);
    R.poke(1, 0) = sin(theta);
    R.poke(1, 1) = cos(theta);
    R.poke(0, 1) = -sin(theta);
	    
    dR.poke(0, 0) = -sin(theta);
    dR.poke(1, 0) = cos(theta);
    dR.poke(1, 1) = -sin(theta);
    dR.poke(0, 1) = -cos(theta);

    r.pokeROI(0, 0, R*p);
    r.pokeROI(0, 1, VISVector(1.0f, 0.0f));
    r.pokeROI(0, 2, VISVector(0.0f, 1.0f));
    r.pokeROI(0, 3, s*dR*p);

    for (int j = 0; j < _num_params; j++)
	r.poke(2, j) = 0.0f;

//    cout << "dS_dB " << r << endl;

    return(r);
}


// this should be faster than calling points one at a time
VISMatrix TwoCircles::dS_dot_dB(const VISMatrix &vects) const
{
    VISMatrix R(2, 2), dR(2, 2);
    VISMatrix tmp(_range_dim + 1, _num_params);
    VISMatrix r(_num_params, _num_points);
    VISVector p;
    int j;

    float theta = _beta.peek(3), s = _beta.peek(0);

    R.poke(0, 0) = cos(theta);
    R.poke(1, 0) = sin(theta);
    R.poke(1, 1) = cos(theta);
    R.poke(0, 1) = -sin(theta);
	    
    dR.poke(0, 0) = -sin(theta);
    dR.poke(1, 0) = cos(theta);
    dR.poke(1, 1) = -sin(theta);
    dR.poke(0, 1) = -cos(theta);


// stuff that does not depend on point
    for (j = 0; j < _num_params; j++)
      tmp.poke(_range_dim, j) = 0.0f;

    tmp.pokeROI(0, 1, VISVector(1.0f, 0.0f));
    tmp.pokeROI(0, 2, VISVector(0.0f, 1.0f));


    for (j = 0; j < _num_points; j++)
	{
	    p = _points_obj.peekROI(0, 1, j, j);
	    tmp.pokeROI(0, 0, R*p);
	    tmp.pokeROI(0, 3, s*dR*p);
//	    cout << "tmp.t()" << tmp.t() << endl;
//	    cout << "grad_g_pts.vec()" << grad_g_pts.vec(j) << endl;
	    r.pokeROI(0, j, (tmp.t()*vects.vec(j)));
	}

//    cout << "dS_dB " << r << endl;

    return(r);
}


// this should be done as a matrix for efficiency
VISVector Model::dS_N_dB(int i) const
{
//    cout << "ds_db " << dS_dB(i) << " normal " <<
//			 _normals.vec(i) << endl;
    return(((dS_dB(i)).t())*_normals.vec(i));
}


//
// not for now.....  work in progress
//
//VISMatrix TwoCircles::dS_N_dB()
//{
//    VISVector p;
//    
//    VISMatrix R, dR;
//    VISMatrix r(4, num_points);
//    VISMatrix point_deriv(4, 2);
//
//    float theta = _beta.peek(3), s = _beta.peek(0);
//
//    R.poke(0, 0) = cos(theta);
//    R.poke(0, 1) = sin(theta);
//    R.poke(1, 1) = cos(theta);
//    R.poke(1, 0) = -sin(theta);
//	    
//    dR.poke(0, 0) = -sin(theta);
//    dR.poke(0, 1) = cos(theta);
//    dR.poke(1, 1) = -sin(theta);
//    dR.poke(1, 0) = -cos(theta);
//
//    
//    for (int i = 0; i < _num_points; i++)
//	{
//	    p = _points_obj.col(i);
//
//	    point_deriv = point_deriv.putroi(0, 0, R*point);
//	    point_deriv = point_deriv.putroi(0, 1, VISVector(1.0f, 0.0f));
//	    point_deriv = point_deriv.putroi(0, 2, VISVector(0.0f, 1.0f));
//	    pointer_deriv = point_deriv.putroi(0, 3, s*dR*p);
//
//	    r.putroi(0, i, (point_deriv*_normals.peek(i)).trans())
//	}
//    return(r);
//}
//


void TwoCircles::updatePoints()
{
    
    VISMatrix T(3, 3);
    float theta = _beta.peek(3);
    float scale = _beta.peek(0);
    
    T.poke(0, 0) = scale*cos(theta);
    T.poke(1, 0) = scale*sin(theta);
    T.poke(1, 1) = scale*cos(theta);
    T.poke(0, 1) = -scale*sin(theta);
    
    T.poke(2, 0) = 0.0f;
    T.poke(2, 1) = 0.0f;
    T.poke(2, 2) = 1.0f;

    T.pokeROI(0, 2, _beta.peekROI(X0, Y0));

    _points = T*_points_obj;

//    for (int i = 0; i < _num_points; i++)
//	{
//	    cout  << R*_points_obj.vec(i) + trans << endl ;
//	    _points.pokeROI(0, i, R*_points_obj.vec(i) + trans);
//	}

//    cout << "_points " << _points ;
//    cout << "_points_obj " << _points_obj ;

}

void TwoCircles::updateNormals()
{

    VISMatrix T = VISIdentity(3);
    float theta = _beta.peek(3);
    
    T.poke(0, 0) = cos(theta);
    T.poke(1, 0) = sin(theta);
    T.poke(1, 1) = cos(theta);
    T.poke(0, 1) = -sin(theta);

    _normals = T*_normals_obj;
}

//void TwoCircles::beta(VISVector beta)
//{
// what to do here!!!!! ?????
//}

// ***********************************************************************************
// ***********************************************************************************
// ****************************** Two Spheres *************************************
// ***********************************************************************************
// ***********************************************************************************

void TwoSpheres::initialize(int num_pts, float separation, float radius)
{

  // take a random selection of points on the surface
  float theta, phi;
  int i;

  //  cout << "radius" << radius << endl;
  
  VISMatrix points(4, num_pts);
  VISMatrix normals(4, num_pts);

  for (i = 0; i < num_pts/2; i++)
    {
      theta = rand1()*2.0*M_PI;
      phi = rand1()*M_PI;
      points.pokeROI(0, i, VISVector(radius*sin(theta) + separation/2.0f, radius*sin(phi)*cos(theta), 
		     radius*cos(phi)*cos(theta), 1.0f));
      normals.pokeROI(0, i, VISVector(sin(theta), sin(phi)*cos(theta), cos(phi)*cos(theta), 0.0f));
    }

  for (i; i < num_pts; i++)
    {
      theta = rand1()*2.0*M_PI;
      phi = rand1()*M_PI;
      points.pokeROI(0, i, VISVector(radius*sin(theta) - separation/2.0f, radius*sin(phi)*cos(theta), 
		     radius*cos(phi)*cos(theta), 1.0f));
      normals.pokeROI(0, i, VISVector(sin(theta), sin(phi)*cos(theta), cos(phi)*cos(theta), 0.0f));
    }
  
  // cout << points.t() << endl;
  // cout << normals.t() << endl;

  _centroid = _Rcentroid = VISVector(0.0f, 0.0f, 0.0f, 0.0f);

  Rigid3DModel::initialize(points, normals);
}


float ModelFit::partition(float v) const
{
    if (v <= 0.0)
	return(0.0);
    else if (v > 1.0f/PART_FACTOR)
	return(1.0);
    else
	return(0.5f - 0.5f*cosFast(v*M_PI*PART_FACTOR));
}


int ModelFit::addScan(Scan* scan, const VISMatrix &trans)
{
    _scans.appendItem(scan);
    _transforms.appendItem(trans);
    _transforms_inv.appendItem(trans.inverseSVD());
//    if (_num_scales > 1)
//    cout << "about to create scales " << _num_scales << endl << flush;
    scan->createScales(_num_scales);
//    cout << "created scales " << scan->numScales() << endl << flush;
    return(_scans.n());
}


VISMatrix exterior(const VISVector &vec1, const VISVector &vec2)
{
    int i, j, n;
    if ((n = vec1.n()) != vec2.n())
	{
	    WARN("exterior: size mismatch\n");
	    return VISMatrix();
	}

    VISMatrix r(n, n);
    
    for (i = 0; i < n; i++)
	for (j = 0; j < n; j++)
	    r.poke(i, j) = vec1.peek(i)*vec2.peek(j);
	    
    return(r);
}


VISVector TwoDOrthoScan::get3DPoint(int i, int j) const
{
    return(VISVector((float)i, 
		     _range_map_orig.peek(i, 0), 
		     1.0f));
}

VISVector TwoDOrthoScan::get3DPoint(int i, int j, int scale) const
{
    return(VISVector((float)i, (_range_maps.peek(scale))
		     .interp((float)i/_scale_factor[scale], 0), 
		     1.0f));
}

VISVector ModelFit::get3DPoint(int i, int j, int scan) const
{
    return((_transforms_inv.peek(scan))
	   *(_scans.peek(scan))->get3DPoint(i, j));
}

VISVector ModelFit::get3DPoint(int i, int j, int scan, int scale) const
{
    return((_transforms_inv.peek(scan))
	   *(_scans.peek(scan))->get3DPoint(i, j, scale));
}


// *****************************************************
// *****************************************************
// ******************* SuperQ **************************
// *****************************************************
// *****************************************************


#define E_THRESH (0.05f)
#define E_MIN ((float)(1.0e-3))
#define E_CONST (1.0e-4)
#define E_TOP (1.0f)
#define E_MAX ((float)(E_TOP - 1.0e-3))
//#define E_MAX ((float)(2.0))

float SuperQ::e_edge(float e) const
{
    if (e > 0.0f)
	if (e >= E_THRESH)
	    if (e <= (E_TOP - E_THRESH))
		{
		    return(0.0f);
		}
	    else if (e < E_MAX)
		return(log((E_TOP - e)/E_THRESH)/E_CONST);
	    else
		return(-HUGE);
	else
	    return(-log(e/E_THRESH)/E_CONST);
    else
	return(HUGE);
}


VISMatrix SuperQ::dS_dB(int i) const
{
    VISVector p = _points_obj.peekROI(0, 1, i, i);
    
    VISMatrix R(2, 2), dR(2, 2);
    VISMatrix r(_range_dim + 1,_num_params);
    VISMatrix total;
    VISVector 
	ex(1.0f, 0.0f), 
	ey(0.0f, 1.0f), 
	vec_tmp(3);

    float 
	theta = _beta.peek(THETA), 
	a = _beta.peek(A), 
	b = _beta.peek(B), 
	e = _beta.peek(E);

    R.poke(0, 0) = cos(theta);
    R.poke(1, 0) = sin(theta);
    R.poke(1, 1) = cos(theta);
    R.poke(0, 1) = -sin(theta);
	    
    dR.poke(0, 0) = -sin(theta);
    dR.poke(1, 0) = cos(theta);
    dR.poke(1, 1) = -sin(theta);
    dR.poke(0, 1) = -cos(theta);

    float 
	sine = _unit_circle.peek(0, i), 
	cosine = _unit_circle.peek(1, i);
    float log_tmp;
    
    r.pokeROI(0, X0, ex);
    r.pokeROI(0, Y0, ey);
    r.pokeROI(0, THETA, dR*p);
    r.pokeROI(0, A, p.peek(0)*R*ex);
    r.pokeROI(0, B, p.peek(1)*R*ey);

    for (int j = 0; j < _num_params; j++)
	r.poke(2, j) = 0.0f;

//    cout << " p " << p << endl;

    if ((sine == 0.0f)||((log_tmp = log(sine)) == -HUGE))
	log_tmp = 0.0f;
    
    vec_tmp.poke(0) = _points_obj.peek(0, i)*log_tmp;

    if ((cosine == 0.0f)||
	((log_tmp = log(cosine)) == -HUGE))
	log_tmp = 0.0f;

    vec_tmp.poke(1) = 
	_points_obj.peek(1, i)*log_tmp;

    vec_tmp.poke(2) = 1.0f; 

    vec_tmp += e_edge(e)*_normals.vec(i);

    r.pokeROI(0, E, vec_tmp);

//    cout << "dS_dB " << r << endl;

    

    return(r);
}


VISMatrix SuperQ::dS_dot_dB(const VISMatrix & vects) const
{
    
    VISMatrix R(2, 2), dR(2, 2);
    VISMatrix r(_range_dim + 1,_num_params);
    VISMatrix total(_num_params, _num_points);
    VISVector 
	ex(1.0f, 0.0f), 
	ey(0.0f, 1.0f), 
	vec_tmp(3);
    VISVector p;

    float 
	theta = _beta.peek(THETA), 
	a = _beta.peek(A), 
	b = _beta.peek(B), 
	e = _beta.peek(E);
    int j;

    float edge_effect = e_edge(e);

    R.poke(0, 0) = cos(theta);
    R.poke(1, 0) = sin(theta);
    R.poke(1, 1) = cos(theta);
    R.poke(0, 1) = -sin(theta);
	    
    dR.poke(0, 0) = -sin(theta);
    dR.poke(1, 0) = cos(theta);
    dR.poke(1, 1) = -sin(theta);
    dR.poke(0, 1) = -cos(theta);

    float sine, cosine;
    float log_tmp;
    
    r.pokeROI(0, X0, ex);
    r.pokeROI(0, Y0, ey);

    for (j = 0; j < _num_params; j++)
	r.poke(2, j) = 0.0f;


    for (j = 0; j < _num_points; j++)
	{
	    p = _points_obj.peekROI(0, 1, j, j);
	    sine = _unit_circle.peek(0, j);
	    cosine = _unit_circle.peek(1, j);

	    r.pokeROI(0, THETA, dR*p);
	    r.pokeROI(0, A, p.peek(0)*R*ex);
	    r.pokeROI(0, B, p.peek(1)*R*ey);

//    cout << " p " << p << endl;

	    if ((sine == 0.0f)||((log_tmp = log(sine)) == -HUGE))
		log_tmp = 0.0f;
    
	    vec_tmp.poke(0) = p.peek(0)*log_tmp;
	    
	    if ((cosine == 0.0f)||
		((log_tmp = log(cosine)) == -HUGE))
		log_tmp = 0.0f;

	    vec_tmp.poke(1) = 
		p.peek(1)*log_tmp;

	    vec_tmp.poke(2) = 1.0f; 
	    
	    vec_tmp += edge_effect*_normals.vec(j);

	    r.pokeROI(0, E, vec_tmp);

//	    cout << "dS_dB " << r << endl;

	    total.pokeROI(0, j, (r.t()*vects.vec(j)));

	}
    
    return(total);
}


void SuperQ::initialize(unsigned num_pts)
{
    
    _num_params = 6;
    _range_dim = 2;
    _num_points = num_pts;
    _beta = VISVector(_num_params);
    
    (_points_obj = VISMatrix(_range_dim + 1, _num_points)) = 1.0f;
    (_normals_obj = VISMatrix(_range_dim + 1, _num_points)) = 0.0f;
    
    (_points = VISMatrix(_range_dim + 1, _num_points)) = 1.0f;
    (_normals = VISMatrix(_range_dim + 1, _num_points)) = 0.0f;
    
    _unit_circle = VISMatrix(_range_dim, _num_points);
    _signs = VISMatrix(_range_dim, _num_points);

    int i;

// this los is assumed to go from the scanner outward, and therefore 
// the surface normals should face the interior of the object.    
    for (i = 0; i < _num_points; i++)
	{
	    _unit_circle.poke(0, i) = fabs(sin((float)i/((float)num_pts)
					  *2.0f*M_PI));
	    _unit_circle.poke(1, i) = fabs(cos((float)i/((float)num_pts)
					   *2.0f*M_PI));

	    if (_unit_circle.peek(0, i) != 0.0f)
		_signs.poke(0, i) = sin((float)i/((float)num_pts)
					*2.0f*M_PI)/_unit_circle.peek(0, i);
	    else 
		_signs.poke(0, i) = 1.0f;

	    if (_unit_circle.peek(1, i) != 0.0f)
		_signs.poke(1, i) = cos((float)i/((float)num_pts)
					*2.0f*M_PI)/_unit_circle.peek(1, i);
	    else 
		_signs.poke(0, i) = 1.0f;
	}
}

void SuperQ::updatePoints()
{
    VISMatrix T(3, 3);
    //    cout << "beta " << _beta; 
    float theta = _beta.peek(THETA);
    float a = _beta.peek(A);
    float b = _beta.peek(B);
    float e = _beta.peek(E);
    
    T.poke(0, 0) = cos(theta);
    T.poke(1, 0) = sin(theta);
    T.poke(1, 1) = cos(theta);
    T.poke(0, 1) = -sin(theta);
    
    T.poke(2, 0) = 0.0f;
    T.poke(2, 1) = 0.0f;
    T.poke(2, 2) = 1.0f;

    T.pokeROI(0, 2, _beta.peekROI(X0, Y0));

    int i;

    for (i = 0; i < _num_points; i++)
	{
//	    cout  << R*_points_obj.vec(i) + trans << endl ;
	    _points_obj.poke(0, i) = _signs.peek(0, i)*
		a*pow(_unit_circle.peek(0, i), e);
	    _points_obj.poke(1, i) = _signs.peek(1, i)*
		b*pow(_unit_circle.peek(1, i), e);
	}


//    for (i = 0; i < _num_points; i++)
//	{
//	    cout  << R*_points_obj.vec(i) + trans << endl ;
//	    _points = _points.pokeROI(0, i, R*_points_obj.vec(i) + trans);
//	}

    _points = T*_points_obj;

//    cout << "_points " << _points ;
//    cout << "_points_obj " << _points_obj ;

}

void SuperQ::updateNormals()
{

    VISMatrix T = VISIdentity(3);
    float theta = _beta.peek(THETA);
    float a = _beta.peek(A);
    float b = _beta.peek(B);
    float e = _beta.peek(E);

    int i;
    
    T.poke(0, 0) = cos(theta);
    T.poke(1, 0) = sin(theta);
    T.poke(1, 1) = cos(theta);
    T.poke(0, 1) = -sin(theta);

    float sine, cosine;
    VISVector normal_raw(2);

    for (i = 0; i < _num_points; i++)    
	{
	    sine = _unit_circle.peek(0, i);
	    cosine = _unit_circle.peek(1, i);

	    if (sine == 1.0f)
		{
		    normal_raw.poke(0) = _signs.peek(0, i);
		    normal_raw.poke(1) = 0.0f;
		}
	    else if (cosine == 1.0f)
		{
		    normal_raw.poke(1) = _signs.peek(1, i);
		    normal_raw.poke(0) = 0.0f;
		}
	    else
		{
		    normal_raw.poke(0) = 
			_signs.peek(0, i)*
			b*pow(cosine, e - 1.0f)
			*sine;
		    
		    normal_raw.poke(1) = 
			_signs.peek(1, i)*
			a*pow(sine, e - 1.0f)
			*cosine;
		}
	    _normals_obj.pokeROI(0, i, normal_raw/normal_raw.norm());
	}

    _normals = T*_normals_obj;
    //    cout << "normals " << _normals;
    //    cout << "points " << _points;
}


// *****************************************************
// *****************************************************
// ****************** OrthoRangeModel*******************
// *****************************************************
// *****************************************************

void OrthoRangeModel::initialize(const VISImage<float> &range_map, 
float normals_scale)
{
     _range_map = range_map;
     _range_map_dx = range_map.derivative(1, 0, normals_scale);
    _num_params = 3;
    _range_dim = 2;
    _num_points = _range_map.width();
    _beta = VISVector(_num_params);
    
    _points_obj = VISMatrix(_range_dim + 1, _num_points);
    _normals_obj = VISMatrix(_range_dim + 1, _num_points);
    
    _points = VISMatrix(_range_dim + 1, _num_points);
    _normals = VISMatrix(_range_dim + 1, _num_points);
    
    int i;

// this appears to be because of a bug in dx_kernel
//    VISImage<float> dx_kernel = gauss_dx_kernel(1, normals_scale);
//    VISImage<float> range_dx = _range_map.convolve(dx_kernel);
//    VISImage<float> range_dx = -1.0f*_range_map.dx();
    float this_dx, this_mag;

    for (i = 0; i < _num_points; i++)
	{
	    _points_obj.poke(0, i) = i;
	    _points_obj.poke(1, i) = _range_map.peek(i, 0);
	    _points_obj.poke(2, i) = 1.0f;

	    this_dx = _range_map_dx.peek(i, 0);
	    this_mag = sqrt(1.0f + this_dx*this_dx);
	    
	    _normals_obj.poke(0, i) = this_dx/this_mag;
	    _normals_obj.poke(1, i) = -1.0f/this_mag;
	    _normals_obj.poke(2, i) = 0.0f;
	    
	}
}


void OrthoRangeModel::updatePoints()
{
    
    VISMatrix T(3, 3);
    float theta = _beta.peek(THETA);
    
    T.poke(0, 0) = cos(theta);
    T.poke(1, 0) = sin(theta);
    T.poke(1, 1) = cos(theta);
    T.poke(0, 1) = -sin(theta);
    
    T.poke(2, 0) = 0.0f;
    T.poke(2, 1) = 0.0f;
    T.poke(2, 2) = 1.0f;

    T.pokeROI(0, 2, _beta.peekROI(X0, Y0));

    _points = T*_points_obj;

}

void OrthoRangeModel::updateNormals()
{

    VISMatrix T = VISIdentity(3);
    float theta = _beta.peek(THETA);
    
    T.poke(0, 0) = cos(theta);
    T.poke(1, 0) = sin(theta);
    T.poke(1, 1) = cos(theta);
    T.poke(0, 1) = -sin(theta);

    _normals = T*_normals_obj;
}

VISMatrix OrthoRangeModel::dS_dB(int i) const
{
    VISVector p = _points_obj.peekROI(0, 1, i, i);
    
    VISMatrix dR(2, 2);
    VISMatrix r(_range_dim + 1, _num_params);

    float theta = _beta.peek(THETA);

    dR.poke(0, 0) = -sin(theta);
    dR.poke(1, 0) = cos(theta);
    dR.poke(1, 1) = -sin(theta);
    dR.poke(0, 1) = -cos(theta);

    r.pokeROI(0, X0, VISVector(1.0f, 0.0f));
    r.pokeROI(0, Y0, VISVector(0.0f, 1.0f));
    r.pokeROI(0, THETA, dR*p);

    for (int j = 0; j < _num_params; j++)
	r.poke(2, j) = 0.0f;

//    cout << "dS_dB " << r << endl;

    return(r);
}


VISMatrix OrthoRangeModel::dS_dot_dB(const VISMatrix &vects) const
{
    VISVector p;    
    VISMatrix dR(2, 2);
    VISMatrix r(_range_dim + 1, _num_params);
    VISMatrix total(_num_params, _num_points);

    float theta = _beta.peek(THETA);
    int j;

    dR.poke(0, 0) = -sin(theta);
    dR.poke(1, 0) = cos(theta);
    dR.poke(1, 1) = -sin(theta);
    dR.poke(0, 1) = -cos(theta);

    r.pokeROI(0, X0, VISVector(1.0f, 0.0f));
    r.pokeROI(0, Y0, VISVector(0.0f, 1.0f));

    for (j = 0; j < _num_params; j++)
	r.poke(2, j) = 0.0f;


    for (j = 0; j < _num_points; j++)
	{
	    p = _points_obj.peekROI(0, 1, j, j);
	    r.pokeROI(0, THETA, dR*p);
	    total.pokeROI(0, j, (r.t()*vects.vec(j)));
	}

//    cout << "dS_dB " << r << endl;

    return(total);
}



//*********************************************************************
//*********************************************************************
//********************** 3D Spherical Scan ****************************
//*********************************************************************
//*********************************************************************


void SphericalScan::setParams(const VISVector &params)
{
  if (params.n() == SphereScan::NUMPARAMS)
    {
      setR0(params.peek(SphereScan::R0)); 
      setC0(params.peek(SphereScan::C0)); 
      setDeltaPhi(params.peek(SphereScan::DPHI)); 
      setDeltaTheta(params.peek(SphereScan::DTHETA)); 
      setRangeDelta(params.peek(SphereScan::RDELTA)); 
      setRangeOffset(params.peek(SphereScan::ROFFSET)); 
    }
}

VISVector SphericalScan::getParams() const
{
  VISVector r(SphereScan::NUMPARAMS);
  r.poke(SphereScan::R0) = _r0; 
  r.poke(SphereScan::C0) = _c0; 
  r.poke(SphereScan::DPHI) = _delta_phi; 
  r.poke(SphereScan::DTHETA) = _delta_theta; 
  r.poke(SphereScan::RDELTA) = _r_delta; 
  r.poke(SphereScan::ROFFSET) = _r_offset; 
  return(r);
}

VISVector SphericalScan::getSurfaceNormal(int i, int j) const
{
    VISVector r(4);
    //    VISImage<float> range_map = _range_map_orig;
    float di = _r_delta*(_range_map_orig_dx.peek(i, j));
    float dj = _r_delta*(_range_map_orig_dy.peek(i, j));
    //    cout << "getSurfaceNormal " << i << " " << j  << " "<< di  << " "<< dj << endl;
    float cos_theta, sin_theta, cos_phi, sin_phi, 
	dR_dPhi, dR_dTheta;

    float range = _range_map_orig.peek(i, j);

    float phi = (j - _r0)*_delta_phi,
      theta = (i - _c0)*_delta_theta;

    //    cout << "_delta " << _delta_phi << " " << _delta_theta  <<  endl;
    //    cout << "getSurfaceNormal " << phi << " " << theta  << " "<< range  <<  endl;

    cos_theta = cos(theta);
    sin_theta = sin(theta);
    cos_phi = cos(phi);
    sin_phi = sin(phi);


    dR_dTheta = di/_delta_theta;
    dR_dPhi = dj/_delta_phi;

    r.poke(0) = cos_theta*(range*sin_theta - dR_dTheta*cos_theta);
    r.poke(1) = cos_theta*sin_phi*(range*cos_theta + dR_dTheta*sin_theta)
	- dR_dPhi*cos_theta;
    r.poke(2) = dR_dPhi*sin_phi
	+ cos_phi*cos_theta*(range*cos_theta + dR_dTheta*sin_theta);
    r.poke(3) = 0.0f;

// should this be positive or negative???
    r /= -1.0f*r.norm();
    //    r /= r.norm();
//    cout << "normal length test " << r.dot(r) << endl;
    return(r);
}

VISVector SphericalScan::getSurfaceNormal(int i, int j, int scale) const
{

    VISVector r(4);
    float this_scale = _scale_factor[scale];

    float di = _r_delta*(_range_maps_dx.peek(scale)).interp(i/this_scale, j/this_scale);
    float dj = _r_delta*(_range_maps_dy.peek(scale)).interp(i/this_scale, j/this_scale);
    float cos_theta, sin_theta, cos_phi, sin_phi, 
	dR_dPhi, dR_dTheta, norm;

    float range = (_range_maps.peek(scale)).interp(i/this_scale, j/this_scale);

    float phi = (j - _r0)*_delta_phi,
      theta = (i - _c0)*_delta_theta;

    cos_theta = cos(theta);
    sin_theta = sin(theta);
    cos_phi = cos(phi);
    sin_phi = sin(phi);


    dR_dTheta = di/_delta_theta;
    dR_dPhi = dj/_delta_phi;

//    norm = 1.0f/sqrt(pow(dR_dTheta*cos_theta, 2) 
//		     + pow(range*cos_theta, 2) 
//		     + pow(dR_dPhi, 2));

    r.poke(0) = cos_theta*(range*sin_theta - dR_dTheta*cos_theta);
    r.poke(1) = cos_theta*sin_phi*(range*cos_theta + dR_dTheta*sin_theta)
      - dR_dPhi*cos_theta;
    r.poke(2) = dR_dPhi*sin_phi
      + cos_phi*cos_theta*(range*cos_theta + dR_dTheta*sin_theta);
    r.poke(3) = 0.0f;

//
// this is what you get if you have (x, y, z) = r*(sin_theta*cos_phi, sin_phi, 
// cos_phi*cos_theta)
//
//    r.poke(0) = cos_phi*sin_theta*(range*cos_phi + dR_dPhi*sin_phi)
//	- dR_dTheta*cos_theta;
//    r.poke(1) = cos_phi*(range*sin_phi - dR_dPhi*cos_phi);
//    r.poke(2) = dR_dTheta*sin_theta
//	+ cos_phi*cos_theta*(range*cos_phi + dR_dPhi*sin_phi);
//    r.poke(3) = 0.0f;

// should this be positive or negative???
    r /= -1.0f*r.norm();
//    r /= r.norm();

//    cout << "normal length test " << r.dot(r) << endl;

    return(r);

}


VISVector SphericalScan::myImageCoord(const VISVector &p) const
{
  float x = p.peek(0), y = p.peek(1), z = p.peek(2);
  float theta, phi;
  float r, c;

  if (z > 0)
    {
      theta = ::atan(x/sqrt(y*y + z*z));
      phi = ::atan(y/z);
      c = (theta/_delta_theta + _c0);
      r = (phi/_delta_phi + _r0);
      //		if ((r >= 0.0)&&(c >= 0.0)&&(c <= (_depth.width() - 1))
      //		    &&(r <= (_depth.height() - 1)))
      return(VISVector(c, r));
    }
  else
    {
      //	    cout << "imageCoord: ERROR " << p << endl;
      return VISVector();
    }
}


VISVector SphericalScan::imageCoord(const VISVector &p, int scale) const
{
    VISVector coord = imageCoord(p);    
    return(coord/_scale_factor[scale]);
}


float SphericalScan::depth(const VISVector &p, int scale) const
{
  VISVector coord = imageCoord(p);
  float c, r;
  VISImage<float> this_map;
    
  if (coord.isValid()&&
      ((this_map = _range_maps.peek(scale))
       .checkBounds(c = (coord.peek(0)/
			 _scale_factor[scale]),
		    r = (coord.peek(1)/
			 _scale_factor[scale]))))
    {
//             cout << "point " << p << endl;
//             (this_map).print();
//             cout << "coord " << coord << endl;
//             cout << "c " << c << " r " << r << endl;
//             cout << "depth " << this_map.VISImage<float>::interp(c, r) << endl;
      return(_r_offset + _r_delta*this_map.VISImage<float>::interp(c, r));
    }
  else 		
    return(errorCode());
}

float SphericalScan::depth(const VISVector &p) const
{
  // geometry of the range map is defined here.....
  VISVector coord = imageCoord(p);
  float c, r;
    
  if (coord.isValid()&&(_range_map_orig
			.checkBounds(c = coord.peek(0),r = coord.peek(1))))
    return(_r_offset + _r_delta*_range_map_orig.VISImage<float>::interp(c, r));
  else 		
    return(errorCode());
}


VISVector SphericalScan::grad_depth(const VISVector &p, int scale) const
{

  VISVector coord = myImageCoord(p);
  float c, r;
  VISImage<float> this_map, im_dx, im_dy;
  float dc, dr;
  float x = p.peek(0), y = p.peek(1), z = p.peek(2);

  if (!(coord.isValid()&&
	((im_dx  = _range_maps_dx.peek(scale))
	 .checkBounds(c = (coord.peek(0)/_scale_factor[scale]),
		      r = (coord.peek(1)/_scale_factor[scale])))))
    return VISVector(0.0f, 0.0f, 0.0f, 0.0f);
  im_dy = _range_maps_dy.peek(scale);

  dc = _r_delta*im_dx.interp(c, r);
  dr = _r_delta*im_dy.interp(c, r);

  //  cout << "dc " << dc << " and dr " << dr << endl;

  float dtheta_dx, dtheta_dy, dtheta_dz, dphi_dx, dphi_dy, dphi_dz;
  float mag1 = y*y + z*z;
  float mag2 = mag1 + x*x;
  float mag1_sqrt = sqrt(mag1);

  dtheta_dx = mag1/(mag2*mag1_sqrt);
  dtheta_dy = -x*y/(mag2*(mag1_sqrt));
  dtheta_dz = -x*z/(mag2*(mag1_sqrt));
  dphi_dx = 0.0f;
  dphi_dy = z/(mag1);
  dphi_dz = -y/(mag1);

  return VISVector
    (
     dtheta_dx*dc/_delta_theta, 
     dtheta_dy*dc/_delta_theta + dphi_dy*dr/_delta_phi,
     dtheta_dz*dc/_delta_theta + dphi_dz*dr/_delta_phi, 
     0.0f
     );
}


float SphericalScan::confidence(const VISVector &p) const
{
  VISVector coord = imageCoord(p);
  float c, r;

  //  cout << "image coord" << coord << endl;
  //  cout << "point " << p << endl;
    
  if (coord.isValid()&&
	(_conf_map.checkBounds(c = coord.peek(0),r = coord.peek(1))))
	return(_conf_map.VISImage<float>::interp(c, r));
    else 		
	{
	  //	    cout << "got conf oob at " << c << " " << r << endl 
	  //		 << flush;
	    return(0.0f);
	}
}

float SphericalScan::confidence(const VISVector &p, int scale) const
{
  //  **************
  //    should we subsample confidence ????
  //  **************
  return(confidence(p));
  
//    VISVector coord = myImageCoord(p);
//    float c, r;
//    VISImage<float> this_map;



//    if (coord.isValid()&&((this_map = _conf_maps.peek(scale))
//   			.checkBounds(c = (coord.peek(0)/_scale_factor[scale]),
//   				     r = (coord.peek(1)/_scale_factor[scale]))))
//      {
//        if (this_map.VISImage<float>::interp(c, r) < 0.5)
//  	cout << this_map.VISImage<float>::interp(c, r)  << 
//  	  "got one less than 0.5 " << p << " " << r << " " << c << endl;
    
//        else
//  	cout << "got one greater ";
//        return(this_map.VISImage<float>::interp(c, r));
//      }
//    else 		
//      return(0.0f);

}

VISVector SphericalScan::lineOfSight(const VISVector &p) const
{
    float 
	mag,
	z = p[2],
	y = p[1],
	x = p[0];

    mag = sqrt(x*x + y*y + z*z);

    //    cout << "los " << p << endl;

    if (depth(p) == errorCode())
      return(VISVector(0.0f, 0.0f, 0.0f, 0.0f));

    if (mag > 0.0)
      mag = -1.0f/mag;
    return(VISVector(x*mag, y*mag, z*mag, 0.0f));
}

VISVector SphericalScan::get3DPoint(int i, int j) const
{
  float 
    phi = (j - _r0)*_delta_phi, 
    theta = (i - _c0)*_delta_theta,
    r = _r_offset + _r_delta*_range_map_orig.peek(i, j);
  return(VISVector(r*sin(theta), r*sin(phi)*cos(theta), 
		   r*cos(phi)*cos(theta), 1.0f)); 
}

VISVector SphericalScan::get3DPoint(int i, int j, int scale) const
// should this assume i, j are scaled or not?
// for now, yes
{
    float this_scale = _scale_factor[scale];
    float phi = (this_scale*j - _r0)*_delta_phi;
    float theta = (this_scale*i - _c0)*_delta_theta;

//    cout << "this scale factor " << this_scale << endl;
    float r = _r_offset + _r_delta*(_range_maps.peek(scale)).peek(i, j);

    return(VISVector(r*sin(theta), r*sin(phi)*cos(theta), 
		     r*cos(phi)*cos(theta), 1.0f)); 
}


VISArray< VISImage<float> > Scan
::get3DPoints(const VISMatrix &transform) const
{
  return get3DPoints(transform, 0);
}

VISArray< VISImage<float> > Scan::get3DPoints(const VISMatrix &transform, int scale) const
{
  VISImage<float> range = _range_maps.peek(scale);
  int w = range.width(), h = range.height();
  VISImage<float>
    x(w - 2*BORDER, h - 2*BORDER), 
    y(w - 2*BORDER, h - 2*BORDER), 
    z(w - 2*BORDER, h - 2*BORDER);
  

  VISVector p(4);
  p.poke(3) = 1.0f;
  int i, j;

  for (i = 0; i < (w - 2*BORDER); i++)
    for (j = 0; j < (h - 2*BORDER); j++)
      {
	//	    cout << "about to do get 3D point, scale = " << scale << endl;
	p.pokeROI(0, get3DPoint(i + BORDER, j + BORDER, scale));
	p = transform*p;
	x.poke(i, j) = p.peek(0);
	y.poke(i, j) = p.peek(1);
	z.poke(i, j) = p.peek(2);
      }

  VISArray< VISImage<float> > r;

  r.poke(0) = x;
  r.poke(1) = y;
  r.poke(2) = z;

  return(r);
}


//*********************************************************************
//*********************************************************************
//************************ 3D Ortho Scan ****************************
//*********************************************************************
//*********************************************************************

VISVector OrthoScan::getSurfaceNormal(int i, int j) const
{
    VISVector r(4);
    r.poke(0) = _range_map_orig_dx.peek(i, j);
    r.poke(1) = _range_map_orig_dy.peek(i, j);
    r.poke(2) = -1.0f;
    r.poke(3) = 0.0f;
    r /= r.norm();
    return(r);
}

VISVector OrthoScan::getSurfaceNormal(int i, int j, int scale) const
{

    VISVector r(4);
    float this_scale = _scale_factor[scale];
    r.poke(0) = (_range_maps_dx.peek(scale)).interp(i/this_scale, j/this_scale);
    r.poke(1) = (_range_maps_dy.peek(scale)).interp(i/this_scale, j/this_scale);
    r.poke(2) = -1.0f;
    r.poke(3) = 0.0f;
    r /= r.norm();

    return(r);

}

VISVector OrthoScan::myImageCoord(const VISVector &p) const
{
  float x = p.peek(0), y = p.peek(1), z = p.peek(2);
  float theta, phi;
  float r, c;

  if (z > 0)
    {
      c = x/_delta_x + _c0;
      r = y/_delta_y + _r0;
      return(VISVector(c, r));
    }
  else
    {
      return VISVector();
    }
}


VISVector OrthoScan::grad_depth(const VISVector &p, int scale) const
{
  VISVector coord = myImageCoord(p);
  float c, r;
  VISImage<float> this_map, im_dx, im_dy;
  float dc, dr;
  float x = p.peek(0), y = p.peek(1), z = p.peek(2);

  if (!(coord.isValid()&&
	((im_dx  = _range_maps_dx.peek(scale))
	 .checkBounds(c = (coord.peek(0)/_scale_factor[scale]),
		      r = (coord.peek(1)/_scale_factor[scale])))))
    {
      //  cout << "got invalid grad_depth" << p << coord <<endl;
      return VISVector(0.0f, 0.0f, 0.0f, 0.0f);
    }
  //  else
  //    cout << "got valid grad_depth" << p << coord <<endl;
  im_dy = _range_maps_dy.peek(scale);

  dc = im_dx.interp(c, r);
  dr = im_dy.interp(c, r);

  return VISVector
    (
     dc/_delta_x, 
     dr/_delta_y, 
     0.0f, 
     0.0f
     );
}


VISVector OrthoScan::get3DPoint(int i, int j) const
{
    float 
      y = (j - _r0)*_delta_y, 
      x = (i - _c0)*_delta_x, 
      z = _range_map_orig.peek(i, j);
    return(VISVector(x, y, z, 1.0f)); 
}

VISVector OrthoScan::get3DPoint(int i, int j, int scale) const
// should this assume i, j are scaled or not?
// for now, yes
{
    float this_scale = _scale_factor[scale];
    float 
      y = (this_scale*j - _r0)*_delta_y, 
      x = (this_scale*i - _c0)*_delta_x, 
      z = (_range_maps.peek(scale)).peek(i, j);
    return(VISVector(x, y, z, 1.0f)); 
}


//*********************************************************************
//*********************************************************************
//********************** Rigid3DModel *********************************
//*********************************************************************
//*********************************************************************

void Rigid3DModel::initialize(const VISMatrix &points, const VISMatrix &normals, 
			      const VISVector &new_beta)
{
  //  beta(new_beta); 
  int num_points = points.c();
  //
  //_points_obj = points;
  //
  // 
  // Build a local coordinate system around centroid????  Not that effective, it seems
  //
  // cout << "points " << points  << endl;
  _centroid = points.mean();
  _centroid.poke(3) = 0.0f;
  //  _centroid = 0.0f;
  //    cout << "cent before" << _centroid << endl;
    //    cout << "new_beta " << new_beta << endl;	
  _points_obj = VISMatrix(points.r(), points.c());
  //_points_obj = points;
	
  for (int i = 0; i < num_points; i++)
    {
      _points_obj.pokeROI(0, i, points.vec(i) - _centroid);
    }

  //  cout << new_beta << endl;
  //  cout << quaternionToTransform(new_beta) << endl;

  //  cout << quaternionToTransform(new_beta.peekROI(0, 6));

  _Rcentroid = quaternionToTransform(new_beta)*_centroid;
  
  //  cout << "cent after" << _Rcentroid << endl;
	
  //	cout << "points obj" << _points_obj << endl;
  //	cout << "points" << points << endl;
  //
  _normals_obj = normals;
  _num_params = 7;
  _range_dim = 3;
  _num_points = _points_obj.c();
  //	_beta = VISVector(_num_params);
  // _beta = beta;
}



// must normalize the quaternions here
int Rigid3DModel::beta(const VISVector &new_beta)
{
  if (new_beta.n() < _num_params)
    {
      cout << "WARNING:  Rigid3DModel::beta -- bad vector size" << endl;
      return(-1);
    }
  
  VISVector q = new_beta.peekROI(Quaternion::W, Quaternion::Z);
//   VISVector q(4);
  
//   q.poke(0) = new_beta.peek(Quaternion::W);
//   q.poke(1) = new_beta.peek(Quaternion::X);
//   q.poke(2) = new_beta.peek(Quaternion::Y);
//   q.poke(3) = new_beta.peek(Quaternion::Z);
  
  _beta = new_beta;
  
  q = q/(q.norm());
  _beta.pokeROI(Quaternion::W, q);

  //  cout << new_beta << endl;
  //  cout << _beta << endl;

  //  _beta = new_beta;
  
  return(0);

}



void Rigid3DModel::updatePoints()
{

  VISMatrix T;
  
  T = quaternionToTransform(_beta);
  T.pokeROI(0, 3, (_beta.peekROI(Quaternion::TX, Quaternion::TZ) + _Rcentroid.peekROI(0, 2)));
  //  T.pokeROI(0, 3, _beta.peekROI(Quaternion::TX, Quaternion::TZ));

  _points = T*_points_obj;
  
  //
  //  cout << "Transform " << T << endl;
  //      cout << "Obj points " << endl << _points_obj.t() << endl;
  //      cout << "updated points " << endl << _points.t() << endl;
}


void Rigid3DModel::updateNormals()
{
  VISMatrix T;
  //
  // this assumes that the quaternion is stored in the beginning of _beta
  // Be careful!!!!!!!!!!!!!!!!!!!!!
  //
    _normals =   (quaternionToTransform(_beta))*_normals_obj;
    //    cout << "normals " << _normals;
}


VISMatrix Rigid3DModel::dS_dot_dB(const VISMatrix &vects) const
{
  int j;
  VISVector p;    
  VISMatrix r(_range_dim + 1, _num_params);
  VISMatrix total(_num_params, _num_points);

        r.pokeROI(0, Quaternion::TX, VISVector(1.0f, 0.0f, 0.0f));
       r.pokeROI(0, Quaternion::TY, VISVector(0.0f, 1.0f, 0.0f));
        r.pokeROI(0, Quaternion::TZ, VISVector(0.0f, 0.0f, 1.0f));

	//        r.pokeROI(0, Quaternion::TX, VISVector(0.0f, 0.0f, 0.0f));
	//          r.pokeROI(0, Quaternion::TY, VISVector(0.0f, 0.0f, 0.0f));
	//          r.pokeROI(0, Quaternion::TZ, VISVector(0.0f, 0.0f, 0.0f));

  for (j = 0; j < _num_params; j++)
    r.poke(3, j) = 0.0f;

  //  cout << "num points" << _num_points << endl;

  for (j = 0; j < _num_points; j++)
    {
            r.pokeROI(0, Quaternion::W, -1.0f*dx_dq(_beta, _points_obj.vec(j)));
      //r.pokeROI(0, Quaternion::W, 0.0f*dx_dq(_beta, _points_obj.vec(j)));
      //       r.pokeROI(0, Quaternion::W, dx_dq(_beta, _points_obj.vec(j)));
      //      cout << "r is " << r << endl;
      //      cout << r.t()*vects.vec(j) <<  endl;
      total.pokeROI(0, j, (r.t()*vects.vec(j)));
    }

  //    cout << "dS_dB " << r << endl;
  return(total);
}



//*********************************************************************
//*********************************************************************
//********************** ThreeDRangeModel *********************************
//*********************************************************************
//*********************************************************************




VISArray< VISImage<float> > ThreeDRangeModel::get3DPoints() const
{
  return(_range_scan->get3DPoints(quaternionToTransform(beta())));
}


VISArray< VISImage<float> > ThreeDRangeModel::get3DPoints(int scale) const
{
    return(_range_scan->get3DPoints(quaternionToTransform(beta()), scale));
}


void ThreeDRangeModel::getObjPointsNormals(int num_pts,
					   int res_level, 
					   VISMatrix &points, 
					   VISMatrix &normals
					   )
{

  int i;

  //    cout << "number of points is " << num_pts << endl;
    
  VISImage<float> depth = _range_scan->Scan::depth(res_level);
  //    VISImage<float> dx_kernel = gauss_dx_kernel(1, normals_scale);
  //    VISImage<float> dy_kernel = gauss_dy_kernel(1, normals_scale);
  //    VISImage<float> range_dx = depth.dx();
  //    VISImage<float> range_dy = depth.dy();

  //    float this_dx, this_dy, this_mag;
  //    float scale_factor = range_scan.scaleFactor(res_level);

  int w = depth.width(), h = depth.height();
  int r, c;

  if (num_pts > 0)
    {
      points = VISMatrix(4, num_pts);
      normals = VISMatrix(4, num_pts);
      for (i = 0; i < num_pts; i++)
      {
	// get random sampling of c's and r's
	c = (int)(rand1()*(w - 2*BORDER - 1) + BORDER);
	r = (int)(rand1()*(h - 2*BORDER - 1) + BORDER);

	// cout << "c and r " << c << " " << r << endl;

	//	points.poke(3, i) = 1.0f;
	points.pokeROI(0, i, _range_scan->get3DPoint(c, r, res_level));
	//	cout << (_range_scan->get3DPoint(c, r, res_level)).t() << endl;
	//	cout << (_range_scan->getSurfaceNormal(c, r, res_level)).t() << endl;
	normals.pokeROI(0, i, _range_scan->getSurfaceNormal(c, r, res_level));
      }
    }
  else
    {
      points = VISMatrix(4, (w - 2*BORDER)*(h - 2*BORDER));
      normals = VISMatrix(4, (w - 2*BORDER)*(h - 2*BORDER));
      i = 0;
      for (c = BORDER; c < (w - BORDER); c++)
	for (r = BORDER; r < (h - BORDER); r++)		
	  {
	    //	    points.poke(3, i) = 1.0f;
	    points.pokeROI(0, i, _range_scan->get3DPoint
			   (c, r, res_level));
	    normals.pokeROI(0, i, _range_scan->getSurfaceNormal
			    (c, r, res_level));
	    i++;
	  }
    }

  // cout << "init points " << points.t() << endl;
  //  cout << "init normals " << normals.t() << endl;

}

void ThreeDRangeModel::initialize(Scan *range_scan, int num_pts, const VISVector &new_beta)
{
  initialize(range_scan, num_pts, 0, new_beta);
}

void ThreeDRangeModel::initialize(Scan *range_scan, 
				  int num_pts)
{
  initialize(range_scan, num_pts, 0);
}

void ThreeDRangeModel::initialize(Scan *range_scan, 
				  int num_pts,
				  int res_level)
{
  initialize(range_scan, num_pts, res_level, VISVector());
}

void ThreeDRangeModel::initialize(Scan *range_scan, 
				  int num_pts,
				  int res_level, 
				  const VISVector &new_beta)
{
  //    cout << "init num points " << num_pts << endl;
  _range_scan = range_scan;
  //    cout << "_range_scan num scales " << _range_scan->numScales();
  //    cout << "_range_scan num scales " << _range_scan->numScales();

  VISMatrix points, normals;
  
  getObjPointsNormals(num_pts, res_level, points, normals);
  //  cout << "points " << points << endl;
  //   cout << "normals " << normals << endl;
  //  cout << "beta in init is  " << new_beta << endl;

  if (new_beta.isValid())
    Rigid3DModel::initialize(points, normals, new_beta);
  else
    Rigid3DModel::initialize(points, normals);

  //  cout << "init points " << points.t() << endl;
  //  cout << "init normals " << normals.t() << endl;
  //    cout << "num pts" << _points.c() << endl;
  //    cout << "num pts obj" << _points_obj.c() << endl;
  //    cout << "done initialize on 3D range model" << endl;

}




// ****************************************************
// ****************************************************
// ******************* RegisterScan *******************
// ****************************************************
// ****************************************************

//
// registers one scan to a set of others
// use parameter files to get the appropriate parameters
//

RegisterScan::RegisterScan(VPF::ParameterFile F)
{
  int num_scans;
  float c0, r0, model_c0, model_r0, c0_default, r0_default;
  float 
    delta_phi, model_delta_phi, delta_phi_default, 
    delta_theta, model_delta_theta, delta_theta_default;

  Scan *s_ptr;
  char scan_keyword[80];
  char filename[80];
  VISMatrix transform(4,4);
  int int_tmp;
  int float_tmp;
  VISImage<float> range_map, conf_im;

  VISImageFile im_file;
  int i, j;
  int num_scales;
    
  int model_num_points;
  float model_normals_scale;
  float force_window;

  float x0, y0, z0, kappa, phi, omega;
  float range_offset, range_scale;

  VISVector params;

  // scanner geometry
#define SCAN_SPHERICAL (0)
#define SCAN_ORTHO (1)

  int scan_type;
  ThreeDRangeModel *model;

  try 
    {
      VPF::set(scan_type, 
	       F["SCAN_TYPE"][0]);
    }
  catch (VPF::Exception e)
    {
      scan_type = SCAN_SPHERICAL;
    }

  try
    {
      VPF::set(range_offset, 
	       F["RANGE_OFFSET"][0]);
    }
  catch (VPF::Exception e)
    {
      cout << "did not get range_offset" << endl;
      range_offset = 0.0f;
    }
    
  try 
    {
      VPF::set(range_scale, F["RANGE_SCALE"][0]);
    }
  catch (VPF::Exception e)
    {
      range_scale = 1.0f;
    }


  try
    {
      VPF::set(r0_default, 
	       F["R0"][0]);
    }
  catch (VPF::Exception e)
    {
      r0_default = -FLT_MAX;
    }

    
  try 
    {
      VPF::set(c0_default, 
	       F["C0"][0]);
    }
  catch (VPF::Exception e)
    {
      c0_default = -FLT_MAX;
    }

  try 
    {
      VPF::set(delta_phi_default, 
	       F["DELTA_PHI"][0]);
    }
  catch (VPF::Exception e)
    {
      delta_phi_default = -FLT_MAX;
    }

  try 
    {
      VPF::set(delta_theta_default, 
	       F["DELTA_THETA"][0]);
    }
  catch (VPF::Exception e)
    {
      delta_theta_default = -FLT_MAX;
    }

  try 
    {
      VPF::set(num_scales, 
	       F["NUM_SCALES"][0]);
    }
  catch (VPF::Exception e)
    {
      num_scales = 1;
    }

  try 
    {
      VPF::set(model_num_points, 
	       F["MODEL_NUM_PTS"][0]);
    }
  catch (VPF::Exception e)
    {
      model_num_points = -1;
    }

  try 
    {
      VPF::set(model_normals_scale, 
	       F["MODEL_NORMALS_SCALE"][0]);
    }
  catch (VPF::Exception e)
    {
      model_normals_scale = 0.0f;
    }

  try 
    {
      VPF::set(force_window, 
	       F["FORCE_WINDOW"][0]);
    }
  catch (VPF::Exception e)
    {
      force_window = DEFAULT_FORCE_WINDOW;
    }
  forceWindow(force_window);	


  // model data and initialization
  try 
    {
      VPF::set(filename, 
	       F["MODEL_FILE"][0]);
	
      range_map = VISImage<float>
	(im_file.read(filename));
      // put in the scale and the offset???
      //	    range_map = range_scale*range_map + range_offset;
      try 
	{
	  VPF::set(x0, 
		   F["X0"][0]);
	}
      catch (VPF::Exception e)
	{
	  x0 = 0.0f;
	}

      try 
	{
	  VPF::set(y0, 
		   F["Y0"][0]);
	}
      catch (VPF::Exception e)
	{
	  y0 = 0.0f;
	}

      try 
	{
	  VPF::set(z0, 
		   F["Z0"][0]);
	}
      catch (VPF::Exception e)
	{
	  z0 = 0.0f;
	}

      try 
	{
	  VPF::set(kappa, 
		   F["KAPPA"][0]);
	}
      catch (VPF::Exception e)
	{
	  kappa = 0.0f;
	}

      try 
	{
	  VPF::set(phi, 
		   F["PHI"][0]);
	}
      catch (VPF::Exception e)
	{
	  phi = 0.0f;
	}

      try 
	{
	  VPF::set(omega, 
		   F["OMEGA"][0]);
	}
      catch (VPF::Exception e)
	{
	  omega = 0.0f;
	}
    }
  catch (VPF::Exception e)
    {
      cout << "WARNING: Did not find model range file --- will ignore." << endl;
    }

  try 
    {
      VPF::set(r0, 
	       F["MODEL_R0"][0]);
    }
  catch (VPF::Exception e)
    {
      if (r0_default != -FLT_MAX)
	r0 = r0_default;
      else
	r0 = range_map.height()/2.0f - 0.5f;
    }

  try 
    {
      VPF::set(c0, 
	       F["MODEL_C0"][0]);
    }
  catch (VPF::Exception e)
    {
      if (c0_default != -FLT_MAX)
	c0 = c0_default;
      else
	c0 = range_map.width()/2.0f - 0.5f;
    }
  try 
    {
      VPF::set(delta_phi, 
	       F["MODEL_DELTA_PHI"][0]);
    }
  catch (VPF::Exception e)
    {
      if (delta_phi_default != -FLT_MAX)
	delta_phi = delta_phi_default;
      else
	delta_phi = DEFAULT_FOV/range_map.height();
    }

  try
    {
      VPF::set(delta_theta, 
	       F["MODEL_DELTA_THETA"][0]);
    }
  catch (VPF::Exception e)
    {
      if (delta_theta_default != -FLT_MAX)
	delta_theta = delta_theta_default;
      else
	delta_theta = DEFAULT_FOV/range_map.width();
    }

  // must accomodate different types of scans

  switch(scan_type)
    {
    case SCAN_SPHERICAL:
      // this smooths the model when it's too noisy
      //      s_ptr = new SphericalScan(range_map.gauss(5.0), model_normals_scale);
      //s_ptr = new SphericalScan(range_map.gauss(model_normals_scale), model_normals_scale);
      for (int kk = 0;  kk < (int)(2.0f*pow(model_normals_scale, 2)); kk++)
	range_map += 0.25f*range_map.anisoDiffuse(2.0);
	//      for (int kk = 0;  kk < (int)(10.0f*model_normals_scale); kk++)
	//	range_map += 0.25f*range_map.anisoDiffuse(2.0);
      s_ptr = new SphericalScan(range_map, model_normals_scale);
      break;
    case SCAN_ORTHO:
      //s_ptr = new OrthoScan(range_map.gauss(model_normals_scale), model_normals_scale);
      //      s_ptr = new OrthoScan(range_map, model_normals_scale);
      break;
    }

  //     s_ptr->setC0(r0);
  //     s_ptr->setR0(c0);
  //     s_ptr->setDeltaTheta(delta_phi);
  //     s_ptr->setDeltaPhi(delta_theta);
  //     s_ptr->setRangeOffset(range_offset);
  //     s_ptr->setRangeDelta(range_scale);
  params = VISVector(SphereScan::NUMPARAMS);
  params.poke(SphereScan::R0) = r0;
  params.poke(SphereScan::C0) = c0;
  params.poke(SphereScan::DPHI) = delta_phi;
  params.poke(SphereScan::DTHETA) = delta_theta;
  params.poke(SphereScan::ROFFSET) = range_offset;
  params.poke(SphereScan::RDELTA) = range_scale;
  s_ptr->setParams(params);
    
  //   

  //    s_ptr->setParams(VISVector(r0, c0, delta_phi, delta_theta));
  //    s_ptr->setParams(VISVector(r0, c0, delta_phi, delta_theta, range_offset, range_scale));

  if (num_scales > 1)
    model = new ThreeDRangeModel(s_ptr, model_num_points, 
				 num_scales);
  else
    // do you need ot set the number of scales here
    model = new ThreeDRangeModel(s_ptr, model_num_points);

  //    scan = SphericalScan(range_map.gauss(model_normals_scale));

  VISVector euler(6);
  euler.poke(Euler::TX) = x0;
  euler.poke(Euler::TY) = y0;
  euler.poke(Euler::TZ) = z0;
  euler.poke(Euler::KAPPA) = kappa;
  euler.poke(Euler::PHI) = phi;
  euler.poke(Euler::OMEGA) = omega;

    
  model->beta(transformToQuaternion(eulerToTransform(euler)));
  //    cout <<     "trans  " <<  eulerToTransform(euler) << endl;   
  //    cout <<     "quat  " <<  transformToQuaternion(eulerToTransform(euler)) << endl;   
  //    cout <<     "quat test " <<  quaternionToTransform(transformToQuaternion(eulerToTransform(euler))) << endl;   
  //    cout <<     "beta " << model->Model::beta() << endl;
  setModel(model, num_scales);


  // the scans
  try {
    VPF::set(num_scans, 
	     F["NUM_SCANS"][0]);
  }
  catch (VPF::Exception e)
    {
      cerr << "Argument 1 of \"NUM_SCANS\" is not valid" << endl;
    }
	    
  for (i = 0; i < num_scans; i++)
    {

      sprintf(scan_keyword, "RANGE_FILE_%d", i);
      try
	{
	  VPF::set(filename, 
		   F[scan_keyword][0]);
	  range_map = VISImage<float>
	    (im_file.read(filename));
	  //		    range_map = range_map.gauss(2.0);

	  sprintf(scan_keyword, "TRANSFORM_%d", i);
	  try
	    {
	      VPF::set(int_tmp, 
		       F[scan_keyword][0]);
	      VPF::set(int_tmp, 
		       F[scan_keyword][1]);
	      for (j = 0; j < 16; j++)
		{
		  VPF::set(
			   transform.poke(j/4, j%4),
			   F[scan_keyword][j]);
		}
	    }
	  catch (VPF::Exception e)
	    {
	      // igore transform
	      transform = VISIdentity(4);
	    }

	  sprintf(scan_keyword, "R0_%d", i);
	  try
	    {
	      VPF::set(r0, F[scan_keyword][0]);
	    }
	  catch (VPF::Exception e)
	    {
	      if (r0_default != -FLT_MAX)
		r0 = r0_default;
	      else
		r0 = range_map.height()/2.0f - 0.5f;
	    }

	  sprintf(scan_keyword, "C0_%d", i);
	  try
	    {
	      VPF::set(c0, 
		       F[scan_keyword][0]);
	    }
	  catch (VPF::Exception e)
	    {
	      if (c0_default != -FLT_MAX)
		c0 = c0_default;
	      else
		c0 = range_map.width()/2.0f - 0.5f;
	    }

	  sprintf(scan_keyword, "DELTA_PHI_%d", i);
	  try
	    {
	      VPF::set(delta_phi, 
		       F[scan_keyword][0]);
	    }
	  catch (VPF::Exception e)
	    {
	      if (delta_phi_default != -FLT_MAX)
		delta_phi = delta_phi_default;
	      else
		delta_phi = DEFAULT_FOV/range_map.height();
	    }

	  sprintf(scan_keyword, "DELTA_THETA_%d", i);
	  try
	    {
	      VPF::set(delta_theta, 
		       F[scan_keyword][0]);
	    }
	  catch (VPF::Exception e)
	    {
	      if (delta_theta_default != -FLT_MAX)
		delta_theta = delta_theta_default;
	      else
		delta_theta = DEFAULT_FOV/range_map.width();
	    }

	  // put in the scale and the offset???
	  //range_map = range_scale*range_map + range_offset;

	  //		    cout << "range min " << range_map.min()
	  //			 << " and max " << range_map.max() << endl;

	  //		    cout << "about to do sphere create " << range_map.width()
	  //			 << range_map.height() << endl;
	  switch(scan_type)
	    {
	    case SCAN_SPHERICAL:
	      s_ptr = new SphericalScan(range_map, model_normals_scale);
	      break;
	    case SCAN_ORTHO:
	      s_ptr = new OrthoScan(range_map, model_normals_scale);
	      break;
	    }
	  //

	  // 		    s_ptr->setC0(r0);
	  // 		    s_ptr->setR0(c0);
	  // 		    s_ptr->setDeltaTheta(delta_phi);
	  // 		    s_ptr->setDeltaPhi(delta_theta);
	  // 		    s_ptr->setRangeOffset(range_offset);
	  // 		    s_ptr->setRangeDelta(range_scale);

	  params = VISVector(SphereScan::NUMPARAMS);
	  params.poke(SphereScan::R0) = r0;
	  params.poke(SphereScan::C0) = c0;
	  params.poke(SphereScan::DPHI) = delta_phi;
	  params.poke(SphereScan::DTHETA) = delta_theta;
	  params.poke(SphereScan::ROFFSET) = range_offset;
	  params.poke(SphereScan::RDELTA) = range_scale;
	  s_ptr->setParams(params);


	  // add confidence information if it exists
	  sprintf(scan_keyword, "CONFIDENCE_FILE_%d", i);
	  try
	    {
	      VPF::set(filename, 
		       F[scan_keyword][0]);
	      conf_im = VISImage<float>
		(im_file.read(filename));
	      if (conf_im.isValid())
		s_ptr->setConfidence(conf_im);
	      //			    cout << "conf min " << conf_im.min()
	      //				 << " and max " << conf_im.max() << endl;
	    }
	  catch (VPF::Exception e)
	    {
	      ;
	    }

	  addScan(s_ptr, transform);
	}
      catch (VPF::Exception e)
	{
	  cout << "WARNING: Did not find range file " << i 
	       << " -- will ignore." << endl;
	}
    }
}

ModelFit::~ModelFit()
{
    Scan* s_ptr;
    for (int i = 0; i < _scans.n(); i++)
      {
	if ((s_ptr = _scans.poke(i)) != NULL)
	  delete s_ptr;
      }

        if (_model != NULL)
	  delete _model;
}



// *******************************************************************
// *******************************************************************
// *************** STAND-ALONE FUNCTIONS  ****************************
// *******************************************************************
// *******************************************************************


// uses the convention of angles given in "Rigid3DModel"
VISMatrix eulerToTransform(const VISVector &euler) 
{
    VISMatrix T(4, 4);
    VISMatrix 	R_Phi, 	R_Omega, R_Kappa;

    R_Phi = VISIdentity(3);
    R_Omega = VISIdentity(3);
    R_Kappa = VISIdentity(3);

    if (euler.n() != 6)
	{
	    cout << "bad euler vector in eulerToTransform" << endl;
	    return VISMatrix();
	}

    float 
	phi = euler.peek(Euler::PHI), 
	kappa = euler.peek(Euler::KAPPA), 
	omega = euler.peek(Euler::OMEGA);

    R_Omega.poke(1, 1) = R_Omega.poke(2, 2) = cos(omega);
    R_Omega.poke(2, 1) = -1.0f*(R_Omega.poke(1, 2) = sin(omega));

    R_Phi.poke(0, 0) = R_Phi.poke(2, 2) = cos(phi);
    R_Phi.poke(0, 2) = -1.0f*(R_Phi.poke(2, 0) = sin(phi));
    
    R_Kappa.poke(1, 1) = R_Kappa.poke(0, 0) = cos(kappa);
    R_Kappa.poke(1, 0) = -1.0f*(R_Kappa.poke(0, 1) = sin(kappa));

    T = 0.0f;
    T.pokeROI(0, 3, euler.peekROI(Euler::TX, Euler::TZ));
    T.pokeROI(0, 0, R_Omega*R_Phi*R_Kappa);
    T.poke(3, 3) = 1.0f;

    return(T);
}

// must check this

VISVector transformToEuler(const VISMatrix &t)
{
    VISVector euler(6);
    float kappa = atan2(t.peek(0, 1), t.peek(0, 0));
    float omega = atan2(t.peek(1, 2), t.peek(2, 2));
//
// must get the sign correct here .....
//
//    float phi = atan2(t.peek(0, 2), 
//		fsqrt(pow(t.peek(0, 0), 2) 
//		      + pow(t.peek(0, 1), 2)));

    float phi = -asin(t.peek(0, 2));

    euler.poke(Euler::OMEGA) = omega;
    euler.poke(Euler::KAPPA) = kappa;
    euler.poke(Euler::PHI) = phi;
    euler.poke(Euler::TX) = t.peek(0, 3);
    euler.poke(Euler::TY) = t.peek(1, 3);
    euler.poke(Euler::TZ) = t.peek(2, 3);
    return(euler);
}


VISVector transformToQuaternion(const VISMatrix &transform)
{
  
  VISVector q;


  if ((transform.r() == 4)&&(transform.c() == 4))
    {
      q = VISVector(7);
      q.pokeROI(Quaternion::TX, transform.peekROI(0, 2, 3, 3));
    }
  else
    q = VISVector(4);
  float S, T = transform.trace();

  if (T > 0)
    {
      //      cout << "got case 0";
      S = 0.5f/sqrt(T);
      q.poke(0) = 0.25/S;
      q.poke(1) = (transform.peek(2, 1) - transform.peek(1, 2))*S;
      q.poke(2) = (transform.peek(0, 2) - transform.peek(2, 0))*S;
      q.poke(3) = (transform.peek(1, 0) - transform.peek(0, 1))*S;
    }
  else
    {
      if (transform.peek(0, 0) > transform.peek(1, 1))
	{
	  if (transform.peek(0, 0) > transform.peek(2, 2))
	    {
	      //	      cout << "got case 1";
	      S  = sqrt(1.0 + transform.peek(0, 0) 
			- transform.peek(1, 1) 
			- transform.peek(2, 2) ) * 2.0f;
	      q.poke(0) = (transform.peek(0, 2) + transform.peek(2, 0))/S;
	      q.poke(1) = 0.5/S;
	      q.poke(2) = (transform.peek(0, 1) + transform.peek(1, 0))/S;
	      q.poke(3) = (transform.peek(2, 1) + transform.peek(1, 2))/S;
	    }
	  else
	    {
	      //	      cout << "got case 2";
	      S  = sqrt(1.0 + transform.peek(2, 2) - transform.peek(0, 0) 
			- transform.peek(1, 1) )* 2.0f;
	  
	      q.poke(0) = (transform.peek(0, 1) + transform.peek(1, 0))/S;
	      q.poke(1) = (transform.peek(0, 2) + transform.peek(2, 0))/S;
	      q.poke(2) = (transform.peek(2, 1) + transform.peek(1, 2))/S;
	      q.poke(3) = 0.5/S;
	    }
	}
      else if (transform.peek(0, 0) > transform.peek(2, 2))
	{
	  //	      cout << "got case 3";
	      S  = sqrt(1.0 + transform.peek(1, 1) - transform.peek(0, 0) 
			- transform.peek(2, 2) )* 2.0f;
	  
	      q.poke(0) = (transform.peek(0, 2) + transform.peek(2, 0))/S;
	      q.poke(1) = (transform.peek(0, 1) + transform.peek(1, 0))/S;
	      q.poke(2) = 0.5/S;
	      q.poke(3) = (transform.peek(2, 1) + transform.peek(1, 2))/S;
	}
      else
	{
	  //	  cout << "got case 4";
	  S  = sqrt(1.0 + transform.peek(2, 2) - transform.peek(0, 0) 
		    - transform.peek(1, 1) )* 2.0f;
	  
	  q.poke(0) = (transform.peek(0, 1) + transform.peek(1, 0))/S;
	  q.poke(1) = (transform.peek(0, 2) + transform.peek(2, 0))/S;
	  q.poke(2) = (transform.peek(2, 1) + transform.peek(1, 2))/S;
	  q.poke(3) = 0.5/S;
	}
    }
  return(q);
}

VISMatrix quaternionToTransform(const VISVector &q)
{
  
  VISMatrix t = VISIdentity(4);

  if (q.n() < 4) 
    {
      cout << "WARNING:  qauternionToTransform - quaternion < 4" << endl;
      return VISMatrix();
    }
  
  t.poke(0, 0) = 1 - 2.0f*(pow(q.peek(2), 2) + pow(q.peek(3), 2));
  t.poke(1, 1) = 1 - 2.0f*(pow(q.peek(1), 2) + pow(q.peek(3), 2));
  t.poke(2, 2) = 1 - 2.0f*(pow(q.peek(1), 2) + pow(q.peek(2), 2));
  
  t.poke(0, 1) = 2.0f*(q.peek(1)*q.peek(2) - q.peek(0)*q.peek(3));
  t.poke(0, 2) = 2.0f*(q.peek(1)*q.peek(3) + q.peek(0)*q.peek(2));

  t.poke(1, 0) = 2.0f*(q.peek(1)*q.peek(2) + q.peek(0)*q.peek(3));
  t.poke(1, 2) = 2.0f*(q.peek(2)*q.peek(3) - q.peek(0)*q.peek(1));

  t.poke(2, 0) = 2.0f*(q.peek(1)*q.peek(3) - q.peek(0)*q.peek(2));
  t.poke(2, 1) = 2.0f*(q.peek(2)*q.peek(3) + q.peek(0)*q.peek(1));

  if (q.n() >= 7)
    t.pokeROI(0, 3, q.peekROI(Quaternion::TX, Quaternion::TZ));
  
  return(t);
}

// from the Terz and Metaxas paper PAMI 91
VISMatrix dx_dq(const VISVector &q, const VISVector &x)
{

  VISMatrix  G(3, 4), P(3, 3);

  G.poke(0, 0) = -1.0f*q.peek(1);
  G.poke(0, 1) = q.peek(0);
  G.poke(0, 2) = q.peek(3);
  G.poke(0, 3) = -1.0f*q.peek(2);

  G.poke(1, 0) = -1.0f*q.peek(2);
  G.poke(1, 1) = -1.0f*q.peek(3);
  G.poke(1, 2) = q.peek(0);
  G.poke(1, 3) = q.peek(1);

  G.poke(2, 0) = -1.0f*q.peek(3);
  G.poke(2, 1) = q.peek(2);
  G.poke(2, 2) = -1.0f*q.peek(1);
  G.poke(2, 3) = q.peek(0);

  P.poke(0, 0) = 0.0f;
  P.poke(0, 1) = -1.0f*x.peek(2);
  P.poke(0, 2) = x.peek(1);

  P.poke(1, 0) = x.peek(2);
  P.poke(1, 1) = 0.0f;
  P.poke(1, 2) = -1.0f*x.peek(0);

  P.poke(2, 0) = -1.0f*x.peek(1);
  P.poke(2, 1) = x.peek(0);
  P.poke(2, 2) = 0.0f;

  return(2.0f*(P*G));

}


void writeIVPoints(const VISMatrix &p, const char* filename)
{
  ofstream out(filename);
  
  if (p.r() < 2)
    return;

  out << "#Inventor V2.0 ascii \n Separator { \n Coordinate3 { \n point [ \n";
  for (int i = 0; i < p.c(); i++)
      out << p.peek(0, i) << " " << p.peek(1, i) << " " << p.peek(2, i) << "," << endl;

  out << "              ]" << endl;
  out << "     }" << endl;
  out << "  PointSet { \n startIndex  0\n numPoints   -1\n } \n } \n" << endl;
}


// ****************************************************
// ****************************************************
// ******************* CalibScan *******************
// ****************************************************
// ****************************************************

//
// registers one scan to a set of others
// use parameter files to get the appropriate parameters
//

CalibScan::CalibScan(VPF::ParameterFile F)
{
    int num_scans;
    float c0, r0, model_c0, model_r0, c0_default, r0_default;
    float 
      delta_phi, model_delta_phi, delta_phi_default, 
      delta_theta, model_delta_theta, delta_theta_default;
    float 
      model_range_offset, model_range_scale, 
      range_offset, range_scale, 
      range_offset_default, range_scale_default;

    Scan *s_ptr;
    char scan_keyword[80];
    char filename[80];
    VISMatrix transform(4,4);
    int int_tmp;
    int float_tmp;
    VISImage<float> range_map, conf_im;

    VISImageFile im_file;
    int i, j;
    int num_scales;
    
    int model_num_points;
    float model_normals_scale;
    float force_window;

    float x0, y0, z0, kappa, phi, omega;

    VISVector params;

// scanner geometry
#define SCAN_SPHERICAL (0)
#define SCAN_ORTHO (1)

    int scan_type;
    RangeModelCalib *model;

    try 
	{
	    VPF::set(scan_type, 
		     F["SCAN_TYPE"][0]);
	}
    catch (VPF::Exception e)
	{
	  scan_type = SCAN_SPHERICAL;
	}

    try
	{
	    VPF::set(range_offset, 
		     F["RANGE_OFFSET"][0]);
	}
    catch (VPF::Exception e)
	{
	    cout << "did not get range_offset" << endl;
	    range_offset = 0.0f;
	}
    
    try 
	{
	    VPF::set(range_scale, F["RANGE_SCALE"][0]);
	}
    catch (VPF::Exception e)
	{
	    range_scale = 1.0f;
	}


    try
	{
	    VPF::set(r0_default, 
		     F["R0"][0]);
	}
    catch (VPF::Exception e)
	{
	    r0_default = -FLT_MAX;
	}

    
    try 
	{
	    VPF::set(c0_default, 
		     F["C0"][0]);
	}
    catch (VPF::Exception e)
	{
	    c0_default = -FLT_MAX;
	}

    try 
	{
	    VPF::set(delta_phi_default, 
		     F["DELTA_PHI"][0]);
	}
    catch (VPF::Exception e)
	{
	    delta_phi_default = -FLT_MAX;
	}

    try 
	{
	    VPF::set(delta_theta_default, 
 		     F["DELTA_THETA"][0]);
	}
    catch (VPF::Exception e)
	{
	    delta_theta_default = -FLT_MAX;
	}

    try 
	{
	    VPF::set(num_scales, 
		     F["NUM_SCALES"][0]);
	}
    catch (VPF::Exception e)
	{
	    num_scales = 1;
	}

    try 
	{
	    VPF::set(model_num_points, 
		     F["MODEL_NUM_PTS"][0]);
	}
    catch (VPF::Exception e)
	{
	    model_num_points = -1;
	}

    try 
	{
	    VPF::set(model_normals_scale, 
		     F["MODEL_NORMALS_SCALE"][0]);
	}
    catch (VPF::Exception e)
	{
	    model_normals_scale = 0.0f;
	}

    try 
	{
	    VPF::set(force_window, 
		     F["FORCE_WINDOW"][0]);
	}
    catch (VPF::Exception e)
	{
	    force_window = DEFAULT_FORCE_WINDOW;
	}
    forceWindow(force_window);	


// model data and initialization
    try 
	{
	    VPF::set(filename, 
		     F["MODEL_FILE"][0]);

	    range_map = VISImage<float>
		(im_file.read(filename));
// put in the scale and the offset???
	    //	    range_map = range_scale*range_map + range_offset;
	    try 
		{
		    VPF::set(x0, 
			     F["X0"][0]);
		}
	    catch (VPF::Exception e)
		{
		    x0 = 0.0f;
		}

	    try 
		{
		    VPF::set(y0, 
			     F["Y0"][0]);
		}
	    catch (VPF::Exception e)
		{
		    y0 = 0.0f;
		}

	    try 
		{
		    VPF::set(z0, 
			     F["Z0"][0]);
		}
	    catch (VPF::Exception e)
		{
		    z0 = 0.0f;
		}

	    try 
		{
		    VPF::set(kappa, 
			     F["KAPPA"][0]);
		}
	    catch (VPF::Exception e)
		{
		    kappa = 0.0f;
		}

	    try 
		{
		    VPF::set(phi, 
			     F["PHI"][0]);
		}
	    catch (VPF::Exception e)
		{
		    phi = 0.0f;
		}

	    try 
		{
		    VPF::set(omega, 
			     F["OMEGA"][0]);
		}
	    catch (VPF::Exception e)
		{
		    omega = 0.0f;
		}
	}
    catch (VPF::Exception e)
	{
	    cout << "WARNING: Did not find model range file --- will ignore." << endl;
	}

    try 
	{
	    VPF::set(model_r0, 
		     F["MODEL_R0"][0]);
	}
    catch (VPF::Exception e)
	{
	    if (r0_default != -FLT_MAX)
		model_r0 = r0_default;
	    else
	      model_r0 = range_map.height()/2.0f - 0.5f;
	}

    try 
	{
	    VPF::set(model_c0, 
		     F["MODEL_C0"][0]);
	}
    catch (VPF::Exception e)
	{
	    if (c0_default != -FLT_MAX)
		model_c0 = c0_default;
	    else
		model_c0 = range_map.width()/2.0f - 0.5f;
	}
	try 
	  {
	    VPF::set(model_delta_phi, 
		     F["MODEL_DELTA_PHI"][0]);
	  }
	catch (VPF::Exception e)
	  {
	    if (delta_phi_default != -FLT_MAX)
	      model_delta_phi = delta_phi_default;
	    else
	      model_delta_phi = DEFAULT_FOV/range_map.height();
	  }

	try 
	  {
	    VPF::set(model_range_scale,
		     F["MODEL_RANGE_SCALE"][0]);
	  }
	catch (VPF::Exception e)
	  {
	      model_range_scale = range_scale;
	  }

	try 
	  {
	    VPF::set(model_range_offset,
		     F["MODEL_RANGE_OFFSET"][0]);
	  }
	catch (VPF::Exception e)
	  {
	      model_range_offset = range_offset;
	  }


	try
	  {
	    VPF::set(model_delta_theta, 
		     F["MODEL_DELTA_THETA"][0]);
	  }
	catch (VPF::Exception e)
	  {
	    if (delta_theta_default != -FLT_MAX)
	      model_delta_theta = delta_theta_default;
	    else
	      model_delta_theta = DEFAULT_FOV/range_map.width();
	  }

    // must accomodate different types of scans

    switch(scan_type)
      {
      case SCAN_SPHERICAL:
	for (int kk = 0;  kk < (int)(2.0f*pow(model_normals_scale, 2)); kk++)
	  range_map += 0.25f*range_map.anisoDiffuse(2.0);
	s_ptr = new SphericalScan(range_map, model_normals_scale);
	break;
      case SCAN_ORTHO:
	s_ptr = new OrthoScan(range_map, model_normals_scale);
	break;
      }
    //    s_ptr->setParams(VISVector(r0, c0, delta_phi, delta_theta));

    //    s_ptr->setParams(VISVector(r0, c0, delta_phi, delta_theta, range_offset, range_scale));
//     s_ptr->setC0(model_r0);
//     s_ptr->setR0(model_c0);
//     s_ptr->setDeltaTheta(model_delta_phi);
//     s_ptr->setDeltaPhi(model_delta_theta);
//     s_ptr->setRangeOffset(model_range_offset);
//     s_ptr->setRangeDelta(model_range_scale);

    params = VISVector(SphereScan::NUMPARAMS);
    params.poke(SphereScan::R0) = model_r0;
    params.poke(SphereScan::C0) = model_c0;
    params.poke(SphereScan::DPHI) = model_delta_phi;
    params.poke(SphereScan::DTHETA) = model_delta_theta;
    params.poke(SphereScan::ROFFSET) = model_range_offset;
    params.poke(SphereScan::RDELTA) = model_range_scale;
    s_ptr->setParams(params);
    //    cout << params;

//     if (num_scales > 1)
//       model = new RangeModelCalib(s_ptr, model_num_points);
//     else
//       // do you need ot set the number of scales here
    model = new RangeModelCalib(s_ptr, model_num_points);

//    scan = SphericalScan(range_map.gauss(model_normals_scale));

    VISVector euler(6);
    euler.poke(Euler::TX) = x0;
    euler.poke(Euler::TY) = y0;
    euler.poke(Euler::TZ) = z0;
    euler.poke(Euler::KAPPA) = kappa;
    euler.poke(Euler::PHI) = phi;
    euler.poke(Euler::OMEGA) = omega;

    model->beta(transformToQuaternion(eulerToTransform(euler))
		.concatRow(params.peekROI(2, 5)));
    //    cout <<     "trans  " <<  eulerToTransform(euler) << endl;   
    //    cout <<     "quat  " <<  transformToQuaternion(eulerToTransform(euler)) << endl;   
    //    cout <<     "quat test " <<  quaternionToTransform(transformToQuaternion(eulerToTransform(euler))) << endl;   
    //    cout <<     "beta " << model->Model::beta() << endl;
    setModel(model, num_scales);


// the scans
    try {
	VPF::set(num_scans, 
		 F["NUM_SCANS"][0]);
	    }
    catch (VPF::Exception e)
	{
	    cerr << "Argument 1 of \"NUM_SCANS\" is not valid" << endl;
	}
	    
    for (i = 0; i < num_scans; i++)
	{

	    sprintf(scan_keyword, "RANGE_FILE_%d", i);
	    try
		{
		    VPF::set(filename, 
			     F[scan_keyword][0]);
		    range_map = VISImage<float>
			(im_file.read(filename));

		    sprintf(scan_keyword, "TRANSFORM_%d", i);
		    try
			{
			    VPF::set(int_tmp, 
				     F[scan_keyword][0]);
			    VPF::set(int_tmp, 
				     F[scan_keyword][1]);
			    for (j = 0; j < 16; j++)
				{
				    VPF::set(
					transform.poke(j/4, j%4),
					F[scan_keyword][j]);
				}
			}
		    catch (VPF::Exception e)
			{
			    // igore transform
			    transform = VISIdentity(4);
			}

		    sprintf(scan_keyword, "R0_%d", i);
		    try
			{
			    VPF::set(r0, F[scan_keyword][0]);
			}
		    catch (VPF::Exception e)
			{
			    if (r0_default != -FLT_MAX)
				r0 = r0_default;
			    else
				r0 = range_map.height()/2.0f - 0.5f;
			}

		    sprintf(scan_keyword, "C0_%d", i);
		    try
			{
			    VPF::set(c0, 
				     F[scan_keyword][0]);
			}
		    catch (VPF::Exception e)
			{
			    if (c0_default != -FLT_MAX)
				c0 = c0_default;
			    else
				c0 = range_map.width()/2.0f - 0.5f;
			}

		    sprintf(scan_keyword, "DELTA_PHI_%d", i);
		    try
			{
			    VPF::set(delta_phi, 
				     F[scan_keyword][0]);
			}
		    catch (VPF::Exception e)
			{
			    if (delta_phi_default != -FLT_MAX)
				delta_phi = delta_phi_default;
			    else
				delta_phi = DEFAULT_FOV/range_map.height();
			}

		    sprintf(scan_keyword, "DELTA_THETA_%d", i);
		    try
			{
			    VPF::set(delta_theta, 
				     F[scan_keyword][0]);
			}
		    catch (VPF::Exception e)
			{
			    if (delta_theta_default != -FLT_MAX)
				delta_theta = delta_theta_default;
			    else
				delta_theta = DEFAULT_FOV/range_map.width();
			}

// put in the scale and the offset???
		    //		    range_map = range_scale*range_map + range_offset;

//		    cout << "range min " << range_map.min()
//			 << " and max " << range_map.max() << endl;

//		    cout << "about to do sphere create " << range_map.width()
//			 << range_map.height() << endl;
		    switch(scan_type)
		      {
		      case SCAN_SPHERICAL:
			s_ptr = new SphericalScan(range_map, model_normals_scale);
			break;
		      case SCAN_ORTHO:
			s_ptr = new OrthoScan(range_map, model_normals_scale);
			break;
		      }

// 		    s_ptr->setC0(r0);
// 		    s_ptr->setR0(c0);
// 		    s_ptr->setDeltaTheta(delta_phi);
// 		    s_ptr->setDeltaPhi(delta_theta);
// 		    s_ptr->setRangeOffset(range_offset);
// 		    s_ptr->setRangeDelta(range_scale);

    params = VISVector(SphereScan::NUMPARAMS);
    params.poke(SphereScan::R0) = r0;
    params.poke(SphereScan::C0) = c0;
    params.poke(SphereScan::DPHI) = delta_phi;
    params.poke(SphereScan::DTHETA) = delta_theta;
    params.poke(SphereScan::ROFFSET) = range_offset;
    params.poke(SphereScan::RDELTA) = range_scale;
    s_ptr->setParams(params);

    // s_ptr->setParams(VISVector(r0, c0, delta_phi, delta_theta));

// add confidence information if it exists
		    sprintf(scan_keyword, "CONFIDENCE_FILE_%d", i);
		    try
			{
			    VPF::set(filename, 
				     F[scan_keyword][0]);
			    conf_im = VISImage<float>
				(im_file.read(filename));
			    if (conf_im.isValid())
				s_ptr->setConfidence(conf_im);
//			    cout << "conf min " << conf_im.min()
//				 << " and max " << conf_im.max() << endl;
			}
		    catch (VPF::Exception e)
			{
			    ;
			}

		    addScan(s_ptr, transform);
		}
	    catch (VPF::Exception e)
		{
		    cout << "WARNING: Did not find range file " << i 
			 << " -- will ignore." << endl;
		}
	}
}


//*********************************************************************
//*********************************************************************
//********************** RangeModelCalib ***************************
//*********************************************************************
//*********************************************************************

// VISVector SphericalScan::rawToPoint(const VISVector &raw)
// {
//   float phi = (raw.peek(1) - _r0)*_delta_phi;
//   float theta = (raw.peek(0) - _c0)*_delta_theta;

// //    cout << "this scale factor " << this_scale << endl;
//   float r = _r_offset + _r_delta*raw.peek(2);

//   return(VISVector(r*sin(theta), r*sin(phi)*cos(theta), 
// 		   r*cos(phi)*cos(theta), 1.0f)); 
// }


// VISVector SphericalScan::rawToNormal(const VISVector &raw)
// {
//   VISVector r(4);
//   float di = _r_delta*raw.peek(3);
//   float dj =  _r_delta*raw.peek(4);
//   float cos_theta, sin_theta, cos_phi, sin_phi, 
//     dR_dPhi, dR_dTheta;

//   float range = raw.peek(2);

//   float phi = (raw.peek(1) - _r0)*_delta_phi,
// 	theta = (raw.peek(0) - _c0)*_delta_theta;
  
//   cos_theta = cos(theta);
//   sin_theta = sin(theta);
//   cos_phi = cos(phi);
//   sin_phi = sin(phi);

//   dR_dTheta = di/_delta_theta;
//   dR_dPhi = dj/_delta_phi;

//   r.poke(0) = cos_theta*(range*sin_theta - dR_dTheta*cos_theta);
//   r.poke(1) = cos_theta*sin_phi*(range*cos_theta + dR_dTheta*sin_theta)
// 	- dR_dPhi*cos_theta;
//   r.poke(2) = dR_dPhi*sin_phi
//     + cos_phi*cos_theta*(range*cos_theta + dR_dTheta*sin_theta);
//   r.poke(3) = 0.0f;
// // should this be positive or negative???
//     //    r /= -1.0f*r.norm();

//     r /= r.norm();
// //    cout << "normal length test " << r.dot(r) << endl;

//     return(r);
// }

void RangeModelCalib::updateObjPoints()
{  
  int w, h;
  int i, j, k;
  VISVector this_raw;
  VISImage<float> depth;
  w = _points_raw.c();
  if (w  > 0)
    {
      _points_obj = VISMatrix(4, w);
    for (i = 0; i < w; i++)
      _points_obj.pokeROI(0, i, _range_scan->get3DPoint
			  ((int)_points_raw.peek(0, i), 
			   (int)_points_raw.peek(1, i)
			   ) - _centroid);
    }
  else
    {
      k = 0;
      depth = _range_scan->depth();
      w = depth.width();
      h = depth.height();
      _points_obj = VISMatrix(4, (w - 2*BORDER)*(h - 2*BORDER));
      for (i = BORDER; i < (w - BORDER); i++)      
	for (j = BORDER; j < (h - BORDER); j++)
	  _points_obj.pokeROI(0, k++, _range_scan->get3DPoint(i, j)
			      - _centroid);
    }
}

void RangeModelCalib::updateObjNormals()
{  
  int w, h;
  int i, j, k;
  VISVector this_raw;
  VISImage<float> depth;
  w = _points_raw.c();
  if (w  > 0)
    {
      _normals_obj = VISMatrix(4, w);
      for (i = 0; i < w; i++)
      _normals_obj.pokeROI(0, i, _range_scan->getSurfaceNormal
			  ((int)_points_raw.peek(0, i), 
			   (int)_points_raw.peek(1, i)
			   ));
    }
  else
    {
      k = 0;
      depth = _range_scan->depth();
      w = depth.width();
      h = depth.height();
      _normals_obj = VISMatrix(4, (w - 2*BORDER)*(h - 2*BORDER));
      for (i = BORDER; i < (w - BORDER); i++)      
	for (j = BORDER; j < (h - BORDER); j++)
	  _normals_obj.pokeROI(0, k++, 
			       _range_scan->getSurfaceNormal(i, j));
    }
  //  cout << "normals obj " << _normals_obj;
}

VISMatrix RangeModelCalib::dS_dot_dB_deform(const VISMatrix &vects) const
{
float phi, theta, sin_theta, cos_theta, sin_phi, cos_phi, range;
//    cout << "this scale factor " << this_scale << endl;
//return(VISVector(r*sin(theta), r*sin(phi)*cos(theta), 
//  r*cos(phi)*cos(theta), 1.0f)); 
 int i, j, k;
 int w, h;
// 4 parameters in the spherical scan --- how to decide?
 VISMatrix total;
 VISVector raw;
 VISImage<float> depth;

// cout << "num_points in ds_dot..." << _num_points << endl;
 
 if (_num_points > 0)
   {
     total = VISMatrix(4, _num_points);
     for (j = 0; j < _num_points; j++)
	 total.pokeROI(0, j, ((_range_scan->dx_dparams
			       ((int)_points_raw.peek(0, j), 
				(int)_points_raw.peek(1, j))).t()*vects.vec(j)));
   }
 else
   {
     k = 0;
     depth = _range_scan->depth();
     w = depth.width();
     h = depth.height();
     total = VISMatrix(4, (w - 2*BORDER)*(h - 2*BORDER));
     for (i = BORDER; i < (w - BORDER); i++)      
       for (j = BORDER; j < (h - BORDER); j++)
	 {
	   total.pokeROI(0, k, (_range_scan->dx_dparams(i, j)).t()*vects.vec(k));
	   k++;
	 }
   }
 // cout << "dS_dB_deform out  " << total << endl;
 return(total);
}


VISMatrix SphericalScan::dx_dparams(int c, int r) const
{
  // 4 coords (homo),  4 free parameters
  VISMatrix ret(4, 4);

  float 
    phi = (r - _r0)*_delta_phi,
    theta = (c - _c0)*_delta_theta,
    sin_theta = sin(theta),
    cos_theta = cos(theta),
    sin_phi = sin(phi),
    cos_phi = cos(phi),
    range_raw, range;
  
  range_raw = _range_map_orig.peek(c, r);
  range = _r_offset + _r_delta*range_raw;

    ret.pokeROI(0, Calib::ROFFSET - 7, 
    VISVector(sin_theta, sin_phi*cos_theta, 
  	    cos_phi*cos_theta, 0.0f));
    //ret.pokeROI(0, Calib::ROFFSET - 7, VISVector(0.0f, 0.0f, 0.0f, 0.0f));
  ret.pokeROI(0, Calib::RDELTA - 7, 
	      range_raw*VISVector(sin_theta, sin_phi*cos_theta, 
				  cos_phi*cos_theta, 0.0f));
  ret.pokeROI(0, Calib::DTHETA - 7, ((c - _c0)*range)
	      *VISVector(cos_theta, -sin_phi*sin_theta, 
			 -cos_phi*sin_theta, 0.0f));
  ret.pokeROI(0, Calib::DPHI - 7, ((r - _r0)*range)
	      *VISVector(sin_theta,  cos_phi*cos_theta, 
			 -sin_phi*cos_theta, 0.0f));
  return(ret);
}

RangeModelCalib::RangeModelCalib(Scan *range_scan, 
				 int num_pts)
{ 
  _num_points = num_pts;
  _range_scan = range_scan;
  _range_scan->createScales(1);
  initialize(_range_scan, num_pts); 
  //  cout << "points raw" << _points_raw << endl;
  //	cout << "done create of 3D range model" << endl;
}


void RangeModelCalib::initialize(Scan *range_scan, 
				 int num_pts,
				 const VISVector &new_beta)
{
  //    cout << "init num points " << num_pts << endl;
  _range_scan = range_scan;
  //    cout << "_range_scan num scales " << _range_scan->numScales();
  //    cout << "_range_scan num scales " << _range_scan->numScales();
  VISVector p;
  VISImage<float> depth = _range_scan->depth();
  //  depth_raw = ((SphericalScan*)_range_scan)->depthRaw();
  int i, c, r, w = depth.width(), h = depth.height();

  if (num_pts > 0)
    {
      _points_raw = VISMatrix(2, num_pts);
      for (i = 0; i < num_pts; i++)
      {
	// get random sampling of c's and r's
	c = (int)(rand1()*(w - 2*BORDER - 1) + BORDER);
	r = (int)(rand1()*(h - 2*BORDER - 1) + BORDER);

	//	cout << "c and r " << c << " " << r << endl;
	//	points.poke(3, i) = 1.0f;
	_points_raw.pokeROI(0, i, VISVector((float)c, (float)r));
      }
    }
  else
    {
      //      _points_raw = VISMatrix(2, num_pts);
      //      for (c = BORDER; c < (w - BORDER); c++)
      //	for (r = BORDER; r < (h - BORDER); r++)
      //	  _points_raw.pokeROI(0, i, VISVector((float)c, (float)r));
      _points_raw = VISMatrix();
    }
      
  // all of the raw data is contained in "points"
  //  cout << "points raw" << _points_raw << endl;


  updateObjPoints();
  //  cout << "about to upd obj normals" << endl;
  updateObjNormals();

  // this is done just to take care of the centroid stuff
  if (new_beta.isValid())
    Rigid3DModel::initialize(_points_obj, _normals_obj, new_beta);
  else
    Rigid3DModel::initialize(_points_obj, _normals_obj);

  //  cout << "init points " << points.t() << endl;
  //  cout << "init normals " << normals.t() << endl;
  //    cout << "num pts obj" << _points_obj.c() << endl;
  //    cout << "done initialize on 3D range model" << endl;

}

int RangeModelCalib::beta(const VISVector &new_beta)
{
  VISVector scan_params = _range_scan->getParams();
  //  cout << new_beta << endl;
  scan_params.poke(SphereScan::RDELTA) = new_beta.peek(Calib::RDELTA);
  scan_params.poke(SphereScan::ROFFSET) = new_beta.peek(Calib::ROFFSET);
  scan_params.poke(SphereScan::DTHETA) = new_beta.peek(Calib::DTHETA);
  scan_params.poke(SphereScan::DPHI) = new_beta.peek(Calib::DPHI);
  //  cout << scan_params << endl;
  _range_scan->setParams(scan_params);
  Rigid3DModel::beta(new_beta);
  return(0);
}
