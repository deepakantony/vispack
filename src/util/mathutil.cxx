// ******************************************************************
// VISPack. Copyright (c) 1994-2000 Ross Whitaker rtw@utk.edu       *
// For conditions of distribution and use, see the file LICENSE.txt *
// accompanying this distribution.                                  *
// ******************************************************************
// 
// $Id: mathutil.cxx,v 1.2 2003/02/14 21:16:29 whitaker Exp $

#include "util/mathutil.h"


// returns a uniformly distributed deviate on the 0-1.0 interval
float rand1()
{
// why is "random" not availible
//    return((float)random()/(float)(LONG_MAX + 1);)
// rand returns values in the (0, 1) range
    return((float)rand()/(float)RAND_MAX);
}

float gasdev()
// returns a normally distributed deviate with zero mean an unit variance, 
// using ran1() as the source of uniform deviates
{
    static int iset=0;
    static float gset;
    float fac, rsq, v1, v2;
    
    if (iset == 0)
	{
	    do 
		{
		    v1 = 2.0f*rand1() - 1;
		    v2 = 2.0f*rand1() - 1;
		    rsq = v1*v1 + v2*v2;
		}
	    // make sure they are within the unit circle
	    while ((rsq >= 1.0f)||(rsq == 0.0f));
	    
	    fac = sqrt(-2.0*log(rsq)/rsq);
// now make the Box-Muller transformation to get two normal deviates.  
// Return one and save the other for next time	    
	    gset = v1*fac;
	    iset = 1;
	    return(v2*fac);
	}
    else
	{
	    iset = 0;
	    return gset;
	}
}


// returns a MEAN=1 gamma dist of order ia                                                
float gamdev(int ia)                                                                      
{                                                                                         
  int j;                                                                                  
  float am,e,s,v1,v2,x,y;                                                                 
                                                                                          
  if (ia < 1)                                                                             
    {                                                                                     
      cout << "Error in routine gamdev: bad int" << endl;                                 
      exit(-1);                                                                           
    }                                                                                     
  if (ia < 6) {                                                                           
    x=1.0;                                                                                
    for (j=1;j<=ia;j++) x *= rand1();                                                     
    x = -log(x);                                                                          
  } else {                                                                                
    do {                                                                                  
      do {                                                                                
        do {                                                                              
          v1=2.0*rand1()-1.0;                                                             
          v2=2.0*rand1()-1.0;                                                             
        } while (v1*v1+v2*v2 > 1.0);                                                      
        y=v2/v1;                                                                          
        am=ia-1;                                                                          
        s=sqrt(2.0*am+1.0);                                                               
        x=s*y+am;                                                                         
      } while (x <= 0.0);                                                                 
      e=(1.0+y*y)*exp(am*log(x/am)-s*y);                                                  
    } while (rand1() > e);                                                                
  }                                                                                       
  return(x/(float)ia);                                                                    
}          


//
// computes exp((-a^2)/2) using a lookuptable
//


float gaussFast(float a)
{
  a = fabs(a)/DELTA_GAUSS;
  if (a > 511.0)
    return(0.0f);
  int a_lo, a_hi;
  float a_in;

  a_lo = (int)floor(a);
  a_hi = (int)ceil(a);
  a_in = a - (float)a_lo;
    
  return(_GAUSS_TABLE[a_lo]*(1.0f - a_in) + 
	 _GAUSS_TABLE[a_hi]*(a_in));
}


float cosFast(float a)
{
    float factor = 1.0f;
    a = fabs(a);
    while (a > 2.0f*M_PI)
	a -= 2.0f*M_PI;

    if (a > M_PI)
	a = 2.0f*M_PI - a;
    if (a > M_PI/2.0f)
	{
	    factor = -1.0f;
	    a = M_PI - a;
	}
    
    a = a/DELTA_COS;
    int a_lo, a_hi;
    float a_in;

    a_lo = (int)floor(a);
    a_hi = (int)ceil(a);
    a_in = a - (float)a_lo;
    
    return(factor*(_COS_TABLE[a_lo]*(1.0f - a_in) + 
	   _COS_TABLE[a_hi]*(a_in)));
}


float acosFast(float a)
{
    float factor = 1.0f;
    if (a < 0.0f)
	factor = -1.0f;
    a = fabs(a);
    if (a > 1.0f)
	a = 1.0f;
    
    a = 511.0 - a/DELTA_ACOS;
    int a_lo, a_hi;
    float a_in;

    a_lo = (int)floor(a);
    a_hi = (int)ceil(a);
    a_in = a - (float)a_lo;
    
    if (factor > 0.0f)
	return((_ACOS_TABLE[a_lo]*(1.0f - a_in) + 
		_ACOS_TABLE[a_hi]*(a_in)));
    else 
	return(M_PI - ((_ACOS_TABLE[a_lo]*(1.0f - a_in) + 
			_ACOS_TABLE[a_hi]*(a_in))));
}


//
// computes int_{-/infty}^a exp((-x^2)/2) dx using a lookuptable
//


float gaussCumFast(float a)
{
    float pos;
    float result;
    pos = fabs(a)/DELTA_GAUSS;
    if (pos > 511.0)
	result = 1.0;
    else
	{
	    int a_lo, a_hi;
	    float a_in;

	    a_lo = (int)floor(pos);
	    a_hi = (int)ceil(pos);
	    a_in = pos - (float)a_lo;

	    result = _GAUSSCUM_TABLE[a_lo]*(1.0f - a_in) + 
		_GAUSSCUM_TABLE[a_hi]*(a_in);
	}

    if (a > 0.0f)
        return(result);
    else 
	return(1.0f - result);
}
