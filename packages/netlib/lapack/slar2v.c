#include "f2c.h"

/* Subroutine */ int slar2v_(integer *n, real *x, real *y, real *z, integer *
	incx, real *c, real *s, integer *incc)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    SLAR2V applies a vector of real plane rotations from both sides to   
    a sequence of 2-by-2 real symmetric matrices, defined by the elements 
  
    of the vectors x, y and z. For i = 1,2,...,n   

       ( x(i)  z(i) ) := (  c(i)  s(i) ) ( x(i)  z(i) ) ( c(i) -s(i) )   
       ( z(i)  y(i) )    ( -s(i)  c(i) ) ( z(i)  y(i) ) ( s(i)  c(i) )   

    Arguments   
    =========   

    N       (input) INTEGER   
            The number of plane rotations to be applied.   

    X       (input/output) REAL array,   
                           dimension (1+(N-1)*INCX)   
            The vector x.   

    Y       (input/output) REAL array,   
                           dimension (1+(N-1)*INCX)   
            The vector y.   

    Z       (input/output) REAL array,   
                           dimension (1+(N-1)*INCX)   
            The vector z.   

    INCX    (input) INTEGER   
            The increment between elements of X, Y and Z. INCX > 0.   

    C       (input) REAL array, dimension (1+(N-1)*INCC)   
            The cosines of the plane rotations.   

    S       (input) REAL array, dimension (1+(N-1)*INCC)   
            The sines of the plane rotations.   

    INCC    (input) INTEGER   
            The increment between elements of C and S. INCC > 0.   

    ===================================================================== 
  


    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    integer i__1;
    /* Local variables */
    static integer i;
    static real t1, t2, t3, t4, t5, t6;
    static integer ic;
    static real ci, si;
    static integer ix;
    static real xi, yi, zi;


#define S(I) s[(I)-1]
#define C(I) c[(I)-1]
#define Z(I) z[(I)-1]
#define Y(I) y[(I)-1]
#define X(I) x[(I)-1]


    ix = 1;
    ic = 1;
    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	xi = X(ix);
	yi = Y(ix);
	zi = Z(ix);
	ci = C(ic);
	si = S(ic);
	t1 = si * zi;
	t2 = ci * zi;
	t3 = t2 - si * xi;
	t4 = t2 + si * yi;
	t5 = ci * xi + t1;
	t6 = ci * yi - t1;
	X(ix) = ci * t5 + si * t4;
	Y(ix) = ci * t6 - si * t3;
	Z(ix) = ci * t4 - si * t5;
	ix += *incx;
	ic += *incc;
/* L10: */
    }

/*     End of SLAR2V */

    return 0;
} /* slar2v_ */

