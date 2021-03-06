#include "f2c.h"

integer icmax1_(integer *n, complex *cx, integer *incx)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    ICMAX1 finds the index of the element whose real part has maximum   
    absolute value.   

    Based on ICAMAX from Level 1 BLAS.   
    The change is to use the 'genuine' absolute value.   

    Contributed by Nick Higham for use with CLACON.   

    Arguments   
    =========   

    N       (input) INTEGER   
            The number of elements in the vector CX.   

    CX      (input) COMPLEX array, dimension (N)   
            The vector whose elements will be summed.   

    INCX    (input) INTEGER   
            The spacing between successive values of CX.  INCX >= 1.   

   ===================================================================== 
  


       NEXT LINE IS THE ONLY MODIFICATION.   

    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    integer ret_val, i__1, i__2;
    real r__1;
    /* Local variables */
    static real smax;
    static integer i, ix;


#define CX(I) cx[(I)-1]


    ret_val = 0;
    if (*n < 1) {
	return ret_val;
    }
    ret_val = 1;
    if (*n == 1) {
	return ret_val;
    }
    if (*incx == 1) {
	goto L30;
    }

/*     CODE FOR INCREMENT NOT EQUAL TO 1 */

    ix = 1;
    smax = (r__1 = CX(1).r, dabs(r__1));
    ix += *incx;
    i__1 = *n;
    for (i = 2; i <= *n; ++i) {
	i__2 = ix;
	if ((r__1 = CX(ix).r, dabs(r__1)) <= smax) {
	    goto L10;
	}
	ret_val = i;
	i__2 = ix;
	smax = (r__1 = CX(ix).r, dabs(r__1));
L10:
	ix += *incx;
/* L20: */
    }
    return ret_val;

/*     CODE FOR INCREMENT EQUAL TO 1 */

L30:
    smax = (r__1 = CX(1).r, dabs(r__1));
    i__1 = *n;
    for (i = 2; i <= *n; ++i) {
	i__2 = i;
	if ((r__1 = CX(i).r, dabs(r__1)) <= smax) {
	    goto L40;
	}
	ret_val = i;
	i__2 = i;
	smax = (r__1 = CX(i).r, dabs(r__1));
L40:
	;
    }
    return ret_val;

/*     End of ICMAX1 */

} /* icmax1_ */

