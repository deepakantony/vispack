#include "f2c.h"

doublereal clangt_(char *norm, integer *n, complex *dl, complex *d, complex *
	du)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    CLANGT  returns the value of the one norm,  or the Frobenius norm, or 
  
    the  infinity norm,  or the  element of  largest absolute value  of a 
  
    complex tridiagonal matrix A.   

    Description   
    ===========   

    CLANGT returns the value   

       CLANGT = ( max(abs(A(i,j))), NORM = 'M' or 'm'   
                (   
                ( norm1(A),         NORM = '1', 'O' or 'o'   
                (   
                ( normI(A),         NORM = 'I' or 'i'   
                (   
                ( normF(A),         NORM = 'F', 'f', 'E' or 'e'   

    where  norm1  denotes the  one norm of a matrix (maximum column sum), 
  
    normI  denotes the  infinity norm  of a matrix  (maximum row sum) and 
  
    normF  denotes the  Frobenius norm of a matrix (square root of sum of 
  
    squares).  Note that  max(abs(A(i,j)))  is not a  matrix norm.   

    Arguments   
    =========   

    NORM    (input) CHARACTER*1   
            Specifies the value to be returned in CLANGT as described   
            above.   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.  When N = 0, CLANGT is   
            set to zero.   

    DL      (input) COMPLEX array, dimension (N-1)   
            The (n-1) sub-diagonal elements of A.   

    D       (input) COMPLEX array, dimension (N)   
            The diagonal elements of A.   

    DU      (input) COMPLEX array, dimension (N-1)   
            The (n-1) super-diagonal elements of A.   

    ===================================================================== 
  


    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__1 = 1;
    
    /* System generated locals */
    integer i__1;
    real ret_val, r__1, r__2;
    /* Builtin functions */
    double c_abs(complex *), sqrt(doublereal);
    /* Local variables */
    static integer i;
    static real scale;
    extern logical lsame_(char *, char *);
    static real anorm;
    extern /* Subroutine */ int classq_(integer *, complex *, integer *, real 
	    *, real *);
    static real sum;



#define DU(I) du[(I)-1]
#define D(I) d[(I)-1]
#define DL(I) dl[(I)-1]


    if (*n <= 0) {
	anorm = 0.f;
    } else if (lsame_(norm, "M")) {

/*        Find max(abs(A(i,j))). */

	anorm = c_abs(&D(*n));
	i__1 = *n - 1;
	for (i = 1; i <= *n-1; ++i) {
/* Computing MAX */
	    r__1 = anorm, r__2 = c_abs(&DL(i));
	    anorm = dmax(r__1,r__2);
/* Computing MAX */
	    r__1 = anorm, r__2 = c_abs(&D(i));
	    anorm = dmax(r__1,r__2);
/* Computing MAX */
	    r__1 = anorm, r__2 = c_abs(&DU(i));
	    anorm = dmax(r__1,r__2);
/* L10: */
	}
    } else if (lsame_(norm, "O") || *(unsigned char *)norm == '1') {

/*        Find norm1(A). */

	if (*n == 1) {
	    anorm = c_abs(&D(1));
	} else {
/* Computing MAX */
	    r__1 = c_abs(&D(1)) + c_abs(&DL(1)), r__2 = c_abs(&D(*n)) + c_abs(
		    &DU(*n - 1));
	    anorm = dmax(r__1,r__2);
	    i__1 = *n - 1;
	    for (i = 2; i <= *n-1; ++i) {
/* Computing MAX */
		r__1 = anorm, r__2 = c_abs(&D(i)) + c_abs(&DL(i)) + c_abs(&DU(
			i - 1));
		anorm = dmax(r__1,r__2);
/* L20: */
	    }
	}
    } else if (lsame_(norm, "I")) {

/*        Find normI(A). */

	if (*n == 1) {
	    anorm = c_abs(&D(1));
	} else {
/* Computing MAX */
	    r__1 = c_abs(&D(1)) + c_abs(&DU(1)), r__2 = c_abs(&D(*n)) + c_abs(
		    &DL(*n - 1));
	    anorm = dmax(r__1,r__2);
	    i__1 = *n - 1;
	    for (i = 2; i <= *n-1; ++i) {
/* Computing MAX */
		r__1 = anorm, r__2 = c_abs(&D(i)) + c_abs(&DU(i)) + c_abs(&DL(
			i - 1));
		anorm = dmax(r__1,r__2);
/* L30: */
	    }
	}
    } else if (lsame_(norm, "F") || lsame_(norm, "E")) {

/*        Find normF(A). */

	scale = 0.f;
	sum = 1.f;
	classq_(n, &D(1), &c__1, &scale, &sum);
	if (*n > 1) {
	    i__1 = *n - 1;
	    classq_(&i__1, &DL(1), &c__1, &scale, &sum);
	    i__1 = *n - 1;
	    classq_(&i__1, &DU(1), &c__1, &scale, &sum);
	}
	anorm = scale * sqrt(sum);
    }

    ret_val = anorm;
    return ret_val;

/*     End of CLANGT */

} /* clangt_ */

