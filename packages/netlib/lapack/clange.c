#include "f2c.h"

doublereal clange_(char *norm, integer *m, integer *n, complex *a, integer *
	lda, real *work)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    CLANGE  returns the value of the one norm,  or the Frobenius norm, or 
  
    the  infinity norm,  or the  element of  largest absolute value  of a 
  
    complex matrix A.   

    Description   
    ===========   

    CLANGE returns the value   

       CLANGE = ( max(abs(A(i,j))), NORM = 'M' or 'm'   
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
            Specifies the value to be returned in CLANGE as described   
            above.   

    M       (input) INTEGER   
            The number of rows of the matrix A.  M >= 0.  When M = 0,   
            CLANGE is set to zero.   

    N       (input) INTEGER   
            The number of columns of the matrix A.  N >= 0.  When N = 0, 
  
            CLANGE is set to zero.   

    A       (input) COMPLEX array, dimension (LDA,N)   
            The m by n matrix A.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(M,1).   

    WORK    (workspace) REAL array, dimension (LWORK),   
            where LWORK >= M when NORM = 'I'; otherwise, WORK is not   
            referenced.   

   ===================================================================== 
  


    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__1 = 1;
    
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    real ret_val, r__1, r__2;
    /* Builtin functions */
    double c_abs(complex *), sqrt(doublereal);
    /* Local variables */
    static integer i, j;
    static real scale;
    extern logical lsame_(char *, char *);
    static real value;
    extern /* Subroutine */ int classq_(integer *, complex *, integer *, real 
	    *, real *);
    static real sum;



#define WORK(I) work[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    if (min(*m,*n) == 0) {
	value = 0.f;
    } else if (lsame_(norm, "M")) {

/*        Find max(abs(A(i,j))). */

	value = 0.f;
	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
	    i__2 = *m;
	    for (i = 1; i <= *m; ++i) {
/* Computing MAX */
		r__1 = value, r__2 = c_abs(&A(i,j));
		value = dmax(r__1,r__2);
/* L10: */
	    }
/* L20: */
	}
    } else if (lsame_(norm, "O") || *(unsigned char *)norm == '1') {

/*        Find norm1(A). */

	value = 0.f;
	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
	    sum = 0.f;
	    i__2 = *m;
	    for (i = 1; i <= *m; ++i) {
		sum += c_abs(&A(i,j));
/* L30: */
	    }
	    value = dmax(value,sum);
/* L40: */
	}
    } else if (lsame_(norm, "I")) {

/*        Find normI(A). */

	i__1 = *m;
	for (i = 1; i <= *m; ++i) {
	    WORK(i) = 0.f;
/* L50: */
	}
	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
	    i__2 = *m;
	    for (i = 1; i <= *m; ++i) {
		WORK(i) += c_abs(&A(i,j));
/* L60: */
	    }
/* L70: */
	}
	value = 0.f;
	i__1 = *m;
	for (i = 1; i <= *m; ++i) {
/* Computing MAX */
	    r__1 = value, r__2 = WORK(i);
	    value = dmax(r__1,r__2);
/* L80: */
	}
    } else if (lsame_(norm, "F") || lsame_(norm, "E")) {

/*        Find normF(A). */

	scale = 0.f;
	sum = 1.f;
	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
	    classq_(m, &A(1,j), &c__1, &scale, &sum);
/* L90: */
	}
	value = scale * sqrt(sum);
    }

    ret_val = value;
    return ret_val;

/*     End of CLANGE */

} /* clange_ */

