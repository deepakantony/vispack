#include "f2c.h"

/* Subroutine */ int cung2l_(integer *m, integer *n, integer *k, complex *a, 
	integer *lda, complex *tau, complex *work, integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    CUNG2L generates an m by n complex matrix Q with orthonormal columns, 
  
    which is defined as the last n columns of a product of k elementary   
    reflectors of order m   

          Q  =  H(k) . . . H(2) H(1)   

    as returned by CGEQLF.   

    Arguments   
    =========   

    M       (input) INTEGER   
            The number of rows of the matrix Q. M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix Q. M >= N >= 0.   

    K       (input) INTEGER   
            The number of elementary reflectors whose product defines the 
  
            matrix Q. N >= K >= 0.   

    A       (input/output) COMPLEX array, dimension (LDA,N)   
            On entry, the (n-k+i)-th column must contain the vector which 
  
            defines the elementary reflector H(i), for i = 1,2,...,k, as 
  
            returned by CGEQLF in the last k columns of its array   
            argument A.   
            On exit, the m-by-n matrix Q.   

    LDA     (input) INTEGER   
            The first dimension of the array A. LDA >= max(1,M).   

    TAU     (input) COMPLEX array, dimension (K)   
            TAU(i) must contain the scalar factor of the elementary   
            reflector H(i), as returned by CGEQLF.   

    WORK    (workspace) COMPLEX array, dimension (N)   

    INFO    (output) INTEGER   
            = 0: successful exit   
            < 0: if INFO = -i, the i-th argument has an illegal value   

    ===================================================================== 
  


       Test the input arguments   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__1 = 1;
    
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    complex q__1;
    /* Local variables */
    static integer i, j, l;
    extern /* Subroutine */ int cscal_(integer *, complex *, complex *, 
	    integer *), clarf_(char *, integer *, integer *, complex *, 
	    integer *, complex *, complex *, integer *, complex *);
    static integer ii;
    extern /* Subroutine */ int xerbla_(char *, integer *);



#define TAU(I) tau[(I)-1]
#define WORK(I) work[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    *info = 0;
    if (*m < 0) {
	*info = -1;
    } else if (*n < 0 || *n > *m) {
	*info = -2;
    } else if (*k < 0 || *k > *n) {
	*info = -3;
    } else if (*lda < max(1,*m)) {
	*info = -5;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("CUNG2L", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*n <= 0) {
	return 0;
    }

/*     Initialise columns 1:n-k to columns of the unit matrix */

    i__1 = *n - *k;
    for (j = 1; j <= *n-*k; ++j) {
	i__2 = *m;
	for (l = 1; l <= *m; ++l) {
	    i__3 = l + j * a_dim1;
	    A(l,j).r = 0.f, A(l,j).i = 0.f;
/* L10: */
	}
	i__2 = *m - *n + j + j * a_dim1;
	A(*m-*n+j,j).r = 1.f, A(*m-*n+j,j).i = 0.f;
/* L20: */
    }

    i__1 = *k;
    for (i = 1; i <= *k; ++i) {
	ii = *n - *k + i;

/*        Apply H(i) to A(1:m-k+i,1:n-k+i) from the left */

	i__2 = *m - *n + ii + ii * a_dim1;
	A(*m-*n+ii,ii).r = 1.f, A(*m-*n+ii,ii).i = 0.f;
	i__2 = *m - *n + ii;
	i__3 = ii - 1;
	clarf_("Left", &i__2, &i__3, &A(1,ii), &c__1, &TAU(i), &A(1,1), lda, &WORK(1));
	i__2 = *m - *n + ii - 1;
	i__3 = i;
	q__1.r = -(doublereal)TAU(i).r, q__1.i = -(doublereal)TAU(i).i;
	cscal_(&i__2, &q__1, &A(1,ii), &c__1);
	i__2 = *m - *n + ii + ii * a_dim1;
	i__3 = i;
	q__1.r = 1.f - TAU(i).r, q__1.i = 0.f - TAU(i).i;
	A(*m-*n+ii,ii).r = q__1.r, A(*m-*n+ii,ii).i = q__1.i;

/*        Set A(m-k+i+1:m,n-k+i) to zero */

	i__2 = *m;
	for (l = *m - *n + ii + 1; l <= *m; ++l) {
	    i__3 = l + ii * a_dim1;
	    A(l,ii).r = 0.f, A(l,ii).i = 0.f;
/* L30: */
	}
/* L40: */
    }
    return 0;

/*     End of CUNG2L */

} /* cung2l_ */

