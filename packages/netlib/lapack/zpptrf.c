#include "f2c.h"

/* Subroutine */ int zpptrf_(char *uplo, integer *n, doublecomplex *ap, 
	integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    ZPPTRF computes the Cholesky factorization of a complex Hermitian   
    positive definite matrix A stored in packed format.   

    The factorization has the form   
       A = U**H * U,  if UPLO = 'U', or   
       A = L  * L**H,  if UPLO = 'L',   
    where U is an upper triangular matrix and L is lower triangular.   

    Arguments   
    =========   

    UPLO    (input) CHARACTER*1   
            = 'U':  Upper triangle of A is stored;   
            = 'L':  Lower triangle of A is stored.   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    AP      (input/output) COMPLEX*16 array, dimension (N*(N+1)/2)   
            On entry, the upper or lower triangle of the Hermitian matrix 
  
            A, packed columnwise in a linear array.  The j-th column of A 
  
            is stored in the array AP as follows:   
            if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;   
            if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.   
            See below for further details.   

            On exit, if INFO = 0, the triangular factor U or L from the   
            Cholesky factorization A = U**H*U or A = L*L**H, in the same 
  
            storage format as A.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   
            > 0:  if INFO = i, the leading minor of order i is not   
                  positive definite, and the factorization could not be   
                  completed.   

    Further Details   
    ===============   

    The packed storage scheme is illustrated by the following example   
    when N = 4, UPLO = 'U':   

    Two-dimensional storage of the Hermitian matrix A:   

       a11 a12 a13 a14   
           a22 a23 a24   
               a33 a34     (aij = conjg(aji))   
                   a44   

    Packed storage of the upper triangle of A:   

    AP = [ a11, a12, a22, a13, a23, a33, a14, a24, a34, a44 ]   

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__1 = 1;
    static doublereal c_b16 = -1.;
    
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1;
    doublecomplex z__1, z__2;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    extern /* Subroutine */ int zhpr_(char *, integer *, doublereal *, 
	    doublecomplex *, integer *, doublecomplex *);
    static integer j;
    extern logical lsame_(char *, char *);
    extern /* Double Complex */ VOID zdotc_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    static logical upper;
    extern /* Subroutine */ int ztpsv_(char *, char *, char *, integer *, 
	    doublecomplex *, doublecomplex *, integer *);
    static integer jc, jj;
    extern /* Subroutine */ int xerbla_(char *, integer *), zdscal_(
	    integer *, doublereal *, doublecomplex *, integer *);
    static doublereal ajj;



#define AP(I) ap[(I)-1]


    *info = 0;
    upper = lsame_(uplo, "U");
    if (! upper && ! lsame_(uplo, "L")) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("ZPPTRF", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return 0;
    }

    if (upper) {

/*        Compute the Cholesky factorization A = U'*U. */

	jj = 0;
	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
	    jc = jj + 1;
	    jj += j;

/*           Compute elements 1:J-1 of column J. */

	    if (j > 1) {
		i__2 = j - 1;
		ztpsv_("Upper", "Conjugate transpose", "Non-unit", &i__2, &AP(
			1), &AP(jc), &c__1);
	    }

/*           Compute U(J,J) and test for non-positive-definiteness
. */

	    i__2 = jj;
	    d__1 = AP(jj).r;
	    i__3 = j - 1;
	    zdotc_(&z__2, &i__3, &AP(jc), &c__1, &AP(jc), &c__1);
	    z__1.r = d__1 - z__2.r, z__1.i = -z__2.i;
	    ajj = z__1.r;
	    if (ajj <= 0.) {
		i__2 = jj;
		AP(jj).r = ajj, AP(jj).i = 0.;
		goto L30;
	    }
	    i__2 = jj;
	    d__1 = sqrt(ajj);
	    AP(jj).r = d__1, AP(jj).i = 0.;
/* L10: */
	}
    } else {

/*        Compute the Cholesky factorization A = L*L'. */

	jj = 1;
	i__1 = *n;
	for (j = 1; j <= *n; ++j) {

/*           Compute L(J,J) and test for non-positive-definiteness
. */

	    i__2 = jj;
	    ajj = AP(jj).r;
	    if (ajj <= 0.) {
		i__2 = jj;
		AP(jj).r = ajj, AP(jj).i = 0.;
		goto L30;
	    }
	    ajj = sqrt(ajj);
	    i__2 = jj;
	    AP(jj).r = ajj, AP(jj).i = 0.;

/*           Compute elements J+1:N of column J and update the tra
iling   
             submatrix. */

	    if (j < *n) {
		i__2 = *n - j;
		d__1 = 1. / ajj;
		zdscal_(&i__2, &d__1, &AP(jj + 1), &c__1);
		i__2 = *n - j;
		zhpr_("Lower", &i__2, &c_b16, &AP(jj + 1), &c__1, &AP(jj + *n 
			- j + 1));
		jj = jj + *n - j + 1;
	    }
/* L20: */
	}
    }
    goto L40;

L30:
    *info = j;

L40:
    return 0;

/*     End of ZPPTRF */

} /* zpptrf_ */
