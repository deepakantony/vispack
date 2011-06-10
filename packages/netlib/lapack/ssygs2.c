#include "f2c.h"

/* Subroutine */ int ssygs2_(integer *itype, char *uplo, integer *n, real *a, 
	integer *lda, real *b, integer *ldb, integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    SSYGS2 reduces a real symmetric-definite generalized eigenproblem   
    to standard form.   

    If ITYPE = 1, the problem is A*x = lambda*B*x,   
    and A is overwritten by inv(U')*A*inv(U) or inv(L)*A*inv(L')   

    If ITYPE = 2 or 3, the problem is A*B*x = lambda*x or   
    B*A*x = lambda*x, and A is overwritten by U*A*U` or L'*A*L.   

    B must have been previously factorized as U'*U or L*L' by SPOTRF.   

    Arguments   
    =========   

    ITYPE   (input) INTEGER   
            = 1: compute inv(U')*A*inv(U) or inv(L)*A*inv(L');   
            = 2 or 3: compute U*A*U' or L'*A*L.   

    UPLO    (input) CHARACTER   
            Specifies whether the upper or lower triangular part of the   
            symmetric matrix A is stored, and how B has been factorized. 
  
            = 'U':  Upper triangular   
            = 'L':  Lower triangular   

    N       (input) INTEGER   
            The order of the matrices A and B.  N >= 0.   

    A       (input/output) REAL array, dimension (LDA,N)   
            On entry, the symmetric matrix A.  If UPLO = 'U', the leading 
  
            n by n upper triangular part of A contains the upper   
            triangular part of the matrix A, and the strictly lower   
            triangular part of A is not referenced.  If UPLO = 'L', the   
            leading n by n lower triangular part of A contains the lower 
  
            triangular part of the matrix A, and the strictly upper   
            triangular part of A is not referenced.   

            On exit, if INFO = 0, the transformed matrix, stored in the   
            same format as A.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,N).   

    B       (input) REAL array, dimension (LDB,N)   
            The triangular factor from the Cholesky factorization of B,   
            as returned by SPOTRF.   

    LDB     (input) INTEGER   
            The leading dimension of the array B.  LDB >= max(1,N).   

    INFO    (output) INTEGER   
            = 0:  successful exit.   
            < 0:  if INFO = -i, the i-th argument had an illegal value.   

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static real c_b6 = -1.f;
    static integer c__1 = 1;
    static real c_b27 = 1.f;
    
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2;
    real r__1;
    /* Local variables */
    static integer k;
    extern /* Subroutine */ int ssyr2_(char *, integer *, real *, real *, 
	    integer *, real *, integer *, real *, integer *);
    extern logical lsame_(char *, char *);
    extern /* Subroutine */ int sscal_(integer *, real *, real *, integer *);
    static logical upper;
    extern /* Subroutine */ int saxpy_(integer *, real *, real *, integer *, 
	    real *, integer *), strmv_(char *, char *, char *, integer *, 
	    real *, integer *, real *, integer *), 
	    strsv_(char *, char *, char *, integer *, real *, integer *, real 
	    *, integer *);
    static real ct;
    extern /* Subroutine */ int xerbla_(char *, integer *);
    static real akk, bkk;




#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]
#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]

    *info = 0;
    upper = lsame_(uplo, "U");
    if (*itype < 1 || *itype > 3) {
	*info = -1;
    } else if (! upper && ! lsame_(uplo, "L")) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*lda < max(1,*n)) {
	*info = -5;
    } else if (*ldb < max(1,*n)) {
	*info = -7;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("SSYGS2", &i__1);
	return 0;
    }

    if (*itype == 1) {
	if (upper) {

/*           Compute inv(U')*A*inv(U) */

	    i__1 = *n;
	    for (k = 1; k <= *n; ++k) {

/*              Update the upper triangle of A(k:n,k:n) */

		akk = A(k,k);
		bkk = B(k,k);
/* Computing 2nd power */
		r__1 = bkk;
		akk /= r__1 * r__1;
		A(k,k) = akk;
		if (k < *n) {
		    i__2 = *n - k;
		    r__1 = 1.f / bkk;
		    sscal_(&i__2, &r__1, &A(k,k+1), lda);
		    ct = akk * -.5f;
		    i__2 = *n - k;
		    saxpy_(&i__2, &ct, &B(k,k+1), ldb, &A(k,k+1), lda);
		    i__2 = *n - k;
		    ssyr2_(uplo, &i__2, &c_b6, &A(k,k+1), lda, 
			    &B(k,k+1), ldb, &A(k+1,k+1), lda);
		    i__2 = *n - k;
		    saxpy_(&i__2, &ct, &B(k,k+1), ldb, &A(k,k+1), lda);
		    i__2 = *n - k;
		    strsv_(uplo, "Transpose", "Non-unit", &i__2, &B(k+1,k+1), ldb, &A(k,k+1), 
			    lda);
		}
/* L10: */
	    }
	} else {

/*           Compute inv(L)*A*inv(L') */

	    i__1 = *n;
	    for (k = 1; k <= *n; ++k) {

/*              Update the lower triangle of A(k:n,k:n) */

		akk = A(k,k);
		bkk = B(k,k);
/* Computing 2nd power */
		r__1 = bkk;
		akk /= r__1 * r__1;
		A(k,k) = akk;
		if (k < *n) {
		    i__2 = *n - k;
		    r__1 = 1.f / bkk;
		    sscal_(&i__2, &r__1, &A(k+1,k), &c__1);
		    ct = akk * -.5f;
		    i__2 = *n - k;
		    saxpy_(&i__2, &ct, &B(k+1,k), &c__1, &A(k+1,k), &c__1);
		    i__2 = *n - k;
		    ssyr2_(uplo, &i__2, &c_b6, &A(k+1,k), &c__1, 
			    &B(k+1,k), &c__1, &A(k+1,k+1), lda);
		    i__2 = *n - k;
		    saxpy_(&i__2, &ct, &B(k+1,k), &c__1, &A(k+1,k), &c__1);
		    i__2 = *n - k;
		    strsv_(uplo, "No transpose", "Non-unit", &i__2, &B(k+1,k+1), ldb, &A(k+1,k), 
			    &c__1);
		}
/* L20: */
	    }
	}
    } else {
	if (upper) {

/*           Compute U*A*U' */

	    i__1 = *n;
	    for (k = 1; k <= *n; ++k) {

/*              Update the upper triangle of A(1:k,1:k) */

		akk = A(k,k);
		bkk = B(k,k);
		i__2 = k - 1;
		strmv_(uplo, "No transpose", "Non-unit", &i__2, &B(1,1), 
			ldb, &A(1,k), &c__1);
		ct = akk * .5f;
		i__2 = k - 1;
		saxpy_(&i__2, &ct, &B(1,k), &c__1, &A(1,k), &c__1);
		i__2 = k - 1;
		ssyr2_(uplo, &i__2, &c_b27, &A(1,k), &c__1, &B(1,k), &c__1, &A(1,1), lda);
		i__2 = k - 1;
		saxpy_(&i__2, &ct, &B(1,k), &c__1, &A(1,k), &c__1);
		i__2 = k - 1;
		sscal_(&i__2, &bkk, &A(1,k), &c__1);
/* Computing 2nd power */
		r__1 = bkk;
		A(k,k) = akk * (r__1 * r__1);
/* L30: */
	    }
	} else {

/*           Compute L'*A*L */

	    i__1 = *n;
	    for (k = 1; k <= *n; ++k) {

/*              Update the lower triangle of A(1:k,1:k) */

		akk = A(k,k);
		bkk = B(k,k);
		i__2 = k - 1;
		strmv_(uplo, "Transpose", "Non-unit", &i__2, &B(1,1), 
			ldb, &A(k,1), lda);
		ct = akk * .5f;
		i__2 = k - 1;
		saxpy_(&i__2, &ct, &B(k,1), ldb, &A(k,1), lda);
		i__2 = k - 1;
		ssyr2_(uplo, &i__2, &c_b27, &A(k,1), lda, &B(k,1), ldb, &A(1,1), lda);
		i__2 = k - 1;
		saxpy_(&i__2, &ct, &B(k,1), ldb, &A(k,1), lda);
		i__2 = k - 1;
		sscal_(&i__2, &bkk, &A(k,1), lda);
/* Computing 2nd power */
		r__1 = bkk;
		A(k,k) = akk * (r__1 * r__1);
/* L40: */
	    }
	}
    }
    return 0;

/*     End of SSYGS2 */

} /* ssygs2_ */
