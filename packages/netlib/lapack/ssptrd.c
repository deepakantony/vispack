#include "f2c.h"

/* Subroutine */ int ssptrd_(char *uplo, integer *n, real *ap, real *d, real *
	e, real *tau, integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       March 31, 1993   


    Purpose   
    =======   

    SSPTRD reduces a real symmetric matrix A stored in packed form to   
    symmetric tridiagonal form T by an orthogonal similarity   
    transformation: Q**T * A * Q = T.   

    Arguments   
    =========   

    UPLO    (input) CHARACTER*1   
            = 'U':  Upper triangle of A is stored;   
            = 'L':  Lower triangle of A is stored.   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    AP      (input/output) REAL array, dimension (N*(N+1)/2)   
            On entry, the upper or lower triangle of the symmetric matrix 
  
            A, packed columnwise in a linear array.  The j-th column of A 
  
            is stored in the array AP as follows:   
            if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;   
            if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n. 
  
            On exit, if UPLO = 'U', the diagonal and first superdiagonal 
  
            of A are overwritten by the corresponding elements of the   
            tridiagonal matrix T, and the elements above the first   
            superdiagonal, with the array TAU, represent the orthogonal   
            matrix Q as a product of elementary reflectors; if UPLO   
            = 'L', the diagonal and first subdiagonal of A are over-   
            written by the corresponding elements of the tridiagonal   
            matrix T, and the elements below the first subdiagonal, with 
  
            the array TAU, represent the orthogonal matrix Q as a product 
  
            of elementary reflectors. See Further Details.   

    D       (output) REAL array, dimension (N)   
            The diagonal elements of the tridiagonal matrix T:   
            D(i) = A(i,i).   

    E       (output) REAL array, dimension (N-1)   
            The off-diagonal elements of the tridiagonal matrix T:   
            E(i) = A(i,i+1) if UPLO = 'U', E(i) = A(i+1,i) if UPLO = 'L'. 
  

    TAU     (output) REAL array, dimension (N-1)   
            The scalar factors of the elementary reflectors (see Further 
  
            Details).   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   

    Further Details   
    ===============   

    If UPLO = 'U', the matrix Q is represented as a product of elementary 
  
    reflectors   

       Q = H(n-1) . . . H(2) H(1).   

    Each H(i) has the form   

       H(i) = I - tau * v * v'   

    where tau is a real scalar, and v is a real vector with   
    v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in AP,   
    overwriting A(1:i-1,i+1), and tau is stored in TAU(i).   

    If UPLO = 'L', the matrix Q is represented as a product of elementary 
  
    reflectors   

       Q = H(1) H(2) . . . H(n-1).   

    Each H(i) has the form   

       H(i) = I - tau * v * v'   

    where tau is a real scalar, and v is a real vector with   
    v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in AP,   
    overwriting A(i+2:n,i), and tau is stored in TAU(i).   

    ===================================================================== 
  


       Test the input parameters   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__1 = 1;
    static real c_b8 = 0.f;
    static real c_b14 = -1.f;
    
    /* System generated locals */
    integer i__1, i__2;
    /* Local variables */
    static real taui;
    extern doublereal sdot_(integer *, real *, integer *, real *, integer *);
    static integer i;
    extern /* Subroutine */ int sspr2_(char *, integer *, real *, real *, 
	    integer *, real *, integer *, real *);
    static real alpha;
    extern logical lsame_(char *, char *);
    static integer i1;
    static logical upper;
    extern /* Subroutine */ int saxpy_(integer *, real *, real *, integer *, 
	    real *, integer *), sspmv_(char *, integer *, real *, real *, 
	    real *, integer *, real *, real *, integer *);
    static integer ii;
    extern /* Subroutine */ int xerbla_(char *, integer *), slarfg_(
	    integer *, real *, real *, integer *, real *);
    static integer i1i1;



#define TAU(I) tau[(I)-1]
#define E(I) e[(I)-1]
#define D(I) d[(I)-1]
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
	xerbla_("SSPTRD", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*n <= 0) {
	return 0;
    }

    if (upper) {

/*        Reduce the upper triangle of A.   
          I1 is the index in AP of A(1,I+1). */

	i1 = *n * (*n - 1) / 2 + 1;
	for (i = *n - 1; i >= 1; --i) {

/*           Generate elementary reflector H(i) = I - tau * v * v'
   
             to annihilate A(1:i-1,i+1) */

	    slarfg_(&i, &AP(i1 + i - 1), &AP(i1), &c__1, &taui);
	    E(i) = AP(i1 + i - 1);

	    if (taui != 0.f) {

/*              Apply H(i) from both sides to A(1:i,1:i) */

		AP(i1 + i - 1) = 1.f;

/*              Compute  y := tau * A * v  storing y in TAU(1:
i) */

		sspmv_(uplo, &i, &taui, &AP(1), &AP(i1), &c__1, &c_b8, &TAU(1)
			, &c__1);

/*              Compute  w := y - 1/2 * tau * (y'*v) * v */

		alpha = taui * -.5f * sdot_(&i, &TAU(1), &c__1, &AP(i1), &
			c__1);
		saxpy_(&i, &alpha, &AP(i1), &c__1, &TAU(1), &c__1);

/*              Apply the transformation as a rank-2 update: 
  
                   A := A - v * w' - w * v' */

		sspr2_(uplo, &i, &c_b14, &AP(i1), &c__1, &TAU(1), &c__1, &AP(
			1));

		AP(i1 + i - 1) = E(i);
	    }
	    D(i + 1) = AP(i1 + i);
	    TAU(i) = taui;
	    i1 -= i;
/* L10: */
	}
	D(1) = AP(1);
    } else {

/*        Reduce the lower triangle of A. II is the index in AP of   
          A(i,i) and I1I1 is the index of A(i+1,i+1). */

	ii = 1;
	i__1 = *n - 1;
	for (i = 1; i <= *n-1; ++i) {
	    i1i1 = ii + *n - i + 1;

/*           Generate elementary reflector H(i) = I - tau * v * v'
   
             to annihilate A(i+2:n,i) */

	    i__2 = *n - i;
	    slarfg_(&i__2, &AP(ii + 1), &AP(ii + 2), &c__1, &taui);
	    E(i) = AP(ii + 1);

	    if (taui != 0.f) {

/*              Apply H(i) from both sides to A(i+1:n,i+1:n) 
*/

		AP(ii + 1) = 1.f;

/*              Compute  y := tau * A * v  storing y in TAU(i:
n-1) */

		i__2 = *n - i;
		sspmv_(uplo, &i__2, &taui, &AP(i1i1), &AP(ii + 1), &c__1, &
			c_b8, &TAU(i), &c__1);

/*              Compute  w := y - 1/2 * tau * (y'*v) * v */

		i__2 = *n - i;
		alpha = taui * -.5f * sdot_(&i__2, &TAU(i), &c__1, &AP(ii + 1)
			, &c__1);
		i__2 = *n - i;
		saxpy_(&i__2, &alpha, &AP(ii + 1), &c__1, &TAU(i), &c__1);

/*              Apply the transformation as a rank-2 update: 
  
                   A := A - v * w' - w * v' */

		i__2 = *n - i;
		sspr2_(uplo, &i__2, &c_b14, &AP(ii + 1), &c__1, &TAU(i), &
			c__1, &AP(i1i1));

		AP(ii + 1) = E(i);
	    }
	    D(i) = AP(ii);
	    TAU(i) = taui;
	    ii = i1i1;
/* L20: */
	}
	D(*n) = AP(ii);
    }

    return 0;

/*     End of SSPTRD */

} /* ssptrd_ */

