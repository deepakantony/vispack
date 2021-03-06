#include "f2c.h"

/* Subroutine */ int sgeqpf_(integer *m, integer *n, real *a, integer *lda, 
	integer *jpvt, real *tau, real *work, integer *info)
{
/*  -- LAPACK test routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       March 31, 1993   


    Purpose   
    =======   

    SGEQPF computes a QR factorization with column pivoting of a   
    real M-by-N matrix A: A*P = Q*R.   

    Arguments   
    =========   

    M       (input) INTEGER   
            The number of rows of the matrix A. M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix A. N >= 0   

    A       (input/output) REAL array, dimension (LDA,N)   
            On entry, the M-by-N matrix A.   
            On exit, the upper triangle of the array contains the   
            min(M,N)-by-N upper triangular matrix R; the elements   
            below the diagonal, together with the array TAU,   
            represent the orthogonal matrix Q as a product of   
            min(m,n) elementary reflectors.   

    LDA     (input) INTEGER   
            The leading dimension of the array A. LDA >= max(1,M).   

    JPVT    (input/output) INTEGER array, dimension (N)   
            On entry, if JPVT(i) .ne. 0, the i-th column of A is permuted 
  
            to the front of A*P (a leading column); if JPVT(i) = 0,   
            the i-th column of A is a free column.   
            On exit, if JPVT(i) = k, then the i-th column of A*P   
            was the k-th column of A.   

    TAU     (output) REAL array, dimension (min(M,N))   
            The scalar factors of the elementary reflectors.   

    WORK    (workspace) REAL array, dimension (3*N)   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   

    Further Details   
    ===============   

    The matrix Q is represented as a product of elementary reflectors   

       Q = H(1) H(2) . . . H(n)   

    Each H(i) has the form   

       H = I - tau * v * v'   

    where tau is a real scalar, and v is a real vector with   
    v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i). 
  

    The matrix P is represented in jpvt as follows: If   
       jpvt(j) = i   
    then the jth column of P is the ith canonical unit vector.   

    ===================================================================== 
  


       Test the input arguments   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__1 = 1;
    
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    real r__1, r__2;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    static real temp, temp2;
    extern doublereal snrm2_(integer *, real *, integer *);
    static integer i, j;
    extern /* Subroutine */ int slarf_(char *, integer *, integer *, real *, 
	    integer *, real *, real *, integer *, real *);
    static integer itemp;
    extern /* Subroutine */ int sswap_(integer *, real *, integer *, real *, 
	    integer *), sgeqr2_(integer *, integer *, real *, integer *, real 
	    *, real *, integer *);
    static integer ma;
    extern /* Subroutine */ int sorm2r_(char *, char *, integer *, integer *, 
	    integer *, real *, integer *, real *, real *, integer *, real *, 
	    integer *);
    static integer mn;
    extern /* Subroutine */ int xerbla_(char *, integer *), slarfg_(
	    integer *, real *, real *, integer *, real *);
    extern integer isamax_(integer *, real *, integer *);
    static real aii;
    static integer pvt;



#define JPVT(I) jpvt[(I)-1]
#define TAU(I) tau[(I)-1]
#define WORK(I) work[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    *info = 0;
    if (*m < 0) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < max(1,*m)) {
	*info = -4;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("SGEQPF", &i__1);
	return 0;
    }

    mn = min(*m,*n);

/*     Move initial columns up front */

    itemp = 1;
    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	if (JPVT(i) != 0) {
	    if (i != itemp) {
		sswap_(m, &A(1,i), &c__1, &A(1,itemp), &
			c__1);
		JPVT(i) = JPVT(itemp);
		JPVT(itemp) = i;
	    } else {
		JPVT(i) = i;
	    }
	    ++itemp;
	} else {
	    JPVT(i) = i;
	}
/* L10: */
    }
    --itemp;

/*     Compute the QR factorization and update remaining columns */

    if (itemp > 0) {
	ma = min(itemp,*m);
	sgeqr2_(m, &ma, &A(1,1), lda, &TAU(1), &WORK(1), info);
	if (ma < *n) {
	    i__1 = *n - ma;
	    sorm2r_("Left", "Transpose", m, &i__1, &ma, &A(1,1), lda, &
		    TAU(1), &A(1,ma+1), lda, &WORK(1), info);
	}
    }

    if (itemp < mn) {

/*        Initialize partial column norms. The first n elements of   
          work store the exact column norms. */

	i__1 = *n;
	for (i = itemp + 1; i <= *n; ++i) {
	    i__2 = *m - itemp;
	    WORK(i) = snrm2_(&i__2, &A(itemp+1,i), &c__1);
	    WORK(*n + i) = WORK(i);
/* L20: */
	}

/*        Compute factorization */

	i__1 = mn;
	for (i = itemp + 1; i <= mn; ++i) {

/*           Determine ith pivot column and swap if necessary */

	    i__2 = *n - i + 1;
	    pvt = i - 1 + isamax_(&i__2, &WORK(i), &c__1);

	    if (pvt != i) {
		sswap_(m, &A(1,pvt), &c__1, &A(1,i), &
			c__1);
		itemp = JPVT(pvt);
		JPVT(pvt) = JPVT(i);
		JPVT(i) = itemp;
		WORK(pvt) = WORK(i);
		WORK(*n + pvt) = WORK(*n + i);
	    }

/*           Generate elementary reflector H(i) */

	    if (i < *m) {
		i__2 = *m - i + 1;
		slarfg_(&i__2, &A(i,i), &A(i+1,i), &
			c__1, &TAU(i));
	    } else {
		slarfg_(&c__1, &A(*m,*m), &A(*m,*m), &
			c__1, &TAU(*m));
	    }

	    if (i < *n) {

/*              Apply H(i) to A(i:m,i+1:n) from the left */

		aii = A(i,i);
		A(i,i) = 1.f;
		i__2 = *m - i + 1;
		i__3 = *n - i;
		slarf_("LEFT", &i__2, &i__3, &A(i,i), &c__1, &TAU(
			i), &A(i,i+1), lda, &WORK((*n << 1) + 
			1));
		A(i,i) = aii;
	    }

/*           Update partial column norms */

	    i__2 = *n;
	    for (j = i + 1; j <= *n; ++j) {
		if (WORK(j) != 0.f) {
/* Computing 2nd power */
		    r__2 = (r__1 = A(i,j), dabs(r__1)) / WORK(j);
		    temp = 1.f - r__2 * r__2;
		    temp = dmax(temp,0.f);
/* Computing 2nd power */
		    r__1 = WORK(j) / WORK(*n + j);
		    temp2 = temp * .05f * (r__1 * r__1) + 1.f;
		    if (temp2 == 1.f) {
			if (*m - i > 0) {
			    i__3 = *m - i;
			    WORK(j) = snrm2_(&i__3, &A(i+1,j), &
				    c__1);
			    WORK(*n + j) = WORK(j);
			} else {
			    WORK(j) = 0.f;
			    WORK(*n + j) = 0.f;
			}
		    } else {
			WORK(j) *= sqrt(temp);
		    }
		}
/* L30: */
	    }

/* L40: */
	}
    }
    return 0;

/*     End of SGEQPF */

} /* sgeqpf_ */

