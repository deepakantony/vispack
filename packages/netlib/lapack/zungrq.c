#include "f2c.h"

/* Subroutine */ int zungrq_(integer *m, integer *n, integer *k, 
	doublecomplex *a, integer *lda, doublecomplex *tau, doublecomplex *
	work, integer *lwork, integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    ZUNGRQ generates an M-by-N complex matrix Q with orthonormal rows,   
    which is defined as the last M rows of a product of K elementary   
    reflectors of order N   

          Q  =  H(1)' H(2)' . . . H(k)'   

    as returned by ZGERQF.   

    Arguments   
    =========   

    M       (input) INTEGER   
            The number of rows of the matrix Q. M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix Q. N >= M.   

    K       (input) INTEGER   
            The number of elementary reflectors whose product defines the 
  
            matrix Q. M >= K >= 0.   

    A       (input/output) COMPLEX*16 array, dimension (LDA,N)   
            On entry, the (m-k+i)-th row must contain the vector which   
            defines the elementary reflector H(i), for i = 1,2,...,k, as 
  
            returned by ZGERQF in the last k rows of its array argument   
            A.   
            On exit, the M-by-N matrix Q.   

    LDA     (input) INTEGER   
            The first dimension of the array A. LDA >= max(1,M).   

    TAU     (input) COMPLEX*16 array, dimension (K)   
            TAU(i) must contain the scalar factor of the elementary   
            reflector H(i), as returned by ZGERQF.   

    WORK    (workspace/output) COMPLEX*16 array, dimension (LWORK)   
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The dimension of the array WORK. LWORK >= max(1,M).   
            For optimum performance LWORK >= M*NB, where NB is the   
            optimal blocksize.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument has an illegal value   

    ===================================================================== 
  


       Test the input arguments   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__1 = 1;
    static integer c_n1 = -1;
    static integer c__3 = 3;
    static integer c__2 = 2;
    
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;
    /* Local variables */
    static integer i, j, l, nbmin, iinfo, ib, nb, ii, kk;
    extern /* Subroutine */ int zungr2_(integer *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *);
    static integer nx;
    extern /* Subroutine */ int xerbla_(char *, integer *);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int zlarfb_(char *, char *, char *, char *, 
	    integer *, integer *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    static integer ldwork;
    extern /* Subroutine */ int zlarft_(char *, char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *);
    static integer iws;



#define TAU(I) tau[(I)-1]
#define WORK(I) work[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    *info = 0;
    if (*m < 0) {
	*info = -1;
    } else if (*n < *m) {
	*info = -2;
    } else if (*k < 0 || *k > *m) {
	*info = -3;
    } else if (*lda < max(1,*m)) {
	*info = -5;
    } else if (*lwork < max(1,*m)) {
	*info = -8;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("ZUNGRQ", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*m <= 0) {
	WORK(1).r = 1., WORK(1).i = 0.;
	return 0;
    }

/*     Determine the block size. */

    nb = ilaenv_(&c__1, "ZUNGRQ", " ", m, n, k, &c_n1, 6L, 1L);
    nbmin = 2;
    nx = 0;
    iws = *m;
    if (nb > 1 && nb < *k) {

/*        Determine when to cross over from blocked to unblocked code.
   

   Computing MAX */
	i__1 = 0, i__2 = ilaenv_(&c__3, "ZUNGRQ", " ", m, n, k, &c_n1, 6L, 1L)
		;
	nx = max(i__1,i__2);
	if (nx < *k) {

/*           Determine if workspace is large enough for blocked co
de. */

	    ldwork = *m;
	    iws = ldwork * nb;
	    if (*lwork < iws) {

/*              Not enough workspace to use optimal NB:  reduc
e NB and   
                determine the minimum value of NB. */

		nb = *lwork / ldwork;
/* Computing MAX */
		i__1 = 2, i__2 = ilaenv_(&c__2, "ZUNGRQ", " ", m, n, k, &c_n1,
			 6L, 1L);
		nbmin = max(i__1,i__2);
	    }
	}
    }

    if (nb >= nbmin && nb < *k && nx < *k) {

/*        Use blocked code after the first block.   
          The last kk rows are handled by the block method.   

   Computing MIN */
	i__1 = *k, i__2 = (*k - nx + nb - 1) / nb * nb;
	kk = min(i__1,i__2);

/*        Set A(1:m-kk,n-kk+1:n) to zero. */

	i__1 = *n;
	for (j = *n - kk + 1; j <= *n; ++j) {
	    i__2 = *m - kk;
	    for (i = 1; i <= *m-kk; ++i) {
		i__3 = i + j * a_dim1;
		A(i,j).r = 0., A(i,j).i = 0.;
/* L10: */
	    }
/* L20: */
	}
    } else {
	kk = 0;
    }

/*     Use unblocked code for the first or only block. */

    i__1 = *m - kk;
    i__2 = *n - kk;
    i__3 = *k - kk;
    zungr2_(&i__1, &i__2, &i__3, &A(1,1), lda, &TAU(1), &WORK(1), &iinfo)
	    ;

    if (kk > 0) {

/*        Use blocked code */

	i__1 = *k;
	i__2 = nb;
	for (i = *k - kk + 1; nb < 0 ? i >= *k : i <= *k; i += nb) {
/* Computing MIN */
	    i__3 = nb, i__4 = *k - i + 1;
	    ib = min(i__3,i__4);
	    ii = *m - *k + i;
	    if (ii > 1) {

/*              Form the triangular factor of the block reflec
tor   
                H = H(i+ib-1) . . . H(i+1) H(i) */

		i__3 = *n - *k + i + ib - 1;
		zlarft_("Backward", "Rowwise", &i__3, &ib, &A(ii,1), 
			lda, &TAU(i), &WORK(1), &ldwork);

/*              Apply H' to A(1:m-k+i-1,1:n-k+i+ib-1) from the
 right */

		i__3 = ii - 1;
		i__4 = *n - *k + i + ib - 1;
		zlarfb_("Right", "Conjugate transpose", "Backward", "Rowwise",
			 &i__3, &i__4, &ib, &A(ii,1), lda, &WORK(1), &
			ldwork, &A(1,1), lda, &WORK(ib + 1), &ldwork);
	    }

/*           Apply H' to columns 1:n-k+i+ib-1 of current block */

	    i__3 = *n - *k + i + ib - 1;
	    zungr2_(&ib, &i__3, &ib, &A(ii,1), lda, &TAU(i), &WORK(1), 
		    &iinfo);

/*           Set columns n-k+i+ib:n of current block to zero */

	    i__3 = *n;
	    for (l = *n - *k + i + ib; l <= *n; ++l) {
		i__4 = ii + ib - 1;
		for (j = ii; j <= ii+ib-1; ++j) {
		    i__5 = j + l * a_dim1;
		    A(j,l).r = 0., A(j,l).i = 0.;
/* L30: */
		}
/* L40: */
	    }
/* L50: */
	}
    }

    WORK(1).r = (doublereal) iws, WORK(1).i = 0.;
    return 0;

/*     End of ZUNGRQ */

} /* zungrq_ */

