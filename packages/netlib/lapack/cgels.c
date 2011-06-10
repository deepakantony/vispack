#include "f2c.h"

/* Subroutine */ int cgels_(char *trans, integer *m, integer *n, integer *
	nrhs, complex *a, integer *lda, complex *b, integer *ldb, complex *
	work, integer *lwork, integer *info)
{
/*  -- LAPACK driver routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    CGELS solves overdetermined or underdetermined complex linear systems 
  
    involving an M-by-N matrix A, or its conjugate-transpose, using a QR 
  
    or LQ factorization of A.  It is assumed that A has full rank.   

    The following options are provided:   

    1. If TRANS = 'N' and m >= n:  find the least squares solution of   
       an overdetermined system, i.e., solve the least squares problem   
                    minimize || B - A*X ||.   

    2. If TRANS = 'N' and m < n:  find the minimum norm solution of   
       an underdetermined system A * X = B.   

    3. If TRANS = 'C' and m >= n:  find the minimum norm solution of   
       an undetermined system A**H * X = B.   

    4. If TRANS = 'C' and m < n:  find the least squares solution of   
       an overdetermined system, i.e., solve the least squares problem   
                    minimize || B - A**H * X ||.   

    Several right hand side vectors b and solution vectors x can be   
    handled in a single call; they are stored as the columns of the   
    M-by-NRHS right hand side matrix B and the N-by-NRHS solution   
    matrix X.   

    Arguments   
    =========   

    TRANS   (input) CHARACTER   
            = 'N': the linear system involves A;   
            = 'C': the linear system involves A**H.   

    M       (input) INTEGER   
            The number of rows of the matrix A.  M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix A.  N >= 0.   

    NRHS    (input) INTEGER   
            The number of right hand sides, i.e., the number of   
            columns of the matrices B and X. NRHS >= 0.   

    A       (input/output) COMPLEX array, dimension (LDA,N)   
            On entry, the M-by-N matrix A.   
              if M >= N, A is overwritten by details of its QR   
                         factorization as returned by CGEQRF;   
              if M <  N, A is overwritten by details of its LQ   
                         factorization as returned by CGELQF.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,M).   

    B       (input/output) COMPLEX array, dimension (LDB,NRHS)   
            On entry, the matrix B of right hand side vectors, stored   
            columnwise; B is M-by-NRHS if TRANS = 'N', or N-by-NRHS   
            if TRANS = 'C'.   
            On exit, B is overwritten by the solution vectors, stored   
            columnwise:   
            if TRANS = 'N' and m >= n, rows 1 to n of B contain the least 
  
            squares solution vectors; the residual sum of squares for the 
  
            solution in each column is given by the sum of squares of   
            elements N+1 to M in that column;   
            if TRANS = 'N' and m < n, rows 1 to N of B contain the   
            minimum norm solution vectors;   
            if TRANS = 'C' and m >= n, rows 1 to M of B contain the   
            minimum norm solution vectors;   
            if TRANS = 'C' and m < n, rows 1 to M of B contain the   
            least squares solution vectors; the residual sum of squares   
            for the solution in each column is given by the sum of   
            squares of elements M+1 to N in that column.   

    LDB     (input) INTEGER   
            The leading dimension of the array B. LDB >= MAX(1,M,N).   

    WORK    (workspace/output) COMPLEX array, dimension (LWORK)   
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The dimension of the array WORK.   
            LWORK >= min(M,N) + MAX(1,M,N,NRHS).   
            For optimal performance,   
            LWORK >= min(M,N) + MAX(1,M,N,NRHS) * NB   
            where NB is the optimum block size.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   

    ===================================================================== 
  


       Test the input arguments.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static complex c_b1 = {0.f,0.f};
    static complex c_b2 = {1.f,0.f};
    static integer c__1 = 1;
    static integer c_n1 = -1;
    static integer c__0 = 0;
    
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3;
    doublereal d__1;
    /* Local variables */
    static real anrm, bnrm;
    static integer brow;
    static logical tpsd;
    static integer i, j, iascl, ibscl;
    extern logical lsame_(char *, char *);
    extern /* Subroutine */ int ctrsm_(char *, char *, char *, char *, 
	    integer *, integer *, complex *, complex *, integer *, complex *, 
	    integer *);
    static integer wsize;
    static real rwork[1];
    static integer nb;
    extern /* Subroutine */ int slabad_(real *, real *);
    extern doublereal clange_(char *, integer *, integer *, complex *, 
	    integer *, real *);
    static integer mn;
    extern /* Subroutine */ int cgelqf_(integer *, integer *, complex *, 
	    integer *, complex *, complex *, integer *, integer *), clascl_(
	    char *, integer *, integer *, real *, real *, integer *, integer *
	    , complex *, integer *, integer *);
    extern doublereal slamch_(char *);
    extern /* Subroutine */ int cgeqrf_(integer *, integer *, complex *, 
	    integer *, complex *, complex *, integer *, integer *), claset_(
	    char *, integer *, integer *, complex *, complex *, complex *, 
	    integer *), xerbla_(char *, integer *);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static integer scllen;
    static real bignum;
    extern /* Subroutine */ int cunmlq_(char *, char *, integer *, integer *, 
	    integer *, complex *, integer *, complex *, complex *, integer *, 
	    complex *, integer *, integer *), cunmqr_(char *, 
	    char *, integer *, integer *, integer *, complex *, integer *, 
	    complex *, complex *, integer *, complex *, integer *, integer *);
    static real smlnum;



#define WORK(I) work[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]
#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]

    *info = 0;
    mn = min(*m,*n);
    if (! (lsame_(trans, "N") || lsame_(trans, "C"))) {
	*info = -1;
    } else if (*m < 0) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*nrhs < 0) {
	*info = -4;
    } else if (*lda < max(1,*m)) {
	*info = -6;
    } else /* if(complicated condition) */ {
/* Computing MAX */
	i__1 = max(1,*m);
	if (*ldb < max(i__1,*n)) {
	    *info = -8;
	} else /* if(complicated condition) */ {
/* Computing MAX   
   Computing MAX */
	    i__3 = max(*m,*n);
	    i__1 = 1, i__2 = mn + max(i__3,*nrhs);
	    if (*lwork < max(i__1,i__2)) {
		*info = -10;
	    }
	}
    }

/*     Figure out optimal block size */

    if (*info == 0 || *info == -10) {

	tpsd = TRUE_;
	if (lsame_(trans, "N")) {
	    tpsd = FALSE_;
	}

	if (*m >= *n) {
	    nb = ilaenv_(&c__1, "CGEQRF", " ", m, n, &c_n1, &c_n1, 6L, 1L);
	    if (tpsd) {
/* Computing MAX */
		i__1 = nb, i__2 = ilaenv_(&c__1, "CUNMQR", "LN", m, nrhs, n, &
			c_n1, 6L, 2L);
		nb = max(i__1,i__2);
	    } else {
/* Computing MAX */
		i__1 = nb, i__2 = ilaenv_(&c__1, "CUNMQR", "LC", m, nrhs, n, &
			c_n1, 6L, 2L);
		nb = max(i__1,i__2);
	    }
	} else {
	    nb = ilaenv_(&c__1, "CGELQF", " ", m, n, &c_n1, &c_n1, 6L, 1L);
	    if (tpsd) {
/* Computing MAX */
		i__1 = nb, i__2 = ilaenv_(&c__1, "CUNMLQ", "LC", n, nrhs, m, &
			c_n1, 6L, 2L);
		nb = max(i__1,i__2);
	    } else {
/* Computing MAX */
		i__1 = nb, i__2 = ilaenv_(&c__1, "CUNMLQ", "LN", n, nrhs, m, &
			c_n1, 6L, 2L);
		nb = max(i__1,i__2);
	    }
	}

/* Computing MAX */
	i__1 = max(*m,*n);
	wsize = mn + max(i__1,*nrhs) * nb;
	d__1 = (real) wsize;
	WORK(1).r = d__1, WORK(1).i = 0.f;

    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("CGELS ", &i__1);
	return 0;
    }

/*     Quick return if possible   

   Computing MIN */
    i__1 = min(*m,*n);
    if (min(i__1,*nrhs) == 0) {
	i__1 = max(*m,*n);
	claset_("Full", &i__1, nrhs, &c_b1, &c_b1, &B(1,1), ldb);
	return 0;
    }

/*     Get machine parameters */

    smlnum = slamch_("S") / slamch_("P");
    bignum = 1.f / smlnum;
    slabad_(&smlnum, &bignum);

/*     Scale A, B if max element outside range [SMLNUM,BIGNUM] */

    anrm = clange_("M", m, n, &A(1,1), lda, rwork);
    iascl = 0;
    if (anrm > 0.f && anrm < smlnum) {

/*        Scale matrix norm up to SMLNUM */

	clascl_("G", &c__0, &c__0, &anrm, &smlnum, m, n, &A(1,1), lda, 
		info);
	iascl = 1;
    } else if (anrm > bignum) {

/*        Scale matrix norm down to BIGNUM */

	clascl_("G", &c__0, &c__0, &anrm, &bignum, m, n, &A(1,1), lda, 
		info);
	iascl = 2;
    } else if (anrm == 0.f) {

/*        VISMatrix all zero. Return zero solution. */

	i__1 = max(*m,*n);
	claset_("F", &i__1, nrhs, &c_b1, &c_b1, &B(1,1), ldb);
	goto L50;
    }

    brow = *m;
    if (tpsd) {
	brow = *n;
    }
    bnrm = clange_("M", &brow, nrhs, &B(1,1), ldb, rwork);
    ibscl = 0;
    if (bnrm > 0.f && bnrm < smlnum) {

/*        Scale matrix norm up to SMLNUM */

	clascl_("G", &c__0, &c__0, &bnrm, &smlnum, &brow, nrhs, &B(1,1), 
		ldb, info);
	ibscl = 1;
    } else if (bnrm > bignum) {

/*        Scale matrix norm down to BIGNUM */

	clascl_("G", &c__0, &c__0, &bnrm, &bignum, &brow, nrhs, &B(1,1), 
		ldb, info);
	ibscl = 2;
    }

    if (*m >= *n) {

/*        compute QR factorization of A */

	i__1 = *lwork - mn;
	cgeqrf_(m, n, &A(1,1), lda, &WORK(1), &WORK(mn + 1), &i__1, info)
		;

/*        workspace at least N, optimally N*NB */

	if (! tpsd) {

/*           Least-Squares Problem min || A * X - B ||   

             B(1:M,1:NRHS) := Q' * B(1:M,1:NRHS) */

	    i__1 = *lwork - mn;
	    cunmqr_("Left", "Conjugate transpose", m, nrhs, n, &A(1,1), 
		    lda, &WORK(1), &B(1,1), ldb, &WORK(mn + 1), &i__1, 
		    info);

/*           workspace at least NRHS, optimally NRHS*NB   

             B(1:N,1:NRHS) := inv(R) * B(1:N,1:NRHS) */

	    ctrsm_("Left", "Upper", "No transpose", "Non-unit", n, nrhs, &
		    c_b2, &A(1,1), lda, &B(1,1), ldb);

	    scllen = *n;

	} else {

/*           Overdetermined system of equations A' * X = B   

             B(1:N,1:NRHS) := inv(R') * B(1:N,1:NRHS) */

	    ctrsm_("Left", "Upper", "Conjugate transpose", "Non-unit", n, 
		    nrhs, &c_b2, &A(1,1), lda, &B(1,1), ldb);

/*           B(N+1:M,1:NRHS) = ZERO */

	    i__1 = *nrhs;
	    for (j = 1; j <= *nrhs; ++j) {
		i__2 = *m;
		for (i = *n + 1; i <= *m; ++i) {
		    i__3 = i + j * b_dim1;
		    B(i,j).r = 0.f, B(i,j).i = 0.f;
/* L10: */
		}
/* L20: */
	    }

/*           B(1:M,1:NRHS) := Q(1:N,:) * B(1:N,1:NRHS) */

	    i__1 = *lwork - mn;
	    cunmqr_("Left", "No transpose", m, nrhs, n, &A(1,1), lda, &
		    WORK(1), &B(1,1), ldb, &WORK(mn + 1), &i__1, info);

/*           workspace at least NRHS, optimally NRHS*NB */

	    scllen = *m;

	}

    } else {

/*        Compute LQ factorization of A */

	i__1 = *lwork - mn;
	cgelqf_(m, n, &A(1,1), lda, &WORK(1), &WORK(mn + 1), &i__1, info)
		;

/*        workspace at least M, optimally M*NB. */

	if (! tpsd) {

/*           underdetermined system of equations A * X = B   

             B(1:M,1:NRHS) := inv(L) * B(1:M,1:NRHS) */

	    ctrsm_("Left", "Lower", "No transpose", "Non-unit", m, nrhs, &
		    c_b2, &A(1,1), lda, &B(1,1), ldb);

/*           B(M+1:N,1:NRHS) = 0 */

	    i__1 = *nrhs;
	    for (j = 1; j <= *nrhs; ++j) {
		i__2 = *n;
		for (i = *m + 1; i <= *n; ++i) {
		    i__3 = i + j * b_dim1;
		    B(i,j).r = 0.f, B(i,j).i = 0.f;
/* L30: */
		}
/* L40: */
	    }

/*           B(1:N,1:NRHS) := Q(1:N,:)' * B(1:M,1:NRHS) */

	    i__1 = *lwork - mn;
	    cunmlq_("Left", "Conjugate transpose", n, nrhs, m, &A(1,1), 
		    lda, &WORK(1), &B(1,1), ldb, &WORK(mn + 1), &i__1, 
		    info);

/*           workspace at least NRHS, optimally NRHS*NB */

	    scllen = *n;

	} else {

/*           overdetermined system min || A' * X - B ||   

             B(1:N,1:NRHS) := Q * B(1:N,1:NRHS) */

	    i__1 = *lwork - mn;
	    cunmlq_("Left", "No transpose", n, nrhs, m, &A(1,1), lda, &
		    WORK(1), &B(1,1), ldb, &WORK(mn + 1), &i__1, info);

/*           workspace at least NRHS, optimally NRHS*NB   

             B(1:M,1:NRHS) := inv(L') * B(1:M,1:NRHS) */

	    ctrsm_("Left", "Lower", "Conjugate transpose", "Non-unit", m, 
		    nrhs, &c_b2, &A(1,1), lda, &B(1,1), ldb);

	    scllen = *m;

	}

    }

/*     Undo scaling */

    if (iascl == 1) {
	clascl_("G", &c__0, &c__0, &anrm, &smlnum, &scllen, nrhs, &B(1,1)
		, ldb, info);
    } else if (iascl == 2) {
	clascl_("G", &c__0, &c__0, &anrm, &bignum, &scllen, nrhs, &B(1,1)
		, ldb, info);
    }
    if (ibscl == 1) {
	clascl_("G", &c__0, &c__0, &smlnum, &bnrm, &scllen, nrhs, &B(1,1)
		, ldb, info);
    } else if (ibscl == 2) {
	clascl_("G", &c__0, &c__0, &bignum, &bnrm, &scllen, nrhs, &B(1,1)
		, ldb, info);
    }

L50:
    d__1 = (real) wsize;
    WORK(1).r = d__1, WORK(1).i = 0.f;

    return 0;

/*     End of CGELS */

} /* cgels_ */
