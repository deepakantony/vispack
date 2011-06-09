#include "f2c.h"

/* Subroutine */ int cheevx_(char *jobz, char *range, char *uplo, integer *n, 
	complex *a, integer *lda, real *vl, real *vu, integer *il, integer *
	iu, real *abstol, integer *m, real *w, complex *z, integer *ldz, 
	complex *work, integer *lwork, real *rwork, integer *iwork, integer *
	ifail, integer *info)
{
/*  -- LAPACK driver routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    CHEEVX computes selected eigenvalues and, optionally, eigenvectors   
    of a complex Hermitian matrix A.  Eigenvalues and eigenvectors can   
    be selected by specifying either a range of values or a range of   
    indices for the desired eigenvalues.   

    Arguments   
    =========   

    JOBZ    (input) CHARACTER*1   
            = 'N':  Compute eigenvalues only;   
            = 'V':  Compute eigenvalues and eigenvectors.   

    RANGE   (input) CHARACTER*1   
            = 'A': all eigenvalues will be found.   
            = 'V': all eigenvalues in the half-open interval (VL,VU]   
                   will be found.   
            = 'I': the IL-th through IU-th eigenvalues will be found.   

    UPLO    (input) CHARACTER*1   
            = 'U':  Upper triangle of A is stored;   
            = 'L':  Lower triangle of A is stored.   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    A       (input/output) COMPLEX array, dimension (LDA, N)   
            On entry, the Hermitian matrix A.  If UPLO = 'U', the   
            leading N-by-N upper triangular part of A contains the   
            upper triangular part of the matrix A.  If UPLO = 'L',   
            the leading N-by-N lower triangular part of A contains   
            the lower triangular part of the matrix A.   
            On exit, the lower triangle (if UPLO='L') or the upper   
            triangle (if UPLO='U') of A, including the diagonal, is   
            destroyed.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,N).   

    VL      (input) REAL   
    VU      (input) REAL   
            If RANGE='V', the lower and upper bounds of the interval to   
            be searched for eigenvalues. VL < VU.   
            Not referenced if RANGE = 'A' or 'I'.   

    IL      (input) INTEGER   
    IU      (input) INTEGER   
            If RANGE='I', the indices (in ascending order) of the   
            smallest and largest eigenvalues to be returned.   
            1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.   
            Not referenced if RANGE = 'A' or 'V'.   

    ABSTOL  (input) REAL   
            The absolute error tolerance for the eigenvalues.   
            An approximate eigenvalue is accepted as converged   
            when it is determined to lie in an interval [a,b]   
            of width less than or equal to   

                    ABSTOL + EPS *   max( |a|,|b| ) ,   

            where EPS is the machine precision.  If ABSTOL is less than   
            or equal to zero, then  EPS*|T|  will be used in its place,   
            where |T| is the 1-norm of the tridiagonal matrix obtained   
            by reducing A to tridiagonal form.   

            Eigenvalues will be computed most accurately when ABSTOL is   
            set to twice the underflow threshold 2*SLAMCH('S'), not zero. 
  
            If this routine returns with INFO>0, indicating that some   
            eigenvectors did not converge, try setting ABSTOL to   
            2*SLAMCH('S').   

            See "Computing Small Singular Values of Bidiagonal Matrices   
            with Guaranteed High Relative Accuracy," by Demmel and   
            Kahan, LAPACK Working Note #3.   

    M       (output) INTEGER   
            The total number of eigenvalues found.  0 <= M <= N.   
            If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1.   

    W       (output) REAL array, dimension (N)   
            On normal exit, the first M elements contain the selected   
            eigenvalues in ascending order.   

    Z       (output) COMPLEX array, dimension (LDZ, max(1,M))   
            If JOBZ = 'V', then if INFO = 0, the first M columns of Z   
            contain the orthonormal eigenvectors of the matrix A   
            corresponding to the selected eigenvalues, with the i-th   
            column of Z holding the eigenvector associated with W(i).   
            If an eigenvector fails to converge, then that column of Z   
            contains the latest approximation to the eigenvector, and the 
  
            index of the eigenvector is returned in IFAIL.   
            If JOBZ = 'N', then Z is not referenced.   
            Note: the user must ensure that at least max(1,M) columns are 
  
            supplied in the array Z; if RANGE = 'V', the exact value of M 
  
            is not known in advance and an upper bound must be used.   

    LDZ     (input) INTEGER   
            The leading dimension of the array Z.  LDZ >= 1, and if   
            JOBZ = 'V', LDZ >= max(1,N).   

    WORK    (workspace/output) COMPLEX array, dimension (LWORK)   
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The length of the array WORK.  LWORK >= max(1,2*N-1).   
            For optimal efficiency, LWORK >= (NB+1)*N,   
            where NB is the blocksize for CHETRD returned by ILAENV.   

    RWORK   (workspace) REAL array, dimension (7*N)   

    IWORK   (workspace) INTEGER array, dimension (5*N)   

    IFAIL   (output) INTEGER array, dimension (N)   
            If JOBZ = 'V', then if INFO = 0, the first M elements of   
            IFAIL are zero.  If INFO > 0, then IFAIL contains the   
            indices of the eigenvectors that failed to converge.   
            If JOBZ = 'N', then IFAIL is not referenced.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   
            > 0:  if INFO = i, then i eigenvectors failed to converge.   
                  Their indices are stored in array IFAIL.   

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__1 = 1;
    
    /* System generated locals */
    integer a_dim1, a_offset, z_dim1, z_offset, i__1, i__2;
    real r__1, r__2;
    doublereal d__1;
    complex q__1;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    static integer indd, inde;
    static real anrm;
    static integer imax;
    static real rmin, rmax;
    static integer lopt, itmp1, i, j, indee;
    static real sigma;
    extern logical lsame_(char *, char *);
    static integer iinfo;
    extern /* Subroutine */ int sscal_(integer *, real *, real *, integer *);
    static char order[1];
    extern /* Subroutine */ int cswap_(integer *, complex *, integer *, 
	    complex *, integer *);
    static logical lower;
    extern /* Subroutine */ int scopy_(integer *, real *, integer *, real *, 
	    integer *);
    static logical wantz;
    static integer jj;
    extern doublereal clanhe_(char *, char *, integer *, complex *, integer *,
	     real *);
    static logical alleig, indeig;
    static integer iscale, indibl;
    static logical valeig;
    extern doublereal slamch_(char *);
    extern /* Subroutine */ int chetrd_(char *, integer *, complex *, integer 
	    *, real *, real *, complex *, complex *, integer *, integer *), csscal_(integer *, real *, complex *, integer *), 
	    clacpy_(char *, integer *, integer *, complex *, integer *, 
	    complex *, integer *);
    static real safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *);
    static real abstll, bignum;
    static integer indiwk, indisp, indtau;
    extern /* Subroutine */ int cstein_(integer *, real *, real *, integer *, 
	    real *, integer *, integer *, complex *, integer *, real *, 
	    integer *, integer *, integer *);
    static integer indrwk, indwrk;
    extern /* Subroutine */ int csteqr_(char *, integer *, real *, real *, 
	    complex *, integer *, real *, integer *), cungtr_(char *, 
	    integer *, complex *, integer *, complex *, complex *, integer *, 
	    integer *), ssterf_(integer *, real *, real *, integer *),
	     cunmtr_(char *, char *, char *, integer *, integer *, complex *, 
	    integer *, complex *, complex *, integer *, complex *, integer *, 
	    integer *);
    static integer nsplit, llwork;
    static real smlnum;
    extern /* Subroutine */ int sstebz_(char *, char *, integer *, real *, 
	    real *, integer *, integer *, real *, real *, real *, integer *, 
	    integer *, real *, integer *, integer *, real *, integer *, 
	    integer *);
    static real eps, vll, vuu, tmp1;



#define W(I) w[(I)-1]
#define WORK(I) work[(I)-1]
#define RWORK(I) rwork[(I)-1]
#define IWORK(I) iwork[(I)-1]
#define IFAIL(I) ifail[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]
#define Z(I,J) z[(I)-1 + ((J)-1)* ( *ldz)]

    lower = lsame_(uplo, "L");
    wantz = lsame_(jobz, "V");
    alleig = lsame_(range, "A");
    valeig = lsame_(range, "V");
    indeig = lsame_(range, "I");

    *info = 0;
    if (! (wantz || lsame_(jobz, "N"))) {
	*info = -1;
    } else if (! (alleig || valeig || indeig)) {
	*info = -2;
    } else if (! (lower || lsame_(uplo, "U"))) {
	*info = -3;
    } else if (*n < 0) {
	*info = -4;
    } else if (*lda < max(1,*n)) {
	*info = -6;
    } else if (valeig && *n > 0 && *vu <= *vl) {
	*info = -8;
    } else if (indeig && *il < 1) {
	*info = -9;
    } else if (indeig && (*iu < min(*n,*il) || *iu > *n)) {
	*info = -10;
    } else if (*ldz < 1 || wantz && *ldz < *n) {
	*info = -15;
    } else /* if(complicated condition) */ {
/* Computing MAX */
	i__1 = 1, i__2 = (*n << 1) - 1;
	if (*lwork < max(i__1,i__2)) {
	    *info = -17;
	}
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("CHEEVX", &i__1);
	return 0;
    }

/*     Quick return if possible */

    *m = 0;
    if (*n == 0) {
	WORK(1).r = 1.f, WORK(1).i = 0.f;
	return 0;
    }

    if (*n == 1) {
	WORK(1).r = 1.f, WORK(1).i = 0.f;
	if (alleig || indeig) {
	    *m = 1;
	    i__1 = a_dim1 + 1;
	    W(1) = A(1,1).r;
	} else if (valeig) {
	    i__1 = a_dim1 + 1;
	    i__2 = a_dim1 + 1;
	    if (*vl < A(1,1).r && *vu >= A(1,1).r) {
		*m = 1;
		i__1 = a_dim1 + 1;
		W(1) = A(1,1).r;
	    }
	}
	if (wantz) {
	    i__1 = z_dim1 + 1;
	    Z(1,1).r = 1.f, Z(1,1).i = 0.f;
	}
	return 0;
    }

/*     Get machine constants. */

    safmin = slamch_("Safe minimum");
    eps = slamch_("Precision");
    smlnum = safmin / eps;
    bignum = 1.f / smlnum;
    rmin = sqrt(smlnum);
/* Computing MIN */
    r__1 = sqrt(bignum), r__2 = 1.f / sqrt(sqrt(safmin));
    rmax = dmin(r__1,r__2);

/*     Scale matrix to allowable range, if necessary. */

    iscale = 0;
    abstll = *abstol;
    if (valeig) {
	vll = *vl;
	vuu = *vu;
    }
    anrm = clanhe_("M", uplo, n, &A(1,1), lda, &RWORK(1));
    if (anrm > 0.f && anrm < rmin) {
	iscale = 1;
	sigma = rmin / anrm;
    } else if (anrm > rmax) {
	iscale = 1;
	sigma = rmax / anrm;
    }
    if (iscale == 1) {
	if (lower) {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		i__2 = *n - j + 1;
		csscal_(&i__2, &sigma, &A(j,j), &c__1);
/* L10: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		csscal_(&j, &sigma, &A(1,j), &c__1);
/* L20: */
	    }
	}
	if (*abstol > 0.f) {
	    abstll = *abstol * sigma;
	}
	if (valeig) {
	    vll = *vl * sigma;
	    vuu = *vu * sigma;
	}
    }

/*     Call CHETRD to reduce Hermitian matrix to tridiagonal form. */

    indd = 1;
    inde = indd + *n;
    indrwk = inde + *n;
    indtau = 1;
    indwrk = indtau + *n;
    llwork = *lwork - indwrk + 1;
    chetrd_(uplo, n, &A(1,1), lda, &RWORK(indd), &RWORK(inde), &WORK(
	    indtau), &WORK(indwrk), &llwork, &iinfo);
    i__1 = indwrk;
    q__1.r = *n + WORK(indwrk).r, q__1.i = WORK(indwrk).i;
    lopt = q__1.r;

/*     If all eigenvalues are desired and ABSTOL is less than or equal to 
  
       zero, then call SSTERF or CUNGTR and CSTEQR.  If this fails for   
       some eigenvalue, then try SSTEBZ. */

    if ((alleig || indeig && *il == 1 && *iu == *n) && *abstol <= 0.f) {
	scopy_(n, &RWORK(indd), &c__1, &W(1), &c__1);
	indee = indrwk + (*n << 1);
	if (! wantz) {
	    i__1 = *n - 1;
	    scopy_(&i__1, &RWORK(inde), &c__1, &RWORK(indee), &c__1);
	    ssterf_(n, &W(1), &RWORK(indee), info);
	} else {
	    clacpy_("A", n, n, &A(1,1), lda, &Z(1,1), ldz);
	    cungtr_(uplo, n, &Z(1,1), ldz, &WORK(indtau), &WORK(indwrk), 
		    &llwork, &iinfo);
	    i__1 = *n - 1;
	    scopy_(&i__1, &RWORK(inde), &c__1, &RWORK(indee), &c__1);
	    csteqr_(jobz, n, &W(1), &RWORK(indee), &Z(1,1), ldz, &RWORK(
		    indrwk), info);
	    if (*info == 0) {
		i__1 = *n;
		for (i = 1; i <= *n; ++i) {
		    IFAIL(i) = 0;
/* L30: */
		}
	    }
	}
	if (*info == 0) {
	    *m = *n;
	    goto L40;
	}
	*info = 0;
    }

/*     Otherwise, call SSTEBZ and, if eigenvectors are desired, CSTEIN. */

    if (wantz) {
	*(unsigned char *)order = 'B';
    } else {
	*(unsigned char *)order = 'E';
    }
    indibl = 1;
    indisp = indibl + *n;
    indiwk = indisp + *n;
    sstebz_(range, order, n, &vll, &vuu, il, iu, &abstll, &RWORK(indd), &
	    RWORK(inde), m, &nsplit, &W(1), &IWORK(indibl), &IWORK(indisp), &
	    RWORK(indrwk), &IWORK(indiwk), info);

    if (wantz) {
	cstein_(n, &RWORK(indd), &RWORK(inde), m, &W(1), &IWORK(indibl), &
		IWORK(indisp), &Z(1,1), ldz, &RWORK(indrwk), &IWORK(
		indiwk), &IFAIL(1), info);

/*        Apply unitary matrix used in reduction to tridiagonal   
          form to eigenvectors returned by CSTEIN. */

	cunmtr_("L", uplo, "N", n, m, &A(1,1), lda, &WORK(indtau), &Z(1,1), ldz, &WORK(indwrk), &llwork, &iinfo);
    }

/*     If matrix was scaled, then rescale eigenvalues appropriately. */

L40:
    if (iscale == 1) {
	if (*info == 0) {
	    imax = *m;
	} else {
	    imax = *info - 1;
	}
	r__1 = 1.f / sigma;
	sscal_(&imax, &r__1, &W(1), &c__1);
    }

/*     If eigenvalues are not in order, then sort them, along with   
       eigenvectors. */

    if (wantz) {
	i__1 = *m - 1;
	for (j = 1; j <= *m-1; ++j) {
	    i = 0;
	    tmp1 = W(j);
	    i__2 = *m;
	    for (jj = j + 1; jj <= *m; ++jj) {
		if (W(jj) < tmp1) {
		    i = jj;
		    tmp1 = W(jj);
		}
/* L50: */
	    }

	    if (i != 0) {
		itmp1 = IWORK(indibl + i - 1);
		W(i) = W(j);
		IWORK(indibl + i - 1) = IWORK(indibl + j - 1);
		W(j) = tmp1;
		IWORK(indibl + j - 1) = itmp1;
		cswap_(n, &Z(1,i), &c__1, &Z(1,j), &
			c__1);
		if (*info != 0) {
		    itmp1 = IFAIL(i);
		    IFAIL(i) = IFAIL(j);
		    IFAIL(j) = itmp1;
		}
	    }
/* L60: */
	}
    }

/*     Set WORK(1) to optimal complex workspace size.   

   Computing MAX */
    i__1 = (*n << 1) - 1;
    d__1 = (doublereal) max(i__1,lopt);
    WORK(1).r = d__1, WORK(1).i = 0.f;

    return 0;

/*     End of CHEEVX */

} /* cheevx_ */
