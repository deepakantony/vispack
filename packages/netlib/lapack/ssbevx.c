#include "f2c.h"

/* Subroutine */ int ssbevx_(char *jobz, char *range, char *uplo, integer *n, 
	integer *kd, real *ab, integer *ldab, real *q, integer *ldq, real *vl,
	 real *vu, integer *il, integer *iu, real *abstol, integer *m, real *
	w, real *z, integer *ldz, real *work, integer *iwork, integer *ifail, 
	integer *info)
{
/*  -- LAPACK driver routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    SSBEVX computes selected eigenvalues and, optionally, eigenvectors   
    of a real symmetric band matrix A.  Eigenvalues and eigenvectors can 
  
    be selected by specifying either a range of values or a range of   
    indices for the desired eigenvalues.   

    Arguments   
    =========   

    JOBZ    (input) CHARACTER*1   
            = 'N':  Compute eigenvalues only;   
            = 'V':  Compute eigenvalues and eigenvectors.   

    RANGE   (input) CHARACTER*1   
            = 'A': all eigenvalues will be found;   
            = 'V': all eigenvalues in the half-open interval (VL,VU]   
                   will be found;   
            = 'I': the IL-th through IU-th eigenvalues will be found.   

    UPLO    (input) CHARACTER*1   
            = 'U':  Upper triangle of A is stored;   
            = 'L':  Lower triangle of A is stored.   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    KD      (input) INTEGER   
            The number of superdiagonals of the matrix A if UPLO = 'U',   
            or the number of subdiagonals if UPLO = 'L'.  KD >= 0.   

    AB      (input/output) REAL array, dimension (LDAB, N)   
            On entry, the upper or lower triangle of the symmetric band   
            matrix A, stored in the first KD+1 rows of the array.  The   
            j-th column of A is stored in the j-th column of the array AB 
  
            as follows:   
            if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j; 
  
            if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd). 
  

            On exit, AB is overwritten by values generated during the   
            reduction to tridiagonal form.  If UPLO = 'U', the first   
            superdiagonal and the diagonal of the tridiagonal matrix T   
            are returned in rows KD and KD+1 of AB, and if UPLO = 'L',   
            the diagonal and first subdiagonal of T are returned in the   
            first two rows of AB.   

    LDAB    (input) INTEGER   
            The leading dimension of the array AB.  LDAB >= KD + 1.   

    Q       (output) REAL array, dimension (LDQ, N)   
            If JOBZ = 'V', the N-by-N orthogonal matrix used in the   
                           reduction to tridiagonal form.   
            If JOBZ = 'N', the array Q is not referenced.   

    LDQ     (input) INTEGER   
            The leading dimension of the array Q.  If JOBZ = 'V', then   
            LDQ >= max(1,N).   

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
            by reducing AB to tridiagonal form.   

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
            The first M elements contain the selected eigenvalues in   
            ascending order.   

    Z       (output) REAL array, dimension (LDZ, max(1,M))   
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

    WORK    (workspace) REAL array, dimension (7*N)   

    IWORK   (workspace) INTEGER array, dimension (5*N)   

    IFAIL   (output) INTEGER array, dimension (N)   
            If JOBZ = 'V', then if INFO = 0, the first M elements of   
            IFAIL are zero.  If INFO > 0, then IFAIL contains the   
            indices of the eigenvectors that failed to converge.   
            If JOBZ = 'N', then IFAIL is not referenced.   

    INFO    (output) INTEGER   
            = 0:  successful exit.   
            < 0:  if INFO = -i, the i-th argument had an illegal value.   
            > 0:  if INFO = i, then i eigenvectors failed to converge.   
                  Their indices are stored in array IFAIL.   

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static real c_b14 = 1.f;
    static integer c__1 = 1;
    static real c_b34 = 0.f;
    
    /* System generated locals */
    integer ab_dim1, ab_offset, q_dim1, q_offset, z_dim1, z_offset, i__1, 
	    i__2;
    real r__1, r__2;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    static integer indd, inde;
    static real anrm;
    static integer imax;
    static real rmin, rmax;
    static integer itmp1, i, j, indee;
    static real sigma;
    extern logical lsame_(char *, char *);
    static integer iinfo;
    extern /* Subroutine */ int sscal_(integer *, real *, real *, integer *);
    static char order[1];
    extern /* Subroutine */ int sgemv_(char *, integer *, integer *, real *, 
	    real *, integer *, real *, integer *, real *, real *, integer *);
    static logical lower;
    extern /* Subroutine */ int scopy_(integer *, real *, integer *, real *, 
	    integer *), sswap_(integer *, real *, integer *, real *, integer *
	    );
    static logical wantz;
    static integer jj;
    static logical alleig, indeig;
    static integer iscale, indibl;
    static logical valeig;
    extern doublereal slamch_(char *);
    static real safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *);
    static real abstll, bignum;
    extern doublereal slansb_(char *, char *, integer *, integer *, real *, 
	    integer *, real *);
    extern /* Subroutine */ int slascl_(char *, integer *, integer *, real *, 
	    real *, integer *, integer *, real *, integer *, integer *);
    static integer indisp, indiwo;
    extern /* Subroutine */ int slacpy_(char *, integer *, integer *, real *, 
	    integer *, real *, integer *);
    static integer indwrk;
    extern /* Subroutine */ int ssbtrd_(char *, char *, integer *, integer *, 
	    real *, integer *, real *, real *, real *, integer *, real *, 
	    integer *), sstein_(integer *, real *, real *, 
	    integer *, real *, integer *, integer *, real *, integer *, real *
	    , integer *, integer *, integer *), ssterf_(integer *, real *, 
	    real *, integer *);
    static integer nsplit;
    static real smlnum;
    extern /* Subroutine */ int sstebz_(char *, char *, integer *, real *, 
	    real *, integer *, integer *, real *, real *, real *, integer *, 
	    integer *, real *, integer *, integer *, real *, integer *, 
	    integer *), ssteqr_(char *, integer *, real *, 
	    real *, real *, integer *, real *, integer *);
    static real eps, vll, vuu, tmp1;



#define W(I) w[(I)-1]
#define WORK(I) work[(I)-1]
#define IWORK(I) iwork[(I)-1]
#define IFAIL(I) ifail[(I)-1]

#define AB(I,J) ab[(I)-1 + ((J)-1)* ( *ldab)]
#define Q(I,J) q[(I)-1 + ((J)-1)* ( *ldq)]
#define Z(I,J) z[(I)-1 + ((J)-1)* ( *ldz)]

    wantz = lsame_(jobz, "V");
    alleig = lsame_(range, "A");
    valeig = lsame_(range, "V");
    indeig = lsame_(range, "I");
    lower = lsame_(uplo, "L");

    *info = 0;
    if (! (wantz || lsame_(jobz, "N"))) {
	*info = -1;
    } else if (! (alleig || valeig || indeig)) {
	*info = -2;
    } else if (! (lower || lsame_(uplo, "U"))) {
	*info = -3;
    } else if (*n < 0) {
	*info = -4;
    } else if (*kd < 0) {
	*info = -5;
    } else if (*ldab < *kd + 1) {
	*info = -7;
    } else if (*ldq < *n) {
	*info = -9;
    } else if (valeig && *n > 0 && *vu <= *vl) {
	*info = -11;
    } else if (indeig && *il < 1) {
	*info = -12;
    } else if (indeig && (*iu < min(*n,*il) || *iu > *n)) {
	*info = -13;
    } else if (*ldz < 1 || wantz && *ldz < *n) {
	*info = -18;
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("SSBEVX", &i__1);
	return 0;
    }

/*     Quick return if possible */

    *m = 0;
    if (*n == 0) {
	return 0;
    }

    if (*n == 1) {
	if (alleig || indeig) {
	    *m = 1;
	    W(1) = AB(1,1);
	} else {
	    if (*vl < AB(1,1) && *vu >= AB(1,1)) {
		*m = 1;
		W(1) = AB(1,1);
	    }
	}
	if (wantz) {
	    Z(1,1) = 1.f;
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
    anrm = slansb_("M", uplo, n, kd, &AB(1,1), ldab, &WORK(1));
    if (anrm > 0.f && anrm < rmin) {
	iscale = 1;
	sigma = rmin / anrm;
    } else if (anrm > rmax) {
	iscale = 1;
	sigma = rmax / anrm;
    }
    if (iscale == 1) {
	if (lower) {
	    slascl_("B", kd, kd, &c_b14, &sigma, n, n, &AB(1,1), ldab, 
		    info);
	} else {
	    slascl_("Q", kd, kd, &c_b14, &sigma, n, n, &AB(1,1), ldab, 
		    info);
	}
	if (*abstol > 0.f) {
	    abstll = *abstol * sigma;
	}
	if (valeig) {
	    vll = *vl * sigma;
	    vuu = *vu * sigma;
	}
    }

/*     Call SSBTRD to reduce symmetric band matrix to tridiagonal form. */

    indd = 1;
    inde = indd + *n;
    indwrk = inde + *n;
    ssbtrd_(jobz, uplo, n, kd, &AB(1,1), ldab, &WORK(indd), &WORK(inde),
	     &Q(1,1), ldq, &WORK(indwrk), &iinfo);

/*     If all eigenvalues are desired and ABSTOL is less than or equal   
       to zero, then call SSTERF or SSTEQR.  If this fails for some   
       eigenvalue, then try SSTEBZ. */

    if ((alleig || indeig && *il == 1 && *iu == *n) && *abstol <= 0.f) {
	scopy_(n, &WORK(indd), &c__1, &W(1), &c__1);
	indee = indwrk + (*n << 1);
	if (! wantz) {
	    i__1 = *n - 1;
	    scopy_(&i__1, &WORK(inde), &c__1, &WORK(indee), &c__1);
	    ssterf_(n, &W(1), &WORK(indee), info);
	} else {
	    slacpy_("A", n, n, &Q(1,1), ldq, &Z(1,1), ldz);
	    i__1 = *n - 1;
	    scopy_(&i__1, &WORK(inde), &c__1, &WORK(indee), &c__1);
	    ssteqr_(jobz, n, &W(1), &WORK(indee), &Z(1,1), ldz, &WORK(
		    indwrk), info);
	    if (*info == 0) {
		i__1 = *n;
		for (i = 1; i <= *n; ++i) {
		    IFAIL(i) = 0;
/* L10: */
		}
	    }
	}
	if (*info == 0) {
	    *m = *n;
	    goto L30;
	}
	*info = 0;
    }

/*     Otherwise, call SSTEBZ and, if eigenvectors are desired, SSTEIN. */

    if (wantz) {
	*(unsigned char *)order = 'B';
    } else {
	*(unsigned char *)order = 'E';
    }
    indibl = 1;
    indisp = indibl + *n;
    indiwo = indisp + *n;
    sstebz_(range, order, n, &vll, &vuu, il, iu, &abstll, &WORK(indd), &WORK(
	    inde), m, &nsplit, &W(1), &IWORK(indibl), &IWORK(indisp), &WORK(
	    indwrk), &IWORK(indiwo), info);

    if (wantz) {
	sstein_(n, &WORK(indd), &WORK(inde), m, &W(1), &IWORK(indibl), &IWORK(
		indisp), &Z(1,1), ldz, &WORK(indwrk), &IWORK(indiwo), &
		IFAIL(1), info);

/*        Apply orthogonal matrix used in reduction to tridiagonal   
          form to eigenvectors returned by SSTEIN. */

	i__1 = *m;
	for (j = 1; j <= *m; ++j) {
	    scopy_(n, &Z(1,j), &c__1, &WORK(1), &c__1);
	    sgemv_("N", n, n, &c_b14, &Q(1,1), ldq, &WORK(1), &c__1, &
		    c_b34, &Z(1,j), &c__1);
/* L20: */
	}
    }

/*     If matrix was scaled, then rescale eigenvalues appropriately. */

L30:
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
/* L40: */
	    }

	    if (i != 0) {
		itmp1 = IWORK(indibl + i - 1);
		W(i) = W(j);
		IWORK(indibl + i - 1) = IWORK(indibl + j - 1);
		W(j) = tmp1;
		IWORK(indibl + j - 1) = itmp1;
		sswap_(n, &Z(1,i), &c__1, &Z(1,j), &
			c__1);
		if (*info != 0) {
		    itmp1 = IFAIL(i);
		    IFAIL(i) = IFAIL(j);
		    IFAIL(j) = itmp1;
		}
	    }
/* L50: */
	}
    }

    return 0;

/*     End of SSBEVX */

} /* ssbevx_ */

