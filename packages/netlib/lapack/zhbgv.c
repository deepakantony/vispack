#include "f2c.h"

/* Subroutine */ int zhbgv_(char *jobz, char *uplo, integer *n, integer *ka, 
	integer *kb, doublecomplex *ab, integer *ldab, doublecomplex *bb, 
	integer *ldbb, doublereal *w, doublecomplex *z, integer *ldz, 
	doublecomplex *work, doublereal *rwork, integer *info)
{
/*  -- LAPACK driver routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    ZHBGV computes all the eigenvalues, and optionally, the eigenvectors 
  
    of a complex generalized Hermitian-definite banded eigenproblem, of   
    the form A*x=(lambda)*B*x. Here A and B are assumed to be Hermitian   
    and banded, and B is also positive definite.   

    Arguments   
    =========   

    JOBZ    (input) CHARACTER*1   
            = 'N':  Compute eigenvalues only;   
            = 'V':  Compute eigenvalues and eigenvectors.   

    UPLO    (input) CHARACTER*1   
            = 'U':  Upper triangles of A and B are stored;   
            = 'L':  Lower triangles of A and B are stored.   

    N       (input) INTEGER   
            The order of the matrices A and B.  N >= 0.   

    KA      (input) INTEGER   
            The number of superdiagonals of the matrix A if UPLO = 'U',   
            or the number of subdiagonals if UPLO = 'L'. KA >= 0.   

    KB      (input) INTEGER   
            The number of superdiagonals of the matrix B if UPLO = 'U',   
            or the number of subdiagonals if UPLO = 'L'. KB >= 0.   

    AB      (input/output) COMPLEX*16 array, dimension (LDAB, N)   
            On entry, the upper or lower triangle of the Hermitian band   
            matrix A, stored in the first ka+1 rows of the array.  The   
            j-th column of A is stored in the j-th column of the array AB 
  
            as follows:   
            if UPLO = 'U', AB(ka+1+i-j,j) = A(i,j) for max(1,j-ka)<=i<=j; 
  
            if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+ka). 
  

            On exit, the contents of AB are destroyed.   

    LDAB    (input) INTEGER   
            The leading dimension of the array AB.  LDAB >= KA+1.   

    BB      (input/output) COMPLEX*16 array, dimension (LDBB, N)   
            On entry, the upper or lower triangle of the Hermitian band   
            matrix B, stored in the first kb+1 rows of the array.  The   
            j-th column of B is stored in the j-th column of the array BB 
  
            as follows:   
            if UPLO = 'U', BB(kb+1+i-j,j) = B(i,j) for max(1,j-kb)<=i<=j; 
  
            if UPLO = 'L', BB(1+i-j,j)    = B(i,j) for j<=i<=min(n,j+kb). 
  

            On exit, the factor S from the split Cholesky factorization   
            B = S**H*S, as returned by ZPBSTF.   

    LDBB    (input) INTEGER   
            The leading dimension of the array BB.  LDBB >= KB+1.   

    W       (output) DOUBLE PRECISION array, dimension (N)   
            If INFO = 0, the eigenvalues in ascending order.   

    Z       (output) COMPLEX*16 array, dimension (LDZ, N)   
            If JOBZ = 'V', then if INFO = 0, Z contains the matrix Z of   
            eigenvectors, with the i-th column of Z holding the   
            eigenvector associated with W(i). The eigenvectors are   
            normalized so that Z**H*B*Z = I.   
            If JOBZ = 'N', then Z is not referenced.   

    LDZ     (input) INTEGER   
            The leading dimension of the array Z.  LDZ >= 1, and if   
            JOBZ = 'V', LDZ >= N.   

    WORK    (workspace) COMPLEX*16 array, dimension (N)   

    RWORK   (workspace) DOUBLE PRECISION array, dimension (3*N)   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   
            > 0:  if INFO = i, and i is:   
               <= N:  the algorithm failed to converge:   
                      i off-diagonal elements of an intermediate   
                      tridiagonal form did not converge to zero;   
               > N:   if INFO = N + i, for 1 <= i <= N, then ZPBSTF   
                      returned INFO = i: B is not positive definite.   
                      The factorization of B could not be completed and   
                      no eigenvalues or eigenvectors were computed.   

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    integer ab_dim1, ab_offset, bb_dim1, bb_offset, z_dim1, z_offset, i__1;
    /* Local variables */
    static integer inde;
    static char vect[1];
    extern logical lsame_(char *, char *);
    static integer iinfo;
    static logical upper, wantz;
    extern /* Subroutine */ int xerbla_(char *, integer *), dsterf_(
	    integer *, doublereal *, doublereal *, integer *), zhbtrd_(char *,
	     char *, integer *, integer *, doublecomplex *, integer *, 
	    doublereal *, doublereal *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    static integer indwrk;
    extern /* Subroutine */ int zhbgst_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *,
	     doublecomplex *, integer *, doublecomplex *, doublereal *, 
	    integer *), zpbstf_(char *, integer *, integer *, 
	    doublecomplex *, integer *, integer *), zsteqr_(char *, 
	    integer *, doublereal *, doublereal *, doublecomplex *, integer *,
	     doublereal *, integer *);


#define W(I) w[(I)-1]
#define WORK(I) work[(I)-1]
#define RWORK(I) rwork[(I)-1]

#define AB(I,J) ab[(I)-1 + ((J)-1)* ( *ldab)]
#define BB(I,J) bb[(I)-1 + ((J)-1)* ( *ldbb)]
#define Z(I,J) z[(I)-1 + ((J)-1)* ( *ldz)]

    wantz = lsame_(jobz, "V");
    upper = lsame_(uplo, "U");

    *info = 0;
    if (! (wantz || lsame_(jobz, "N"))) {
	*info = -1;
    } else if (! (upper || lsame_(uplo, "L"))) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*ka < 0) {
	*info = -4;
    } else if (*kb < 0 || *kb > *ka) {
	*info = -5;
    } else if (*ldab < *ka + 1) {
	*info = -7;
    } else if (*ldbb < *kb + 1) {
	*info = -9;
    } else if (*ldz < 1 || wantz && *ldz < *n) {
	*info = -12;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("ZHBGV ", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return 0;
    }

/*     Form a split Cholesky factorization of B. */

    zpbstf_(uplo, n, kb, &BB(1,1), ldbb, info);
    if (*info != 0) {
	*info = *n + *info;
	return 0;
    }

/*     Transform problem to standard eigenvalue problem. */

    inde = 1;
    indwrk = inde + *n;
    zhbgst_(jobz, uplo, n, ka, kb, &AB(1,1), ldab, &BB(1,1), ldbb,
	     &Z(1,1), ldz, &WORK(1), &RWORK(indwrk), &iinfo);

/*     Reduce to tridiagonal form. */

    if (wantz) {
	*(unsigned char *)vect = 'U';
    } else {
	*(unsigned char *)vect = 'N';
    }
    zhbtrd_(vect, uplo, n, ka, &AB(1,1), ldab, &W(1), &RWORK(inde), &Z(1,1), ldz, &WORK(1), &iinfo);

/*     For eigenvalues only, call DSTERF.  For eigenvectors, call ZSTEQR. 
*/

    if (! wantz) {
	dsterf_(n, &W(1), &RWORK(inde), info);
    } else {
	zsteqr_(jobz, n, &W(1), &RWORK(inde), &Z(1,1), ldz, &RWORK(
		indwrk), info);
    }
    return 0;

/*     End of ZHBGV */

} /* zhbgv_ */
