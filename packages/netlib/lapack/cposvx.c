#include "f2c.h"

/* Subroutine */ int cposvx_(char *fact, char *uplo, integer *n, integer *
	nrhs, complex *a, integer *lda, complex *af, integer *ldaf, char *
	equed, real *s, complex *b, integer *ldb, complex *x, integer *ldx, 
	real *rcond, real *ferr, real *berr, complex *work, real *rwork, 
	integer *info)
{
/*  -- LAPACK driver routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    CPOSVX uses the Cholesky factorization A = U**H*U or A = L*L**H to   
    compute the solution to a complex system of linear equations   
       A * X = B,   
    where A is an N-by-N Hermitian positive definite matrix and X and B   
    are N-by-NRHS matrices.   

    Error bounds on the solution and a condition estimate are also   
    provided.   

    Description   
    ===========   

    The following steps are performed:   

    1. If FACT = 'E', real scaling factors are computed to equilibrate   
       the system:   
          diag(S) * A * diag(S) * inv(diag(S)) * X = diag(S) * B   
       Whether or not the system will be equilibrated depends on the   
       scaling of the matrix A, but if equilibration is used, A is   
       overwritten by diag(S)*A*diag(S) and B by diag(S)*B.   

    2. If FACT = 'N' or 'E', the Cholesky decomposition is used to   
       factor the matrix A (after equilibration if FACT = 'E') as   
          A = U**H* U,  if UPLO = 'U', or   
          A = L * L**H,  if UPLO = 'L',   
       where U is an upper triangular matrix and L is a lower triangular 
  
       matrix.   

    3. The factored form of A is used to estimate the condition number   
       of the matrix A.  If the reciprocal of the condition number is   
       less than machine precision, steps 4-6 are skipped.   

    4. The system of equations is solved for X using the factored form   
       of A.   

    5. Iterative refinement is applied to improve the computed solution   
       matrix and calculate error bounds and backward error estimates   
       for it.   

    6. If equilibration was used, the matrix X is premultiplied by   
       diag(S) so that it solves the original system before   
       equilibration.   

    Arguments   
    =========   

    FACT    (input) CHARACTER*1   
            Specifies whether or not the factored form of the matrix A is 
  
            supplied on entry, and if not, whether the matrix A should be 
  
            equilibrated before it is factored.   
            = 'F':  On entry, AF contains the factored form of A.   
                    If EQUED = 'Y', the matrix A has been equilibrated   
                    with scaling factors given by S.  A and AF will not   
                    be modified.   
            = 'N':  The matrix A will be copied to AF and factored.   
            = 'E':  The matrix A will be equilibrated if necessary, then 
  
                    copied to AF and factored.   

    UPLO    (input) CHARACTER*1   
            = 'U':  Upper triangle of A is stored;   
            = 'L':  Lower triangle of A is stored.   

    N       (input) INTEGER   
            The number of linear equations, i.e., the order of the   
            matrix A.  N >= 0.   

    NRHS    (input) INTEGER   
            The number of right hand sides, i.e., the number of columns   
            of the matrices B and X.  NRHS >= 0.   

    A       (input/output) COMPLEX array, dimension (LDA,N)   
            On entry, the Hermitian matrix A, except if FACT = 'F' and   
            EQUED = 'Y', then A must contain the equilibrated matrix   
            diag(S)*A*diag(S).  If UPLO = 'U', the leading   
            N-by-N upper triangular part of A contains the upper   
            triangular part of the matrix A, and the strictly lower   
            triangular part of A is not referenced.  If UPLO = 'L', the   
            leading N-by-N lower triangular part of A contains the lower 
  
            triangular part of the matrix A, and the strictly upper   
            triangular part of A is not referenced.  A is not modified if 
  
            FACT = 'F' or 'N', or if FACT = 'E' and EQUED = 'N' on exit. 
  

            On exit, if FACT = 'E' and EQUED = 'Y', A is overwritten by   
            diag(S)*A*diag(S).   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,N).   

    AF      (input or output) COMPLEX array, dimension (LDAF,N)   
            If FACT = 'F', then AF is an input argument and on entry   
            contains the triangular factor U or L from the Cholesky   
            factorization A = U**H*U or A = L*L**H, in the same storage   
            format as A.  If EQUED .ne. 'N', then AF is the factored form 
  
            of the equilibrated matrix diag(S)*A*diag(S).   

            If FACT = 'N', then AF is an output argument and on exit   
            returns the triangular factor U or L from the Cholesky   
            factorization A = U**H*U or A = L*L**H of the original   
            matrix A.   

            If FACT = 'E', then AF is an output argument and on exit   
            returns the triangular factor U or L from the Cholesky   
            factorization A = U**H*U or A = L*L**H of the equilibrated   
            matrix A (see the description of A for the form of the   
            equilibrated matrix).   

    LDAF    (input) INTEGER   
            The leading dimension of the array AF.  LDAF >= max(1,N).   

    EQUED   (input or output) CHARACTER*1   
            Specifies the form of equilibration that was done.   
            = 'N':  No equilibration (always true if FACT = 'N').   
            = 'Y':  Equilibration was done, i.e., A has been replaced by 
  
                    diag(S) * A * diag(S).   
            EQUED is an input argument if FACT = 'F'; otherwise, it is an 
  
            output argument.   

    S       (input or output) REAL array, dimension (N)   
            The scale factors for A; not accessed if EQUED = 'N'.  S is   
            an input argument if FACT = 'F'; otherwise, S is an output   
            argument.  If FACT = 'F' and EQUED = 'Y', each element of S   
            must be positive.   

    B       (input/output) COMPLEX array, dimension (LDB,NRHS)   
            On entry, the N-by-NRHS righthand side matrix B.   
            On exit, if EQUED = 'N', B is not modified; if EQUED = 'Y',   
            B is overwritten by diag(S) * B.   

    LDB     (input) INTEGER   
            The leading dimension of the array B.  LDB >= max(1,N).   

    X       (output) COMPLEX array, dimension (LDX,NRHS)   
            If INFO = 0, the N-by-NRHS solution matrix X to the original 
  
            system of equations.  Note that if EQUED = 'Y', A and B are   
            modified on exit, and the solution to the equilibrated system 
  
            is inv(diag(S))*X.   

    LDX     (input) INTEGER   
            The leading dimension of the array X.  LDX >= max(1,N).   

    RCOND   (output) REAL   
            The estimate of the reciprocal condition number of the matrix 
  
            A after equilibration (if done).  If RCOND is less than the   
            machine precision (in particular, if RCOND = 0), the matrix   
            is singular to working precision.  This condition is   
            indicated by a return code of INFO > 0, and the solution and 
  
            error bounds are not computed.   

    FERR    (output) REAL array, dimension (NRHS)   
            The estimated forward error bound for each solution vector   
            X(j) (the j-th column of the solution matrix X).   
            If XTRUE is the true solution corresponding to X(j), FERR(j) 
  
            is an estimated upper bound for the magnitude of the largest 
  
            element in (X(j) - XTRUE) divided by the magnitude of the   
            largest element in X(j).  The estimate is as reliable as   
            the estimate for RCOND, and is almost always a slight   
            overestimate of the true error.   

    BERR    (output) REAL array, dimension (NRHS)   
            The componentwise relative backward error of each solution   
            vector X(j) (i.e., the smallest relative change in   
            any element of A or B that makes X(j) an exact solution).   

    WORK    (workspace) COMPLEX array, dimension (2*N)   

    RWORK   (workspace) REAL array, dimension (N)   

    INFO    (output) INTEGER   
            = 0: successful exit   
            < 0: if INFO = -i, the i-th argument had an illegal value   
            > 0: if INFO = i, and i is   
                 <= N: the leading minor of order i of A   
                       is not positive definite, so the factorization   
                       could not be completed, and the solution and error 
  
                       bounds could not be computed.   
                 = N+1: RCOND is less than machine precision.  The   
                       factorization has been completed, but the matrix   
                       is singular to working precision, and the solution 
  
                       and error bounds have not been computed.   

    ===================================================================== 
  


    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    integer a_dim1, a_offset, af_dim1, af_offset, b_dim1, b_offset, x_dim1, 
	    x_offset, i__1, i__2, i__3, i__4, i__5;
    real r__1, r__2;
    complex q__1;
    /* Local variables */
    static real amax, smin, smax;
    static integer i, j;
    extern logical lsame_(char *, char *);
    static real scond, anorm;
    static logical equil, rcequ;
    extern doublereal clanhe_(char *, char *, integer *, complex *, integer *,
	     real *);
    extern /* Subroutine */ int claqhe_(char *, integer *, complex *, integer 
	    *, real *, real *, real *, char *);
    extern doublereal slamch_(char *);
    static logical nofact;
    extern /* Subroutine */ int clacpy_(char *, integer *, integer *, complex 
	    *, integer *, complex *, integer *), xerbla_(char *, 
	    integer *);
    static real bignum;
    extern /* Subroutine */ int cpocon_(char *, integer *, complex *, integer 
	    *, real *, real *, complex *, real *, integer *);
    static integer infequ;
    extern /* Subroutine */ int cpoequ_(integer *, complex *, integer *, real 
	    *, real *, real *, integer *), cporfs_(char *, integer *, integer 
	    *, complex *, integer *, complex *, integer *, complex *, integer 
	    *, complex *, integer *, real *, real *, complex *, real *, 
	    integer *), cpotrf_(char *, integer *, complex *, integer 
	    *, integer *), cpotrs_(char *, integer *, integer *, 
	    complex *, integer *, complex *, integer *, integer *);
    static real smlnum;


#define S(I) s[(I)-1]
#define FERR(I) ferr[(I)-1]
#define BERR(I) berr[(I)-1]
#define WORK(I) work[(I)-1]
#define RWORK(I) rwork[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]
#define AF(I,J) af[(I)-1 + ((J)-1)* ( *ldaf)]
#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]
#define X(I,J) x[(I)-1 + ((J)-1)* ( *ldx)]

    *info = 0;
    nofact = lsame_(fact, "N");
    equil = lsame_(fact, "E");
    if (nofact || equil) {
	*(unsigned char *)equed = 'N';
	rcequ = FALSE_;
    } else {
	rcequ = lsame_(equed, "Y");
	smlnum = slamch_("Safe minimum");
	bignum = 1.f / smlnum;
    }

/*     Test the input parameters. */

    if (! nofact && ! equil && ! lsame_(fact, "F")) {
	*info = -1;
    } else if (! lsame_(uplo, "U") && ! lsame_(uplo, "L")) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*nrhs < 0) {
	*info = -4;
    } else if (*lda < max(1,*n)) {
	*info = -6;
    } else if (*ldaf < max(1,*n)) {
	*info = -8;
    } else if (lsame_(fact, "F") && ! (rcequ || lsame_(equed, "N"))) {
	*info = -9;
    } else {
	if (rcequ) {
	    smin = bignum;
	    smax = 0.f;
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
/* Computing MIN */
		r__1 = smin, r__2 = S(j);
		smin = dmin(r__1,r__2);
/* Computing MAX */
		r__1 = smax, r__2 = S(j);
		smax = dmax(r__1,r__2);
/* L10: */
	    }
	    if (smin <= 0.f) {
		*info = -10;
	    } else if (*n > 0) {
		scond = dmax(smin,smlnum) / dmin(smax,bignum);
	    } else {
		scond = 1.f;
	    }
	}
	if (*info == 0) {
	    if (*ldb < max(1,*n)) {
		*info = -12;
	    } else if (*ldx < max(1,*n)) {
		*info = -14;
	    }
	}
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("CPOSVX", &i__1);
	return 0;
    }

    if (equil) {

/*        Compute row and column scalings to equilibrate the matrix A.
 */

	cpoequ_(n, &A(1,1), lda, &S(1), &scond, &amax, &infequ);
	if (infequ == 0) {

/*           Equilibrate the matrix. */

	    claqhe_(uplo, n, &A(1,1), lda, &S(1), &scond, &amax, equed);
	    rcequ = lsame_(equed, "Y");
	}
    }

/*     Scale the right hand side. */

    if (rcequ) {
	i__1 = *nrhs;
	for (j = 1; j <= *nrhs; ++j) {
	    i__2 = *n;
	    for (i = 1; i <= *n; ++i) {
		i__3 = i + j * b_dim1;
		i__4 = i;
		i__5 = i + j * b_dim1;
		q__1.r = S(i) * B(i,j).r, q__1.i = S(i) * B(i,j).i;
		B(i,j).r = q__1.r, B(i,j).i = q__1.i;
/* L20: */
	    }
/* L30: */
	}
    }

    if (nofact || equil) {

/*        Compute the Cholesky factorization A = U'*U or A = L*L'. */

	clacpy_(uplo, n, n, &A(1,1), lda, &AF(1,1), ldaf);
	cpotrf_(uplo, n, &AF(1,1), ldaf, info);

/*        Return if INFO is non-zero. */

	if (*info != 0) {
	    if (*info > 0) {
		*rcond = 0.f;
	    }
	    return 0;
	}
    }

/*     Compute the norm of the matrix A. */

    anorm = clanhe_("1", uplo, n, &A(1,1), lda, &RWORK(1));

/*     Compute the reciprocal of the condition number of A. */

    cpocon_(uplo, n, &AF(1,1), ldaf, &anorm, rcond, &WORK(1), &RWORK(1),
	     info);

/*     Return if the matrix is singular to working precision. */

    if (*rcond < slamch_("Epsilon")) {
	*info = *n + 1;
	return 0;
    }

/*     Compute the solution matrix X. */

    clacpy_("Full", n, nrhs, &B(1,1), ldb, &X(1,1), ldx);
    cpotrs_(uplo, n, nrhs, &AF(1,1), ldaf, &X(1,1), ldx, info);

/*     Use iterative refinement to improve the computed solution and   
       compute error bounds and backward error estimates for it. */

    cporfs_(uplo, n, nrhs, &A(1,1), lda, &AF(1,1), ldaf, &B(1,1), ldb, &X(1,1), ldx, &FERR(1), &BERR(1), &WORK(1), &
	    RWORK(1), info);

/*     Transform the solution matrix X to a solution of the original   
       system. */

    if (rcequ) {
	i__1 = *nrhs;
	for (j = 1; j <= *nrhs; ++j) {
	    i__2 = *n;
	    for (i = 1; i <= *n; ++i) {
		i__3 = i + j * x_dim1;
		i__4 = i;
		i__5 = i + j * x_dim1;
		q__1.r = S(i) * X(i,j).r, q__1.i = S(i) * X(i,j).i;
		X(i,j).r = q__1.r, X(i,j).i = q__1.i;
/* L40: */
	    }
/* L50: */
	}
	i__1 = *nrhs;
	for (j = 1; j <= *nrhs; ++j) {
	    FERR(j) /= scond;
/* L60: */
	}
    }

    return 0;

/*     End of CPOSVX */

} /* cposvx_ */

