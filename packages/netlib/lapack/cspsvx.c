#include "f2c.h"

/* Subroutine */ int cspsvx_(char *fact, char *uplo, integer *n, integer *
	nrhs, complex *ap, complex *afp, integer *ipiv, complex *b, integer *
	ldb, complex *x, integer *ldx, real *rcond, real *ferr, real *berr, 
	complex *work, real *rwork, integer *info)
{
/*  -- LAPACK driver routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    CSPSVX uses the diagonal pivoting factorization A = U*D*U**T or   
    A = L*D*L**T to compute the solution to a complex system of linear   
    equations A * X = B, where A is an N-by-N symmetric matrix stored   
    in packed format and X and B are N-by-NRHS matrices.   

    Error bounds on the solution and a condition estimate are also   
    provided.   

    Description   
    ===========   

    The following steps are performed:   

    1. If FACT = 'N', the diagonal pivoting method is used to factor A as 
  
          A = U * D * U**T,  if UPLO = 'U', or   
          A = L * D * L**T,  if UPLO = 'L',   
       where U (or L) is a product of permutation and unit upper (lower) 
  
       triangular matrices and D is symmetric and block diagonal with   
       1-by-1 and 2-by-2 diagonal blocks.   

    2. The factored form of A is used to estimate the condition number   
       of the matrix A.  If the reciprocal of the condition number is   
       less than machine precision, steps 3 and 4 are skipped.   

    3. The system of equations is solved for X using the factored form   
       of A.   

    4. Iterative refinement is applied to improve the computed solution   
       matrix and calculate error bounds and backward error estimates   
       for it.   

    Arguments   
    =========   

    FACT    (input) CHARACTER*1   
            Specifies whether or not the factored form of A has been   
            supplied on entry.   
            = 'F':  On entry, AFP and IPIV contain the factored form   
                    of A.  AP, AFP and IPIV will not be modified.   
            = 'N':  The matrix A will be copied to AFP and factored.   

    UPLO    (input) CHARACTER*1   
            = 'U':  Upper triangle of A is stored;   
            = 'L':  Lower triangle of A is stored.   

    N       (input) INTEGER   
            The number of linear equations, i.e., the order of the   
            matrix A.  N >= 0.   

    NRHS    (input) INTEGER   
            The number of right hand sides, i.e., the number of columns   
            of the matrices B and X.  NRHS >= 0.   

    AP      (input) COMPLEX array, dimension (N*(N+1)/2)   
            The upper or lower triangle of the symmetric matrix A, packed 
  
            columnwise in a linear array.  The j-th column of A is stored 
  
            in the array AP as follows:   
            if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;   
            if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n. 
  
            See below for further details.   

    AFP     (input or output) COMPLEX array, dimension (N*(N+1)/2)   
            If FACT = 'F', then AFP is an input argument and on entry   
            contains the block diagonal matrix D and the multipliers used 
  
            to obtain the factor U or L from the factorization   
            A = U*D*U**T or A = L*D*L**T as computed by CSPTRF, stored as 
  
            a packed triangular matrix in the same storage format as A.   

            If FACT = 'N', then AFP is an output argument and on exit   
            contains the block diagonal matrix D and the multipliers used 
  
            to obtain the factor U or L from the factorization   
            A = U*D*U**T or A = L*D*L**T as computed by CSPTRF, stored as 
  
            a packed triangular matrix in the same storage format as A.   

    IPIV    (input or output) INTEGER array, dimension (N)   
            If FACT = 'F', then IPIV is an input argument and on entry   
            contains details of the interchanges and the block structure 
  
            of D, as determined by CSPTRF.   
            If IPIV(k) > 0, then rows and columns k and IPIV(k) were   
            interchanged and D(k,k) is a 1-by-1 diagonal block.   
            If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0, then rows and   
            columns k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k) 
  
            is a 2-by-2 diagonal block.  If UPLO = 'L' and IPIV(k) =   
            IPIV(k+1) < 0, then rows and columns k+1 and -IPIV(k) were   
            interchanged and D(k:k+1,k:k+1) is a 2-by-2 diagonal block.   

            If FACT = 'N', then IPIV is an output argument and on exit   
            contains details of the interchanges and the block structure 
  
            of D, as determined by CSPTRF.   

    B       (input) COMPLEX array, dimension (LDB,NRHS)   
            The N-by-NRHS right hand side matrix B.   

    LDB     (input) INTEGER   
            The leading dimension of the array B.  LDB >= max(1,N).   

    X       (output) COMPLEX array, dimension (LDX,NRHS)   
            If INFO = 0, the N-by-NRHS solution matrix X.   

    LDX     (input) INTEGER   
            The leading dimension of the array X.  LDX >= max(1,N).   

    RCOND   (output) REAL   
            The estimate of the reciprocal condition number of the matrix 
  
            A.  If RCOND is less than the machine precision (in   
            particular, if RCOND = 0), the matrix is singular to working 
  
            precision.  This condition is indicated by a return code of   
            INFO > 0, and the solution and error bounds are not computed. 
  

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
            > 0 and <= N: if INFO = i, D(i,i) is exactly zero.  The   
                 factorization has been completed, but the block diagonal 
  
                 matrix D is exactly singular, so the solution and error 
  
                 bounds could not be computed.   
            = N+1: the block diagonal matrix D is nonsingular, but RCOND 
  
                 is less than machine precision.  The factorization has   
                 been completed, but the matrix is singular to working   
                 precision, so the solution and error bounds have not   
                 been computed.   

    Further Details   
    ===============   

    The packed storage scheme is illustrated by the following example   
    when N = 4, UPLO = 'U':   

    Two-dimensional storage of the symmetric matrix A:   

       a11 a12 a13 a14   
           a22 a23 a24   
               a33 a34     (aij = aji)   
                   a44   

    Packed storage of the upper triangle of A:   

    AP = [ a11, a12, a22, a13, a23, a33, a14, a24, a34, a44 ]   

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__1 = 1;
    
    /* System generated locals */
    integer b_dim1, b_offset, x_dim1, x_offset, i__1;
    /* Local variables */
    extern logical lsame_(char *, char *);
    static real anorm;
    extern /* Subroutine */ int ccopy_(integer *, complex *, integer *, 
	    complex *, integer *);
    extern doublereal slamch_(char *);
    static logical nofact;
    extern /* Subroutine */ int clacpy_(char *, integer *, integer *, complex 
	    *, integer *, complex *, integer *), xerbla_(char *, 
	    integer *);
    extern doublereal clansp_(char *, char *, integer *, complex *, real *);
    extern /* Subroutine */ int cspcon_(char *, integer *, complex *, integer 
	    *, real *, real *, complex *, integer *), csprfs_(char *, 
	    integer *, integer *, complex *, complex *, integer *, complex *, 
	    integer *, complex *, integer *, real *, real *, complex *, real *
	    , integer *), csptrf_(char *, integer *, complex *, 
	    integer *, integer *), csptrs_(char *, integer *, integer 
	    *, complex *, integer *, complex *, integer *, integer *);



#define AP(I) ap[(I)-1]
#define AFP(I) afp[(I)-1]
#define IPIV(I) ipiv[(I)-1]
#define FERR(I) ferr[(I)-1]
#define BERR(I) berr[(I)-1]
#define WORK(I) work[(I)-1]
#define RWORK(I) rwork[(I)-1]

#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]
#define X(I,J) x[(I)-1 + ((J)-1)* ( *ldx)]

    *info = 0;
    nofact = lsame_(fact, "N");
    if (! nofact && ! lsame_(fact, "F")) {
	*info = -1;
    } else if (! lsame_(uplo, "U") && ! lsame_(uplo, "L")) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*nrhs < 0) {
	*info = -4;
    } else if (*ldb < max(1,*n)) {
	*info = -9;
    } else if (*ldx < max(1,*n)) {
	*info = -11;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("CSPSVX", &i__1);
	return 0;
    }

    if (nofact) {

/*        Compute the factorization A = U*D*U' or A = L*D*L'. */

	i__1 = *n * (*n + 1) / 2;
	ccopy_(&i__1, &AP(1), &c__1, &AFP(1), &c__1);
	csptrf_(uplo, n, &AFP(1), &IPIV(1), info);

/*        Return if INFO is non-zero. */

	if (*info != 0) {
	    if (*info > 0) {
		*rcond = 0.f;
	    }
	    return 0;
	}
    }

/*     Compute the norm of the matrix A. */

    anorm = clansp_("I", uplo, n, &AP(1), &RWORK(1));

/*     Compute the reciprocal of the condition number of A. */

    cspcon_(uplo, n, &AFP(1), &IPIV(1), &anorm, rcond, &WORK(1), info);

/*     Return if the matrix is singular to working precision. */

    if (*rcond < slamch_("Epsilon")) {
	*info = *n + 1;
	return 0;
    }

/*     Compute the solution vectors X. */

    clacpy_("Full", n, nrhs, &B(1,1), ldb, &X(1,1), ldx);
    csptrs_(uplo, n, nrhs, &AFP(1), &IPIV(1), &X(1,1), ldx, info);

/*     Use iterative refinement to improve the computed solutions and   
       compute error bounds and backward error estimates for them. */

    csprfs_(uplo, n, nrhs, &AP(1), &AFP(1), &IPIV(1), &B(1,1), ldb, &X(1,1), ldx, &FERR(1), &BERR(1), &WORK(1), &RWORK(1), info)
	    ;

    return 0;

/*     End of CSPSVX */

} /* cspsvx_ */
