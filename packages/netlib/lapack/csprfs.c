#include "f2c.h"

/* Subroutine */ int csprfs_(char *uplo, integer *n, integer *nrhs, complex *
	ap, complex *afp, integer *ipiv, complex *b, integer *ldb, complex *x,
	 integer *ldx, real *ferr, real *berr, complex *work, real *rwork, 
	integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    CSPRFS improves the computed solution to a system of linear   
    equations when the coefficient matrix is symmetric indefinite   
    and packed, and provides error bounds and backward error estimates   
    for the solution.   

    Arguments   
    =========   

    UPLO    (input) CHARACTER*1   
            = 'U':  Upper triangle of A is stored;   
            = 'L':  Lower triangle of A is stored.   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    NRHS    (input) INTEGER   
            The number of right hand sides, i.e., the number of columns   
            of the matrices B and X.  NRHS >= 0.   

    AP      (input) COMPLEX array, dimension (N*(N+1)/2)   
            The upper or lower triangle of the symmetric matrix A, packed 
  
            columnwise in a linear array.  The j-th column of A is stored 
  
            in the array AP as follows:   
            if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;   
            if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n. 
  

    AFP     (input) COMPLEX array, dimension (N*(N+1)/2)   
            The factored form of the matrix A.  AFP contains the block   
            diagonal matrix D and the multipliers used to obtain the   
            factor U or L from the factorization A = U*D*U**T or   
            A = L*D*L**T as computed by CSPTRF, stored as a packed   
            triangular matrix.   

    IPIV    (input) INTEGER array, dimension (N)   
            Details of the interchanges and the block structure of D   
            as determined by CSPTRF.   

    B       (input) COMPLEX array, dimension (LDB,NRHS)   
            The right hand side matrix B.   

    LDB     (input) INTEGER   
            The leading dimension of the array B.  LDB >= max(1,N).   

    X       (input/output) COMPLEX array, dimension (LDX,NRHS)   
            On entry, the solution matrix X, as computed by CSPTRS.   
            On exit, the improved solution matrix X.   

    LDX     (input) INTEGER   
            The leading dimension of the array X.  LDX >= max(1,N).   

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
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   

    Internal Parameters   
    ===================   

    ITMAX is the maximum number of steps of iterative refinement.   

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static complex c_b1 = {1.f,0.f};
    static integer c__1 = 1;
    
    /* System generated locals */
    integer b_dim1, b_offset, x_dim1, x_offset, i__1, i__2, i__3, i__4, i__5;
    real r__1, r__2, r__3, r__4;
    complex q__1;
    /* Builtin functions */
    double r_imag(complex *);
    /* Local variables */
    static integer kase;
    static real safe1, safe2;
    static integer i, j, k;
    static real s;
    extern logical lsame_(char *, char *);
    extern /* Subroutine */ int ccopy_(integer *, complex *, integer *, 
	    complex *, integer *), caxpy_(integer *, complex *, complex *, 
	    integer *, complex *, integer *);
    static integer count;
    extern /* Subroutine */ int cspmv_(char *, integer *, complex *, complex *
	    , complex *, integer *, complex *, complex *, integer *);
    static logical upper;
    static integer ik, kk;
    extern /* Subroutine */ int clacon_(integer *, complex *, complex *, real 
	    *, integer *);
    static real xk;
    extern doublereal slamch_(char *);
    static integer nz;
    static real safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *);
    static real lstres;
    extern /* Subroutine */ int csptrs_(char *, integer *, integer *, complex 
	    *, integer *, complex *, integer *, integer *);
    static real eps;



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
    upper = lsame_(uplo, "U");
    if (! upper && ! lsame_(uplo, "L")) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*nrhs < 0) {
	*info = -3;
    } else if (*ldb < max(1,*n)) {
	*info = -8;
    } else if (*ldx < max(1,*n)) {
	*info = -10;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("CSPRFS", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0 || *nrhs == 0) {
	i__1 = *nrhs;
	for (j = 1; j <= *nrhs; ++j) {
	    FERR(j) = 0.f;
	    BERR(j) = 0.f;
/* L10: */
	}
	return 0;
    }

/*     NZ = maximum number of nonzero elements in each row of A, plus 1 */

    nz = *n + 1;
    eps = slamch_("Epsilon");
    safmin = slamch_("Safe minimum");
    safe1 = nz * safmin;
    safe2 = safe1 / eps;

/*     Do for each right hand side */

    i__1 = *nrhs;
    for (j = 1; j <= *nrhs; ++j) {

	count = 1;
	lstres = 3.f;
L20:

/*        Loop until stopping criterion is satisfied.   

          Compute residual R = B - A * X */

	ccopy_(n, &B(1,j), &c__1, &WORK(1), &c__1);
	q__1.r = -1.f, q__1.i = 0.f;
	cspmv_(uplo, n, &q__1, &AP(1), &X(1,j), &c__1, &c_b1, &
		WORK(1), &c__1);

/*        Compute componentwise relative backward error from formula 
  

          max(i) ( abs(R(i)) / ( abs(A)*abs(X) + abs(B) )(i) )   

          where abs(Z) is the componentwise absolute value of the matr
ix   
          or vector Z.  If the i-th component of the denominator is le
ss   
          than SAFE2, then SAFE1 is added to the i-th components of th
e   
          numerator and denominator before dividing. */

	i__2 = *n;
	for (i = 1; i <= *n; ++i) {
	    i__3 = i + j * b_dim1;
	    RWORK(i) = (r__1 = B(i,j).r, dabs(r__1)) + (r__2 = r_imag(&B(i,j)), dabs(r__2));
/* L30: */
	}

/*        Compute abs(A)*abs(X) + abs(B). */

	kk = 1;
	if (upper) {
	    i__2 = *n;
	    for (k = 1; k <= *n; ++k) {
		s = 0.f;
		i__3 = k + j * x_dim1;
		xk = (r__1 = X(k,j).r, dabs(r__1)) + (r__2 = r_imag(&X(k,j)), dabs(r__2));
		ik = kk;
		i__3 = k - 1;
		for (i = 1; i <= k-1; ++i) {
		    i__4 = ik;
		    RWORK(i) += ((r__1 = AP(ik).r, dabs(r__1)) + (r__2 = 
			    r_imag(&AP(ik)), dabs(r__2))) * xk;
		    i__4 = ik;
		    i__5 = i + j * x_dim1;
		    s += ((r__1 = AP(ik).r, dabs(r__1)) + (r__2 = r_imag(&
			    AP(ik)), dabs(r__2))) * ((r__3 = X(i,j).r, dabs(
			    r__3)) + (r__4 = r_imag(&X(i,j)), dabs(
			    r__4)));
		    ++ik;
/* L40: */
		}
		i__3 = kk + k - 1;
		RWORK(k) = RWORK(k) + ((r__1 = AP(kk+k-1).r, dabs(r__1)) + (
			r__2 = r_imag(&AP(kk + k - 1)), dabs(r__2))) * xk + s;
		kk += k;
/* L50: */
	    }
	} else {
	    i__2 = *n;
	    for (k = 1; k <= *n; ++k) {
		s = 0.f;
		i__3 = k + j * x_dim1;
		xk = (r__1 = X(k,j).r, dabs(r__1)) + (r__2 = r_imag(&X(k,j)), dabs(r__2));
		i__3 = kk;
		RWORK(k) += ((r__1 = AP(kk).r, dabs(r__1)) + (r__2 = r_imag(
			&AP(kk)), dabs(r__2))) * xk;
		ik = kk + 1;
		i__3 = *n;
		for (i = k + 1; i <= *n; ++i) {
		    i__4 = ik;
		    RWORK(i) += ((r__1 = AP(ik).r, dabs(r__1)) + (r__2 = 
			    r_imag(&AP(ik)), dabs(r__2))) * xk;
		    i__4 = ik;
		    i__5 = i + j * x_dim1;
		    s += ((r__1 = AP(ik).r, dabs(r__1)) + (r__2 = r_imag(&
			    AP(ik)), dabs(r__2))) * ((r__3 = X(i,j).r, dabs(
			    r__3)) + (r__4 = r_imag(&X(i,j)), dabs(
			    r__4)));
		    ++ik;
/* L60: */
		}
		RWORK(k) += s;
		kk += *n - k + 1;
/* L70: */
	    }
	}
	s = 0.f;
	i__2 = *n;
	for (i = 1; i <= *n; ++i) {
	    if (RWORK(i) > safe2) {
/* Computing MAX */
		i__3 = i;
		r__3 = s, r__4 = ((r__1 = WORK(i).r, dabs(r__1)) + (r__2 = 
			r_imag(&WORK(i)), dabs(r__2))) / RWORK(i);
		s = dmax(r__3,r__4);
	    } else {
/* Computing MAX */
		i__3 = i;
		r__3 = s, r__4 = ((r__1 = WORK(i).r, dabs(r__1)) + (r__2 = 
			r_imag(&WORK(i)), dabs(r__2)) + safe1) / (RWORK(i) + 
			safe1);
		s = dmax(r__3,r__4);
	    }
/* L80: */
	}
	BERR(j) = s;

/*        Test stopping criterion. Continue iterating if   
             1) The residual BERR(J) is larger than machine epsilon, a
nd   
             2) BERR(J) decreased by at least a factor of 2 during the
   
                last iteration, and   
             3) At most ITMAX iterations tried. */

	if (BERR(j) > eps && BERR(j) * 2.f <= lstres && count <= 5) {

/*           Update solution and try again. */

	    csptrs_(uplo, n, &c__1, &AFP(1), &IPIV(1), &WORK(1), n, info);
	    caxpy_(n, &c_b1, &WORK(1), &c__1, &X(1,j), &c__1);
	    lstres = BERR(j);
	    ++count;
	    goto L20;
	}

/*        Bound error from formula   

          norm(X - XTRUE) / norm(X) .le. FERR =   
          norm( abs(inv(A))*   
             ( abs(R) + NZ*EPS*( abs(A)*abs(X)+abs(B) ))) / norm(X)   

          where   
            norm(Z) is the magnitude of the largest component of Z   
            inv(A) is the inverse of A   
            abs(Z) is the componentwise absolute value of the matrix o
r   
               vector Z   
            NZ is the maximum number of nonzeros in any row of A, plus
 1   
            EPS is machine epsilon   

          The i-th component of abs(R)+NZ*EPS*(abs(A)*abs(X)+abs(B)) 
  
          is incremented by SAFE1 if the i-th component of   
          abs(A)*abs(X) + abs(B) is less than SAFE2.   

          Use CLACON to estimate the infinity-norm of the matrix   
             inv(A) * diag(W),   
          where W = abs(R) + NZ*EPS*( abs(A)*abs(X)+abs(B) ))) */

	i__2 = *n;
	for (i = 1; i <= *n; ++i) {
	    if (RWORK(i) > safe2) {
		i__3 = i;
		RWORK(i) = (r__1 = WORK(i).r, dabs(r__1)) + (r__2 = r_imag(
			&WORK(i)), dabs(r__2)) + nz * eps * RWORK(i);
	    } else {
		i__3 = i;
		RWORK(i) = (r__1 = WORK(i).r, dabs(r__1)) + (r__2 = r_imag(
			&WORK(i)), dabs(r__2)) + nz * eps * RWORK(i) + safe1;
	    }
/* L90: */
	}

	kase = 0;
L100:
	clacon_(n, &WORK(*n + 1), &WORK(1), &FERR(j), &kase);
	if (kase != 0) {
	    if (kase == 1) {

/*              Multiply by diag(W)*inv(A'). */

		csptrs_(uplo, n, &c__1, &AFP(1), &IPIV(1), &WORK(1), n, info);
		i__2 = *n;
		for (i = 1; i <= *n; ++i) {
		    i__3 = i;
		    i__4 = i;
		    i__5 = i;
		    q__1.r = RWORK(i) * WORK(i).r, q__1.i = RWORK(i) 
			    * WORK(i).i;
		    WORK(i).r = q__1.r, WORK(i).i = q__1.i;
/* L110: */
		}
	    } else if (kase == 2) {

/*              Multiply by inv(A)*diag(W). */

		i__2 = *n;
		for (i = 1; i <= *n; ++i) {
		    i__3 = i;
		    i__4 = i;
		    i__5 = i;
		    q__1.r = RWORK(i) * WORK(i).r, q__1.i = RWORK(i) 
			    * WORK(i).i;
		    WORK(i).r = q__1.r, WORK(i).i = q__1.i;
/* L120: */
		}
		csptrs_(uplo, n, &c__1, &AFP(1), &IPIV(1), &WORK(1), n, info);
	    }
	    goto L100;
	}

/*        Normalize error. */

	lstres = 0.f;
	i__2 = *n;
	for (i = 1; i <= *n; ++i) {
/* Computing MAX */
	    i__3 = i + j * x_dim1;
	    r__3 = lstres, r__4 = (r__1 = X(i,j).r, dabs(r__1)) + (r__2 = 
		    r_imag(&X(i,j)), dabs(r__2));
	    lstres = dmax(r__3,r__4);
/* L130: */
	}
	if (lstres != 0.f) {
	    FERR(j) /= lstres;
	}

/* L140: */
    }

    return 0;

/*     End of CSPRFS */

} /* csprfs_ */

