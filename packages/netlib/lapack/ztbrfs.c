#include "f2c.h"

/* Subroutine */ int ztbrfs_(char *uplo, char *trans, char *diag, integer *n, 
	integer *kd, integer *nrhs, doublecomplex *ab, integer *ldab, 
	doublecomplex *b, integer *ldb, doublecomplex *x, integer *ldx, 
	doublereal *ferr, doublereal *berr, doublecomplex *work, doublereal *
	rwork, integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    ZTBRFS provides error bounds and backward error estimates for the   
    solution to a system of linear equations with a triangular band   
    coefficient matrix.   

    The solution matrix X must be computed by ZTBTRS or some other   
    means before entering this routine.  ZTBRFS does not do iterative   
    refinement because doing so cannot improve the backward error.   

    Arguments   
    =========   

    UPLO    (input) CHARACTER*1   
            = 'U':  A is upper triangular;   
            = 'L':  A is lower triangular.   

    TRANS   (input) CHARACTER*1   
            Specifies the form of the system of equations:   
            = 'N':  A * X = B     (No transpose)   
            = 'T':  A**T * X = B  (Transpose)   
            = 'C':  A**H * X = B  (Conjugate transpose)   

    DIAG    (input) CHARACTER*1   
            = 'N':  A is non-unit triangular;   
            = 'U':  A is unit triangular.   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    KD      (input) INTEGER   
            The number of superdiagonals or subdiagonals of the   
            triangular band matrix A.  KD >= 0.   

    NRHS    (input) INTEGER   
            The number of right hand sides, i.e., the number of columns   
            of the matrices B and X.  NRHS >= 0.   

    AB      (input) COMPLEX*16 array, dimension (LDAB,N)   
            The upper or lower triangular band matrix A, stored in the   
            first kd+1 rows of the array. The j-th column of A is stored 
  
            in the j-th column of the array AB as follows:   
            if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j; 
  
            if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd). 
  
            If DIAG = 'U', the diagonal elements of A are not referenced 
  
            and are assumed to be 1.   

    LDAB    (input) INTEGER   
            The leading dimension of the array AB.  LDAB >= KD+1.   

    B       (input) COMPLEX*16 array, dimension (LDB,NRHS)   
            The right hand side matrix B.   

    LDB     (input) INTEGER   
            The leading dimension of the array B.  LDB >= max(1,N).   

    X       (input) COMPLEX*16 array, dimension (LDX,NRHS)   
            The solution matrix X.   

    LDX     (input) INTEGER   
            The leading dimension of the array X.  LDX >= max(1,N).   

    FERR    (output) DOUBLE PRECISION array, dimension (NRHS)   
            The estimated forward error bound for each solution vector   
            X(j) (the j-th column of the solution matrix X).   
            If XTRUE is the true solution corresponding to X(j), FERR(j) 
  
            is an estimated upper bound for the magnitude of the largest 
  
            element in (X(j) - XTRUE) divided by the magnitude of the   
            largest element in X(j).  The estimate is as reliable as   
            the estimate for RCOND, and is almost always a slight   
            overestimate of the true error.   

    BERR    (output) DOUBLE PRECISION array, dimension (NRHS)   
            The componentwise relative backward error of each solution   
            vector X(j) (i.e., the smallest relative change in   
            any element of A or B that makes X(j) an exact solution).   

    WORK    (workspace) COMPLEX*16 array, dimension (2*N)   

    RWORK   (workspace) DOUBLE PRECISION array, dimension (N)   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__1 = 1;
    
    /* System generated locals */
    integer ab_dim1, ab_offset, b_dim1, b_offset, x_dim1, x_offset, i__1, 
	    i__2, i__3, i__4, i__5;
    doublereal d__1, d__2, d__3, d__4;
    doublecomplex z__1;
    /* Builtin functions */
    double d_imag(doublecomplex *);
    /* Local variables */
    static integer kase;
    static doublereal safe1, safe2;
    static integer i, j, k;
    static doublereal s;
    extern logical lsame_(char *, char *);
    static logical upper;
    extern /* Subroutine */ int ztbmv_(char *, char *, char *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *), zcopy_(integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *), ztbsv_(char *, char *, 
	    char *, integer *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), zaxpy_(
	    integer *, doublecomplex *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    extern doublereal dlamch_(char *);
    static doublereal xk;
    static integer nz;
    static doublereal safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *), zlacon_(
	    integer *, doublecomplex *, doublecomplex *, doublereal *, 
	    integer *);
    static logical notran;
    static char transn[1], transt[1];
    static logical nounit;
    static doublereal lstres, eps;



#define FERR(I) ferr[(I)-1]
#define BERR(I) berr[(I)-1]
#define WORK(I) work[(I)-1]
#define RWORK(I) rwork[(I)-1]

#define AB(I,J) ab[(I)-1 + ((J)-1)* ( *ldab)]
#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]
#define X(I,J) x[(I)-1 + ((J)-1)* ( *ldx)]

    *info = 0;
    upper = lsame_(uplo, "U");
    notran = lsame_(trans, "N");
    nounit = lsame_(diag, "N");

    if (! upper && ! lsame_(uplo, "L")) {
	*info = -1;
    } else if (! notran && ! lsame_(trans, "T") && ! lsame_(trans, 
	    "C")) {
	*info = -2;
    } else if (! nounit && ! lsame_(diag, "U")) {
	*info = -3;
    } else if (*n < 0) {
	*info = -4;
    } else if (*kd < 0) {
	*info = -5;
    } else if (*nrhs < 0) {
	*info = -6;
    } else if (*ldab < *kd + 1) {
	*info = -8;
    } else if (*ldb < max(1,*n)) {
	*info = -10;
    } else if (*ldx < max(1,*n)) {
	*info = -12;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("ZTBRFS", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0 || *nrhs == 0) {
	i__1 = *nrhs;
	for (j = 1; j <= *nrhs; ++j) {
	    FERR(j) = 0.;
	    BERR(j) = 0.;
/* L10: */
	}
	return 0;
    }

    if (notran) {
	*(unsigned char *)transn = 'N';
	*(unsigned char *)transt = 'C';
    } else {
	*(unsigned char *)transn = 'C';
	*(unsigned char *)transt = 'N';
    }

/*     NZ = maximum number of nonzero elements in each row of A, plus 1 */

    nz = *kd + 2;
    eps = dlamch_("Epsilon");
    safmin = dlamch_("Safe minimum");
    safe1 = nz * safmin;
    safe2 = safe1 / eps;

/*     Do for each right hand side */

    i__1 = *nrhs;
    for (j = 1; j <= *nrhs; ++j) {

/*        Compute residual R = B - op(A) * X,   
          where op(A) = A, A**T, or A**H, depending on TRANS. */

	zcopy_(n, &X(1,j), &c__1, &WORK(1), &c__1);
	ztbmv_(uplo, trans, diag, n, kd, &AB(1,1), ldab, &WORK(1), &
		c__1);
	z__1.r = -1., z__1.i = 0.;
	zaxpy_(n, &z__1, &B(1,j), &c__1, &WORK(1), &c__1);

/*        Compute componentwise relative backward error from formula 
  

          max(i) ( abs(R(i)) / ( abs(op(A))*abs(X) + abs(B) )(i) )   

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
	    RWORK(i) = (d__1 = B(i,j).r, abs(d__1)) + (d__2 = d_imag(&B(i,j)), abs(d__2));
/* L20: */
	}

	if (notran) {

/*           Compute abs(A)*abs(X) + abs(B). */

	    if (upper) {
		if (nounit) {
		    i__2 = *n;
		    for (k = 1; k <= *n; ++k) {
			i__3 = k + j * x_dim1;
			xk = (d__1 = X(k,j).r, abs(d__1)) + (d__2 = d_imag(&
				X(k,j)), abs(d__2));
/* Computing MAX */
			i__3 = 1, i__4 = k - *kd;
			i__5 = k;
			for (i = max(1,k-*kd); i <= k; ++i) {
			    i__3 = *kd + 1 + i - k + k * ab_dim1;
			    RWORK(i) += ((d__1 = AB(*kd+1+i-k,k).r, abs(d__1)) + (
				    d__2 = d_imag(&AB(*kd+1+i-k,k)), abs(d__2))) * xk;
/* L30: */
			}
/* L40: */
		    }
		} else {
		    i__2 = *n;
		    for (k = 1; k <= *n; ++k) {
			i__5 = k + j * x_dim1;
			xk = (d__1 = X(k,j).r, abs(d__1)) + (d__2 = d_imag(&
				X(k,j)), abs(d__2));
/* Computing MAX */
			i__5 = 1, i__3 = k - *kd;
			i__4 = k - 1;
			for (i = max(1,k-*kd); i <= k-1; ++i) {
			    i__5 = *kd + 1 + i - k + k * ab_dim1;
			    RWORK(i) += ((d__1 = AB(*kd+1+i-k,k).r, abs(d__1)) + (
				    d__2 = d_imag(&AB(*kd+1+i-k,k)), abs(d__2))) * xk;
/* L50: */
			}
			RWORK(k) += xk;
/* L60: */
		    }
		}
	    } else {
		if (nounit) {
		    i__2 = *n;
		    for (k = 1; k <= *n; ++k) {
			i__4 = k + j * x_dim1;
			xk = (d__1 = X(k,j).r, abs(d__1)) + (d__2 = d_imag(&
				X(k,j)), abs(d__2));
/* Computing MIN */
			i__5 = *n, i__3 = k + *kd;
			i__4 = min(i__5,i__3);
			for (i = k; i <= min(*n,k+*kd); ++i) {
			    i__5 = i + 1 - k + k * ab_dim1;
			    RWORK(i) += ((d__1 = AB(i+1-k,k).r, abs(d__1)) + (
				    d__2 = d_imag(&AB(i+1-k,k)
				    ), abs(d__2))) * xk;
/* L70: */
			}
/* L80: */
		    }
		} else {
		    i__2 = *n;
		    for (k = 1; k <= *n; ++k) {
			i__4 = k + j * x_dim1;
			xk = (d__1 = X(k,j).r, abs(d__1)) + (d__2 = d_imag(&
				X(k,j)), abs(d__2));
/* Computing MIN */
			i__5 = *n, i__3 = k + *kd;
			i__4 = min(i__5,i__3);
			for (i = k + 1; i <= min(*n,k+*kd); ++i) {
			    i__5 = i + 1 - k + k * ab_dim1;
			    RWORK(i) += ((d__1 = AB(i+1-k,k).r, abs(d__1)) + (
				    d__2 = d_imag(&AB(i+1-k,k)
				    ), abs(d__2))) * xk;
/* L90: */
			}
			RWORK(k) += xk;
/* L100: */
		    }
		}
	    }
	} else {

/*           Compute abs(A**H)*abs(X) + abs(B). */

	    if (upper) {
		if (nounit) {
		    i__2 = *n;
		    for (k = 1; k <= *n; ++k) {
			s = 0.;
/* Computing MAX */
			i__4 = 1, i__5 = k - *kd;
			i__3 = k;
			for (i = max(1,k-*kd); i <= k; ++i) {
			    i__4 = *kd + 1 + i - k + k * ab_dim1;
			    i__5 = i + j * x_dim1;
			    s += ((d__1 = AB(*kd+1+i-k,k).r, abs(d__1)) + (d__2 = 
				    d_imag(&AB(*kd+1+i-k,k))
				    , abs(d__2))) * ((d__3 = X(i,j).r, abs(
				    d__3)) + (d__4 = d_imag(&X(i,j)
				    ), abs(d__4)));
/* L110: */
			}
			RWORK(k) += s;
/* L120: */
		    }
		} else {
		    i__2 = *n;
		    for (k = 1; k <= *n; ++k) {
			i__3 = k + j * x_dim1;
			s = (d__1 = X(k,j).r, abs(d__1)) + (d__2 = d_imag(&X(k,j)), abs(d__2));
/* Computing MAX */
			i__3 = 1, i__4 = k - *kd;
			i__5 = k - 1;
			for (i = max(1,k-*kd); i <= k-1; ++i) {
			    i__3 = *kd + 1 + i - k + k * ab_dim1;
			    i__4 = i + j * x_dim1;
			    s += ((d__1 = AB(*kd+1+i-k,k).r, abs(d__1)) + (d__2 = 
				    d_imag(&AB(*kd+1+i-k,k))
				    , abs(d__2))) * ((d__3 = X(i,j).r, abs(
				    d__3)) + (d__4 = d_imag(&X(i,j)
				    ), abs(d__4)));
/* L130: */
			}
			RWORK(k) += s;
/* L140: */
		    }
		}
	    } else {
		if (nounit) {
		    i__2 = *n;
		    for (k = 1; k <= *n; ++k) {
			s = 0.;
/* Computing MIN */
			i__3 = *n, i__4 = k + *kd;
			i__5 = min(i__3,i__4);
			for (i = k; i <= min(*n,k+*kd); ++i) {
			    i__3 = i + 1 - k + k * ab_dim1;
			    i__4 = i + j * x_dim1;
			    s += ((d__1 = AB(i+1-k,k).r, abs(d__1)) + (d__2 = 
				    d_imag(&AB(i+1-k,k)), abs(
				    d__2))) * ((d__3 = X(i,j).r, abs(d__3)) 
				    + (d__4 = d_imag(&X(i,j)), abs(
				    d__4)));
/* L150: */
			}
			RWORK(k) += s;
/* L160: */
		    }
		} else {
		    i__2 = *n;
		    for (k = 1; k <= *n; ++k) {
			i__5 = k + j * x_dim1;
			s = (d__1 = X(k,j).r, abs(d__1)) + (d__2 = d_imag(&X(k,j)), abs(d__2));
/* Computing MIN */
			i__3 = *n, i__4 = k + *kd;
			i__5 = min(i__3,i__4);
			for (i = k + 1; i <= min(*n,k+*kd); ++i) {
			    i__3 = i + 1 - k + k * ab_dim1;
			    i__4 = i + j * x_dim1;
			    s += ((d__1 = AB(i+1-k,k).r, abs(d__1)) + (d__2 = 
				    d_imag(&AB(i+1-k,k)), abs(
				    d__2))) * ((d__3 = X(i,j).r, abs(d__3)) 
				    + (d__4 = d_imag(&X(i,j)), abs(
				    d__4)));
/* L170: */
			}
			RWORK(k) += s;
/* L180: */
		    }
		}
	    }
	}
	s = 0.;
	i__2 = *n;
	for (i = 1; i <= *n; ++i) {
	    if (RWORK(i) > safe2) {
/* Computing MAX */
		i__5 = i;
		d__3 = s, d__4 = ((d__1 = WORK(i).r, abs(d__1)) + (d__2 = 
			d_imag(&WORK(i)), abs(d__2))) / RWORK(i);
		s = max(d__3,d__4);
	    } else {
/* Computing MAX */
		i__5 = i;
		d__3 = s, d__4 = ((d__1 = WORK(i).r, abs(d__1)) + (d__2 = 
			d_imag(&WORK(i)), abs(d__2)) + safe1) / (RWORK(i) + 
			safe1);
		s = max(d__3,d__4);
	    }
/* L190: */
	}
	BERR(j) = s;

/*        Bound error from formula   

          norm(X - XTRUE) / norm(X) .le. FERR =   
          norm( abs(inv(op(A)))*   
             ( abs(R) + NZ*EPS*( abs(op(A))*abs(X)+abs(B) ))) / norm(X
)   

          where   
            norm(Z) is the magnitude of the largest component of Z   
            inv(op(A)) is the inverse of op(A)   
            abs(Z) is the componentwise absolute value of the matrix o
r   
               vector Z   
            NZ is the maximum number of nonzeros in any row of A, plus
 1   
            EPS is machine epsilon   

          The i-th component of abs(R)+NZ*EPS*(abs(op(A))*abs(X)+abs(B
))   
          is incremented by SAFE1 if the i-th component of   
          abs(op(A))*abs(X) + abs(B) is less than SAFE2.   

          Use ZLACON to estimate the infinity-norm of the matrix   
             inv(op(A)) * diag(W),   
          where W = abs(R) + NZ*EPS*( abs(op(A))*abs(X)+abs(B) ))) */

	i__2 = *n;
	for (i = 1; i <= *n; ++i) {
	    if (RWORK(i) > safe2) {
		i__5 = i;
		RWORK(i) = (d__1 = WORK(i).r, abs(d__1)) + (d__2 = d_imag(&
			WORK(i)), abs(d__2)) + nz * eps * RWORK(i);
	    } else {
		i__5 = i;
		RWORK(i) = (d__1 = WORK(i).r, abs(d__1)) + (d__2 = d_imag(&
			WORK(i)), abs(d__2)) + nz * eps * RWORK(i) + safe1;
	    }
/* L200: */
	}

	kase = 0;
L210:
	zlacon_(n, &WORK(*n + 1), &WORK(1), &FERR(j), &kase);
	if (kase != 0) {
	    if (kase == 1) {

/*              Multiply by diag(W)*inv(op(A)**H). */

		ztbsv_(uplo, transt, diag, n, kd, &AB(1,1), ldab, &WORK(
			1), &c__1);
		i__2 = *n;
		for (i = 1; i <= *n; ++i) {
		    i__5 = i;
		    i__3 = i;
		    i__4 = i;
		    z__1.r = RWORK(i) * WORK(i).r, z__1.i = RWORK(i) 
			    * WORK(i).i;
		    WORK(i).r = z__1.r, WORK(i).i = z__1.i;
/* L220: */
		}
	    } else {

/*              Multiply by inv(op(A))*diag(W). */

		i__2 = *n;
		for (i = 1; i <= *n; ++i) {
		    i__5 = i;
		    i__3 = i;
		    i__4 = i;
		    z__1.r = RWORK(i) * WORK(i).r, z__1.i = RWORK(i) 
			    * WORK(i).i;
		    WORK(i).r = z__1.r, WORK(i).i = z__1.i;
/* L230: */
		}
		ztbsv_(uplo, transn, diag, n, kd, &AB(1,1), ldab, &WORK(
			1), &c__1);
	    }
	    goto L210;
	}

/*        Normalize error. */

	lstres = 0.;
	i__2 = *n;
	for (i = 1; i <= *n; ++i) {
/* Computing MAX */
	    i__5 = i + j * x_dim1;
	    d__3 = lstres, d__4 = (d__1 = X(i,j).r, abs(d__1)) + (d__2 = 
		    d_imag(&X(i,j)), abs(d__2));
	    lstres = max(d__3,d__4);
/* L240: */
	}
	if (lstres != 0.) {
	    FERR(j) /= lstres;
	}

/* L250: */
    }

    return 0;

/*     End of ZTBRFS */

} /* ztbrfs_ */

