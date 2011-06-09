#include "f2c.h"

/* Subroutine */ int clatps_(char *uplo, char *trans, char *diag, char *
	normin, integer *n, complex *ap, complex *x, real *scale, real *cnorm,
	 integer *info)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    CLATPS solves one of the triangular systems   

       A * x = s*b,  A**T * x = s*b,  or  A**H * x = s*b,   

    with scaling to prevent overflow, where A is an upper or lower   
    triangular matrix stored in packed form.  Here A**T denotes the   
    transpose of A, A**H denotes the conjugate transpose of A, x and b   
    are n-element vectors, and s is a scaling factor, usually less than   
    or equal to 1, chosen so that the components of x will be less than   
    the overflow threshold.  If the unscaled problem will not cause   
    overflow, the Level 2 BLAS routine CTPSV is called. If the matrix A   
    is singular (A(j,j) = 0 for some j), then s is set to 0 and a   
    non-trivial solution to A*x = 0 is returned.   

    Arguments   
    =========   

    UPLO    (input) CHARACTER*1   
            Specifies whether the matrix A is upper or lower triangular. 
  
            = 'U':  Upper triangular   
            = 'L':  Lower triangular   

    TRANS   (input) CHARACTER*1   
            Specifies the operation applied to A.   
            = 'N':  Solve A * x = s*b     (No transpose)   
            = 'T':  Solve A**T * x = s*b  (Transpose)   
            = 'C':  Solve A**H * x = s*b  (Conjugate transpose)   

    DIAG    (input) CHARACTER*1   
            Specifies whether or not the matrix A is unit triangular.   
            = 'N':  Non-unit triangular   
            = 'U':  Unit triangular   

    NORMIN  (input) CHARACTER*1   
            Specifies whether CNORM has been set or not.   
            = 'Y':  CNORM contains the column norms on entry   
            = 'N':  CNORM is not set on entry.  On exit, the norms will   
                    be computed and stored in CNORM.   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    AP      (input) COMPLEX array, dimension (N*(N+1)/2)   
            The upper or lower triangular matrix A, packed columnwise in 
  
            a linear array.  The j-th column of A is stored in the array 
  
            AP as follows:   
            if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;   
            if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.   

    X       (input/output) COMPLEX array, dimension (N)   
            On entry, the right hand side b of the triangular system.   
            On exit, X is overwritten by the solution vector x.   

    SCALE   (output) REAL   
            The scaling factor s for the triangular system   
               A * x = s*b,  A**T * x = s*b,  or  A**H * x = s*b.   
            If SCALE = 0, the matrix A is singular or badly scaled, and   
            the vector x is an exact or approximate solution to A*x = 0. 
  

    CNORM   (input or output) REAL array, dimension (N)   

            If NORMIN = 'Y', CNORM is an input argument and CNORM(j)   
            contains the norm of the off-diagonal part of the j-th column 
  
            of A.  If TRANS = 'N', CNORM(j) must be greater than or equal 
  
            to the infinity-norm, and if TRANS = 'T' or 'C', CNORM(j)   
            must be greater than or equal to the 1-norm.   

            If NORMIN = 'N', CNORM is an output argument and CNORM(j)   
            returns the 1-norm of the offdiagonal part of the j-th column 
  
            of A.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -k, the k-th argument had an illegal value   

    Further Details   
    ======= =======   

    A rough bound on x is computed; if that is less than overflow, CTPSV 
  
    is called, otherwise, specific code is used which checks for possible 
  
    overflow or divide-by-zero at every operation.   

    A columnwise scheme is used for solving A*x = b.  The basic algorithm 
  
    if A is lower triangular is   

         x[1:n] := b[1:n]   
         for j = 1, ..., n   
              x(j) := x(j) / A(j,j)   
              x[j+1:n] := x[j+1:n] - x(j) * A[j+1:n,j]   
         end   

    Define bounds on the components of x after j iterations of the loop: 
  
       M(j) = bound on x[1:j]   
       G(j) = bound on x[j+1:n]   
    Initially, let M(0) = 0 and G(0) = max{x(i), i=1,...,n}.   

    Then for iteration j+1 we have   
       M(j+1) <= G(j) / | A(j+1,j+1) |   
       G(j+1) <= G(j) + M(j+1) * | A[j+2:n,j+1] |   
              <= G(j) ( 1 + CNORM(j+1) / | A(j+1,j+1) | )   

    where CNORM(j+1) is greater than or equal to the infinity-norm of   
    column j+1 of A, not counting the diagonal.  Hence   

       G(j) <= G(0) product ( 1 + CNORM(i) / | A(i,i) | )   
                    1<=i<=j   
    and   

       |x(j)| <= ( G(0) / |A(j,j)| ) product ( 1 + CNORM(i) / |A(i,i)| ) 
  
                                     1<=i< j   

    Since |x(j)| <= M(j), we use the Level 2 BLAS routine CTPSV if the   
    reciprocal of the largest M(j), j=1,..,n, is larger than   
    max(underflow, 1/overflow).   

    The bound on x(j) is also used to determine when a step in the   
    columnwise method can be performed without fear of overflow.  If   
    the computed bound is greater than a large constant, x is scaled to   
    prevent overflow, but if the bound overflows, x is set to 0, x(j) to 
  
    1, and scale to 0, and a non-trivial solution to A*x = 0 is found.   

    Similarly, a row-wise scheme is used to solve A**T *x = b  or   
    A**H *x = b.  The basic algorithm for A upper triangular is   

         for j = 1, ..., n   
              x(j) := ( b(j) - A[1:j-1,j]' * x[1:j-1] ) / A(j,j)   
         end   

    We simultaneously compute two bounds   
         G(j) = bound on ( b(i) - A[1:i-1,i]' * x[1:i-1] ), 1<=i<=j   
         M(j) = bound on x(i), 1<=i<=j   

    The initial values are G(0) = 0, M(0) = max{b(i), i=1,..,n}, and we   
    add the constraint G(j) >= G(j-1) and M(j) >= M(j-1) for j >= 1.   
    Then the bound on x(j) is   

         M(j) <= M(j-1) * ( 1 + CNORM(j) ) / | A(j,j) |   

              <= M(0) * product ( ( 1 + CNORM(i) ) / |A(i,i)| )   
                        1<=i<=j   

    and we can safely call CTPSV if 1/M(n) and 1/G(n) are both greater   
    than max(underflow, 1/overflow).   

    ===================================================================== 
  


    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__1 = 1;
    static real c_b36 = .5f;
    
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5;
    real r__1, r__2, r__3, r__4;
    complex q__1, q__2, q__3, q__4;
    /* Builtin functions */
    double r_imag(complex *);
    void r_cnjg(complex *, complex *);
    /* Local variables */
    static integer jinc, jlen;
    static real xbnd;
    static integer imax;
    static real tmax;
    static complex tjjs;
    static real xmax, grow;
    static integer i, j;
    extern /* Complex */ VOID cdotc_(complex *, integer *, complex *, integer 
	    *, complex *, integer *);
    extern logical lsame_(char *, char *);
    extern /* Subroutine */ int sscal_(integer *, real *, real *, integer *);
    static real tscal;
    static complex uscal;
    static integer jlast;
    extern /* Complex */ VOID cdotu_(complex *, integer *, complex *, integer 
	    *, complex *, integer *);
    static complex csumj;
    extern /* Subroutine */ int caxpy_(integer *, complex *, complex *, 
	    integer *, complex *, integer *);
    static logical upper;
    extern /* Subroutine */ int ctpsv_(char *, char *, char *, integer *, 
	    complex *, complex *, integer *), slabad_(
	    real *, real *);
    static integer ip;
    static real xj;
    extern integer icamax_(integer *, complex *, integer *);
    extern /* Complex */ VOID cladiv_(complex *, complex *, complex *);
    extern doublereal slamch_(char *);
    extern /* Subroutine */ int csscal_(integer *, real *, complex *, integer 
	    *), xerbla_(char *, integer *);
    static real bignum;
    extern integer isamax_(integer *, real *, integer *);
    extern doublereal scasum_(integer *, complex *, integer *);
    static logical notran;
    static integer jfirst;
    static real smlnum;
    static logical nounit;
    static real rec, tjj;



#define CNORM(I) cnorm[(I)-1]
#define X(I) x[(I)-1]
#define AP(I) ap[(I)-1]


    *info = 0;
    upper = lsame_(uplo, "U");
    notran = lsame_(trans, "N");
    nounit = lsame_(diag, "N");

/*     Test the input parameters. */

    if (! upper && ! lsame_(uplo, "L")) {
	*info = -1;
    } else if (! notran && ! lsame_(trans, "T") && ! lsame_(trans, 
	    "C")) {
	*info = -2;
    } else if (! nounit && ! lsame_(diag, "U")) {
	*info = -3;
    } else if (! lsame_(normin, "Y") && ! lsame_(normin, "N"))
	     {
	*info = -4;
    } else if (*n < 0) {
	*info = -5;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("CLATPS", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return 0;
    }

/*     Determine machine dependent parameters to control overflow. */

    smlnum = slamch_("Safe minimum");
    bignum = 1.f / smlnum;
    slabad_(&smlnum, &bignum);
    smlnum /= slamch_("Precision");
    bignum = 1.f / smlnum;
    *scale = 1.f;

    if (lsame_(normin, "N")) {

/*        Compute the 1-norm of each column, not including the diagona
l. */

	if (upper) {

/*           A is upper triangular. */

	    ip = 1;
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		i__2 = j - 1;
		CNORM(j) = scasum_(&i__2, &AP(ip), &c__1);
		ip += j;
/* L10: */
	    }
	} else {

/*           A is lower triangular. */

	    ip = 1;
	    i__1 = *n - 1;
	    for (j = 1; j <= *n-1; ++j) {
		i__2 = *n - j;
		CNORM(j) = scasum_(&i__2, &AP(ip + 1), &c__1);
		ip = ip + *n - j + 1;
/* L20: */
	    }
	    CNORM(*n) = 0.f;
	}
    }

/*     Scale the column norms by TSCAL if the maximum element in CNORM is 
  
       greater than BIGNUM/2. */

    imax = isamax_(n, &CNORM(1), &c__1);
    tmax = CNORM(imax);
    if (tmax <= bignum * .5f) {
	tscal = 1.f;
    } else {
	tscal = .5f / (smlnum * tmax);
	sscal_(n, &tscal, &CNORM(1), &c__1);
    }

/*     Compute a bound on the computed solution vector to see if the   
       Level 2 BLAS routine CTPSV can be used. */

    xmax = 0.f;
    i__1 = *n;
    for (j = 1; j <= *n; ++j) {
/* Computing MAX */
	i__2 = j;
	r__3 = xmax, r__4 = (r__1 = X(j).r / 2.f, dabs(r__1)) + (r__2 = 
		r_imag(&X(j)) / 2.f, dabs(r__2));
	xmax = dmax(r__3,r__4);
/* L30: */
    }
    xbnd = xmax;
    if (notran) {

/*        Compute the growth in A * x = b. */

	if (upper) {
	    jfirst = *n;
	    jlast = 1;
	    jinc = -1;
	} else {
	    jfirst = 1;
	    jlast = *n;
	    jinc = 1;
	}

	if (tscal != 1.f) {
	    grow = 0.f;
	    goto L60;
	}

	if (nounit) {

/*           A is non-unit triangular.   

             Compute GROW = 1/G(j) and XBND = 1/M(j).   
             Initially, G(0) = max{x(i), i=1,...,n}. */

	    grow = .5f / dmax(xbnd,smlnum);
	    xbnd = grow;
	    ip = jfirst * (jfirst + 1) / 2;
	    jlen = *n;
	    i__1 = jlast;
	    i__2 = jinc;
	    for (j = jfirst; jinc < 0 ? j >= jlast : j <= jlast; j += jinc) {

/*              Exit the loop if the growth factor is too smal
l. */

		if (grow <= smlnum) {
		    goto L60;
		}

		i__3 = ip;
		tjjs.r = AP(ip).r, tjjs.i = AP(ip).i;
		tjj = (r__1 = tjjs.r, dabs(r__1)) + (r__2 = r_imag(&tjjs), 
			dabs(r__2));

		if (tjj >= smlnum) {

/*                 M(j) = G(j-1) / abs(A(j,j))   

   Computing MIN */
		    r__1 = xbnd, r__2 = dmin(1.f,tjj) * grow;
		    xbnd = dmin(r__1,r__2);
		} else {

/*                 M(j) could overflow, set XBND to 0. */

		    xbnd = 0.f;
		}

		if (tjj + CNORM(j) >= smlnum) {

/*                 G(j) = G(j-1)*( 1 + CNORM(j) / abs(A(j,
j)) ) */

		    grow *= tjj / (tjj + CNORM(j));
		} else {

/*                 G(j) could overflow, set GROW to 0. */

		    grow = 0.f;
		}
		ip += jinc * jlen;
		--jlen;
/* L40: */
	    }
	    grow = xbnd;
	} else {

/*           A is unit triangular.   

             Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...
,n}.   

   Computing MIN */
	    r__1 = 1.f, r__2 = .5f / dmax(xbnd,smlnum);
	    grow = dmin(r__1,r__2);
	    i__2 = jlast;
	    i__1 = jinc;
	    for (j = jfirst; jinc < 0 ? j >= jlast : j <= jlast; j += jinc) {

/*              Exit the loop if the growth factor is too smal
l. */

		if (grow <= smlnum) {
		    goto L60;
		}

/*              G(j) = G(j-1)*( 1 + CNORM(j) ) */

		grow *= 1.f / (CNORM(j) + 1.f);
/* L50: */
	    }
	}
L60:

	;
    } else {

/*        Compute the growth in A**T * x = b  or  A**H * x = b. */

	if (upper) {
	    jfirst = 1;
	    jlast = *n;
	    jinc = 1;
	} else {
	    jfirst = *n;
	    jlast = 1;
	    jinc = -1;
	}

	if (tscal != 1.f) {
	    grow = 0.f;
	    goto L90;
	}

	if (nounit) {

/*           A is non-unit triangular.   

             Compute GROW = 1/G(j) and XBND = 1/M(j).   
             Initially, M(0) = max{x(i), i=1,...,n}. */

	    grow = .5f / dmax(xbnd,smlnum);
	    xbnd = grow;
	    ip = jfirst * (jfirst + 1) / 2;
	    jlen = 1;
	    i__1 = jlast;
	    i__2 = jinc;
	    for (j = jfirst; jinc < 0 ? j >= jlast : j <= jlast; j += jinc) {

/*              Exit the loop if the growth factor is too smal
l. */

		if (grow <= smlnum) {
		    goto L90;
		}

/*              G(j) = max( G(j-1), M(j-1)*( 1 + CNORM(j) ) ) 
*/

		xj = CNORM(j) + 1.f;
/* Computing MIN */
		r__1 = grow, r__2 = xbnd / xj;
		grow = dmin(r__1,r__2);

		i__3 = ip;
		tjjs.r = AP(ip).r, tjjs.i = AP(ip).i;
		tjj = (r__1 = tjjs.r, dabs(r__1)) + (r__2 = r_imag(&tjjs), 
			dabs(r__2));

		if (tjj >= smlnum) {

/*                 M(j) = M(j-1)*( 1 + CNORM(j) ) / abs(A(
j,j)) */

		    if (xj > tjj) {
			xbnd *= tjj / xj;
		    }
		} else {

/*                 M(j) could overflow, set XBND to 0. */

		    xbnd = 0.f;
		}
		++jlen;
		ip += jinc * jlen;
/* L70: */
	    }
	    grow = dmin(grow,xbnd);
	} else {

/*           A is unit triangular.   

             Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...
,n}.   

   Computing MIN */
	    r__1 = 1.f, r__2 = .5f / dmax(xbnd,smlnum);
	    grow = dmin(r__1,r__2);
	    i__2 = jlast;
	    i__1 = jinc;
	    for (j = jfirst; jinc < 0 ? j >= jlast : j <= jlast; j += jinc) {

/*              Exit the loop if the growth factor is too smal
l. */

		if (grow <= smlnum) {
		    goto L90;
		}

/*              G(j) = ( 1 + CNORM(j) )*G(j-1) */

		xj = CNORM(j) + 1.f;
		grow /= xj;
/* L80: */
	    }
	}
L90:
	;
    }

    if (grow * tscal > smlnum) {

/*        Use the Level 2 BLAS solve if the reciprocal of the bound on
   
          elements of X is not too small. */

	ctpsv_(uplo, trans, diag, n, &AP(1), &X(1), &c__1);
    } else {

/*        Use a Level 1 BLAS solve, scaling intermediate results. */

	if (xmax > bignum * .5f) {

/*           Scale X so that its components are less than or equal
 to   
             BIGNUM in absolute value. */

	    *scale = bignum * .5f / xmax;
	    csscal_(n, scale, &X(1), &c__1);
	    xmax = bignum;
	} else {
	    xmax *= 2.f;
	}

	if (notran) {

/*           Solve A * x = b */

	    ip = jfirst * (jfirst + 1) / 2;
	    i__1 = jlast;
	    i__2 = jinc;
	    for (j = jfirst; jinc < 0 ? j >= jlast : j <= jlast; j += jinc) {

/*              Compute x(j) = b(j) / A(j,j), scaling x if nec
essary. */

		i__3 = j;
		xj = (r__1 = X(j).r, dabs(r__1)) + (r__2 = r_imag(&X(j)), 
			dabs(r__2));
		if (nounit) {
		    i__3 = ip;
		    q__1.r = tscal * AP(ip).r, q__1.i = tscal * AP(ip).i;
		    tjjs.r = q__1.r, tjjs.i = q__1.i;
		} else {
		    tjjs.r = tscal, tjjs.i = 0.f;
		    if (tscal == 1.f) {
			goto L105;
		    }
		}
		tjj = (r__1 = tjjs.r, dabs(r__1)) + (r__2 = r_imag(&tjjs), 
			dabs(r__2));
		if (tjj > smlnum) {

/*                    abs(A(j,j)) > SMLNUM: */

		    if (tjj < 1.f) {
			if (xj > tjj * bignum) {

/*                          Scale x by 1/b(j). */

			    rec = 1.f / xj;
			    csscal_(n, &rec, &X(1), &c__1);
			    *scale *= rec;
			    xmax *= rec;
			}
		    }
		    i__3 = j;
		    cladiv_(&q__1, &X(j), &tjjs);
		    X(j).r = q__1.r, X(j).i = q__1.i;
		    i__3 = j;
		    xj = (r__1 = X(j).r, dabs(r__1)) + (r__2 = r_imag(&X(j)
			    ), dabs(r__2));
		} else if (tjj > 0.f) {

/*                    0 < abs(A(j,j)) <= SMLNUM: */

		    if (xj > tjj * bignum) {

/*                       Scale x by (1/abs(x(j)))*abs(
A(j,j))*BIGNUM   
                         to avoid overflow when dividi
ng by A(j,j). */

			rec = tjj * bignum / xj;
			if (CNORM(j) > 1.f) {

/*                          Scale by 1/CNORM(j) to
 avoid overflow when   
                            multiplying x(j) times
 column j. */

			    rec /= CNORM(j);
			}
			csscal_(n, &rec, &X(1), &c__1);
			*scale *= rec;
			xmax *= rec;
		    }
		    i__3 = j;
		    cladiv_(&q__1, &X(j), &tjjs);
		    X(j).r = q__1.r, X(j).i = q__1.i;
		    i__3 = j;
		    xj = (r__1 = X(j).r, dabs(r__1)) + (r__2 = r_imag(&X(j)
			    ), dabs(r__2));
		} else {

/*                    A(j,j) = 0:  Set x(1:n) = 0, x(j) = 
1, and   
                      scale = 0, and compute a solution to
 A*x = 0. */

		    i__3 = *n;
		    for (i = 1; i <= *n; ++i) {
			i__4 = i;
			X(i).r = 0.f, X(i).i = 0.f;
/* L100: */
		    }
		    i__3 = j;
		    X(j).r = 1.f, X(j).i = 0.f;
		    xj = 1.f;
		    *scale = 0.f;
		    xmax = 0.f;
		}
L105:

/*              Scale x if necessary to avoid overflow when ad
ding a   
                multiple of column j of A. */

		if (xj > 1.f) {
		    rec = 1.f / xj;
		    if (CNORM(j) > (bignum - xmax) * rec) {

/*                    Scale x by 1/(2*abs(x(j))). */

			rec *= .5f;
			csscal_(n, &rec, &X(1), &c__1);
			*scale *= rec;
		    }
		} else if (xj * CNORM(j) > bignum - xmax) {

/*                 Scale x by 1/2. */

		    csscal_(n, &c_b36, &X(1), &c__1);
		    *scale *= .5f;
		}

		if (upper) {
		    if (j > 1) {

/*                    Compute the update   
                         x(1:j-1) := x(1:j-1) - x(j) *
 A(1:j-1,j) */

			i__3 = j - 1;
			i__4 = j;
			q__2.r = -(doublereal)X(j).r, q__2.i = -(
				doublereal)X(j).i;
			q__1.r = tscal * q__2.r, q__1.i = tscal * q__2.i;
			caxpy_(&i__3, &q__1, &AP(ip - j + 1), &c__1, &X(1), &
				c__1);
			i__3 = j - 1;
			i = icamax_(&i__3, &X(1), &c__1);
			i__3 = i;
			xmax = (r__1 = X(i).r, dabs(r__1)) + (r__2 = 
				r_imag(&X(i)), dabs(r__2));
		    }
		    ip -= j;
		} else {
		    if (j < *n) {

/*                    Compute the update   
                         x(j+1:n) := x(j+1:n) - x(j) *
 A(j+1:n,j) */

			i__3 = *n - j;
			i__4 = j;
			q__2.r = -(doublereal)X(j).r, q__2.i = -(
				doublereal)X(j).i;
			q__1.r = tscal * q__2.r, q__1.i = tscal * q__2.i;
			caxpy_(&i__3, &q__1, &AP(ip + 1), &c__1, &X(j + 1), &
				c__1);
			i__3 = *n - j;
			i = j + icamax_(&i__3, &X(j + 1), &c__1);
			i__3 = i;
			xmax = (r__1 = X(i).r, dabs(r__1)) + (r__2 = 
				r_imag(&X(i)), dabs(r__2));
		    }
		    ip = ip + *n - j + 1;
		}
/* L110: */
	    }

	} else if (lsame_(trans, "T")) {

/*           Solve A**T * x = b */

	    ip = jfirst * (jfirst + 1) / 2;
	    jlen = 1;
	    i__2 = jlast;
	    i__1 = jinc;
	    for (j = jfirst; jinc < 0 ? j >= jlast : j <= jlast; j += jinc) {

/*              Compute x(j) = b(j) - sum A(k,j)*x(k).   
                                      k<>j */

		i__3 = j;
		xj = (r__1 = X(j).r, dabs(r__1)) + (r__2 = r_imag(&X(j)), 
			dabs(r__2));
		uscal.r = tscal, uscal.i = 0.f;
		rec = 1.f / dmax(xmax,1.f);
		if (CNORM(j) > (bignum - xj) * rec) {

/*                 If x(j) could overflow, scale x by 1/(2
*XMAX). */

		    rec *= .5f;
		    if (nounit) {
			i__3 = ip;
			q__1.r = tscal * AP(ip).r, q__1.i = tscal * AP(ip)
				.i;
			tjjs.r = q__1.r, tjjs.i = q__1.i;
		    } else {
			tjjs.r = tscal, tjjs.i = 0.f;
		    }
		    tjj = (r__1 = tjjs.r, dabs(r__1)) + (r__2 = r_imag(&tjjs),
			     dabs(r__2));
		    if (tjj > 1.f) {

/*                       Divide by A(j,j) when scaling
 x if A(j,j) > 1.   

   Computing MIN */
			r__1 = 1.f, r__2 = rec * tjj;
			rec = dmin(r__1,r__2);
			cladiv_(&q__1, &uscal, &tjjs);
			uscal.r = q__1.r, uscal.i = q__1.i;
		    }
		    if (rec < 1.f) {
			csscal_(n, &rec, &X(1), &c__1);
			*scale *= rec;
			xmax *= rec;
		    }
		}

		csumj.r = 0.f, csumj.i = 0.f;
		if (uscal.r == 1.f && uscal.i == 0.f) {

/*                 If the scaling needed for A in the dot 
product is 1,   
                   call CDOTU to perform the dot product. 
*/

		    if (upper) {
			i__3 = j - 1;
			cdotu_(&q__1, &i__3, &AP(ip - j + 1), &c__1, &X(1), &
				c__1);
			csumj.r = q__1.r, csumj.i = q__1.i;
		    } else if (j < *n) {
			i__3 = *n - j;
			cdotu_(&q__1, &i__3, &AP(ip + 1), &c__1, &X(j + 1), &
				c__1);
			csumj.r = q__1.r, csumj.i = q__1.i;
		    }
		} else {

/*                 Otherwise, use in-line code for the dot
 product. */

		    if (upper) {
			i__3 = j - 1;
			for (i = 1; i <= j-1; ++i) {
			    i__4 = ip - j + i;
			    q__3.r = AP(ip-j+i).r * uscal.r - AP(ip-j+i).i * 
				    uscal.i, q__3.i = AP(ip-j+i).r * uscal.i + 
				    AP(ip-j+i).i * uscal.r;
			    i__5 = i;
			    q__2.r = q__3.r * X(i).r - q__3.i * X(i).i, 
				    q__2.i = q__3.r * X(i).i + q__3.i * X(
				    i).r;
			    q__1.r = csumj.r + q__2.r, q__1.i = csumj.i + 
				    q__2.i;
			    csumj.r = q__1.r, csumj.i = q__1.i;
/* L120: */
			}
		    } else if (j < *n) {
			i__3 = *n - j;
			for (i = 1; i <= *n-j; ++i) {
			    i__4 = ip + i;
			    q__3.r = AP(ip+i).r * uscal.r - AP(ip+i).i * 
				    uscal.i, q__3.i = AP(ip+i).r * uscal.i + 
				    AP(ip+i).i * uscal.r;
			    i__5 = j + i;
			    q__2.r = q__3.r * X(j+i).r - q__3.i * X(j+i).i, 
				    q__2.i = q__3.r * X(j+i).i + q__3.i * X(
				    j+i).r;
			    q__1.r = csumj.r + q__2.r, q__1.i = csumj.i + 
				    q__2.i;
			    csumj.r = q__1.r, csumj.i = q__1.i;
/* L130: */
			}
		    }
		}

		q__1.r = tscal, q__1.i = 0.f;
		if (uscal.r == q__1.r && uscal.i == q__1.i) {

/*                 Compute x(j) := ( x(j) - CSUMJ ) / A(j,
j) if 1/A(j,j)   
                   was not used to scale the dotproduct. 
*/

		    i__3 = j;
		    i__4 = j;
		    q__1.r = X(j).r - csumj.r, q__1.i = X(j).i - 
			    csumj.i;
		    X(j).r = q__1.r, X(j).i = q__1.i;
		    i__3 = j;
		    xj = (r__1 = X(j).r, dabs(r__1)) + (r__2 = r_imag(&X(j)
			    ), dabs(r__2));
		    if (nounit) {

/*                    Compute x(j) = x(j) / A(j,j), sc
aling if necessary. */

			i__3 = ip;
			q__1.r = tscal * AP(ip).r, q__1.i = tscal * AP(ip)
				.i;
			tjjs.r = q__1.r, tjjs.i = q__1.i;
		    } else {
			tjjs.r = tscal, tjjs.i = 0.f;
			if (tscal == 1.f) {
			    goto L145;
			}
		    }
		    tjj = (r__1 = tjjs.r, dabs(r__1)) + (r__2 = r_imag(&tjjs),
			     dabs(r__2));
		    if (tjj > smlnum) {

/*                       abs(A(j,j)) > SMLNUM: */

			if (tjj < 1.f) {
			    if (xj > tjj * bignum) {

/*                             Scale X by 1/ab
s(x(j)). */

				rec = 1.f / xj;
				csscal_(n, &rec, &X(1), &c__1);
				*scale *= rec;
				xmax *= rec;
			    }
			}
			i__3 = j;
			cladiv_(&q__1, &X(j), &tjjs);
			X(j).r = q__1.r, X(j).i = q__1.i;
		    } else if (tjj > 0.f) {

/*                       0 < abs(A(j,j)) <= SMLNUM: */

			if (xj > tjj * bignum) {

/*                          Scale x by (1/abs(x(j)
))*abs(A(j,j))*BIGNUM. */

			    rec = tjj * bignum / xj;
			    csscal_(n, &rec, &X(1), &c__1);
			    *scale *= rec;
			    xmax *= rec;
			}
			i__3 = j;
			cladiv_(&q__1, &X(j), &tjjs);
			X(j).r = q__1.r, X(j).i = q__1.i;
		    } else {

/*                       A(j,j) = 0:  Set x(1:n) = 0, 
x(j) = 1, and   
                         scale = 0 and compute a solut
ion to A**T *x = 0. */

			i__3 = *n;
			for (i = 1; i <= *n; ++i) {
			    i__4 = i;
			    X(i).r = 0.f, X(i).i = 0.f;
/* L140: */
			}
			i__3 = j;
			X(j).r = 1.f, X(j).i = 0.f;
			*scale = 0.f;
			xmax = 0.f;
		    }
L145:
		    ;
		} else {

/*                 Compute x(j) := x(j) / A(j,j) - CSUMJ i
f the dot   
                   product has already been divided by 1/A
(j,j). */

		    i__3 = j;
		    cladiv_(&q__2, &X(j), &tjjs);
		    q__1.r = q__2.r - csumj.r, q__1.i = q__2.i - csumj.i;
		    X(j).r = q__1.r, X(j).i = q__1.i;
		}
/* Computing MAX */
		i__3 = j;
		r__3 = xmax, r__4 = (r__1 = X(j).r, dabs(r__1)) + (r__2 = 
			r_imag(&X(j)), dabs(r__2));
		xmax = dmax(r__3,r__4);
		++jlen;
		ip += jinc * jlen;
/* L150: */
	    }

	} else {

/*           Solve A**H * x = b */

	    ip = jfirst * (jfirst + 1) / 2;
	    jlen = 1;
	    i__1 = jlast;
	    i__2 = jinc;
	    for (j = jfirst; jinc < 0 ? j >= jlast : j <= jlast; j += jinc) {

/*              Compute x(j) = b(j) - sum A(k,j)*x(k).   
                                      k<>j */

		i__3 = j;
		xj = (r__1 = X(j).r, dabs(r__1)) + (r__2 = r_imag(&X(j)), 
			dabs(r__2));
		uscal.r = tscal, uscal.i = 0.f;
		rec = 1.f / dmax(xmax,1.f);
		if (CNORM(j) > (bignum - xj) * rec) {

/*                 If x(j) could overflow, scale x by 1/(2
*XMAX). */

		    rec *= .5f;
		    if (nounit) {
			r_cnjg(&q__2, &AP(ip));
			q__1.r = tscal * q__2.r, q__1.i = tscal * q__2.i;
			tjjs.r = q__1.r, tjjs.i = q__1.i;
		    } else {
			tjjs.r = tscal, tjjs.i = 0.f;
		    }
		    tjj = (r__1 = tjjs.r, dabs(r__1)) + (r__2 = r_imag(&tjjs),
			     dabs(r__2));
		    if (tjj > 1.f) {

/*                       Divide by A(j,j) when scaling
 x if A(j,j) > 1.   

   Computing MIN */
			r__1 = 1.f, r__2 = rec * tjj;
			rec = dmin(r__1,r__2);
			cladiv_(&q__1, &uscal, &tjjs);
			uscal.r = q__1.r, uscal.i = q__1.i;
		    }
		    if (rec < 1.f) {
			csscal_(n, &rec, &X(1), &c__1);
			*scale *= rec;
			xmax *= rec;
		    }
		}

		csumj.r = 0.f, csumj.i = 0.f;
		if (uscal.r == 1.f && uscal.i == 0.f) {

/*                 If the scaling needed for A in the dot 
product is 1,   
                   call CDOTC to perform the dot product. 
*/

		    if (upper) {
			i__3 = j - 1;
			cdotc_(&q__1, &i__3, &AP(ip - j + 1), &c__1, &X(1), &
				c__1);
			csumj.r = q__1.r, csumj.i = q__1.i;
		    } else if (j < *n) {
			i__3 = *n - j;
			cdotc_(&q__1, &i__3, &AP(ip + 1), &c__1, &X(j + 1), &
				c__1);
			csumj.r = q__1.r, csumj.i = q__1.i;
		    }
		} else {

/*                 Otherwise, use in-line code for the dot
 product. */

		    if (upper) {
			i__3 = j - 1;
			for (i = 1; i <= j-1; ++i) {
			    r_cnjg(&q__4, &AP(ip - j + i));
			    q__3.r = q__4.r * uscal.r - q__4.i * uscal.i, 
				    q__3.i = q__4.r * uscal.i + q__4.i * 
				    uscal.r;
			    i__4 = i;
			    q__2.r = q__3.r * X(i).r - q__3.i * X(i).i, 
				    q__2.i = q__3.r * X(i).i + q__3.i * X(
				    i).r;
			    q__1.r = csumj.r + q__2.r, q__1.i = csumj.i + 
				    q__2.i;
			    csumj.r = q__1.r, csumj.i = q__1.i;
/* L160: */
			}
		    } else if (j < *n) {
			i__3 = *n - j;
			for (i = 1; i <= *n-j; ++i) {
			    r_cnjg(&q__4, &AP(ip + i));
			    q__3.r = q__4.r * uscal.r - q__4.i * uscal.i, 
				    q__3.i = q__4.r * uscal.i + q__4.i * 
				    uscal.r;
			    i__4 = j + i;
			    q__2.r = q__3.r * X(j+i).r - q__3.i * X(j+i).i, 
				    q__2.i = q__3.r * X(j+i).i + q__3.i * X(
				    j+i).r;
			    q__1.r = csumj.r + q__2.r, q__1.i = csumj.i + 
				    q__2.i;
			    csumj.r = q__1.r, csumj.i = q__1.i;
/* L170: */
			}
		    }
		}

		q__1.r = tscal, q__1.i = 0.f;
		if (uscal.r == q__1.r && uscal.i == q__1.i) {

/*                 Compute x(j) := ( x(j) - CSUMJ ) / A(j,
j) if 1/A(j,j)   
                   was not used to scale the dotproduct. 
*/

		    i__3 = j;
		    i__4 = j;
		    q__1.r = X(j).r - csumj.r, q__1.i = X(j).i - 
			    csumj.i;
		    X(j).r = q__1.r, X(j).i = q__1.i;
		    i__3 = j;
		    xj = (r__1 = X(j).r, dabs(r__1)) + (r__2 = r_imag(&X(j)
			    ), dabs(r__2));
		    if (nounit) {

/*                    Compute x(j) = x(j) / A(j,j), sc
aling if necessary. */

			r_cnjg(&q__2, &AP(ip));
			q__1.r = tscal * q__2.r, q__1.i = tscal * q__2.i;
			tjjs.r = q__1.r, tjjs.i = q__1.i;
		    } else {
			tjjs.r = tscal, tjjs.i = 0.f;
			if (tscal == 1.f) {
			    goto L185;
			}
		    }
		    tjj = (r__1 = tjjs.r, dabs(r__1)) + (r__2 = r_imag(&tjjs),
			     dabs(r__2));
		    if (tjj > smlnum) {

/*                       abs(A(j,j)) > SMLNUM: */

			if (tjj < 1.f) {
			    if (xj > tjj * bignum) {

/*                             Scale X by 1/ab
s(x(j)). */

				rec = 1.f / xj;
				csscal_(n, &rec, &X(1), &c__1);
				*scale *= rec;
				xmax *= rec;
			    }
			}
			i__3 = j;
			cladiv_(&q__1, &X(j), &tjjs);
			X(j).r = q__1.r, X(j).i = q__1.i;
		    } else if (tjj > 0.f) {

/*                       0 < abs(A(j,j)) <= SMLNUM: */

			if (xj > tjj * bignum) {

/*                          Scale x by (1/abs(x(j)
))*abs(A(j,j))*BIGNUM. */

			    rec = tjj * bignum / xj;
			    csscal_(n, &rec, &X(1), &c__1);
			    *scale *= rec;
			    xmax *= rec;
			}
			i__3 = j;
			cladiv_(&q__1, &X(j), &tjjs);
			X(j).r = q__1.r, X(j).i = q__1.i;
		    } else {

/*                       A(j,j) = 0:  Set x(1:n) = 0, 
x(j) = 1, and   
                         scale = 0 and compute a solut
ion to A**H *x = 0. */

			i__3 = *n;
			for (i = 1; i <= *n; ++i) {
			    i__4 = i;
			    X(i).r = 0.f, X(i).i = 0.f;
/* L180: */
			}
			i__3 = j;
			X(j).r = 1.f, X(j).i = 0.f;
			*scale = 0.f;
			xmax = 0.f;
		    }
L185:
		    ;
		} else {

/*                 Compute x(j) := x(j) / A(j,j) - CSUMJ i
f the dot   
                   product has already been divided by 1/A
(j,j). */

		    i__3 = j;
		    cladiv_(&q__2, &X(j), &tjjs);
		    q__1.r = q__2.r - csumj.r, q__1.i = q__2.i - csumj.i;
		    X(j).r = q__1.r, X(j).i = q__1.i;
		}
/* Computing MAX */
		i__3 = j;
		r__3 = xmax, r__4 = (r__1 = X(j).r, dabs(r__1)) + (r__2 = 
			r_imag(&X(j)), dabs(r__2));
		xmax = dmax(r__3,r__4);
		++jlen;
		ip += jinc * jlen;
/* L190: */
	    }
	}
	*scale /= tscal;
    }

/*     Scale the column norms by 1/TSCAL for return. */

    if (tscal != 1.f) {
	r__1 = 1.f / tscal;
	sscal_(n, &r__1, &CNORM(1), &c__1);
    }

    return 0;

/*     End of CLATPS */

} /* clatps_ */
