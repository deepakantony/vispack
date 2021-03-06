#include "f2c.h"

/* Subroutine */ int dlaqtr_(logical *ltran, logical *lreal, integer *n, 
	doublereal *t, integer *ldt, doublereal *b, doublereal *w, doublereal 
	*scale, doublereal *x, doublereal *work, integer *info)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       March 31, 1993   


    Purpose   
    =======   

    DLAQTR solves the real quasi-triangular system   

                 op(T)*p = scale*c,               if LREAL = .TRUE.   

    or the complex quasi-triangular systems   

               op(T + iB)*(p+iq) = scale*(c+id),  if LREAL = .FALSE.   

    in real arithmetic, where T is upper quasi-triangular.   
    If LREAL = .FALSE., then the first diagonal block of T must be   
    1 by 1, B is the specially structured matrix   

                   B = [ b(1) b(2) ... b(n) ]   
                       [       w            ]   
                       [           w        ]   
                       [              .     ]   
                       [                 w  ]   

    op(A) = A or A', A' denotes the conjugate transpose of   
    matrix A.   

    On input, X = [ c ].  On output, X = [ p ].   
                  [ d ]                  [ q ]   

    This subroutine is designed for the condition number estimation   
    in routine DTRSNA.   

    Arguments   
    =========   

    LTRAN   (input) LOGICAL   
            On entry, LTRAN specifies the option of conjugate transpose: 
  
               = .FALSE.,    op(T+i*B) = T+i*B,   
               = .TRUE.,     op(T+i*B) = (T+i*B)'.   

    LREAL   (input) LOGICAL   
            On entry, LREAL specifies the input matrix structure:   
               = .FALSE.,    the input is complex   
               = .TRUE.,     the input is real   

    N       (input) INTEGER   
            On entry, N specifies the order of T+i*B. N >= 0.   

    T       (input) DOUBLE PRECISION array, dimension (LDT,N)   
            On entry, T contains a matrix in Schur canonical form.   
            If LREAL = .FALSE., then the first diagonal block of T mu   
            be 1 by 1.   

    LDT     (input) INTEGER   
            The leading dimension of the matrix T. LDT >= max(1,N).   

    B       (input) DOUBLE PRECISION array, dimension (N)   
            On entry, B contains the elements to form the matrix   
            B as described above.   
            If LREAL = .TRUE., B is not referenced.   

    W       (input) DOUBLE PRECISION   
            On entry, W is the diagonal element of the matrix B.   
            If LREAL = .TRUE., W is not referenced.   

    SCALE   (output) DOUBLE PRECISION   
            On exit, SCALE is the scale factor.   

    X       (input/output) DOUBLE PRECISION array, dimension (2*N)   
            On entry, X contains the right hand side of the system.   
            On exit, X is overwritten by the solution.   

    WORK    (workspace) DOUBLE PRECISION array, dimension (N)   

    INFO    (output) INTEGER   
            On exit, INFO is set to   
               0: successful exit.   
                 1: the some diagonal 1 by 1 block has been perturbed by 
  
                    a small number SMIN to keep nonsingularity.   
                 2: the some diagonal 2 by 2 block has been perturbed by 
  
                    a small number in DLALN2 to keep nonsingularity.   
            NOTE: In the interests of speed, this routine does not   
                  check the inputs for errors.   

    ===================================================================== 
  


       Do not test the input parameters for errors   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__1 = 1;
    static logical c_false = FALSE_;
    static integer c__2 = 2;
    static doublereal c_b21 = 1.;
    static doublereal c_b25 = 0.;
    static logical c_true = TRUE_;
    
    /* System generated locals */
    integer t_dim1, t_offset, i__1, i__2;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6;
    /* Local variables */
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer ierr;
    static doublereal smin, xmax, d[4]	/* was [2][2] */;
    static integer i, j, k;
    static doublereal v[4]	/* was [2][2] */, z;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    extern doublereal dasum_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    static integer jnext, j1, j2;
    static doublereal sminw;
    static integer n1, n2;
    static doublereal xnorm;
    extern /* Subroutine */ int dlaln2_(logical *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     doublereal *, doublereal *, integer *, doublereal *, doublereal *
	    , doublereal *, integer *, doublereal *, doublereal *, integer *);
    extern doublereal dlamch_(char *), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *);
    static doublereal si, xj;
    extern integer idamax_(integer *, doublereal *, integer *);
    static doublereal scaloc, sr;
    extern /* Subroutine */ int dladiv_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);
    static doublereal bignum;
    static logical notran;
    static doublereal smlnum, rec, eps, tjj, tmp;



#define V(I) v[(I)]
#define WAS(I) was[(I)]
#define B(I) b[(I)-1]
#define X(I) x[(I)-1]
#define WORK(I) work[(I)-1]

#define T(I,J) t[(I)-1 + ((J)-1)* ( *ldt)]

    notran = ! (*ltran);
    *info = 0;

/*     Quick return if possible */

    if (*n == 0) {
	return 0;
    }

/*     Set constants to control overflow */

    eps = dlamch_("P");
    smlnum = dlamch_("S") / eps;
    bignum = 1. / smlnum;

    xnorm = dlange_("M", n, n, &T(1,1), ldt, d);
    if (! (*lreal)) {
/* Computing MAX */
	d__1 = xnorm, d__2 = abs(*w), d__1 = max(d__1,d__2), d__2 = dlange_(
		"M", n, &c__1, &B(1), n, d);
	xnorm = max(d__1,d__2);
    }
/* Computing MAX */
    d__1 = smlnum, d__2 = eps * xnorm;
    smin = max(d__1,d__2);

/*     Compute 1-norm of each column of strictly upper triangular   
       part of T to control overflow in triangular solver. */

    WORK(1) = 0.;
    i__1 = *n;
    for (j = 2; j <= *n; ++j) {
	i__2 = j - 1;
	WORK(j) = dasum_(&i__2, &T(1,j), &c__1);
/* L10: */
    }

    if (! (*lreal)) {
	i__1 = *n;
	for (i = 2; i <= *n; ++i) {
	    WORK(i) += (d__1 = B(i), abs(d__1));
/* L20: */
	}
    }

    n2 = *n << 1;
    n1 = *n;
    if (! (*lreal)) {
	n1 = n2;
    }
    k = idamax_(&n1, &X(1), &c__1);
    xmax = (d__1 = X(k), abs(d__1));
    *scale = 1.;

    if (xmax > bignum) {
	*scale = bignum / xmax;
	dscal_(&n1, scale, &X(1), &c__1);
	xmax = bignum;
    }

    if (*lreal) {

	if (notran) {

/*           Solve T*p = scale*c */

	    jnext = *n;
	    for (j = *n; j >= 1; --j) {
		if (j > jnext) {
		    goto L30;
		}
		j1 = j;
		j2 = j;
		jnext = j - 1;
		if (j > 1 && T(j,j-1) != 0.) {
		    j1 = j - 1;
		    j2 = j;
		    jnext = j - 2;
		}

		if (j1 == j2) {

/*                 Meet 1 by 1 diagonal block   

                   Scale to avoid overflow when computing 
  
                       x(j) = b(j)/T(j,j) */

		    xj = (d__1 = X(j1), abs(d__1));
		    tjj = (d__1 = T(j1,j1), abs(d__1));
		    tmp = T(j1,j1);
		    if (tjj < smin) {
			tmp = smin;
			tjj = smin;
			*info = 1;
		    }

		    if (xj == 0.) {
			goto L30;
		    }

		    if (tjj < 1.) {
			if (xj > bignum * tjj) {
			    rec = 1. / xj;
			    dscal_(n, &rec, &X(1), &c__1);
			    *scale *= rec;
			    xmax *= rec;
			}
		    }
		    X(j1) /= tmp;
		    xj = (d__1 = X(j1), abs(d__1));

/*                 Scale x if necessary to avoid overflow 
when adding a   
                   multiple of column j1 of T. */

		    if (xj > 1.) {
			rec = 1. / xj;
			if (WORK(j1) > (bignum - xmax) * rec) {
			    dscal_(n, &rec, &X(1), &c__1);
			    *scale *= rec;
			}
		    }
		    if (j1 > 1) {
			i__1 = j1 - 1;
			d__1 = -X(j1);
			daxpy_(&i__1, &d__1, &T(1,j1), &c__1, &X(1)
				, &c__1);
			i__1 = j1 - 1;
			k = idamax_(&i__1, &X(1), &c__1);
			xmax = (d__1 = X(k), abs(d__1));
		    }

		} else {

/*                 Meet 2 by 2 diagonal block   

                   Call 2 by 2 linear system solve, to tak
e   
                   care of possible overflow by scaling fa
ctor. */

		    d[0] = X(j1);
		    d[1] = X(j2);
		    dlaln2_(&c_false, &c__2, &c__1, &smin, &c_b21, &T(j1,j1), ldt, &c_b21, &c_b21, d, &c__2, &c_b25, 
			    &c_b25, v, &c__2, &scaloc, &xnorm, &ierr);
		    if (ierr != 0) {
			*info = 2;
		    }

		    if (scaloc != 1.) {
			dscal_(n, &scaloc, &X(1), &c__1);
			*scale *= scaloc;
		    }
		    X(j1) = V(0);
		    X(j2) = V(1);

/*                 Scale V(1,1) (= X(J1)) and/or V(2,1) (=
X(J2))   
                   to avoid overflow in updating right-han
d side.   

   Computing MAX */
		    d__1 = abs(V(0)), d__2 = abs(V(1));
		    xj = max(d__1,d__2);
		    if (xj > 1.) {
			rec = 1. / xj;
/* Computing MAX */
			d__1 = WORK(j1), d__2 = WORK(j2);
			if (max(d__1,d__2) > (bignum - xmax) * rec) {
			    dscal_(n, &rec, &X(1), &c__1);
			    *scale *= rec;
			}
		    }

/*                 Update right-hand side */

		    if (j1 > 1) {
			i__1 = j1 - 1;
			d__1 = -X(j1);
			daxpy_(&i__1, &d__1, &T(1,j1), &c__1, &X(1)
				, &c__1);
			i__1 = j1 - 1;
			d__1 = -X(j2);
			daxpy_(&i__1, &d__1, &T(1,j2), &c__1, &X(1)
				, &c__1);
			i__1 = j1 - 1;
			k = idamax_(&i__1, &X(1), &c__1);
			xmax = (d__1 = X(k), abs(d__1));
		    }

		}

L30:
		;
	    }

	} else {

/*           Solve T'*p = scale*c */

	    jnext = 1;
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		if (j < jnext) {
		    goto L40;
		}
		j1 = j;
		j2 = j;
		jnext = j + 1;
		if (j < *n && T(j+1,j) != 0.) {
		    j1 = j;
		    j2 = j + 1;
		    jnext = j + 2;
		}

		if (j1 == j2) {

/*                 1 by 1 diagonal block   

                   Scale if necessary to avoid overflow in
 forming the   
                   right-hand side element by inner produc
t. */

		    xj = (d__1 = X(j1), abs(d__1));
		    if (xmax > 1.) {
			rec = 1. / xmax;
			if (WORK(j1) > (bignum - xj) * rec) {
			    dscal_(n, &rec, &X(1), &c__1);
			    *scale *= rec;
			    xmax *= rec;
			}
		    }

		    i__2 = j1 - 1;
		    X(j1) -= ddot_(&i__2, &T(1,j1), &c__1, &X(1), &
			    c__1);

		    xj = (d__1 = X(j1), abs(d__1));
		    tjj = (d__1 = T(j1,j1), abs(d__1));
		    tmp = T(j1,j1);
		    if (tjj < smin) {
			tmp = smin;
			tjj = smin;
			*info = 1;
		    }

		    if (tjj < 1.) {
			if (xj > bignum * tjj) {
			    rec = 1. / xj;
			    dscal_(n, &rec, &X(1), &c__1);
			    *scale *= rec;
			    xmax *= rec;
			}
		    }
		    X(j1) /= tmp;
/* Computing MAX */
		    d__2 = xmax, d__3 = (d__1 = X(j1), abs(d__1));
		    xmax = max(d__2,d__3);

		} else {

/*                 2 by 2 diagonal block   

                   Scale if necessary to avoid overflow in
 forming the   
                   right-hand side elements by inner produ
ct.   

   Computing MAX */
		    d__3 = (d__1 = X(j1), abs(d__1)), d__4 = (d__2 = X(j2), 
			    abs(d__2));
		    xj = max(d__3,d__4);
		    if (xmax > 1.) {
			rec = 1. / xmax;
/* Computing MAX */
			d__1 = WORK(j2), d__2 = WORK(j1);
			if (max(d__1,d__2) > (bignum - xj) * rec) {
			    dscal_(n, &rec, &X(1), &c__1);
			    *scale *= rec;
			    xmax *= rec;
			}
		    }

		    i__2 = j1 - 1;
		    d[0] = X(j1) - ddot_(&i__2, &T(1,j1), &c__1, &
			    X(1), &c__1);
		    i__2 = j1 - 1;
		    d[1] = X(j2) - ddot_(&i__2, &T(1,j2), &c__1, &
			    X(1), &c__1);

		    dlaln2_(&c_true, &c__2, &c__1, &smin, &c_b21, &T(j1,j1), ldt, &c_b21, &c_b21, d, &c__2, &c_b25, &
			    c_b25, v, &c__2, &scaloc, &xnorm, &ierr);
		    if (ierr != 0) {
			*info = 2;
		    }

		    if (scaloc != 1.) {
			dscal_(n, &scaloc, &X(1), &c__1);
			*scale *= scaloc;
		    }
		    X(j1) = V(0);
		    X(j2) = V(1);
/* Computing MAX */
		    d__3 = (d__1 = X(j1), abs(d__1)), d__4 = (d__2 = X(j2), 
			    abs(d__2)), d__3 = max(d__3,d__4);
		    xmax = max(d__3,xmax);

		}
L40:
		;
	    }
	}

    } else {

/* Computing MAX */
	d__1 = eps * abs(*w);
	sminw = max(d__1,smin);
	if (notran) {

/*           Solve (T + iB)*(p+iq) = c+id */

	    jnext = *n;
	    for (j = *n; j >= 1; --j) {
		if (j > jnext) {
		    goto L70;
		}
		j1 = j;
		j2 = j;
		jnext = j - 1;
		if (j > 1 && T(j,j-1) != 0.) {
		    j1 = j - 1;
		    j2 = j;
		    jnext = j - 2;
		}

		if (j1 == j2) {

/*                 1 by 1 diagonal block   

                   Scale if necessary to avoid overflow in
 division */

		    z = *w;
		    if (j1 == 1) {
			z = B(1);
		    }
		    xj = (d__1 = X(j1), abs(d__1)) + (d__2 = X(*n + j1), abs(
			    d__2));
		    tjj = (d__1 = T(j1,j1), abs(d__1)) + abs(z);
		    tmp = T(j1,j1);
		    if (tjj < sminw) {
			tmp = sminw;
			tjj = sminw;
			*info = 1;
		    }

		    if (xj == 0.) {
			goto L70;
		    }

		    if (tjj < 1.) {
			if (xj > bignum * tjj) {
			    rec = 1. / xj;
			    dscal_(&n2, &rec, &X(1), &c__1);
			    *scale *= rec;
			    xmax *= rec;
			}
		    }
		    dladiv_(&X(j1), &X(*n + j1), &tmp, &z, &sr, &si);
		    X(j1) = sr;
		    X(*n + j1) = si;
		    xj = (d__1 = X(j1), abs(d__1)) + (d__2 = X(*n + j1), abs(
			    d__2));

/*                 Scale x if necessary to avoid overflow 
when adding a   
                   multiple of column j1 of T. */

		    if (xj > 1.) {
			rec = 1. / xj;
			if (WORK(j1) > (bignum - xmax) * rec) {
			    dscal_(&n2, &rec, &X(1), &c__1);
			    *scale *= rec;
			}
		    }

		    if (j1 > 1) {
			i__1 = j1 - 1;
			d__1 = -X(j1);
			daxpy_(&i__1, &d__1, &T(1,j1), &c__1, &X(1)
				, &c__1);
			i__1 = j1 - 1;
			d__1 = -X(*n + j1);
			daxpy_(&i__1, &d__1, &T(1,j1), &c__1, &X(*
				n + 1), &c__1);

			X(1) += B(j1) * X(*n + j1);
			X(*n + 1) -= B(j1) * X(j1);

			xmax = 0.;
			i__1 = j1 - 1;
			for (k = 1; k <= j1-1; ++k) {
/* Computing MAX */
			    d__3 = xmax, d__4 = (d__1 = X(k), abs(d__1)) + (
				    d__2 = X(k + *n), abs(d__2));
			    xmax = max(d__3,d__4);
/* L50: */
			}
		    }

		} else {

/*                 Meet 2 by 2 diagonal block */

		    d[0] = X(j1);
		    d[1] = X(j2);
		    d[2] = X(*n + j1);
		    d[3] = X(*n + j2);
		    d__1 = -(*w);
		    dlaln2_(&c_false, &c__2, &c__2, &sminw, &c_b21, &T(j1,j1), ldt, &c_b21, &c_b21, d, &c__2, &
			    c_b25, &d__1, v, &c__2, &scaloc, &xnorm, &ierr);
		    if (ierr != 0) {
			*info = 2;
		    }

		    if (scaloc != 1.) {
			i__1 = *n << 1;
			dscal_(&i__1, &scaloc, &X(1), &c__1);
			*scale = scaloc * *scale;
		    }
		    X(j1) = V(0);
		    X(j2) = V(1);
		    X(*n + j1) = V(2);
		    X(*n + j2) = V(3);

/*                 Scale X(J1), .... to avoid overflow in 
  
                   updating right hand side.   

   Computing MAX */
		    d__1 = abs(V(0)) + abs(V(2)), d__2 = abs(V(1)) + abs(V(3))
			    ;
		    xj = max(d__1,d__2);
		    if (xj > 1.) {
			rec = 1. / xj;
/* Computing MAX */
			d__1 = WORK(j1), d__2 = WORK(j2);
			if (max(d__1,d__2) > (bignum - xmax) * rec) {
			    dscal_(&n2, &rec, &X(1), &c__1);
			    *scale *= rec;
			}
		    }

/*                 Update the right-hand side. */

		    if (j1 > 1) {
			i__1 = j1 - 1;
			d__1 = -X(j1);
			daxpy_(&i__1, &d__1, &T(1,j1), &c__1, &X(1)
				, &c__1);
			i__1 = j1 - 1;
			d__1 = -X(j2);
			daxpy_(&i__1, &d__1, &T(1,j2), &c__1, &X(1)
				, &c__1);

			i__1 = j1 - 1;
			d__1 = -X(*n + j1);
			daxpy_(&i__1, &d__1, &T(1,j1), &c__1, &X(*
				n + 1), &c__1);
			i__1 = j1 - 1;
			d__1 = -X(*n + j2);
			daxpy_(&i__1, &d__1, &T(1,j2), &c__1, &X(*
				n + 1), &c__1);

			X(1) = X(1) + B(j1) * X(*n + j1) + B(j2) * X(*n + j2);
			X(*n + 1) = X(*n + 1) - B(j1) * X(j1) - B(j2) * X(j2);

			xmax = 0.;
			i__1 = j1 - 1;
			for (k = 1; k <= j1-1; ++k) {
/* Computing MAX */
			    d__3 = (d__1 = X(k), abs(d__1)) + (d__2 = X(k + *
				    n), abs(d__2));
			    xmax = max(d__3,xmax);
/* L60: */
			}
		    }

		}
L70:
		;
	    }

	} else {

/*           Solve (T + iB)'*(p+iq) = c+id */

	    jnext = 1;
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		if (j < jnext) {
		    goto L80;
		}
		j1 = j;
		j2 = j;
		jnext = j + 1;
		if (j < *n && T(j+1,j) != 0.) {
		    j1 = j;
		    j2 = j + 1;
		    jnext = j + 2;
		}

		if (j1 == j2) {

/*                 1 by 1 diagonal block   

                   Scale if necessary to avoid overflow in
 forming the   
                   right-hand side element by inner produc
t. */

		    xj = (d__1 = X(j1), abs(d__1)) + (d__2 = X(j1 + *n), abs(
			    d__2));
		    if (xmax > 1.) {
			rec = 1. / xmax;
			if (WORK(j1) > (bignum - xj) * rec) {
			    dscal_(&n2, &rec, &X(1), &c__1);
			    *scale *= rec;
			    xmax *= rec;
			}
		    }

		    i__2 = j1 - 1;
		    X(j1) -= ddot_(&i__2, &T(1,j1), &c__1, &X(1), &
			    c__1);
		    i__2 = j1 - 1;
		    X(*n + j1) -= ddot_(&i__2, &T(1,j1), &c__1, &X(
			    *n + 1), &c__1);
		    if (j1 > 1) {
			X(j1) -= B(j1) * X(*n + 1);
			X(*n + j1) += B(j1) * X(1);
		    }
		    xj = (d__1 = X(j1), abs(d__1)) + (d__2 = X(j1 + *n), abs(
			    d__2));

		    z = *w;
		    if (j1 == 1) {
			z = B(1);
		    }

/*                 Scale if necessary to avoid overflow in
   
                   complex division */

		    tjj = (d__1 = T(j1,j1), abs(d__1)) + abs(z);
		    tmp = T(j1,j1);
		    if (tjj < sminw) {
			tmp = sminw;
			tjj = sminw;
			*info = 1;
		    }

		    if (tjj < 1.) {
			if (xj > bignum * tjj) {
			    rec = 1. / xj;
			    dscal_(&n2, &rec, &X(1), &c__1);
			    *scale *= rec;
			    xmax *= rec;
			}
		    }
		    d__1 = -z;
		    dladiv_(&X(j1), &X(*n + j1), &tmp, &d__1, &sr, &si);
		    X(j1) = sr;
		    X(j1 + *n) = si;
/* Computing MAX */
		    d__3 = (d__1 = X(j1), abs(d__1)) + (d__2 = X(j1 + *n), 
			    abs(d__2));
		    xmax = max(d__3,xmax);

		} else {

/*                 2 by 2 diagonal block   

                   Scale if necessary to avoid overflow in
 forming the   
                   right-hand side element by inner produc
t.   

   Computing MAX */
		    d__5 = (d__1 = X(j1), abs(d__1)) + (d__2 = X(*n + j1), 
			    abs(d__2)), d__6 = (d__3 = X(j2), abs(d__3)) + (
			    d__4 = X(*n + j2), abs(d__4));
		    xj = max(d__5,d__6);
		    if (xmax > 1.) {
			rec = 1. / xmax;
/* Computing MAX */
			d__1 = WORK(j1), d__2 = WORK(j2);
			if (max(d__1,d__2) > (bignum - xj) / xmax) {
			    dscal_(&n2, &rec, &X(1), &c__1);
			    *scale *= rec;
			    xmax *= rec;
			}
		    }

		    i__2 = j1 - 1;
		    d[0] = X(j1) - ddot_(&i__2, &T(1,j1), &c__1, &
			    X(1), &c__1);
		    i__2 = j1 - 1;
		    d[1] = X(j2) - ddot_(&i__2, &T(1,j2), &c__1, &
			    X(1), &c__1);
		    i__2 = j1 - 1;
		    d[2] = X(*n + j1) - ddot_(&i__2, &T(1,j1), &
			    c__1, &X(*n + 1), &c__1);
		    i__2 = j1 - 1;
		    d[3] = X(*n + j2) - ddot_(&i__2, &T(1,j2), &
			    c__1, &X(*n + 1), &c__1);
		    d[0] -= B(j1) * X(*n + 1);
		    d[1] -= B(j2) * X(*n + 1);
		    d[2] += B(j1) * X(1);
		    d[3] += B(j2) * X(1);

		    dlaln2_(&c_true, &c__2, &c__2, &sminw, &c_b21, &T(j1,j1), ldt, &c_b21, &c_b21, d, &c__2, &c_b25, 
			    w, v, &c__2, &scaloc, &xnorm, &ierr);
		    if (ierr != 0) {
			*info = 2;
		    }

		    if (scaloc != 1.) {
			dscal_(&n2, &scaloc, &X(1), &c__1);
			*scale = scaloc * *scale;
		    }
		    X(j1) = V(0);
		    X(j2) = V(1);
		    X(*n + j1) = V(2);
		    X(*n + j2) = V(3);
/* Computing MAX */
		    d__5 = (d__1 = X(j1), abs(d__1)) + (d__2 = X(*n + j1), 
			    abs(d__2)), d__6 = (d__3 = X(j2), abs(d__3)) + (
			    d__4 = X(*n + j2), abs(d__4)), d__5 = max(d__5,
			    d__6);
		    xmax = max(d__5,xmax);

		}

L80:
		;
	    }

	}

    }

    return 0;

/*     End of DLAQTR */

} /* dlaqtr_ */

