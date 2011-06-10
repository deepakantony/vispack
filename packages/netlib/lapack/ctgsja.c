#include "f2c.h"

/* Subroutine */ int ctgsja_(char *jobu, char *jobv, char *jobq, integer *m, 
	integer *p, integer *n, integer *k, integer *l, complex *a, integer *
	lda, complex *b, integer *ldb, real *tola, real *tolb, real *alpha, 
	real *beta, complex *u, integer *ldu, complex *v, integer *ldv, 
	complex *q, integer *ldq, complex *work, integer *ncycle, integer *
	info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    CTGSJA computes the generalized singular value decomposition (GSVD)   
    of two complex upper triangular (or trapezoidal) matrices A and B.   

    On entry, it is assumed that matrices A and B have the following   
    forms, which may be obtained by the preprocessing subroutine CGGSVP   
    from a general M-by-N matrix A and P-by-N matrix B:   

                 N-K-L  K    L   
       A =    K ( 0    A12  A13 ) if M-K-L >= 0;   
              L ( 0     0   A23 )   
          M-K-L ( 0     0    0  )   

               N-K-L  K    L   
       A =  K ( 0    A12  A13 ) if M-K-L < 0;   
          M-K ( 0     0   A23 )   

               N-K-L  K    L   
       B =  L ( 0     0   B13 )   
          P-L ( 0     0    0  )   

    where the K-by-K matrix A12 and L-by-L matrix B13 are nonsingular   
    upper triangular; A23 is L-by-L upper triangular if M-K-L >= 0,   
    otherwise A23 is (M-K)-by-L upper trapezoidal.   

    On exit,   

           U'*A*Q = D1*( 0 R ),    V'*B*Q = D2*( 0 R ),   

    where U, V and Q are unitary matrices, Z' denotes the conjugate   
    transpose of Z, R is a nonsingular upper triangular matrix, and D1   
    and D2 are ``diagonal'' matrices, which are of the following   
    structures:   

    If M-K-L >= 0,   

                        K  L   
           D1 =     K ( I  0 )   
                    L ( 0  C )   
                M-K-L ( 0  0 )   

                       K  L   
           D2 = L   ( 0  S )   
                P-L ( 0  0 )   

                   N-K-L  K    L   
      ( 0 R ) = K (  0   R11  R12 ) K   
                L (  0    0   R22 ) L   

    where   

      C = diag( ALPHA(K+1), ... , ALPHA(K+L) ),   
      S = diag( BETA(K+1),  ... , BETA(K+L) ),   
      C**2 + S**2 = I.   

      R is stored in A(1:K+L,N-K-L+1:N) on exit.   

    If M-K-L < 0,   

                   K M-K K+L-M   
        D1 =   K ( I  0    0   )   
             M-K ( 0  C    0   )   

                     K M-K K+L-M   
        D2 =   M-K ( 0  S    0   )   
             K+L-M ( 0  0    I   )   
               P-L ( 0  0    0   )   

                   N-K-L  K   M-K  K+L-M   
   ( 0 R ) =    K ( 0    R11  R12  R13  )   
              M-K ( 0     0   R22  R23  )   
            K+L-M ( 0     0    0   R33  )   

    where   
    C = diag( ALPHA(K+1), ... , ALPHA(M) ),   
    S = diag( BETA(K+1),  ... , BETA(M) ),   
    C**2 + S**2 = I.   

    R = ( R11 R12 R13 ) is stored in A(1:M, N-K-L+1:N) and R33 is stored 
  
        (  0  R22 R23 )   
    in B(M-K+1:L,N+M-K-L+1:N) on exit.   

    The computation of the unitary transformation matrices U, V or Q   
    is optional.  These matrices may either be formed explicitly, or they 
  
    may be postmultiplied into input matrices U1, V1, or Q1.   

    Arguments   
    =========   

    JOBU    (input) CHARACTER*1   
            = 'U':  U must contain a unitary matrix U1 on entry, and   
                    the product U1*U is returned;   
            = 'I':  U is initialized to the unit matrix, and the   
                    unitary matrix U is returned;   
            = 'N':  U is not computed.   

    JOBV    (input) CHARACTER*1   
            = 'V':  V must contain a unitary matrix V1 on entry, and   
                    the product V1*V is returned;   
            = 'I':  V is initialized to the unit matrix, and the   
                    unitary matrix V is returned;   
            = 'N':  V is not computed.   

    JOBQ    (input) CHARACTER*1   
            = 'Q':  Q must contain a unitary matrix Q1 on entry, and   
                    the product Q1*Q is returned;   
            = 'I':  Q is initialized to the unit matrix, and the   
                    unitary matrix Q is returned;   
            = 'N':  Q is not computed.   

    M       (input) INTEGER   
            The number of rows of the matrix A.  M >= 0.   

    P       (input) INTEGER   
            The number of rows of the matrix B.  P >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrices A and B.  N >= 0.   

    K       (input) INTEGER   
    L       (input) INTEGER   
            K and L specify the subblocks in the input matrices A and B: 
  
            A23 = A(K+1:MIN(K+L,M),N-L+1:N) and B13 = B(1:L,,N-L+1:N)   
            of A and B, whose GSVD is going to be computed by CTGSJA.   
            See Further details.   

    A       (input/output) COMPLEX array, dimension (LDA,N)   
            On entry, the M-by-N matrix A.   
            On exit, A(N-K+1:N,1:MIN(K+L,M) ) contains the triangular   
            matrix R or part of R.  See Purpose for details.   

    LDA     (input) INTEGER   
            The leading dimension of the array A. LDA >= max(1,M).   

    B       (input/output) COMPLEX array, dimension (LDB,N)   
            On entry, the P-by-N matrix B.   
            On exit, if necessary, B(M-K+1:L,N+M-K-L+1:N) contains   
            a part of R.  See Purpose for details.   

    LDB     (input) INTEGER   
            The leading dimension of the array B. LDB >= max(1,P).   

    TOLA    (input) REAL   
    TOLB    (input) REAL   
            TOLA and TOLB are the convergence criteria for the Jacobi-   
            Kogbetliantz iteration procedure. Generally, they are the   
            same as used in the preprocessing step, say   
                TOLA = MAX(M,N)*norm(A)*MACHEPS,   
                TOLB = MAX(P,N)*norm(B)*MACHEPS.   

    ALPHA   (output) REAL array, dimension (N)   
    BETA    (output) REAL array, dimension (N)   
            On exit, ALPHA and BETA contain the generalized singular   
            value pairs of A and B;   
              ALPHA(1:K) = 1,   
              BETA(1:K)  = 0,   
            and if M-K-L >= 0,   
              ALPHA(K+1:K+L) = diag(C),   
              BETA(K+1:K+L)  = diag(S),   
            or if M-K-L < 0,   
              ALPHA(K+1:M)= C, ALPHA(M+1:K+L)= 0   
              BETA(K+1:M) = S, BETA(M+1:K+L) = 1.   
            Furthermore, if K+L < N,   
              ALPHA(K+L+1:N) = 0   
              BETA(K+L+1:N)  = 0.   

    U       (input/output) COMPLEX array, dimension (LDU,M)   
            On entry, if JOBU = 'U', U must contain a matrix U1 (usually 
  
            the unitary matrix returned by CGGSVP).   
            On exit,   
            if JOBU = 'I', U contains the unitary matrix U;   
            if JOBU = 'U', U contains the product U1*U.   
            If JOBU = 'N', U is not referenced.   

    LDU     (input) INTEGER   
            The leading dimension of the array U. LDU >= max(1,M) if   
            JOBU = 'U'; LDU >= 1 otherwise.   

    V       (input/output) COMPLEX array, dimension (LDV,P)   
            On entry, if JOBV = 'V', V must contain a matrix V1 (usually 
  
            the unitary matrix returned by CGGSVP).   
            On exit,   
            if JOBV = 'I', V contains the unitary matrix V;   
            if JOBV = 'V', V contains the product V1*V.   
            If JOBV = 'N', V is not referenced.   

    LDV     (input) INTEGER   
            The leading dimension of the array V. LDV >= max(1,P) if   
            JOBV = 'V'; LDV >= 1 otherwise.   

    Q       (input/output) COMPLEX array, dimension (LDQ,N)   
            On entry, if JOBQ = 'Q', Q must contain a matrix Q1 (usually 
  
            the unitary matrix returned by CGGSVP).   
            On exit,   
            if JOBQ = 'I', Q contains the unitary matrix Q;   
            if JOBQ = 'Q', Q contains the product Q1*Q.   
            If JOBQ = 'N', Q is not referenced.   

    LDQ     (input) INTEGER   
            The leading dimension of the array Q. LDQ >= max(1,N) if   
            JOBQ = 'Q'; LDQ >= 1 otherwise.   

    WORK    (workspace) COMPLEX array, dimension (2*N)   

    NCYCLE  (output) INTEGER   
            The number of cycles required for convergence.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value.   
            = 1:  the procedure does not converge after MAXIT cycles.   

    Internal Parameters   
    ===================   

    MAXIT   INTEGER   
            MAXIT specifies the total loops that the iterative procedure 
  
            may take. If after MAXIT cycles, the routine fails to   
            converge, we return INFO = 1.   

    Further Details   
    ===============   

    CTGSJA essentially uses a variant of Kogbetliantz algorithm to reduce 
  
    min(L,M-K)-by-L triangular (or trapezoidal) matrix A23 and L-by-L   
    matrix B13 to the form:   

             U1'*A13*Q1 = C1*R1; V1'*B13*Q1 = S1*R1,   

    where U1, V1 and Q1 are unitary matrix, and Z' is the conjugate   
    transpose of Z.  C1 and S1 are diagonal matrices satisfying   

                  C1**2 + S1**2 = I,   

    and R1 is an L-by-L nonsingular upper triangular matrix.   

    ===================================================================== 
  



       Decode and test the input parameters   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static complex c_b1 = {0.f,0.f};
    static complex c_b2 = {1.f,0.f};
    static integer c__1 = 1;
    static real c_b39 = -1.f;
    static real c_b42 = 1.f;
    
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, q_dim1, q_offset, u_dim1, 
	    u_offset, v_dim1, v_offset, i__1, i__2, i__3, i__4;
    real r__1;
    doublereal d__1;
    complex q__1;
    /* Builtin functions */
    void r_cnjg(complex *, complex *);
    /* Local variables */
    extern /* Subroutine */ int crot_(integer *, complex *, integer *, 
	    complex *, integer *, real *, complex *);
    static integer i, j;
    static real gamma;
    extern logical lsame_(char *, char *);
    extern /* Subroutine */ int ccopy_(integer *, complex *, integer *, 
	    complex *, integer *);
    static logical initq;
    static real a1, a3, b1;
    static logical initu, initv, wantq, upper;
    static real b3, error;
    static logical wantu, wantv;
    static real ssmin;
    static complex a2, b2;
    extern /* Subroutine */ int clags2_(logical *, real *, complex *, real *, 
	    real *, complex *, real *, real *, complex *, real *, complex *, 
	    real *, complex *), clapll_(integer *, complex *, integer *, 
	    complex *, integer *, real *), csscal_(integer *, real *, complex 
	    *, integer *);
    static integer kcycle;
    extern /* Subroutine */ int claset_(char *, integer *, integer *, complex 
	    *, complex *, complex *, integer *), xerbla_(char *, 
	    integer *), slartg_(real *, real *, real *, real *, real *
	    );
    static real csq, csu, csv;
    static complex snq;
    static real rwk;
    static complex snu, snv;



#define ALPHA(I) alpha[(I)-1]
#define BETA(I) beta[(I)-1]
#define WORK(I) work[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]
#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]
#define U(I,J) u[(I)-1 + ((J)-1)* ( *ldu)]
#define V(I,J) v[(I)-1 + ((J)-1)* ( *ldv)]
#define Q(I,J) q[(I)-1 + ((J)-1)* ( *ldq)]

    initu = lsame_(jobu, "I");
    wantu = initu || lsame_(jobu, "U");

    initv = lsame_(jobv, "I");
    wantv = initv || lsame_(jobv, "V");

    initq = lsame_(jobq, "I");
    wantq = initq || lsame_(jobq, "Q");

    *info = 0;
    if (! (initu || wantu || lsame_(jobu, "N"))) {
	*info = -1;
    } else if (! (initv || wantv || lsame_(jobv, "N"))) {
	*info = -2;
    } else if (! (initq || wantq || lsame_(jobq, "N"))) {
	*info = -3;
    } else if (*m < 0) {
	*info = -4;
    } else if (*p < 0) {
	*info = -5;
    } else if (*n < 0) {
	*info = -6;
    } else if (*lda < max(1,*m)) {
	*info = -10;
    } else if (*ldb < max(1,*p)) {
	*info = -12;
    } else if (*ldu < 1 || wantu && *ldu < *m) {
	*info = -18;
    } else if (*ldv < 1 || wantv && *ldv < *p) {
	*info = -20;
    } else if (*ldq < 1 || wantq && *ldq < *n) {
	*info = -22;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("CTGSJA", &i__1);
	return 0;
    }

/*     Initialize U, V and Q, if necessary */

    if (initu) {
	claset_("Full", m, m, &c_b1, &c_b2, &U(1,1), ldu);
    }
    if (initv) {
	claset_("Full", p, p, &c_b1, &c_b2, &V(1,1), ldv);
    }
    if (initq) {
	claset_("Full", n, n, &c_b1, &c_b2, &Q(1,1), ldq);
    }

/*     Loop until convergence */

    upper = FALSE_;
    for (kcycle = 1; kcycle <= 40; ++kcycle) {

	upper = ! upper;

	i__1 = *l - 1;
	for (i = 1; i <= *l-1; ++i) {
	    i__2 = *l;
	    for (j = i + 1; j <= *l; ++j) {

		a1 = 0.f;
		a2.r = 0.f, a2.i = 0.f;
		a3 = 0.f;
		if (*k + i <= *m) {
		    i__3 = *k + i + (*n - *l + i) * a_dim1;
		    a1 = A(*k+i,*n-*l+i).r;
		}
		if (*k + j <= *m) {
		    i__3 = *k + j + (*n - *l + j) * a_dim1;
		    a3 = A(*k+j,*n-*l+j).r;
		}

		i__3 = i + (*n - *l + i) * b_dim1;
		b1 = B(i,*n-*l+i).r;
		i__3 = j + (*n - *l + j) * b_dim1;
		b3 = B(j,*n-*l+j).r;

		if (upper) {
		    if (*k + i <= *m) {
			i__3 = *k + i + (*n - *l + j) * a_dim1;
			a2.r = A(*k+i,*n-*l+j).r, a2.i = A(*k+i,*n-*l+j).i;
		    }
		    i__3 = i + (*n - *l + j) * b_dim1;
		    b2.r = B(i,*n-*l+j).r, b2.i = B(i,*n-*l+j).i;
		} else {
		    if (*k + j <= *m) {
			i__3 = *k + j + (*n - *l + i) * a_dim1;
			a2.r = A(*k+j,*n-*l+i).r, a2.i = A(*k+j,*n-*l+i).i;
		    }
		    i__3 = j + (*n - *l + i) * b_dim1;
		    b2.r = B(j,*n-*l+i).r, b2.i = B(j,*n-*l+i).i;
		}

		clags2_(&upper, &a1, &a2, &a3, &b1, &b2, &b3, &csu, &snu, &
			csv, &snv, &csq, &snq);

/*              Update (K+I)-th and (K+J)-th rows of matrix A:
 U'*A */

		if (*k + j <= *m) {
		    r_cnjg(&q__1, &snu);
		    crot_(l, &A(*k+j,*n-*l+1), lda, &A(*k+i,*n-*l+1), lda, &csu, &q__1);
		}

/*              Update I-th and J-th rows of matrix B: V'*B */

		r_cnjg(&q__1, &snv);
		crot_(l, &B(j,*n-*l+1), ldb, &B(i,*n-*l+1), ldb, &csv, &q__1);

/*              Update (N-L+I)-th and (N-L+J)-th columns of ma
trices   
                A and B: A*Q and B*Q   

   Computing MIN */
		i__4 = *k + *l;
		i__3 = min(i__4,*m);
		crot_(&i__3, &A(1,*n-*l+j), &c__1, &A(1,*n-*l+i), &c__1, &csq, &snq);

		crot_(l, &B(1,*n-*l+j), &c__1, &B(1,*n-*l+i), &c__1, &csq, &snq);

		if (upper) {
		    if (*k + i <= *m) {
			i__3 = *k + i + (*n - *l + j) * a_dim1;
			A(*k+i,*n-*l+j).r = 0.f, A(*k+i,*n-*l+j).i = 0.f;
		    }
		    i__3 = i + (*n - *l + j) * b_dim1;
		    B(i,*n-*l+j).r = 0.f, B(i,*n-*l+j).i = 0.f;
		} else {
		    if (*k + j <= *m) {
			i__3 = *k + j + (*n - *l + i) * a_dim1;
			A(*k+j,*n-*l+i).r = 0.f, A(*k+j,*n-*l+i).i = 0.f;
		    }
		    i__3 = j + (*n - *l + i) * b_dim1;
		    B(j,*n-*l+i).r = 0.f, B(j,*n-*l+i).i = 0.f;
		}

/*              Ensure that the diagonal elements of A and B a
re real. */

		if (*k + i <= *m) {
		    i__3 = *k + i + (*n - *l + i) * a_dim1;
		    i__4 = *k + i + (*n - *l + i) * a_dim1;
		    d__1 = A(*k+i,*n-*l+i).r;
		    A(*k+i,*n-*l+i).r = d__1, A(*k+i,*n-*l+i).i = 0.f;
		}
		if (*k + j <= *m) {
		    i__3 = *k + j + (*n - *l + j) * a_dim1;
		    i__4 = *k + j + (*n - *l + j) * a_dim1;
		    d__1 = A(*k+j,*n-*l+j).r;
		    A(*k+j,*n-*l+j).r = d__1, A(*k+j,*n-*l+j).i = 0.f;
		}
		i__3 = i + (*n - *l + i) * b_dim1;
		i__4 = i + (*n - *l + i) * b_dim1;
		d__1 = B(i,*n-*l+i).r;
		B(i,*n-*l+i).r = d__1, B(i,*n-*l+i).i = 0.f;
		i__3 = j + (*n - *l + j) * b_dim1;
		i__4 = j + (*n - *l + j) * b_dim1;
		d__1 = B(j,*n-*l+j).r;
		B(j,*n-*l+j).r = d__1, B(j,*n-*l+j).i = 0.f;

/*              Update unitary matrices U, V, Q, if desired. 
*/

		if (wantu && *k + j <= *m) {
		    crot_(m, &U(1,*k+j), &c__1, &U(1,*k+i), &c__1, &csu, &snu);
		}

		if (wantv) {
		    crot_(p, &V(1,j), &c__1, &V(1,i), &
			    c__1, &csv, &snv);
		}

		if (wantq) {
		    crot_(n, &Q(1,*n-*l+j), &c__1, &Q(1,*n-*l+i), &c__1, &csq, &snq);
		}

/* L10: */
	    }
/* L20: */
	}

	if (! upper) {

/*           The matrices A13 and B13 were lower triangular at the
 start   
             of the cycle, and are now upper triangular.   

             Convergence test: test the parallelism of the corresp
onding   
             rows of A and B. */

	    error = 0.f;
/* Computing MIN */
	    i__2 = *l, i__3 = *m - *k;
	    i__1 = min(i__2,i__3);
	    for (i = 1; i <= min(*l,*m-*k); ++i) {
		i__2 = *l - i + 1;
		ccopy_(&i__2, &A(*k+i,*n-*l+i), lda, &WORK(
			1), &c__1);
		i__2 = *l - i + 1;
		ccopy_(&i__2, &B(i,*n-*l+i), ldb, &WORK(*l + 
			1), &c__1);
		i__2 = *l - i + 1;
		clapll_(&i__2, &WORK(1), &c__1, &WORK(*l + 1), &c__1, &ssmin);
		error = dmax(error,ssmin);
/* L30: */
	    }

	    if (dabs(error) <= (real) (*n) * dmin(*tola,*tolb)) {
		goto L50;
	    }
	}

/*        End of cycle loop   

   L40: */
    }

/*     The algorithm has not converged after MAXIT cycles. */

    *info = 1;
    goto L90;

L50:

/*     If ERROR <= N*MIN(TOLA,TOLB), then the algorithm has converged.   
       Compute the generalized singular value pairs (ALPHA, BETA), and   
       set the triangular matrix R to array A. */

    i__1 = *k;
    for (i = 1; i <= *k; ++i) {
	ALPHA(i) = 1.f;
	BETA(i) = 0.f;
/* L60: */
    }

/* Computing MIN */
    i__2 = *l, i__3 = *m - *k;
    i__1 = min(i__2,i__3);
    for (i = 1; i <= min(*l,*m-*k); ++i) {

	i__2 = *k + i + (*n - *l + i) * a_dim1;
	a1 = A(*k+i,*n-*l+i).r;
	i__2 = i + (*n - *l + i) * b_dim1;
	b1 = B(i,*n-*l+i).r;

	if (a1 != 0.f) {
	    gamma = b1 / a1;

	    if (gamma < 0.f) {
		i__2 = *l - i + 1;
		csscal_(&i__2, &c_b39, &B(i,*n-*l+i), ldb);
		if (wantv) {
		    csscal_(p, &c_b39, &V(1,i), &c__1);
		}
	    }

	    r__1 = dabs(gamma);
	    slartg_(&r__1, &c_b42, &BETA(*k + i), &ALPHA(*k + i), &rwk);

	    if (ALPHA(*k + i) >= BETA(*k + i)) {
		i__2 = *l - i + 1;
		r__1 = 1.f / ALPHA(*k + i);
		csscal_(&i__2, &r__1, &A(*k+i,*n-*l+i), 
			lda);
	    } else {
		i__2 = *l - i + 1;
		r__1 = 1.f / BETA(*k + i);
		csscal_(&i__2, &r__1, &B(i,*n-*l+i), ldb);
		i__2 = *l - i + 1;
		ccopy_(&i__2, &B(i,*n-*l+i), ldb, &A(*k+i,*n-*l+i), lda);
	    }

	} else {
	    ALPHA(*k + i) = 0.f;
	    BETA(*k + i) = 1.f;
	    i__2 = *l - i + 1;
	    ccopy_(&i__2, &B(i,*n-*l+i), ldb, &A(*k+i,*n-*l+i), lda);
	}
/* L70: */
    }

/*     Post-assignment */

    i__1 = *k + *l;
    for (i = *m + 1; i <= *k+*l; ++i) {
	ALPHA(i) = 0.f;
	BETA(i) = 1.f;
/* L80: */
    }

    if (*k + *l < *n) {
	i__1 = *n;
	for (i = *k + *l + 1; i <= *n; ++i) {
	    ALPHA(i) = 0.f;
	    BETA(i) = 0.f;
/* L85: */
	}
    }

L90:
    *ncycle = kcycle;

    return 0;

/*     End of CTGSJA */

} /* ctgsja_ */

