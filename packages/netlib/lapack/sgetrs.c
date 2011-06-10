#include "f2c.h"

/* Subroutine */ int sgetrs_(char *trans, integer *n, integer *nrhs, real *a, 
	integer *lda, integer *ipiv, real *b, integer *ldb, integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       March 31, 1993   


    Purpose   
    =======   

    SGETRS solves a system of linear equations   
       A * X = B  or  A' * X = B   
    with a general N-by-N matrix A using the LU factorization computed   
    by SGETRF.   

    Arguments   
    =========   

    TRANS   (input) CHARACTER*1   
            Specifies the form of the system of equations:   
            = 'N':  A * X = B  (No transpose)   
            = 'T':  A'* X = B  (Transpose)   
            = 'C':  A'* X = B  (Conjugate transpose = Transpose)   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    NRHS    (input) INTEGER   
            The number of right hand sides, i.e., the number of columns   
            of the matrix B.  NRHS >= 0.   

    A       (input) REAL array, dimension (LDA,N)   
            The factors L and U from the factorization A = P*L*U   
            as computed by SGETRF.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,N).   

    IPIV    (input) INTEGER array, dimension (N)   
            The pivot indices from SGETRF; for 1<=i<=N, row i of the   
            matrix was interchanged with row IPIV(i).   

    B       (input/output) REAL array, dimension (LDB,NRHS)   
            On entry, the right hand side matrix B.   
            On exit, the solution matrix X.   

    LDB     (input) INTEGER   
            The leading dimension of the array B.  LDB >= max(1,N).   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__1 = 1;
    static real c_b12 = 1.f;
    static integer c_n1 = -1;
    
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1;
    /* Local variables */
    extern logical lsame_(char *, char *);
    extern /* Subroutine */ int strsm_(char *, char *, char *, char *, 
	    integer *, integer *, real *, real *, integer *, real *, integer *
	    ), xerbla_(char *, integer *);
    static logical notran;
    extern /* Subroutine */ int slaswp_(integer *, real *, integer *, integer 
	    *, integer *, integer *, integer *);



#define IPIV(I) ipiv[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]
#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]

    *info = 0;
    notran = lsame_(trans, "N");
    if (! notran && ! lsame_(trans, "T") && ! lsame_(trans, "C")) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*nrhs < 0) {
	*info = -3;
    } else if (*lda < max(1,*n)) {
	*info = -5;
    } else if (*ldb < max(1,*n)) {
	*info = -8;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("SGETRS", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0 || *nrhs == 0) {
	return 0;
    }

    if (notran) {

/*        Solve A * X = B.   

          Apply row interchanges to the right hand sides. */

	slaswp_(nrhs, &B(1,1), ldb, &c__1, n, &IPIV(1), &c__1);

/*        Solve L*X = B, overwriting B with X. */

	strsm_("Left", "Lower", "No transpose", "Unit", n, nrhs, &c_b12, &A(1,1), lda, &B(1,1), ldb);

/*        Solve U*X = B, overwriting B with X. */

	strsm_("Left", "Upper", "No transpose", "Non-unit", n, nrhs, &c_b12, &
		A(1,1), lda, &B(1,1), ldb);
    } else {

/*        Solve A' * X = B.   

          Solve U'*X = B, overwriting B with X. */

	strsm_("Left", "Upper", "Transpose", "Non-unit", n, nrhs, &c_b12, &A(1,1), lda, &B(1,1), ldb);

/*        Solve L'*X = B, overwriting B with X. */

	strsm_("Left", "Lower", "Transpose", "Unit", n, nrhs, &c_b12, &A(1,1), lda, &B(1,1), ldb);

/*        Apply row interchanges to the solution vectors. */

	slaswp_(nrhs, &B(1,1), ldb, &c__1, n, &IPIV(1), &c_n1);
    }

    return 0;

/*     End of SGETRS */

} /* sgetrs_ */
