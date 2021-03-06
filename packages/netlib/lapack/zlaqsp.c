#include "f2c.h"

/* Subroutine */ int zlaqsp_(char *uplo, integer *n, doublecomplex *ap, 
	doublereal *s, doublereal *scond, doublereal *amax, char *equed)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    ZLAQSP equilibrates a symmetric matrix A using the scaling factors   
    in the vector S.   

    Arguments   
    =========   

    UPLO    (input) CHARACTER*1   
            Specifies whether the upper or lower triangular part of the   
            symmetric matrix A is stored.   
            = 'U':  Upper triangular   
            = 'L':  Lower triangular   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    AP      (input/output) COMPLEX*16 array, dimension (N*(N+1)/2)   
            On entry, the upper or lower triangle of the symmetric matrix 
  
            A, packed columnwise in a linear array.  The j-th column of A 
  
            is stored in the array AP as follows:   
            if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;   
            if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.   

            On exit, the equilibrated matrix:  diag(S) * A * diag(S), in 
  
            the same storage format as A.   

    S       (input) DOUBLE PRECISION array, dimension (N)   
            The scale factors for A.   

    SCOND   (input) DOUBLE PRECISION   
            Ratio of the smallest S(i) to the largest S(i).   

    AMAX    (input) DOUBLE PRECISION   
            Absolute value of largest matrix entry.   

    EQUED   (output) CHARACTER*1   
            Specifies whether or not equilibration was done.   
            = 'N':  No equilibration.   
            = 'Y':  Equilibration was done, i.e., A has been replaced by 
  
                    diag(S) * A * diag(S).   

    Internal Parameters   
    ===================   

    THRESH is a threshold value used to decide if scaling should be done 
  
    based on the ratio of the scaling factors.  If SCOND < THRESH,   
    scaling is done.   

    LARGE and SMALL are threshold values used to decide if scaling should 
  
    be done based on the absolute size of the largest matrix element.   
    If AMAX > LARGE or AMAX < SMALL, scaling is done.   

    ===================================================================== 
  


       Quick return if possible   

    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    doublereal d__1;
    doublecomplex z__1;
    /* Local variables */
    static integer i, j;
    static doublereal large;
    extern logical lsame_(char *, char *);
    static doublereal small;
    static integer jc;
    static doublereal cj;
    extern doublereal dlamch_(char *);


#define S(I) s[(I)-1]
#define AP(I) ap[(I)-1]


    if (*n <= 0) {
	*(unsigned char *)equed = 'N';
	return 0;
    }

/*     Initialize LARGE and SMALL. */

    small = dlamch_("Safe minimum") / dlamch_("Precision");
    large = 1. / small;

    if (*scond >= .1 && *amax >= small && *amax <= large) {

/*        No equilibration */

	*(unsigned char *)equed = 'N';
    } else {

/*        Replace A by diag(S) * A * diag(S). */

	if (lsame_(uplo, "U")) {

/*           Upper triangle of A is stored. */

	    jc = 1;
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		cj = S(j);
		i__2 = j;
		for (i = 1; i <= j; ++i) {
		    i__3 = jc + i - 1;
		    d__1 = cj * S(i);
		    i__4 = jc + i - 1;
		    z__1.r = d__1 * AP(jc+i-1).r, z__1.i = d__1 * AP(jc+i-1).i;
		    AP(jc+i-1).r = z__1.r, AP(jc+i-1).i = z__1.i;
/* L10: */
		}
		jc += j;
/* L20: */
	    }
	} else {

/*           Lower triangle of A is stored. */

	    jc = 1;
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		cj = S(j);
		i__2 = *n;
		for (i = j; i <= *n; ++i) {
		    i__3 = jc + i - j;
		    d__1 = cj * S(i);
		    i__4 = jc + i - j;
		    z__1.r = d__1 * AP(jc+i-j).r, z__1.i = d__1 * AP(jc+i-j).i;
		    AP(jc+i-j).r = z__1.r, AP(jc+i-j).i = z__1.i;
/* L30: */
		}
		jc = jc + *n - j + 1;
/* L40: */
	    }
	}
	*(unsigned char *)equed = 'Y';
    }

    return 0;

/*     End of ZLAQSP */

} /* zlaqsp_ */

