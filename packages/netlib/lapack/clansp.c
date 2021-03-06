#include "f2c.h"

doublereal clansp_(char *norm, char *uplo, integer *n, complex *ap, real *
	work)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    CLANSP  returns the value of the one norm,  or the Frobenius norm, or 
  
    the  infinity norm,  or the  element of  largest absolute value  of a 
  
    complex symmetric matrix A,  supplied in packed form.   

    Description   
    ===========   

    CLANSP returns the value   

       CLANSP = ( max(abs(A(i,j))), NORM = 'M' or 'm'   
                (   
                ( norm1(A),         NORM = '1', 'O' or 'o'   
                (   
                ( normI(A),         NORM = 'I' or 'i'   
                (   
                ( normF(A),         NORM = 'F', 'f', 'E' or 'e'   

    where  norm1  denotes the  one norm of a matrix (maximum column sum), 
  
    normI  denotes the  infinity norm  of a matrix  (maximum row sum) and 
  
    normF  denotes the  Frobenius norm of a matrix (square root of sum of 
  
    squares).  Note that  max(abs(A(i,j)))  is not a  matrix norm.   

    Arguments   
    =========   

    NORM    (input) CHARACTER*1   
            Specifies the value to be returned in CLANSP as described   
            above.   

    UPLO    (input) CHARACTER*1   
            Specifies whether the upper or lower triangular part of the   
            symmetric matrix A is supplied.   
            = 'U':  Upper triangular part of A is supplied   
            = 'L':  Lower triangular part of A is supplied   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.  When N = 0, CLANSP is   
            set to zero.   

    AP      (input) COMPLEX array, dimension (N*(N+1)/2)   
            The upper or lower triangle of the symmetric matrix A, packed 
  
            columnwise in a linear array.  The j-th column of A is stored 
  
            in the array AP as follows:   
            if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;   
            if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.   

    WORK    (workspace) REAL array, dimension (LWORK),   
            where LWORK >= N when NORM = 'I' or '1' or 'O'; otherwise,   
            WORK is not referenced.   

   ===================================================================== 
  


    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__1 = 1;
    
    /* System generated locals */
    integer i__1, i__2;
    real ret_val, r__1, r__2;
    /* Builtin functions */
    double c_abs(complex *), r_imag(complex *), sqrt(doublereal);
    /* Local variables */
    static real absa;
    static integer i, j, k;
    static real scale;
    extern logical lsame_(char *, char *);
    static real value;
    extern /* Subroutine */ int classq_(integer *, complex *, integer *, real 
	    *, real *);
    static real sum;



#define WORK(I) work[(I)-1]
#define AP(I) ap[(I)-1]


    if (*n == 0) {
	value = 0.f;
    } else if (lsame_(norm, "M")) {

/*        Find max(abs(A(i,j))). */

	value = 0.f;
	if (lsame_(uplo, "U")) {
	    k = 1;
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		i__2 = k + j - 1;
		for (i = k; i <= k+j-1; ++i) {
/* Computing MAX */
		    r__1 = value, r__2 = c_abs(&AP(i));
		    value = dmax(r__1,r__2);
/* L10: */
		}
		k += j;
/* L20: */
	    }
	} else {
	    k = 1;
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		i__2 = k + *n - j;
		for (i = k; i <= k+*n-j; ++i) {
/* Computing MAX */
		    r__1 = value, r__2 = c_abs(&AP(i));
		    value = dmax(r__1,r__2);
/* L30: */
		}
		k = k + *n - j + 1;
/* L40: */
	    }
	}
    } else if (lsame_(norm, "I") || lsame_(norm, "O") || *(
	    unsigned char *)norm == '1') {

/*        Find normI(A) ( = norm1(A), since A is symmetric). */

	value = 0.f;
	k = 1;
	if (lsame_(uplo, "U")) {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		sum = 0.f;
		i__2 = j - 1;
		for (i = 1; i <= j-1; ++i) {
		    absa = c_abs(&AP(k));
		    sum += absa;
		    WORK(i) += absa;
		    ++k;
/* L50: */
		}
		WORK(j) = sum + c_abs(&AP(k));
		++k;
/* L60: */
	    }
	    i__1 = *n;
	    for (i = 1; i <= *n; ++i) {
/* Computing MAX */
		r__1 = value, r__2 = WORK(i);
		value = dmax(r__1,r__2);
/* L70: */
	    }
	} else {
	    i__1 = *n;
	    for (i = 1; i <= *n; ++i) {
		WORK(i) = 0.f;
/* L80: */
	    }
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		sum = WORK(j) + c_abs(&AP(k));
		++k;
		i__2 = *n;
		for (i = j + 1; i <= *n; ++i) {
		    absa = c_abs(&AP(k));
		    sum += absa;
		    WORK(i) += absa;
		    ++k;
/* L90: */
		}
		value = dmax(value,sum);
/* L100: */
	    }
	}
    } else if (lsame_(norm, "F") || lsame_(norm, "E")) {

/*        Find normF(A). */

	scale = 0.f;
	sum = 1.f;
	k = 2;
	if (lsame_(uplo, "U")) {
	    i__1 = *n;
	    for (j = 2; j <= *n; ++j) {
		i__2 = j - 1;
		classq_(&i__2, &AP(k), &c__1, &scale, &sum);
		k += j;
/* L110: */
	    }
	} else {
	    i__1 = *n - 1;
	    for (j = 1; j <= *n-1; ++j) {
		i__2 = *n - j;
		classq_(&i__2, &AP(k), &c__1, &scale, &sum);
		k = k + *n - j + 1;
/* L120: */
	    }
	}
	sum *= 2;
	k = 1;
	i__1 = *n;
	for (i = 1; i <= *n; ++i) {
	    i__2 = k;
	    if (AP(k).r != 0.f) {
		i__2 = k;
		absa = (r__1 = AP(k).r, dabs(r__1));
		if (scale < absa) {
/* Computing 2nd power */
		    r__1 = scale / absa;
		    sum = sum * (r__1 * r__1) + 1.f;
		    scale = absa;
		} else {
/* Computing 2nd power */
		    r__1 = absa / scale;
		    sum += r__1 * r__1;
		}
	    }
	    if (r_imag(&AP(k)) != 0.f) {
		absa = (r__1 = r_imag(&AP(k)), dabs(r__1));
		if (scale < absa) {
/* Computing 2nd power */
		    r__1 = scale / absa;
		    sum = sum * (r__1 * r__1) + 1.f;
		    scale = absa;
		} else {
/* Computing 2nd power */
		    r__1 = absa / scale;
		    sum += r__1 * r__1;
		}
	    }
	    if (lsame_(uplo, "U")) {
		k = k + i + 1;
	    } else {
		k = k + *n - i + 1;
	    }
/* L130: */
	}
	value = scale * sqrt(sum);
    }

    ret_val = value;
    return ret_val;

/*     End of CLANSP */

} /* clansp_ */

