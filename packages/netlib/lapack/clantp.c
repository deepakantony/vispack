#include "f2c.h"

doublereal clantp_(char *norm, char *uplo, char *diag, integer *n, complex *
	ap, real *work)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    CLANTP  returns the value of the one norm,  or the Frobenius norm, or 
  
    the  infinity norm,  or the  element of  largest absolute value  of a 
  
    triangular matrix A, supplied in packed form.   

    Description   
    ===========   

    CLANTP returns the value   

       CLANTP = ( max(abs(A(i,j))), NORM = 'M' or 'm'   
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
            Specifies the value to be returned in CLANTP as described   
            above.   

    UPLO    (input) CHARACTER*1   
            Specifies whether the matrix A is upper or lower triangular. 
  
            = 'U':  Upper triangular   
            = 'L':  Lower triangular   

    DIAG    (input) CHARACTER*1   
            Specifies whether or not the matrix A is unit triangular.   
            = 'N':  Non-unit triangular   
            = 'U':  Unit triangular   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.  When N = 0, CLANTP is   
            set to zero.   

    AP      (input) COMPLEX array, dimension (N*(N+1)/2)   
            The upper or lower triangular matrix A, packed columnwise in 
  
            a linear array.  The j-th column of A is stored in the array 
  
            AP as follows:   
            if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;   
            if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.   
            Note that when DIAG = 'U', the elements of the array AP   
            corresponding to the diagonal elements of the matrix A are   
            not referenced, but are assumed to be one.   

    WORK    (workspace) REAL array, dimension (LWORK),   
            where LWORK >= N when NORM = 'I'; otherwise, WORK is not   
            referenced.   

   ===================================================================== 
  


    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__1 = 1;
    
    /* System generated locals */
    integer i__1, i__2;
    real ret_val, r__1, r__2;
    /* Builtin functions */
    double c_abs(complex *), sqrt(doublereal);
    /* Local variables */
    static integer i, j, k;
    static real scale;
    static logical udiag;
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

	k = 1;
	if (lsame_(diag, "U")) {
	    value = 1.f;
	    if (lsame_(uplo, "U")) {
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
		    i__2 = k + j - 2;
		    for (i = k; i <= k+j-2; ++i) {
/* Computing MAX */
			r__1 = value, r__2 = c_abs(&AP(i));
			value = dmax(r__1,r__2);
/* L10: */
		    }
		    k += j;
/* L20: */
		}
	    } else {
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
		    i__2 = k + *n - j;
		    for (i = k + 1; i <= k+*n-j; ++i) {
/* Computing MAX */
			r__1 = value, r__2 = c_abs(&AP(i));
			value = dmax(r__1,r__2);
/* L30: */
		    }
		    k = k + *n - j + 1;
/* L40: */
		}
	    }
	} else {
	    value = 0.f;
	    if (lsame_(uplo, "U")) {
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
		    i__2 = k + j - 1;
		    for (i = k; i <= k+j-1; ++i) {
/* Computing MAX */
			r__1 = value, r__2 = c_abs(&AP(i));
			value = dmax(r__1,r__2);
/* L50: */
		    }
		    k += j;
/* L60: */
		}
	    } else {
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
		    i__2 = k + *n - j;
		    for (i = k; i <= k+*n-j; ++i) {
/* Computing MAX */
			r__1 = value, r__2 = c_abs(&AP(i));
			value = dmax(r__1,r__2);
/* L70: */
		    }
		    k = k + *n - j + 1;
/* L80: */
		}
	    }
	}
    } else if (lsame_(norm, "O") || *(unsigned char *)norm == '1') {

/*        Find norm1(A). */

	value = 0.f;
	k = 1;
	udiag = lsame_(diag, "U");
	if (lsame_(uplo, "U")) {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		if (udiag) {
		    sum = 1.f;
		    i__2 = k + j - 2;
		    for (i = k; i <= k+j-2; ++i) {
			sum += c_abs(&AP(i));
/* L90: */
		    }
		} else {
		    sum = 0.f;
		    i__2 = k + j - 1;
		    for (i = k; i <= k+j-1; ++i) {
			sum += c_abs(&AP(i));
/* L100: */
		    }
		}
		k += j;
		value = dmax(value,sum);
/* L110: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		if (udiag) {
		    sum = 1.f;
		    i__2 = k + *n - j;
		    for (i = k + 1; i <= k+*n-j; ++i) {
			sum += c_abs(&AP(i));
/* L120: */
		    }
		} else {
		    sum = 0.f;
		    i__2 = k + *n - j;
		    for (i = k; i <= k+*n-j; ++i) {
			sum += c_abs(&AP(i));
/* L130: */
		    }
		}
		k = k + *n - j + 1;
		value = dmax(value,sum);
/* L140: */
	    }
	}
    } else if (lsame_(norm, "I")) {

/*        Find normI(A). */

	k = 1;
	if (lsame_(uplo, "U")) {
	    if (lsame_(diag, "U")) {
		i__1 = *n;
		for (i = 1; i <= *n; ++i) {
		    WORK(i) = 1.f;
/* L150: */
		}
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
		    i__2 = j - 1;
		    for (i = 1; i <= j-1; ++i) {
			WORK(i) += c_abs(&AP(k));
			++k;
/* L160: */
		    }
		    ++k;
/* L170: */
		}
	    } else {
		i__1 = *n;
		for (i = 1; i <= *n; ++i) {
		    WORK(i) = 0.f;
/* L180: */
		}
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
		    i__2 = j;
		    for (i = 1; i <= j; ++i) {
			WORK(i) += c_abs(&AP(k));
			++k;
/* L190: */
		    }
/* L200: */
		}
	    }
	} else {
	    if (lsame_(diag, "U")) {
		i__1 = *n;
		for (i = 1; i <= *n; ++i) {
		    WORK(i) = 1.f;
/* L210: */
		}
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
		    ++k;
		    i__2 = *n;
		    for (i = j + 1; i <= *n; ++i) {
			WORK(i) += c_abs(&AP(k));
			++k;
/* L220: */
		    }
/* L230: */
		}
	    } else {
		i__1 = *n;
		for (i = 1; i <= *n; ++i) {
		    WORK(i) = 0.f;
/* L240: */
		}
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
		    i__2 = *n;
		    for (i = j; i <= *n; ++i) {
			WORK(i) += c_abs(&AP(k));
			++k;
/* L250: */
		    }
/* L260: */
		}
	    }
	}
	value = 0.f;
	i__1 = *n;
	for (i = 1; i <= *n; ++i) {
/* Computing MAX */
	    r__1 = value, r__2 = WORK(i);
	    value = dmax(r__1,r__2);
/* L270: */
	}
    } else if (lsame_(norm, "F") || lsame_(norm, "E")) {

/*        Find normF(A). */

	if (lsame_(uplo, "U")) {
	    if (lsame_(diag, "U")) {
		scale = 1.f;
		sum = (real) (*n);
		k = 2;
		i__1 = *n;
		for (j = 2; j <= *n; ++j) {
		    i__2 = j - 1;
		    classq_(&i__2, &AP(k), &c__1, &scale, &sum);
		    k += j;
/* L280: */
		}
	    } else {
		scale = 0.f;
		sum = 1.f;
		k = 1;
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
		    classq_(&j, &AP(k), &c__1, &scale, &sum);
		    k += j;
/* L290: */
		}
	    }
	} else {
	    if (lsame_(diag, "U")) {
		scale = 1.f;
		sum = (real) (*n);
		k = 2;
		i__1 = *n - 1;
		for (j = 1; j <= *n-1; ++j) {
		    i__2 = *n - j;
		    classq_(&i__2, &AP(k), &c__1, &scale, &sum);
		    k = k + *n - j + 1;
/* L300: */
		}
	    } else {
		scale = 0.f;
		sum = 1.f;
		k = 1;
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
		    i__2 = *n - j + 1;
		    classq_(&i__2, &AP(k), &c__1, &scale, &sum);
		    k = k + *n - j + 1;
/* L310: */
		}
	    }
	}
	value = scale * sqrt(sum);
    }

    ret_val = value;
    return ret_val;

/*     End of CLANTP */

} /* clantp_ */

