#include "f2c.h"

doublereal zlantb_(char *norm, char *uplo, char *diag, integer *n, integer *k,
	 doublecomplex *ab, integer *ldab, doublereal *work)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    ZLANTB  returns the value of the one norm,  or the Frobenius norm, or 
  
    the  infinity norm,  or the element of  largest absolute value  of an 
  
    n by n triangular band matrix A,  with ( k + 1 ) diagonals.   

    Description   
    ===========   

    ZLANTB returns the value   

       ZLANTB = ( max(abs(A(i,j))), NORM = 'M' or 'm'   
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
            Specifies the value to be returned in ZLANTB as described   
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
            The order of the matrix A.  N >= 0.  When N = 0, ZLANTB is   
            set to zero.   

    K       (input) INTEGER   
            The number of super-diagonals of the matrix A if UPLO = 'U', 
  
            or the number of sub-diagonals of the matrix A if UPLO = 'L'. 
  
            K >= 0.   

    AB      (input) COMPLEX*16 array, dimension (LDAB,N)   
            The upper or lower triangular band matrix A, stored in the   
            first k+1 rows of AB.  The j-th column of A is stored   
            in the j-th column of the array AB as follows:   
            if UPLO = 'U', AB(k+1+i-j,j) = A(i,j) for max(1,j-k)<=i<=j;   
            if UPLO = 'L', AB(1+i-j,j)   = A(i,j) for j<=i<=min(n,j+k).   
            Note that when DIAG = 'U', the elements of the array AB   
            corresponding to the diagonal elements of the matrix A are   
            not referenced, but are assumed to be one.   

    LDAB    (input) INTEGER   
            The leading dimension of the array AB.  LDAB >= K+1.   

    WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK),   
            where LWORK >= N when NORM = 'I'; otherwise, WORK is not   
            referenced.   

   ===================================================================== 
  


    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__1 = 1;
    
    /* System generated locals */
    integer ab_dim1, ab_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal ret_val, d__1, d__2;
    /* Builtin functions */
    double z_abs(doublecomplex *), sqrt(doublereal);
    /* Local variables */
    static integer i, j, l;
    static doublereal scale;
    static logical udiag;
    extern logical lsame_(char *, char *);
    static doublereal value;
    extern /* Subroutine */ int zlassq_(integer *, doublecomplex *, integer *,
	     doublereal *, doublereal *);
    static doublereal sum;



#define WORK(I) work[(I)-1]

#define AB(I,J) ab[(I)-1 + ((J)-1)* ( *ldab)]

    if (*n == 0) {
	value = 0.;
    } else if (lsame_(norm, "M")) {

/*        Find max(abs(A(i,j))). */

	if (lsame_(diag, "U")) {
	    value = 1.;
	    if (lsame_(uplo, "U")) {
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
/* Computing MAX */
		    i__2 = *k + 2 - j;
		    i__3 = *k;
		    for (i = max(*k+2-j,1); i <= *k; ++i) {
/* Computing MAX */
			d__1 = value, d__2 = z_abs(&AB(i,j));
			value = max(d__1,d__2);
/* L10: */
		    }
/* L20: */
		}
	    } else {
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
/* Computing MIN */
		    i__2 = *n + 1 - j, i__4 = *k + 1;
		    i__3 = min(i__2,i__4);
		    for (i = 2; i <= min(*n+1-j,*k+1); ++i) {
/* Computing MAX */
			d__1 = value, d__2 = z_abs(&AB(i,j));
			value = max(d__1,d__2);
/* L30: */
		    }
/* L40: */
		}
	    }
	} else {
	    value = 0.;
	    if (lsame_(uplo, "U")) {
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
/* Computing MAX */
		    i__3 = *k + 2 - j;
		    i__2 = *k + 1;
		    for (i = max(*k+2-j,1); i <= *k+1; ++i) {
/* Computing MAX */
			d__1 = value, d__2 = z_abs(&AB(i,j));
			value = max(d__1,d__2);
/* L50: */
		    }
/* L60: */
		}
	    } else {
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
/* Computing MIN */
		    i__3 = *n + 1 - j, i__4 = *k + 1;
		    i__2 = min(i__3,i__4);
		    for (i = 1; i <= min(*n+1-j,*k+1); ++i) {
/* Computing MAX */
			d__1 = value, d__2 = z_abs(&AB(i,j));
			value = max(d__1,d__2);
/* L70: */
		    }
/* L80: */
		}
	    }
	}
    } else if (lsame_(norm, "O") || *(unsigned char *)norm == '1') {

/*        Find norm1(A). */

	value = 0.;
	udiag = lsame_(diag, "U");
	if (lsame_(uplo, "U")) {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		if (udiag) {
		    sum = 1.;
/* Computing MAX */
		    i__2 = *k + 2 - j;
		    i__3 = *k;
		    for (i = max(*k+2-j,1); i <= *k; ++i) {
			sum += z_abs(&AB(i,j));
/* L90: */
		    }
		} else {
		    sum = 0.;
/* Computing MAX */
		    i__3 = *k + 2 - j;
		    i__2 = *k + 1;
		    for (i = max(*k+2-j,1); i <= *k+1; ++i) {
			sum += z_abs(&AB(i,j));
/* L100: */
		    }
		}
		value = max(value,sum);
/* L110: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		if (udiag) {
		    sum = 1.;
/* Computing MIN */
		    i__3 = *n + 1 - j, i__4 = *k + 1;
		    i__2 = min(i__3,i__4);
		    for (i = 2; i <= min(*n+1-j,*k+1); ++i) {
			sum += z_abs(&AB(i,j));
/* L120: */
		    }
		} else {
		    sum = 0.;
/* Computing MIN */
		    i__3 = *n + 1 - j, i__4 = *k + 1;
		    i__2 = min(i__3,i__4);
		    for (i = 1; i <= min(*n+1-j,*k+1); ++i) {
			sum += z_abs(&AB(i,j));
/* L130: */
		    }
		}
		value = max(value,sum);
/* L140: */
	    }
	}
    } else if (lsame_(norm, "I")) {

/*        Find normI(A). */

	value = 0.;
	if (lsame_(uplo, "U")) {
	    if (lsame_(diag, "U")) {
		i__1 = *n;
		for (i = 1; i <= *n; ++i) {
		    WORK(i) = 1.;
/* L150: */
		}
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
		    l = *k + 1 - j;
/* Computing MAX */
		    i__2 = 1, i__3 = j - *k;
		    i__4 = j - 1;
		    for (i = max(1,j-*k); i <= j-1; ++i) {
			WORK(i) += z_abs(&AB(l+i,j));
/* L160: */
		    }
/* L170: */
		}
	    } else {
		i__1 = *n;
		for (i = 1; i <= *n; ++i) {
		    WORK(i) = 0.;
/* L180: */
		}
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
		    l = *k + 1 - j;
/* Computing MAX */
		    i__4 = 1, i__2 = j - *k;
		    i__3 = j;
		    for (i = max(1,j-*k); i <= j; ++i) {
			WORK(i) += z_abs(&AB(l+i,j));
/* L190: */
		    }
/* L200: */
		}
	    }
	} else {
	    if (lsame_(diag, "U")) {
		i__1 = *n;
		for (i = 1; i <= *n; ++i) {
		    WORK(i) = 1.;
/* L210: */
		}
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
		    l = 1 - j;
/* Computing MIN */
		    i__4 = *n, i__2 = j + *k;
		    i__3 = min(i__4,i__2);
		    for (i = j + 1; i <= min(*n,j+*k); ++i) {
			WORK(i) += z_abs(&AB(l+i,j));
/* L220: */
		    }
/* L230: */
		}
	    } else {
		i__1 = *n;
		for (i = 1; i <= *n; ++i) {
		    WORK(i) = 0.;
/* L240: */
		}
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
		    l = 1 - j;
/* Computing MIN */
		    i__4 = *n, i__2 = j + *k;
		    i__3 = min(i__4,i__2);
		    for (i = j; i <= min(*n,j+*k); ++i) {
			WORK(i) += z_abs(&AB(l+i,j));
/* L250: */
		    }
/* L260: */
		}
	    }
	}
	i__1 = *n;
	for (i = 1; i <= *n; ++i) {
/* Computing MAX */
	    d__1 = value, d__2 = WORK(i);
	    value = max(d__1,d__2);
/* L270: */
	}
    } else if (lsame_(norm, "F") || lsame_(norm, "E")) {

/*        Find normF(A). */

	if (lsame_(uplo, "U")) {
	    if (lsame_(diag, "U")) {
		scale = 1.;
		sum = (doublereal) (*n);
		if (*k > 0) {
		    i__1 = *n;
		    for (j = 2; j <= *n; ++j) {
/* Computing MIN */
			i__4 = j - 1;
			i__3 = min(i__4,*k);
/* Computing MAX */
			i__2 = *k + 2 - j;
			zlassq_(&i__3, &AB(max(*k+2-j,1),j), &c__1, 
				&scale, &sum);
/* L280: */
		    }
		}
	    } else {
		scale = 0.;
		sum = 1.;
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
/* Computing MIN */
		    i__4 = j, i__2 = *k + 1;
		    i__3 = min(i__4,i__2);
/* Computing MAX */
		    i__5 = *k + 2 - j;
		    zlassq_(&i__3, &AB(max(*k+2-j,1),j), &c__1, &
			    scale, &sum);
/* L290: */
		}
	    }
	} else {
	    if (lsame_(diag, "U")) {
		scale = 1.;
		sum = (doublereal) (*n);
		if (*k > 0) {
		    i__1 = *n - 1;
		    for (j = 1; j <= *n-1; ++j) {
/* Computing MIN */
			i__4 = *n - j;
			i__3 = min(i__4,*k);
			zlassq_(&i__3, &AB(2,j), &c__1, &scale, &
				sum);
/* L300: */
		    }
		}
	    } else {
		scale = 0.;
		sum = 1.;
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
/* Computing MIN */
		    i__4 = *n - j + 1, i__2 = *k + 1;
		    i__3 = min(i__4,i__2);
		    zlassq_(&i__3, &AB(1,j), &c__1, &scale, &sum);
/* L310: */
		}
	    }
	}
	value = scale * sqrt(sum);
    }

    ret_val = value;
    return ret_val;

/*     End of ZLANTB */

} /* zlantb_ */

