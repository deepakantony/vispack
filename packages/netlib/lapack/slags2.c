#include "f2c.h"

/* Subroutine */ int slags2_(logical *upper, real *a1, real *a2, real *a3, 
	real *b1, real *b2, real *b3, real *csu, real *snu, real *csv, real *
	snv, real *csq, real *snq)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    SLAGS2 computes 2-by-2 orthogonal matrices U, V and Q, such   
    that if ( UPPER ) then   

              U'*A*Q = U'*( A1 A2 )*Q = ( x  0  )   
                          ( 0  A3 )     ( x  x  )   
    and   
              V'*B*Q = V'*( B1 B2 )*Q = ( x  0  )   
                          ( 0  B3 )     ( x  x  )   

    or if ( .NOT.UPPER ) then   

              U'*A*Q = U'*( A1 0  )*Q = ( x  x  )   
                          ( A2 A3 )     ( 0  x  )   
    and   
              V'*B*Q = V'*( B1 0  )*Q = ( x  x  )   
                          ( B2 B3 )     ( 0  x  )   

    The rows of the transformed A and B are parallel, where   

      U = (  CSU  SNU ), V = (  CSV SNV ), Q = (  CSQ   SNQ )   
          ( -SNU  CSU )      ( -SNV CSV )      ( -SNQ   CSQ )   

    Z' denotes the transpose of Z.   


    Arguments   
    =========   

    UPPER   (input) LOGICAL   
            = .TRUE.: the input matrices A and B are upper triangular.   
            = .FALSE.: the input matrices A and B are lower triangular.   

    A1      (input) REAL   
    A2      (input) REAL   
    A3      (input) REAL   
            On entry, A1, A2 and A3 are elements of the input 2-by-2   
            upper (lower) triangular matrix A.   

    B1      (input) REAL   
    B2      (input) REAL   
    B3      (input) REAL   
            On entry, B1, B2 and B3 are elements of the input 2-by-2   
            upper (lower) triangular matrix B.   

    CSU     (output) REAL   
    SNU     (output) REAL   
            The desired orthogonal matrix U.   

    CSV     (output) REAL   
    SNV     (output) REAL   
            The desired orthogonal matrix V.   

    CSQ     (output) REAL   
    SNQ     (output) REAL   
            The desired orthogonal matrix Q.   

    ===================================================================== 
*/
    /* System generated locals */
    real r__1;
    /* Local variables */
    static real aua11, aua12, aua21, aua22, avb11, avb12, avb21, avb22, ua11r,
	     ua22r, vb11r, vb22r, a, b, c, d, r, s1, s2;
    extern /* Subroutine */ int slasv2_(real *, real *, real *, real *, real *
	    , real *, real *, real *, real *), slartg_(real *, real *, real *,
	     real *, real *);
    static real ua11, ua12, ua21, ua22, vb11, vb12, vb21, vb22, csl, csr, snl,
	     snr;



    if (*upper) {

/*        Input matrices A and B are upper triangular matrices   

          Form matrix C = A*adj(B) = ( a b )   
                                     ( 0 d ) */

	a = *a1 * *b3;
	d = *a3 * *b1;
	b = *a2 * *b1 - *a1 * *b2;

/*        The SVD of real 2-by-2 triangular C   

           ( CSL -SNL )*( A B )*(  CSR  SNR ) = ( R 0 )   
           ( SNL  CSL ) ( 0 D ) ( -SNR  CSR )   ( 0 T ) */

	slasv2_(&a, &b, &d, &s1, &s2, &snr, &csr, &snl, &csl);

	if (dabs(csl) >= dabs(snl) || dabs(csr) >= dabs(snr)) {

/*           Compute the (1,1) and (1,2) elements of U'*A and V'*B
,   
             and (1,2) element of |U|'*|A| and |V|'*|B|. */

	    ua11r = csl * *a1;
	    ua12 = csl * *a2 + snl * *a3;

	    vb11r = csr * *b1;
	    vb12 = csr * *b2 + snr * *b3;

	    aua12 = dabs(csl) * dabs(*a2) + dabs(snl) * dabs(*a3);
	    avb12 = dabs(csr) * dabs(*b2) + dabs(snr) * dabs(*b3);

/*           zero (1,2) elements of U'*A and V'*B */

	    if (dabs(ua11r) + dabs(ua12) != 0.f) {
		if (aua12 / (dabs(ua11r) + dabs(ua12)) <= avb12 / (dabs(vb11r)
			 + dabs(vb12))) {
		    r__1 = -(doublereal)ua11r;
		    slartg_(&r__1, &ua12, csq, snq, &r);
		} else {
		    r__1 = -(doublereal)vb11r;
		    slartg_(&r__1, &vb12, csq, snq, &r);
		}
	    } else {
		r__1 = -(doublereal)vb11r;
		slartg_(&r__1, &vb12, csq, snq, &r);
	    }

	    *csu = csl;
	    *snu = -(doublereal)snl;
	    *csv = csr;
	    *snv = -(doublereal)snr;

	} else {

/*           Compute the (2,1) and (2,2) elements of U'*A and V'*B
,   
             and (2,2) element of |U|'*|A| and |V|'*|B|. */

	    ua21 = -(doublereal)snl * *a1;
	    ua22 = -(doublereal)snl * *a2 + csl * *a3;

	    vb21 = -(doublereal)snr * *b1;
	    vb22 = -(doublereal)snr * *b2 + csr * *b3;

	    aua22 = dabs(snl) * dabs(*a2) + dabs(csl) * dabs(*a3);
	    avb22 = dabs(snr) * dabs(*b2) + dabs(csr) * dabs(*b3);

/*           zero (2,2) elements of U'*A and V'*B, and then swap. 
*/

	    if (dabs(ua21) + dabs(ua22) != 0.f) {
		if (aua22 / (dabs(ua21) + dabs(ua22)) <= avb22 / (dabs(vb21) 
			+ dabs(vb22))) {
		    r__1 = -(doublereal)ua21;
		    slartg_(&r__1, &ua22, csq, snq, &r);
		} else {
		    r__1 = -(doublereal)vb21;
		    slartg_(&r__1, &vb22, csq, snq, &r);
		}
	    } else {
		r__1 = -(doublereal)vb21;
		slartg_(&r__1, &vb22, csq, snq, &r);
	    }

	    *csu = snl;
	    *snu = csl;
	    *csv = snr;
	    *snv = csr;

	}

    } else {

/*        Input matrices A and B are lower triangular matrices   

          Form matrix C = A*adj(B) = ( a 0 )   
                                     ( c d ) */

	a = *a1 * *b3;
	d = *a3 * *b1;
	c = *a2 * *b3 - *a3 * *b2;

/*        The SVD of real 2-by-2 triangular C   

           ( CSL -SNL )*( A 0 )*(  CSR  SNR ) = ( R 0 )   
           ( SNL  CSL ) ( C D ) ( -SNR  CSR )   ( 0 T ) */

	slasv2_(&a, &c, &d, &s1, &s2, &snr, &csr, &snl, &csl);

	if (dabs(csr) >= dabs(snr) || dabs(csl) >= dabs(snl)) {

/*           Compute the (2,1) and (2,2) elements of U'*A and V'*B
,   
             and (2,1) element of |U|'*|A| and |V|'*|B|. */

	    ua21 = -(doublereal)snr * *a1 + csr * *a2;
	    ua22r = csr * *a3;

	    vb21 = -(doublereal)snl * *b1 + csl * *b2;
	    vb22r = csl * *b3;

	    aua21 = dabs(snr) * dabs(*a1) + dabs(csr) * dabs(*a2);
	    avb21 = dabs(snl) * dabs(*b1) + dabs(csl) * dabs(*b2);

/*           zero (2,1) elements of U'*A and V'*B. */

	    if (dabs(ua21) + dabs(ua22r) != 0.f) {
		if (aua21 / (dabs(ua21) + dabs(ua22r)) <= avb21 / (dabs(vb21) 
			+ dabs(vb22r))) {
		    slartg_(&ua22r, &ua21, csq, snq, &r);
		} else {
		    slartg_(&vb22r, &vb21, csq, snq, &r);
		}
	    } else {
		slartg_(&vb22r, &vb21, csq, snq, &r);
	    }

	    *csu = csr;
	    *snu = -(doublereal)snr;
	    *csv = csl;
	    *snv = -(doublereal)snl;

	} else {

/*           Compute the (1,1) and (1,2) elements of U'*A and V'*B
,   
             and (1,1) element of |U|'*|A| and |V|'*|B|. */

	    ua11 = csr * *a1 + snr * *a2;
	    ua12 = snr * *a3;

	    vb11 = csl * *b1 + snl * *b2;
	    vb12 = snl * *b3;

	    aua11 = dabs(csr) * dabs(*a1) + dabs(snr) * dabs(*a2);
	    avb11 = dabs(csl) * dabs(*b1) + dabs(snl) * dabs(*b2);

/*           zero (1,1) elements of U'*A and V'*B, and then swap. 
*/

	    if (dabs(ua11) + dabs(ua12) != 0.f) {
		if (aua11 / (dabs(ua11) + dabs(ua12)) <= avb11 / (dabs(vb11) 
			+ dabs(vb12))) {
		    slartg_(&ua12, &ua11, csq, snq, &r);
		} else {
		    slartg_(&vb12, &vb11, csq, snq, &r);
		}
	    } else {
		slartg_(&vb12, &vb11, csq, snq, &r);
	    }

	    *csu = snr;
	    *snu = csr;
	    *csv = snl;
	    *snv = csl;

	}

    }

    return 0;

/*     End of SLAGS2 */

} /* slags2_ */

