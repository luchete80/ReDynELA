#include "blaswrap.h"
#include "f2c.h"

/* Subroutine */ int dtgsyl_(char *trans, integer *ijob, integer *m, integer *
	n, doublereal *a, integer *lda, doublereal *b, integer *ldb, 
	doublereal *c__, integer *ldc, doublereal *d__, integer *ldd, 
	doublereal *e, integer *lde, doublereal *f, integer *ldf, doublereal *
	scale, doublereal *dif, doublereal *work, integer *lwork, integer *
	iwork, integer *info)
{
/*  -- LAPACK routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       June 30, 1999   


    Purpose   
    =======   

    DTGSYL solves the generalized Sylvester equation:   

                A * R - L * B = scale * C                 (1)   
                D * R - L * E = scale * F   

    where R and L are unknown m-by-n matrices, (A, D), (B, E) and   
    (C, F) are given matrix pairs of size m-by-m, n-by-n and m-by-n,   
    respectively, with real entries. (A, D) and (B, E) must be in   
    generalized (real) Schur canonical form, i.e. A, B are upper quasi   
    triangular and D, E are upper triangular.   

    The solution (R, L) overwrites (C, F). 0 <= SCALE <= 1 is an output   
    scaling factor chosen to avoid overflow.   

    In matrix notation (1) is equivalent to solve  Zx = scale b, where   
    Z is defined as   

               Z = [ kron(In, A)  -kron(B', Im) ]         (2)   
                   [ kron(In, D)  -kron(E', Im) ].   

    Here Ik is the identity matrix of size k and X' is the transpose of   
    X. kron(X, Y) is the Kronecker product between the matrices X and Y.   

    If TRANS = 'T', DTGSYL solves the transposed system Z'*y = scale*b,   
    which is equivalent to solve for R and L in   

                A' * R  + D' * L   = scale *  C           (3)   
                R  * B' + L  * E'  = scale * (-F)   

    This case (TRANS = 'T') is used to compute an one-norm-based estimate   
    of Dif[(A,D), (B,E)], the separation between the matrix pairs (A,D)   
    and (B,E), using DLACON.   

    If IJOB >= 1, DTGSYL computes a Frobenius norm-based estimate   
    of Dif[(A,D),(B,E)]. That is, the reciprocal of a lower bound on the   
    reciprocal of the smallest singular value of Z. See [1-2] for more   
    information.   

    This is a level 3 BLAS algorithm.   

    Arguments   
    =========   

    TRANS   (input) CHARACTER*1   
            = 'N', solve the generalized Sylvester equation (1).   
            = 'T', solve the 'transposed' system (3).   

    IJOB    (input) INTEGER   
            Specifies what kind of functionality to be performed.   
             =0: solve (1) only.   
             =1: The functionality of 0 and 3.   
             =2: The functionality of 0 and 4.   
             =3: Only an estimate of Dif[(A,D), (B,E)] is computed.   
                 (look ahead strategy IJOB  = 1 is used).   
             =4: Only an estimate of Dif[(A,D), (B,E)] is computed.   
                 ( DGECON on sub-systems is used ).   
            Not referenced if TRANS = 'T'.   

    M       (input) INTEGER   
            The order of the matrices A and D, and the row dimension of   
            the matrices C, F, R and L.   

    N       (input) INTEGER   
            The order of the matrices B and E, and the column dimension   
            of the matrices C, F, R and L.   

    A       (input) DOUBLE PRECISION array, dimension (LDA, M)   
            The upper quasi triangular matrix A.   

    LDA     (input) INTEGER   
            The leading dimension of the array A. LDA >= max(1, M).   

    B       (input) DOUBLE PRECISION array, dimension (LDB, N)   
            The upper quasi triangular matrix B.   

    LDB     (input) INTEGER   
            The leading dimension of the array B. LDB >= max(1, N).   

    C       (input/output) DOUBLE PRECISION array, dimension (LDC, N)   
            On entry, C contains the right-hand-side of the first matrix   
            equation in (1) or (3).   
            On exit, if IJOB = 0, 1 or 2, C has been overwritten by   
            the solution R. If IJOB = 3 or 4 and TRANS = 'N', C holds R,   
            the solution achieved during the computation of the   
            Dif-estimate.   

    LDC     (input) INTEGER   
            The leading dimension of the array C. LDC >= max(1, M).   

    D       (input) DOUBLE PRECISION array, dimension (LDD, M)   
            The upper triangular matrix D.   

    LDD     (input) INTEGER   
            The leading dimension of the array D. LDD >= max(1, M).   

    E       (input) DOUBLE PRECISION array, dimension (LDE, N)   
            The upper triangular matrix E.   

    LDE     (input) INTEGER   
            The leading dimension of the array E. LDE >= max(1, N).   

    F       (input/output) DOUBLE PRECISION array, dimension (LDF, N)   
            On entry, F contains the right-hand-side of the second matrix   
            equation in (1) or (3).   
            On exit, if IJOB = 0, 1 or 2, F has been overwritten by   
            the solution L. If IJOB = 3 or 4 and TRANS = 'N', F holds L,   
            the solution achieved during the computation of the   
            Dif-estimate.   

    LDF     (input) INTEGER   
            The leading dimension of the array F. LDF >= max(1, M).   

    DIF     (output) DOUBLE PRECISION   
            On exit DIF is the reciprocal of a lower bound of the   
            reciprocal of the Dif-function, i.e. DIF is an upper bound of   
            Dif[(A,D), (B,E)] = sigma_min(Z), where Z as in (2).   
            IF IJOB = 0 or TRANS = 'T', DIF is not touched.   

    SCALE   (output) DOUBLE PRECISION   
            On exit SCALE is the scaling factor in (1) or (3).   
            If 0 < SCALE < 1, C and F hold the solutions R and L, resp.,   
            to a slightly perturbed system but the input matrices A, B, D   
            and E have not been changed. If SCALE = 0, C and F hold the   
            solutions R and L, respectively, to the homogeneous system   
            with C = F = 0. Normally, SCALE = 1.   

    WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)   
            If IJOB = 0, WORK is not referenced.  Otherwise,   
            on exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The dimension of the array WORK. LWORK > = 1.   
            If IJOB = 1 or 2 and TRANS = 'N', LWORK >= 2*M*N.   

            If LWORK = -1, then a workspace query is assumed; the routine   
            only calculates the optimal size of the WORK array, returns   
            this value as the first entry of the WORK array, and no error   
            message related to LWORK is issued by XERBLA.   

    IWORK   (workspace) INTEGER array, dimension (M+N+6)   

    INFO    (output) INTEGER   
              =0: successful exit   
              <0: If INFO = -i, the i-th argument had an illegal value.   
              >0: (A, D) and (B, E) have common or close eigenvalues.   

    Further Details   
    ===============   

    Based on contributions by   
       Bo Kagstrom and Peter Poromaa, Department of Computing Science,   
       Umea University, S-901 87 Umea, Sweden.   

    [1] B. Kagstrom and P. Poromaa, LAPACK-Style Algorithms and Software   
        for Solving the Generalized Sylvester Equation and Estimating the   
        Separation between Regular Matrix Pairs, Report UMINF - 93.23,   
        Department of Computing Science, Umea University, S-901 87 Umea,   
        Sweden, December 1993, Revised April 1994, Also as LAPACK Working   
        Note 75.  To appear in ACM Trans. on Math. Software, Vol 22,   
        No 1, 1996.   

    [2] B. Kagstrom, A Perturbation Analysis of the Generalized Sylvester   
        Equation (AR - LB, DR - LE ) = (C, F), SIAM J. Matrix Anal.   
        Appl., 15(4):1045-1060, 1994   

    [3] B. Kagstrom and L. Westin, Generalized Schur Methods with   
        Condition Estimators for Solving the Generalized Sylvester   
        Equation, IEEE Transactions on Automatic Control, Vol. 34, No. 7,   
        July 1989, pp 745-751.   

    =====================================================================   


       Decode and test input parameters   

       Parameter adjustments */
    /* Table of constant values */
    static integer c__2 = 2;
    static integer c_n1 = -1;
    static integer c__5 = 5;
    static doublereal c_b14 = 0.;
    static integer c__0 = 0;
    static integer c__1 = 1;
    static doublereal c_b53 = -1.;
    static doublereal c_b54 = 1.;
    
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, d_dim1, 
	    d_offset, e_dim1, e_offset, f_dim1, f_offset, i__1, i__2, i__3, 
	    i__4;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    static doublereal dsum;
    static integer ppqq, i__, j, k, p, q;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), dgemm_(char *, char *, integer *, integer *, integer *
	    , doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *);
    extern logical lsame_(char *, char *);
    static integer ifunc, linfo;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer lwmin;
    static doublereal scale2;
    extern /* Subroutine */ int dtgsy2_(char *, integer *, integer *, integer 
	    *, doublereal *, integer *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     integer *, integer *, integer *);
    static integer ie, je, mb, nb;
    static doublereal dscale;
    static integer is, js, pq;
    static doublereal scaloc;
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *), 
	    xerbla_(char *, integer *);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static integer iround;
    static logical notran;
    static integer isolve;
    static logical lquery;
#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]
#define b_ref(a_1,a_2) b[(a_2)*b_dim1 + a_1]
#define c___ref(a_1,a_2) c__[(a_2)*c_dim1 + a_1]
#define d___ref(a_1,a_2) d__[(a_2)*d_dim1 + a_1]
#define e_ref(a_1,a_2) e[(a_2)*e_dim1 + a_1]
#define f_ref(a_1,a_2) f[(a_2)*f_dim1 + a_1]


    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1 * 1;
    b -= b_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1 * 1;
    c__ -= c_offset;
    d_dim1 = *ldd;
    d_offset = 1 + d_dim1 * 1;
    d__ -= d_offset;
    e_dim1 = *lde;
    e_offset = 1 + e_dim1 * 1;
    e -= e_offset;
    f_dim1 = *ldf;
    f_offset = 1 + f_dim1 * 1;
    f -= f_offset;
    --work;
    --iwork;

    /* Function Body */
    *info = 0;
    notran = lsame_(trans, "N");
    lquery = *lwork == -1;

    if ((*ijob == 1 || *ijob == 2) && notran) {
/* Computing MAX */
	i__1 = 1, i__2 = (*m << 1) * *n;
	lwmin = max(i__1,i__2);
    } else {
	lwmin = 1;
    }

    if (! notran && ! lsame_(trans, "T")) {
	*info = -1;
    } else if (*ijob < 0 || *ijob > 4) {
	*info = -2;
    } else if (*m <= 0) {
	*info = -3;
    } else if (*n <= 0) {
	*info = -4;
    } else if (*lda < max(1,*m)) {
	*info = -6;
    } else if (*ldb < max(1,*n)) {
	*info = -8;
    } else if (*ldc < max(1,*m)) {
	*info = -10;
    } else if (*ldd < max(1,*m)) {
	*info = -12;
    } else if (*lde < max(1,*n)) {
	*info = -14;
    } else if (*ldf < max(1,*m)) {
	*info = -16;
    } else if (*lwork < lwmin && ! lquery) {
	*info = -20;
    }

    if (*info == 0) {
	work[1] = (doublereal) lwmin;
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DTGSYL", &i__1);
	return 0;
    } else if (lquery) {
	return 0;
    }

/*     Determine optimal block sizes MB and NB */

    mb = ilaenv_(&c__2, "DTGSYL", trans, m, n, &c_n1, &c_n1, (ftnlen)6, (
	    ftnlen)1);
    nb = ilaenv_(&c__5, "DTGSYL", trans, m, n, &c_n1, &c_n1, (ftnlen)6, (
	    ftnlen)1);

    isolve = 1;
    ifunc = 0;
    if (*ijob >= 3 && notran) {
	ifunc = *ijob - 2;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    dcopy_(m, &c_b14, &c__0, &c___ref(1, j), &c__1);
	    dcopy_(m, &c_b14, &c__0, &f_ref(1, j), &c__1);
/* L10: */
	}
    } else if (*ijob >= 1 && notran) {
	isolve = 2;
    }

    if (mb <= 1 && nb <= 1 || mb >= *m && nb >= *n) {

	i__1 = isolve;
	for (iround = 1; iround <= i__1; ++iround) {

/*           Use unblocked Level 2 solver */

	    dscale = 0.;
	    dsum = 1.;
	    pq = 0;
	    dtgsy2_(trans, &ifunc, m, n, &a[a_offset], lda, &b[b_offset], ldb,
		     &c__[c_offset], ldc, &d__[d_offset], ldd, &e[e_offset], 
		    lde, &f[f_offset], ldf, scale, &dsum, &dscale, &iwork[1], 
		    &pq, info);
	    if (dscale != 0.) {
		if (*ijob == 1 || *ijob == 3) {
		    *dif = sqrt((doublereal) ((*m << 1) * *n)) / (dscale * 
			    sqrt(dsum));
		} else {
		    *dif = sqrt((doublereal) pq) / (dscale * sqrt(dsum));
		}
	    }

	    if (isolve == 2 && iround == 1) {
		ifunc = *ijob;
		scale2 = *scale;
		dlacpy_("F", m, n, &c__[c_offset], ldc, &work[1], m);
		dlacpy_("F", m, n, &f[f_offset], ldf, &work[*m * *n + 1], m);
		i__2 = *n;
		for (j = 1; j <= i__2; ++j) {
		    dcopy_(m, &c_b14, &c__0, &c___ref(1, j), &c__1);
		    dcopy_(m, &c_b14, &c__0, &f_ref(1, j), &c__1);
/* L20: */
		}
	    } else if (isolve == 2 && iround == 2) {
		dlacpy_("F", m, n, &work[1], m, &c__[c_offset], ldc);
		dlacpy_("F", m, n, &work[*m * *n + 1], m, &f[f_offset], ldf);
		*scale = scale2;
	    }
/* L30: */
	}

	return 0;
    }

/*     Determine block structure of A */

    p = 0;
    i__ = 1;
L40:
    if (i__ > *m) {
	goto L50;
    }
    ++p;
    iwork[p] = i__;
    i__ += mb;
    if (i__ >= *m) {
	goto L50;
    }
    if (a_ref(i__, i__ - 1) != 0.) {
	++i__;
    }
    goto L40;
L50:

    iwork[p + 1] = *m + 1;
    if (iwork[p] == iwork[p + 1]) {
	--p;
    }

/*     Determine block structure of B */

    q = p + 1;
    j = 1;
L60:
    if (j > *n) {
	goto L70;
    }
    ++q;
    iwork[q] = j;
    j += nb;
    if (j >= *n) {
	goto L70;
    }
    if (b_ref(j, j - 1) != 0.) {
	++j;
    }
    goto L60;
L70:

    iwork[q + 1] = *n + 1;
    if (iwork[q] == iwork[q + 1]) {
	--q;
    }

    if (notran) {

	i__1 = isolve;
	for (iround = 1; iround <= i__1; ++iround) {

/*           Solve (I, J)-subsystem   
                 A(I, I) * R(I, J) - L(I, J) * B(J, J) = C(I, J)   
                 D(I, I) * R(I, J) - L(I, J) * E(J, J) = F(I, J)   
             for I = P, P - 1,..., 1; J = 1, 2,..., Q */

	    dscale = 0.;
	    dsum = 1.;
	    pq = 0;
	    *scale = 1.;
	    i__2 = q;
	    for (j = p + 2; j <= i__2; ++j) {
		js = iwork[j];
		je = iwork[j + 1] - 1;
		nb = je - js + 1;
		for (i__ = p; i__ >= 1; --i__) {
		    is = iwork[i__];
		    ie = iwork[i__ + 1] - 1;
		    mb = ie - is + 1;
		    ppqq = 0;
		    dtgsy2_(trans, &ifunc, &mb, &nb, &a_ref(is, is), lda, &
			    b_ref(js, js), ldb, &c___ref(is, js), ldc, &
			    d___ref(is, is), ldd, &e_ref(js, js), lde, &f_ref(
			    is, js), ldf, &scaloc, &dsum, &dscale, &iwork[q + 
			    2], &ppqq, &linfo);
		    if (linfo > 0) {
			*info = linfo;
		    }

		    pq += ppqq;
		    if (scaloc != 1.) {
			i__3 = js - 1;
			for (k = 1; k <= i__3; ++k) {
			    dscal_(m, &scaloc, &c___ref(1, k), &c__1);
			    dscal_(m, &scaloc, &f_ref(1, k), &c__1);
/* L80: */
			}
			i__3 = je;
			for (k = js; k <= i__3; ++k) {
			    i__4 = is - 1;
			    dscal_(&i__4, &scaloc, &c___ref(1, k), &c__1);
			    i__4 = is - 1;
			    dscal_(&i__4, &scaloc, &f_ref(1, k), &c__1);
/* L90: */
			}
			i__3 = je;
			for (k = js; k <= i__3; ++k) {
			    i__4 = *m - ie;
			    dscal_(&i__4, &scaloc, &c___ref(ie + 1, k), &c__1)
				    ;
			    i__4 = *m - ie;
			    dscal_(&i__4, &scaloc, &f_ref(ie + 1, k), &c__1);
/* L100: */
			}
			i__3 = *n;
			for (k = je + 1; k <= i__3; ++k) {
			    dscal_(m, &scaloc, &c___ref(1, k), &c__1);
			    dscal_(m, &scaloc, &f_ref(1, k), &c__1);
/* L110: */
			}
			*scale *= scaloc;
		    }

/*                 Substitute R(I, J) and L(I, J) into remaining   
                   equation. */

		    if (i__ > 1) {
			i__3 = is - 1;
			dgemm_("N", "N", &i__3, &nb, &mb, &c_b53, &a_ref(1, 
				is), lda, &c___ref(is, js), ldc, &c_b54, &
				c___ref(1, js), ldc);
			i__3 = is - 1;
			dgemm_("N", "N", &i__3, &nb, &mb, &c_b53, &d___ref(1, 
				is), ldd, &c___ref(is, js), ldc, &c_b54, &
				f_ref(1, js), ldf);
		    }
		    if (j < q) {
			i__3 = *n - je;
			dgemm_("N", "N", &mb, &i__3, &nb, &c_b54, &f_ref(is, 
				js), ldf, &b_ref(js, je + 1), ldb, &c_b54, &
				c___ref(is, je + 1), ldc);
			i__3 = *n - je;
			dgemm_("N", "N", &mb, &i__3, &nb, &c_b54, &f_ref(is, 
				js), ldf, &e_ref(js, je + 1), lde, &c_b54, &
				f_ref(is, je + 1), ldf);
		    }
/* L120: */
		}
/* L130: */
	    }
	    if (dscale != 0.) {
		if (*ijob == 1 || *ijob == 3) {
		    *dif = sqrt((doublereal) ((*m << 1) * *n)) / (dscale * 
			    sqrt(dsum));
		} else {
		    *dif = sqrt((doublereal) pq) / (dscale * sqrt(dsum));
		}
	    }
	    if (isolve == 2 && iround == 1) {
		ifunc = *ijob;
		scale2 = *scale;
		dlacpy_("F", m, n, &c__[c_offset], ldc, &work[1], m);
		dlacpy_("F", m, n, &f[f_offset], ldf, &work[*m * *n + 1], m);
		i__2 = *n;
		for (j = 1; j <= i__2; ++j) {
		    dcopy_(m, &c_b14, &c__0, &c___ref(1, j), &c__1);
		    dcopy_(m, &c_b14, &c__0, &f_ref(1, j), &c__1);
/* L140: */
		}
	    } else if (isolve == 2 && iround == 2) {
		dlacpy_("F", m, n, &work[1], m, &c__[c_offset], ldc);
		dlacpy_("F", m, n, &work[*m * *n + 1], m, &f[f_offset], ldf);
		*scale = scale2;
	    }
/* L150: */
	}

    } else {

/*        Solve transposed (I, J)-subsystem   
               A(I, I)' * R(I, J)  + D(I, I)' * L(I, J)  =  C(I, J)   
               R(I, J)  * B(J, J)' + L(I, J)  * E(J, J)' = -F(I, J)   
          for I = 1,2,..., P; J = Q, Q-1,..., 1 */

	*scale = 1.;
	i__1 = p;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    is = iwork[i__];
	    ie = iwork[i__ + 1] - 1;
	    mb = ie - is + 1;
	    i__2 = p + 2;
	    for (j = q; j >= i__2; --j) {
		js = iwork[j];
		je = iwork[j + 1] - 1;
		nb = je - js + 1;
		dtgsy2_(trans, &ifunc, &mb, &nb, &a_ref(is, is), lda, &b_ref(
			js, js), ldb, &c___ref(is, js), ldc, &d___ref(is, is),
			 ldd, &e_ref(js, js), lde, &f_ref(is, js), ldf, &
			scaloc, &dsum, &dscale, &iwork[q + 2], &ppqq, &linfo);
		if (linfo > 0) {
		    *info = linfo;
		}
		if (scaloc != 1.) {
		    i__3 = js - 1;
		    for (k = 1; k <= i__3; ++k) {
			dscal_(m, &scaloc, &c___ref(1, k), &c__1);
			dscal_(m, &scaloc, &f_ref(1, k), &c__1);
/* L160: */
		    }
		    i__3 = je;
		    for (k = js; k <= i__3; ++k) {
			i__4 = is - 1;
			dscal_(&i__4, &scaloc, &c___ref(1, k), &c__1);
			i__4 = is - 1;
			dscal_(&i__4, &scaloc, &f_ref(1, k), &c__1);
/* L170: */
		    }
		    i__3 = je;
		    for (k = js; k <= i__3; ++k) {
			i__4 = *m - ie;
			dscal_(&i__4, &scaloc, &c___ref(ie + 1, k), &c__1);
			i__4 = *m - ie;
			dscal_(&i__4, &scaloc, &f_ref(ie + 1, k), &c__1);
/* L180: */
		    }
		    i__3 = *n;
		    for (k = je + 1; k <= i__3; ++k) {
			dscal_(m, &scaloc, &c___ref(1, k), &c__1);
			dscal_(m, &scaloc, &f_ref(1, k), &c__1);
/* L190: */
		    }
		    *scale *= scaloc;
		}

/*              Substitute R(I, J) and L(I, J) into remaining equation. */

		if (j > p + 2) {
		    i__3 = js - 1;
		    dgemm_("N", "T", &mb, &i__3, &nb, &c_b54, &c___ref(is, js)
			    , ldc, &b_ref(1, js), ldb, &c_b54, &f_ref(is, 1), 
			    ldf);
		    i__3 = js - 1;
		    dgemm_("N", "T", &mb, &i__3, &nb, &c_b54, &f_ref(is, js), 
			    ldf, &e_ref(1, js), lde, &c_b54, &f_ref(is, 1), 
			    ldf);
		}
		if (i__ < p) {
		    i__3 = *m - ie;
		    dgemm_("T", "N", &i__3, &nb, &mb, &c_b53, &a_ref(is, ie + 
			    1), lda, &c___ref(is, js), ldc, &c_b54, &c___ref(
			    ie + 1, js), ldc);
		    i__3 = *m - ie;
		    dgemm_("T", "N", &i__3, &nb, &mb, &c_b53, &d___ref(is, ie 
			    + 1), ldd, &f_ref(is, js), ldf, &c_b54, &c___ref(
			    ie + 1, js), ldc);
		}
/* L200: */
	    }
/* L210: */
	}

    }

    work[1] = (doublereal) lwmin;

    return 0;

/*     End of DTGSYL */

} /* dtgsyl_ */

#undef f_ref
#undef e_ref
#undef d___ref
#undef c___ref
#undef b_ref
#undef a_ref


